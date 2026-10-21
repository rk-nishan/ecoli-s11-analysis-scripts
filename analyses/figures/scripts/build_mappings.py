"""Rebuild auditable class/mechanism mappings from frozen primary-source metadata.

No network access, fuzzy gene matching, keyword-based mechanism assignment, or
hard-coded mechanism counts are used. Original Galaxy tables are read-only.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
import csv
import io
import re
import tarfile
from pathlib import Path

from common import ROOT, read_csv, sha256, verify_frozen_sources, write_csv, write_json

ARCHIVE = "references/ontology-v4.0.2.tar.bz2"
CLASSES = "references/resfinder_antibiotic_classes_65c9ba3.txt"
ALLOWED_RELATIONSHIPS = {"part_of", "participates_in"}


def normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.casefold())


def ontology_data():
    with tarfile.open(ROOT / ARCHIVE, "r:bz2") as archive:
        obo = archive.extractfile("./aro.obo").read().decode("utf-8")
        index = list(csv.DictReader(io.StringIO(archive.extractfile("./aro.tsv").read().decode("utf-8")), delimiter="\t"))
    terms = {}
    for block in obo.split("\n[Term]\n")[1:]:
        attributes = defaultdict(list)
        for line in block.splitlines():
            if ": " in line:
                key, value = line.split(": ", 1)
                attributes[key].append(value)
        if attributes["id"]:
            terms[attributes["id"][0]] = dict(attributes)
    return terms, index, obo.split("\n[Term]")[0]


def edges(term):
    for value in term.get("is_a", []):
        yield "is_a", value.split()[0]
    for value in term.get("relationship", []):
        relation, destination = value.split()[:2]
        if relation in ALLOWED_RELATIONSHIPS:
            yield relation, destination


def classify(accession, terms):
    mechanisms = {aid for aid, term in terms.items()
                  if "classified_as_resistance_mechanism" in term.get("property_value", [])}
    queue = deque([(accession, [])])
    seen, found = set(), {}
    efflux_regulator = False
    while queue:
        aid, path = queue.popleft()
        if aid in seen:
            continue
        seen.add(aid)
        term = terms.get(aid, {})
        if "classified_as_efflux_regulator" in term.get("property_value", []):
            efflux_regulator = True
        if aid in mechanisms:
            found[aid] = path
        for relation, destination in sorted(edges(term)):
            queue.append((destination, path + [(aid, relation, destination)]))
    return found, efflux_regulator


def build(outdir: Path):
    verify_frozen_sources()
    terms, index, header = ontology_data()
    card_rows = read_csv(ROOT / "data/raw/card_abricate.tsv", "\t")
    mappings = []
    for gene in sorted({row["GENE"] for row in card_rows}, key=str.casefold):
        matches = [row for row in index if normalize(gene) in
                   {normalize(row["Name"]), normalize(row["CARD Short Name"])}]
        if len(matches) != 1:
            raise ValueError(f"Expected one exact normalized ARO name match for {gene}; found {len(matches)}")
        match = matches[0]
        found, regulator = classify(match["Accession"], terms)
        if not found:
            raise ValueError(f"No supported mechanism path for {gene}; do not silently omit it")
        for mechanism_id, path in sorted(found.items()):
            mappings.append({
                "abricate_gene": gene, "aro_accession": match["Accession"], "aro_name": match["Name"],
                "card_short_name": match["CARD Short Name"], "mechanism_aro_accession": mechanism_id,
                "mechanism": terms[mechanism_id]["name"][0], "efflux_regulator": str(regulator).lower(),
                "ontology_path": " | ".join(f"{a} --{r}--> {b}" for a, r, b in path),
                "source_url": f"https://card.mcmaster.ca/ontology/{match['ID']}",
                "ontology_version": "4.0.2", "ontology_date": "2026-08-11",
            })
    write_csv(outdir / "data/mappings/card_mechanism_mapping.csv", mappings)

    drug_index = defaultdict(set)
    with (ROOT / CLASSES).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            group, *drugs = line.rstrip("\r\n").split("\t")
            for drug in drugs:
                drug_index[normalize(drug)].add(group)
    ref_rows = read_csv(ROOT / "data/raw/resfinder_abricate.tsv", "\t")
    drug_mappings = []
    for drug in sorted({d for row in ref_rows for d in row["RESISTANCE"].split(";") if d}, key=str.casefold):
        matches = drug_index[normalize(drug)]
        if len(matches) != 1:
            raise ValueError(f"Missing/ambiguous official drug class for {drug}: {matches}")
        drug_mappings.append({"resfinder_resistance_term": drug, "antimicrobial_class": next(iter(matches)),
                              "reference_file": CLASSES,
                              "source_url": "https://bitbucket.org/genomicepidemiology/resfinder_db/src/65c9ba3/antibiotic_classes.txt"})
    write_csv(outdir / "data/mappings/resfinder_drug_class_mapping.csv", drug_mappings)
    write_json(outdir / "provenance/mapping_method.json", {
        "card_ontology_version": "4.0.2", "obo_header": header,
        "ontology_archive_sha256": sha256(ROOT / ARCHIVE),
        "gene_match": "Unique exact alphanumeric/case-normalized match to ARO Name or CARD Short Name; no fuzzy matching",
        "mechanism_traversal": ["is_a", "part_of", "participates_in"],
        "excluded_edges": "Drug-resistance, regulates, has_part, and other relationships are not traversed",
        "mechanism_targets": "Terms explicitly tagged classified_as_resistance_mechanism",
        "path_policy": "One shortest supported path per gene/mechanism; all distinct mechanisms retained",
        "gene_count": len({r["abricate_gene"] for r in mappings}), "gene_mechanism_pairs": len(mappings),
        "resfinder_class_source": CLASSES, "resfinder_drug_terms": len(drug_mappings),
        "no_amr_detection_rerun": True,
    })
    print(f"Mapped {len(card_rows)} CARD hits to {len(mappings)} gene-mechanism pairs; {len(drug_mappings)} ResFinder drug terms")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--outdir", type=Path, default=ROOT)
    build(cli.parse_args().outdir)
