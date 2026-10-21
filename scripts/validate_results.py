"""Check saved QUAST, ANI and PGAP outputs against the exact S11 assembly."""

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from compare_fasta import read_fasta

ROOT = Path(__file__).resolve().parents[1]
ASSEMBLY_SHA256 = "406a7003fa816dd990fae93d2a9148659b35d2b56c21667c2322e3f65478372f"


def require(value, message):
    if not value:
        raise ValueError(message)


def normalize(identifier):
    return identifier.removeprefix("lcl|")


def protein_identifier(identifier):
    # NCBI serializes a general ID as gnl|database|id in FASTA and database:id in GenBank.
    fields = identifier.split("|", 2)
    return fields[1] + ":" + fields[2] if len(fields) == 3 and fields[0] == "gnl" else identifier


def metrics(records, minimum):
    lengths = sorted((len(s) for _, s in records if len(s) >= minimum), reverse=True)
    total = sum(lengths)
    result = {"contigs": len(lengths), "total": total, "largest": lengths[0]}
    for fraction in [50, 90]:
        cumulative = 0
        for number, length in enumerate(lengths, 1):
            cumulative += length
            if cumulative * 100 >= total * fraction:
                result[f"N{fraction}"] = length
                result[f"L{fraction}"] = number
                break
    return result


def verify(assembly, data, quast_report=None, pgap_directory=None):
    require(hashlib.sha256(assembly.read_bytes()).hexdigest() == ASSEMBLY_SHA256,
            "Assembly differs from the frozen S11 input")
    records = read_fasta(assembly)
    genome = {normalize(h.split()[0]): s for h, s in records}
    require(len(genome) == len(records), "Duplicate assembly identifiers")
    all_stats, selected = metrics(records, 0), metrics(records, 500)
    report_path = quast_report or data / "quast/reports/report.tsv"
    with report_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    observed = {r[0]: r[1] for r in rows[1:] if len(r) >= 2}
    expected = {"# contigs (>= 0 bp)": all_stats["contigs"], "Total length (>= 0 bp)": all_stats["total"],
                "# contigs": selected["contigs"], "Total length": selected["total"], "Largest contig": selected["largest"],
                **{k: selected[k] for k in ["N50", "L50", "N90", "L90"]}}
    for name, value in expected.items():
        require(int(observed[name].replace(",", "")) == value, f"QUAST metric mismatch: {name}")

    pgap = pgap_directory or data / "pgap/results/S11_PGAP"
    annotated_records = read_fasta(pgap / "annot.fna")
    annotated = {normalize(h.split()[0]): s for h, s in annotated_records}
    require(len(annotated) == len(annotated_records), "Duplicate annotated genomic identifiers")
    require(annotated == genome, "PGAP genomic FASTA differs from the source assembly")
    feature_counts, genes, pseudogenes = Counter(), set(), set()
    with (pgap / "annot.gff").open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == 9, "Malformed GFF row")
            contig, start, end = normalize(fields[0]), int(fields[3]), int(fields[4])
            require(contig in genome and 1 <= start <= end <= len(genome[contig]), "GFF coordinates outside source contig")
            feature_counts[fields[2]] += 1
            attrs = dict(p.split("=", 1) for p in fields[8].split(";") if "=" in p)
            if fields[2] in {"gene", "pseudogene"}:
                require("ID" in attrs, "Gene feature lacks an ID")
                target = genes if fields[2] == "gene" else pseudogenes
                require(attrs["ID"] not in target, "Duplicate gene feature")
                target.add(attrs["ID"])

    from Bio import SeqIO
    protein_records = list(SeqIO.parse(pgap / "annot.faa", "fasta"))
    proteins = {protein_identifier(r.id): str(r.seq) for r in protein_records}
    require(len(proteins) == len(protein_records), "Duplicate protein identifiers after namespace normalization")
    gb_genome, translations = {}, {}
    for record in SeqIO.parse(pgap / "annot.gbk", "genbank"):
        require(normalize(record.id) not in gb_genome, "Duplicate GenBank genomic identifier")
        gb_genome[normalize(record.id)] = str(record.seq).upper()
        require(record.annotations.get("organism") == "Escherichia coli", "Unexpected GenBank organism")
        for feature in record.features:
            if feature.type == "source":
                require(feature.qualifiers.get("strain") == ["S11"], "Unexpected strain metadata")
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                protein_id = protein_identifier(feature.qualifiers["protein_id"][0])
                require(protein_id not in translations, "Duplicate translated protein ID")
                translations[protein_id] = feature.qualifiers["translation"][0]
    require(gb_genome == genome, "GenBank sequence differs from the source assembly")
    require(translations == proteins, "GenBank translations differ from the protein FASTA")
    rna = sum(feature_counts[k] for k in ["RNase_P_RNA", "SRP_RNA", "antisense_RNA", "ncRNA", "rRNA", "tRNA", "tmRNA"])
    require(len(genes) == len(proteins) + rna, "Coding/RNA gene totals do not reconcile")
    require((len(genes), len(pseudogenes), len(proteins), rna) == (4809, 182, 4701, 108),
            "Annotation counts differ from the S11 release")

    ani = ET.parse(pgap / "ani-tax-report.xml").getroot()
    result = ani.find("results")
    require(result is not None and result.attrib.get("status") == "CONFIRMED"
            and result.attrib.get("confidence") == "HIGH", "ANI taxonomy was not confirmed with high confidence")
    submitted = result.find("submitted-taxid")
    require(submitted is not None and submitted.attrib.get("taxid") == "562", "Unexpected ANI input taxon")
    comparisons = [dict(x.attrib) for x in ani.findall("./data/subject")]
    require(len(comparisons) == 40, "ANI reference count differs from the saved run")
    return {"passed": True, "assembly_sha256": ASSEMBLY_SHA256,
            "all_contigs": all_stats, "quast_minimum_500bp": selected, "quast_metrics_checked": len(expected),
            "pgap_genomic_fasta_and_genbank_match_input": True, "gff_coordinates_valid": True,
            "protein_fasta_matches_genbank_translations": True,
            "protein_identifier_handling": "Bijective normalization of NCBI general-ID notation: gnl|database|id to database:id; sequence files unchanged",
            "annotation": {"total_genes": len(genes) + len(pseudogenes), "coding_genes": len(proteins),
                           "rna_genes": rna, "pseudogenes": len(pseudogenes), "gff_feature_rows": dict(feature_counts)},
            "ani": {"result": dict(result.attrib), "reference_comparisons": len(comparisons),
                    "highest_ani_reference": max(comparisons, key=lambda x: float(x["ANI"]))}}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--assembly", type=Path, required=True)
    cli.add_argument("--data-dir", type=Path, default=ROOT / "data")
    cli.add_argument("--quast-report", type=Path, help="Existing Galaxy QUAST report.tsv")
    cli.add_argument("--pgap-dir", type=Path, help="Existing directory containing annot.gbk, annot.fna, annot.faa, annot.gff and the ANI report")
    cli.add_argument("--report", type=Path, default=ROOT / "outputs/validation/results.json")
    args = cli.parse_args()
    result = verify(args.assembly, args.data_dir, args.quast_report, args.pgap_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
