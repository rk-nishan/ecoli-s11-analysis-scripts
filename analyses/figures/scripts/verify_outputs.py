"""Audit source integrity, ontology evidence, counts, and delivered image resolution."""

from __future__ import annotations

import argparse
from collections import Counter
import datetime as dt
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from PIL import Image

from common import ROOT, read_csv, sha256, verify_frozen_sources, write_json
from build_mappings import edges, normalize, ontology_data


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def verify(outdir: Path, dpi: int):
    verify_frozen_sources()
    mapping_dir = outdir / "data/mappings"
    if not mapping_dir.exists():
        mapping_dir = ROOT / "data/mappings"
    tables = outdir / "tables"
    fig2 = read_csv(tables / "figure2_gene_class_associations.csv")
    fig3 = read_csv(tables / "figure3_gene_mechanism_assignments.csv")
    summary = read_csv(tables / "figure3_mechanism_counts.csv")
    meta2 = json.loads((outdir / "provenance/Figure_2_S11_resfinder_sankey_render.json").read_text(encoding="utf-8"))
    meta3 = json.loads((outdir / "provenance/Figure_3_S11_card_mechanisms_render.json").read_text(encoding="utf-8"))

    expected2 = set()
    drug_mapping = {r["resfinder_resistance_term"]: r["antimicrobial_class"]
                    for r in read_csv(mapping_dir / "resfinder_drug_class_mapping.csv")}
    selected_ref = [r for r in read_csv(ROOT / "data/raw/resfinder_abricate.tsv", "\t")
                    if float(r["%COVERAGE"]) >= meta2["filters"]["minimum_coverage_percent"]
                    and float(r["%IDENTITY"]) >= meta2["filters"]["minimum_identity_percent"]]
    for row in selected_ref:
        for drug in row["RESISTANCE"].split(";"):
            expected2.add((row["PRODUCT"], drug_mapping[drug]))
    actual2 = {(r["gene"], r["antimicrobial_class"]) for r in fig2}
    require(actual2 == expected2 and len(fig2) == len(actual2), "Figure 2 associations differ from the raw data")
    require(Counter(r["antimicrobial_class"] for r in fig2) == meta2["class_counts"], "Figure 2 plotted class counts differ")
    require(all(r["association_weight"] == "1" for r in fig2), "Sankey links must have unit association weight")

    terms, index, _ = ontology_data()
    maps = read_csv(mapping_dir / "card_mechanism_mapping.csv")
    by_gene = {}
    for mapping in maps:
        aid = mapping["aro_accession"]
        leaf = [r for r in index if r["Accession"] == aid]
        require(len(leaf) == 1, "ARO accession is absent or ambiguous")
        require(normalize(mapping["abricate_gene"]) in
                {normalize(leaf[0]["Name"]), normalize(leaf[0]["CARD Short Name"])}, "Gene name does not match its ARO term")
        current = aid
        for step in mapping["ontology_path"].split(" | "):
            match = re.fullmatch(r"(ARO:\d+) --([a-z_]+)--> (ARO:\d+)", step)
            require(match is not None, f"Malformed ontology path: {step}")
            origin, relation, target = match.groups()
            require(origin == current and (relation, target) in set(edges(terms[origin])), "Ontology path contains an unsupported edge")
            current = target
        require(current == mapping["mechanism_aro_accession"], "Ontology path ends at the wrong mechanism")
        require("classified_as_resistance_mechanism" in terms[current].get("property_value", []), "Path target is not a tagged mechanism")
        require(terms[current]["name"][0] == mapping["mechanism"], "Mechanism label differs from the ontology")
        by_gene.setdefault(mapping["abricate_gene"], set()).add(mapping["mechanism"])

    raw_card = [r for r in read_csv(ROOT / "data/raw/card_abricate.tsv", "\t")
                if float(r["%COVERAGE"]) >= meta3["filters"]["minimum_coverage_percent"]
                and float(r["%IDENTITY"]) >= meta3["filters"]["minimum_identity_percent"]]
    expected3 = {(r["GENE"], r["SEQUENCE"], r["START"], r["END"], mechanism)
                 for r in raw_card for mechanism in by_gene[r["GENE"]]}
    actual3 = {(r["gene"], r["contig"], r["start"], r["end"], r["mechanism"]) for r in fig3}
    require(expected3 == actual3 and len(actual3) == len(fig3), "Figure 3 omits or duplicates gene-locus/mechanism pairs")
    counts = Counter(r["mechanism"] for r in fig3)
    require(counts == meta3["mechanism_counts"], "Figure 3 rendered counts differ from its data table")
    require({r["mechanism"] for r in summary} == set(counts), "Unexpected category omitted from summary")
    for row in summary:
        n = counts[row["mechanism"]]
        require(int(row["count_gene_locus_hits"]) == n, "Incorrect mechanism count")
        require(int(row["denominator_hits"]) == len(raw_card), "Incorrect denominator")
        require(row["percent_of_hits"] == f"{100*n/len(raw_card):.1f}", "Incorrect percentage")
    actual_loci = {(r["gene"], r["contig"], r["start"], r["end"]) for r in fig3}
    require(len(actual_loci) == len(raw_card), "At least one CARD locus lacks a classification")

    regenerated = []
    for name in ["card_mechanism_mapping.csv", "resfinder_drug_class_mapping.csv"]:
        candidate = outdir / "data/mappings" / name
        if candidate.exists():
            regenerated.append(name)

    exports = []
    for stem, expected_height in [("Figure_2_S11_resfinder_sankey", 100), ("Figure_3_S11_card_mechanisms", 90)]:
        base = outdir / "figures"
        png = base / f"{stem}_{dpi}dpi.png"
        tiff = base / f"{stem}_{dpi}dpi.tiff"
        with Image.open(png) as image:
            image.load()
            raster_size = image.size
            png_dpi = image.info.get("dpi", ())
            require(len(png_dpi) == 2 and min(png_dpi) >= 299.9, "PNG dpi metadata below 300")
            require(abs(png_dpi[0] - dpi) < 0.1, "PNG dpi differs from requested export setting")
            require(abs(raster_size[0] / dpi * 25.4 - 178) < 0.05, "PNG physical width differs")
            require(abs(raster_size[1] / dpi * 25.4 - expected_height) < 0.05, "PNG physical height differs")
            pixels = image.convert("RGB")
            require(pixels.getpixel((0, 0)) == (255, 255, 255), "Unexpected nonwhite figure background")
        with Image.open(tiff) as image:
            image.load()
            require(image.size == raster_size and image.mode == "RGB", "TIFF pixels/mode differ from RGB PNG")
            require(min(image.info.get("dpi", (0, 0))) >= 300, "TIFF dpi metadata below 300")
            require(image.tobytes() == pixels.tobytes(), "TIFF pixel values differ from PNG")
        svg = ET.parse(base / f"{stem}.svg").getroot()
        require(not svg.findall(".//{http://www.w3.org/2000/svg}image"), "SVG unexpectedly embeds raster data")
        require(bool(svg.findall(".//{http://www.w3.org/2000/svg}text")), "SVG lacks editable text")
        require((base / f"{stem}.pdf").read_bytes().startswith(b"%PDF-"), "Invalid PDF header")
        exports.append({"figure": stem, "pixels": list(raster_size), "png_dpi": list(png_dpi),
                        "tiff_matches_png_pixels": True, "vector_svg_has_editable_text": True,
                        "svg_embedded_raster_images": 0, "pdf_present": True})
    for metadata in [meta2, meta3]:
        require(metadata["geometry_qa"]["all_text_within_canvas"], "Text clipping check failed")
        require(metadata["geometry_qa"]["minimum_font_pt"] >= 8, "Figure text is too small")

    record = {"passed": True, "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "figure2_genes": len({r["gene"] for r in fig2}), "figure2_pairs": len(fig2),
              "figure2_class_counts": dict(Counter(r["antimicrobial_class"] for r in fig2)),
              "figure3_hits": len(raw_card), "figure3_gene_mechanism_pairs": len(fig3),
              "figure3_mechanism_counts": dict(counts), "all_ontology_paths_valid": True,
              "regenerated_mappings_checked": regenerated, "exports": exports}
    write_json(outdir / "qa/validation.json", record)
    print(f"PASS: {len(fig2)} Figure 2 associations; {len(raw_card)} CARD hits / {len(fig3)} associations; raster and vector exports checked")
    return record


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--outdir", type=Path, default=ROOT)
    cli.add_argument("--dpi", type=int, default=600)
    args = cli.parse_args()
    verify(args.outdir, args.dpi)
