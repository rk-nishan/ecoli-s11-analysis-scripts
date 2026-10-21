"""Shared input validation and publication exports for the S11 figures."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import platform
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
from PIL import Image

CODE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("S11_FIGURE_DATA_DIR", CODE_ROOT)).resolve()
FROZEN_SOURCE_FILES = [{'path': 'data/raw/resfinder_abricate.tsv',
  'sha256': '05f981400872d696676f56c1ec030b6821f0e09e5623ba8aa2fd839ff2165440'},
 {'path': 'data/raw/card_abricate.tsv',
  'sha256': 'af481abd2ec9d660af789616fa3a9cfaa4991695a405ced78359de0da1e47c10'},
 {'path': 'references/ontology-v4.0.2.tar.bz2',
  'sha256': 'b6c48b73d80dbe58eb9aefec9bf133b1c552a3993d363fd104aa161d45b294e2'},
 {'path': 'references/resfinder_antibiotic_classes_65c9ba3.txt',
  'sha256': 'f69f5ff72f8aaacdd3c3ef28ce9a21f7f1451c742c35a2e0d98914a9c4cf65fd'}]
INK = "#171717"
SECONDARY = "#4B5350"
TEAL = "#1F6F5F"
GRAY = "#8D9590"
LIGHT_GRAY = "#B2B8B4"
RULE = "#D8DDD9"
FIXED_DATE = dt.datetime(2026, 10, 22, tzinfo=dt.timezone.utc)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        if not rows:
            raise ValueError(f"No rows or field names supplied for {path}")
        fields = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def verify_frozen_sources() -> dict:
    manifest = {"files": FROZEN_SOURCE_FILES}
    for entry in manifest["files"]:
        actual = sha256(ROOT / entry["path"])
        if actual != entry["sha256"]:
            raise ValueError(f"Frozen input checksum differs: {entry['path']}")
    return manifest


def read_hits(path: Path, database: str, min_coverage: float, min_identity: float):
    rows = read_csv(path, "\t")
    fields = {"GENE", "DATABASE", "SEQUENCE", "START", "END", "%COVERAGE", "%IDENTITY", "PRODUCT", "RESISTANCE"}
    if not rows or not fields.issubset(rows[0]):
        raise ValueError(f"Missing ABRicate columns or empty input: {path}")
    passing, audit, seen = [], [], set()
    for number, row in enumerate(rows, start=2):
        if row["DATABASE"].lower() != database:
            raise ValueError(f"Wrong database in {path}, line {number}")
        coverage, identity = float(row["%COVERAGE"]), float(row["%IDENTITY"])
        if not (0 <= coverage <= 100 and 0 <= identity <= 100):
            raise ValueError(f"Invalid percentages at line {number}")
        if int(row["START"]) > int(row["END"]):
            raise ValueError(f"Invalid coordinates at line {number}")
        key = (row["GENE"], row["SEQUENCE"], row["START"], row["END"], row.get("STRAND", ""))
        if key in seen:
            raise ValueError(f"Duplicate gene-locus row at line {number}; resolve explicitly")
        seen.add(key)
        included = coverage >= min_coverage and identity >= min_identity
        audit.append({"source_line": number, **row, "included_in_figure": included})
        if included:
            passing.append(row)
    if not passing:
        raise ValueError("No ABRicate hits pass the selected thresholds")
    return passing, audit


def parser(description: str) -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=description)
    result.add_argument("--outdir", type=Path, default=ROOT,
                        help="Output root; only its figures/, tables/, and provenance/ files are generated")
    result.add_argument("--mapping-dir", type=Path, default=ROOT / "data/mappings")
    result.add_argument("--dpi", type=int, default=600)
    result.add_argument("--min-coverage", type=float, default=90.0)
    result.add_argument("--min-identity", type=float, default=80.0)
    return result


def configure(args) -> None:
    if args.dpi < 300:
        raise ValueError("Raster output must be at least 300 dpi")
    if not all(0 <= x <= 100 for x in (args.min_coverage, args.min_identity)):
        raise ValueError("Coverage and identity thresholds must be between 0 and 100")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": SECONDARY, "ytick.color": INK,
        "axes.linewidth": 0.65, "axes.edgecolor": RULE,
        "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "svg.hashsalt": "S11-AMR-figures-v1",
        "pdf.fonttype": 42, "pdf.compression": 9,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "mathtext.fontset": "dejavusans", "axes.unicode_minus": True,
    })


def assert_text_inside_canvas(fig) -> dict:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    frame = fig.bbox
    measurements, outside = [], []
    for artist in fig.findobj(Text):
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        box = artist.get_window_extent(renderer)
        record = {"text": artist.get_text(), "font_pt": artist.get_fontsize(),
                  "bounds_px": [float(v) for v in box.bounds]}
        measurements.append(record)
        if (box.x0 < frame.x0 - 1 or box.y0 < frame.y0 - 1 or
                box.x1 > frame.x1 + 1 or box.y1 > frame.y1 + 1):
            outside.append(record)
        if artist.get_fontsize() < 8:
            raise ValueError(f"Text below 8 pt: {artist.get_text()}")
    if outside:
        raise ValueError(f"Clipped figure text: {outside}")
    return {"all_text_within_canvas": True, "minimum_font_pt": min(x["font_pt"] for x in measurements),
            "text_measurements": measurements}


def export_figure(fig, stem: str, args, details: dict) -> None:
    """Save true 600 dpi raster images and vector PDF/SVG at fixed physical size."""
    qa = assert_text_inside_canvas(fig)
    folder = args.outdir / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    png = folder / f"{stem}_{args.dpi}dpi.png"
    fig.savefig(png, dpi=args.dpi, metadata={"Software": "S11 reproducible figures / Matplotlib",
                "Description": details["description"]}, bbox_inches=None)
    with Image.open(png) as source:
        rgb = source.convert("RGB")
        rgb.save(folder / f"{stem}_{args.dpi}dpi.tiff", compression="tiff_lzw", dpi=(args.dpi, args.dpi))
        dimensions = list(source.size)
        embedded_dpi = list(source.info["dpi"])
    fig.savefig(folder / f"{stem}.pdf", metadata={"Title": details["description"],
                "Author": "S11 study", "Creator": "Matplotlib",
                "CreationDate": FIXED_DATE, "ModDate": FIXED_DATE}, bbox_inches=None)
    fig.savefig(folder / f"{stem}.svg", metadata={"Title": details["description"],
                "Date": FIXED_DATE.isoformat(), "Creator": "S11 reproducible figures"}, bbox_inches=None)
    write_json(args.outdir / "provenance" / f"{stem}_render.json", {
        **details, "dpi": args.dpi, "raster_pixels": dimensions, "embedded_png_dpi": embedded_dpi,
        "physical_size_mm": [float(x) * 25.4 for x in fig.get_size_inches()],
        "filters": {"minimum_coverage_percent": args.min_coverage,
                    "minimum_identity_percent": args.min_identity},
        "runtime": {"python": platform.python_version(), **{x: importlib.metadata.version(x)
                    for x in ["matplotlib", "numpy", "Pillow"]}},
        "randomness": "None", "font": "DejaVu Sans from Matplotlib; no font files redistributed",
        "input_sha256": {p: sha256(args.mapping_dir / Path(p).name
                                  if p.startswith("data/mappings/") else ROOT / p)
                         for p in details["input_files"]},
        "geometry_qa": qa,
        "vector_date_metadata": "Fixed release date; actual execution time is recorded separately",
    })
    plt.close(fig)
