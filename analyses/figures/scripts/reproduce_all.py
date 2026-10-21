"""Rebuild both figures offline into a fresh directory, preserving supplied outputs."""

import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--outdir", type=Path, default=ROOT.parents[1] / "outputs/figures")
    cli.add_argument("--input-dir", type=Path, required=True, help="External input folder prepared by prepare_figure_inputs.py and fetch_references.py")
    cli.add_argument("--dpi", type=int, default=600)
    args = cli.parse_args()
    destination = args.outdir.resolve()
    env = dict(os.environ, S11_FIGURE_DATA_DIR=str(args.input_dir.resolve()))
    if destination.exists():
        raise SystemExit(f"Output directory already exists: {destination}. Select a fresh --outdir.")
    destination.mkdir(parents=True)
    # build_mappings is checked against the bundled mapping in verify_outputs.
    subprocess.run([sys.executable, "-B", str(ROOT / "scripts/build_mappings.py"),
                    "--outdir", str(destination)], check=True, env=env)
    for script in ["figure2_sankey.py", "figure3_mechanisms.py"]:
        subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / script),
                        "--outdir", str(destination), "--dpi", str(args.dpi),
                        "--mapping-dir", str(destination / "data/mappings")], check=True, env=env)
    subprocess.run([sys.executable, "-B", str(ROOT / "scripts/verify_outputs.py"),
                    "--outdir", str(destination), "--dpi", str(args.dpi)], check=True, env=env)
    print(f"Both figures reproduced and verified: {destination}")


if __name__ == "__main__":
    main()
