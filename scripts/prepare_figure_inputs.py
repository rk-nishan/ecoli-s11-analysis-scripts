"""Copy externally supplied S11 ABRicate reports into a new input workspace."""

import argparse
import hashlib
from pathlib import Path
import shutil

EXPECTED = {
    "card": "af481abd2ec9d660af789616fa3a9cfaa4991695a405ced78359de0da1e47c10",
    "resfinder": "05f981400872d696676f56c1ec030b6821f0e09e5623ba8aa2fd839ff2165440",
}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--card", type=Path, required=True)
    cli.add_argument("--resfinder", type=Path, required=True)
    cli.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parents[1] / "inputs")
    args = cli.parse_args()
    for name in ["card", "resfinder"]:
        source = getattr(args, name)
        if hashlib.sha256(source.read_bytes()).hexdigest() != EXPECTED[name]:
            raise ValueError(f"{name} report differs from the frozen S11 figure input")
    for name in ["card", "resfinder"]:
        target = args.outdir / "data/raw" / f"{name}_abricate.tsv"
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest() != EXPECTED[name]:
                raise FileExistsError(target)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with getattr(args, name).open("rb") as source, target.open("xb") as output:
            shutil.copyfileobj(source, output)
    print(f"Verified external figure reports are ready in {args.outdir}")


if __name__ == "__main__":
    main()
