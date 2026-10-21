"""Fetch pinned public classification references into an external input workspace."""

import argparse
import hashlib
from pathlib import Path
import shutil
import urllib.request

REFERENCES = [{'path': 'references/ontology-v4.0.2.tar.bz2',
  'url': 'https://card.mcmaster.ca/download/5/ontology-v4.0.2.tar.bz2',
  'sha256': 'b6c48b73d80dbe58eb9aefec9bf133b1c552a3993d363fd104aa161d45b294e2'},
 {'path': 'references/resfinder_antibiotic_classes_65c9ba3.txt',
  'url': 'https://bitbucket.org/genomicepidemiology/resfinder_db/raw/65c9ba3/antibiotic_classes.txt',
  'sha256': 'f69f5ff72f8aaacdd3c3ef28ce9a21f7f1451c742c35a2e0d98914a9c4cf65fd'},
 {'path': 'clermon/primers.fasta',
  'url': 'https://raw.githubusercontent.com/A-BN/ClermonTyping/69da4db3cd98ab40528bacc289529323ff1407a0/data/primers.fasta',
  'sha256': 'c28c71a205f107baa25b513ca4872d057e55ccdb9a5ba582a48d6a2b4b99d73e'}]


def fetch(reference, root):
    target = root / reference["path"]
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest() != reference["sha256"]:
            raise FileExistsError(f"Existing reference differs: {target}")
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".download-part")
    request = urllib.request.Request(reference["url"], headers={"User-Agent": "S11-analysis-scripts/1.0"})
    with urllib.request.urlopen(request, timeout=120) as source, partial.open("xb") as output:
        shutil.copyfileobj(source, output)
    if hashlib.sha256(partial.read_bytes()).hexdigest() != reference["sha256"]:
        raise ValueError(f"Reference failed SHA-256 verification: {target.name}")
    if target.exists():
        raise FileExistsError(target)
    partial.rename(target)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parents[1] / "inputs")
    args = cli.parse_args()
    for reference in REFERENCES:
        fetch(reference, args.outdir)
        print(f"Verified reference: {reference['path']}")


if __name__ == "__main__":
    main()
