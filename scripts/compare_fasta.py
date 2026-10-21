"""Compare FASTA bytes, ordered sequences, and sequence multisets exactly."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def read_fasta(path):
    records, header, fragments = [], None, []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(fragments).upper()))
            header, fragments = line[1:], []
        else:
            if header is None:
                raise ValueError("Sequence data precedes the first FASTA header")
            fragments.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(fragments).upper()))
    if not records or any(not h or not s for h, s in records):
        raise ValueError("Empty FASTA record or header")
    if set("".join(s for _, s in records)) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("Unexpected nucleotide characters")
    return records


def compare(first, second):
    a, b = read_fasta(first), read_fasta(second)
    sequences_a, sequences_b = [s for _, s in a], [s for _, s in b]
    counts_a, counts_b = Counter(sequences_a), Counter(sequences_b)
    return {
        "files": [{"name": p.name, "bytes": p.stat().st_size,
                   "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                   "contigs": len(r), "bases": sum(len(s) for _, s in r)}
                  for p, r in [(first, a), (second, b)]],
        "byte_for_byte_identical": first.read_bytes() == second.read_bytes(),
        "headers_identical_in_order": [h for h, _ in a] == [h for h, _ in b],
        "sequences_identical_in_order": sequences_a == sequences_b,
        "sequences_identical_ignoring_order": counts_a == counts_b,
        "matching_contigs_including_duplicates": sum((counts_a & counts_b).values()),
        "unmatched_first": sum((counts_a - counts_b).values()),
        "unmatched_second": sum((counts_b - counts_a).values()),
        "comparison": "Exact sequence equality after removing whitespace and uppercasing; reverse complements remain distinct",
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("first", type=Path)
    cli.add_argument("second", type=Path)
    cli.add_argument("--report", type=Path)
    args = cli.parse_args()
    report = compare(args.first, args.second)
    text = json.dumps(report, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("x", encoding="utf-8") as handle:
            handle.write(text)
    print(text, end="")
    return 0 if report["sequences_identical_ignoring_order"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
