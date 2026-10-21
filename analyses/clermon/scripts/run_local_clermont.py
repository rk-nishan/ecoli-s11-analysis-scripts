"""Run the author's standalone ClermonTyping PCR workflow using local BLAST+.

The original upstream classifier and primer file are executed unmodified.
This follows the upstream README's 'Clermont Typing without mash and R' mode,
with the word size and effective classifier defaults from its current launcher.
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ASSEMBLY_SHA256 = '406a7003fa816dd990fae93d2a9148659b35d2b56c21667c2322e3f65478372f'
CLASSIFIER_SHA256 = 'c51a48282e784ada174d684286cc04476becb4a3aa73504f4f7f95030aaee829'
PRIMER_SHA256 = 'c28c71a205f107baa25b513ca4872d057e55ccdb9a5ba582a48d6a2b4b99d73e'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT.parents[1] / "outputs/clermon")
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--primers", type=Path, default=ROOT.parents[1] / "inputs/clermon/primers.fasta")
    parser.add_argument("--blast-bin", type=Path, help="Directory containing installed blastn and makeblastdb")
    args = parser.parse_args()
    target = args.outdir.resolve()
    if args.blast_bin:
        os.environ["PATH"] = str(args.blast_bin.resolve()) + os.pathsep + os.environ.get("PATH", "")
    if target.exists():
        raise SystemExit("An output folder already exists here. Use a fresh --outdir to preserve it.")
    # NCBI's native Windows database library mishandles a working path containing spaces.
    # Run in a fresh space-free local workspace and return all records to the chosen folder.
    use_temporary = os.name == "nt" and " " in str(target)
    output = Path(tempfile.mkdtemp(prefix="S11_ClermonTyping_")) if use_temporary else target
    if " " in str(output) and os.name == "nt":
        raise RuntimeError("A Windows BLAST working directory without spaces is required")
    for name in ["input", "results", "reports", "logs", "provenance", "db"]:
        (output / name).mkdir(parents=True, exist_ok=True)

    source_manifest = {"repository": "https://github.com/A-BN/ClermonTyping", "commit": "69da4db3cd98ab40528bacc289529323ff1407a0"}
    source_core = ROOT / "software/ClermonTyping/bin/clermont.py"
    frozen_primers = args.primers.resolve()
    if digest(source_core) != CLASSIFIER_SHA256 or digest(frozen_primers) != PRIMER_SHA256:
        raise ValueError("Pinned classifier or primer checksum changed")
    for source, relative in [(source_core, "bin/clermont.py"), (frozen_primers, "data/primers.fasta")]:
        staged = output / "software/ClermonTyping" / relative
        staged.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, staged)
        if digest(source) != digest(staged):
            raise RuntimeError("Staged upstream source checksum differs")
    frozen_input = args.assembly.resolve()
    expected = ASSEMBLY_SHA256
    if digest(frozen_input) != expected:
        raise ValueError("Frozen S11 input checksum changed")
    destination = output / "input/S11_contigs.fasta"
    if destination != frozen_input:
        if destination.exists():
            raise FileExistsError(destination)
        shutil.copyfile(frozen_input, destination)
    save_json(output / "provenance/execution_staging.json", {
        "requested_result_folder": str(target), "actual_working_folder": str(output),
        "temporary_local_folder": use_temporary,
        "reason": "Native Windows BLAST database paths must not contain spaces" if use_temporary else "No path workaround needed",
        "original_inputs_modified": False, "staged_files_verified_by_sha256": True,
    })

    blastn = shutil.which("blastn")
    makeblastdb = shutil.which("makeblastdb")
    if not blastn or not makeblastdb:
        raise RuntimeError("Existing local BLAST+ blastn and makeblastdb must be on PATH")
    import Bio
    commands = []

    def execute(name, argv, stdout_file=None):
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        result = subprocess.run(argv, cwd=output, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=600)
        stdout_path = output / (stdout_file or f"logs/{name}.stdout.txt")
        stderr_path = output / f"logs/{name}.stderr.txt"
        stdout_path.write_text(result.stdout, encoding="utf-8")
        stderr_path.write_text(result.stderr, encoding="utf-8")
        record = {"stage": name, "argv": [str(a) for a in argv], "cwd": str(output),
                  "started_utc": started, "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                  "exit_code": result.returncode, "stdout": stdout_path.relative_to(output).as_posix(),
                  "stderr": stderr_path.relative_to(output).as_posix()}
        commands.append(record)
        save_json(output / "provenance/commands.json", commands)
        if result.returncode != 0:
            raise RuntimeError(f"{name} exited {result.returncode}: {result.stderr}")
        return result.stdout

    blast_version = execute("blastn_version", [blastn, "-version"])
    make_version = execute("makeblastdb_version", [makeblastdb, "-version"])
    execute("blastn_help", [blastn, "-help"])
    execute("makeblastdb_help", [makeblastdb, "-help"])
    execute("classifier_help", [sys.executable, "-B", "software/ClermonTyping/bin/clermont.py", "--help"])

    parameters = {
        "sample": "S11", "mode": "Upstream-supported standalone PCR typing without Mash/R",
        "software_repository": source_manifest["repository"], "software_commit": source_manifest["commit"],
        "launcher_release": "24.02 (February 2024)", "python_classifier_internal_version": "1.4.0 (February 2019)",
        "assembly_sha256": expected, "primer_file": "software/ClermonTyping/data/primers.fasta",
        "primer_sha256": digest(frozen_primers),
        "blast": {"task": "blastn", "minimum_identity_percent": 90, "word_size": 6,
                  "outfmt": 5, "num_threads": 1, "other_options": "Installed BLAST+ defaults; actual version, complete help and XML parameters preserved"},
        "classifier": {"max_query_length_minus_identities": 6, "three_prime_seed_length": 3,
                       "mismatches_allowed_in_seed": 0, "minimum_contig_size_bp": 0,
                       "expected_primer_start_distance_tolerance_percent": 10},
        "defaults_note": "The actual pinned executable code uses mismatch=6 and seed length=3; the upstream help/README states 2 and 5. Explicit 6/3 values reproduce the actual launcher behavior. The help-documented settings are checked separately by audit_clermont.py.",
        "primer_pair_rules": "Same contig, opposite strands, primer-start distance within +/-10% of the configured marker distance; raw pairing evidence retained separately",
        "mash_performed": False, "r_reporting_performed": False, "automatic_package_installation": False,
    }
    save_json(output / "provenance/parameters.json", parameters)
    save_json(output / "provenance/runtime.json", {
        "python": sys.version, "python_executable": sys.executable, "platform": platform.platform(),
        "biopython_version": Bio.__version__, "blastn_version": blast_version.strip(),
        "makeblastdb_version": make_version.strip(),
        "executables": {"blastn": {"path": blastn, "sha256": digest(blastn)},
                        "makeblastdb": {"path": makeblastdb, "sha256": digest(makeblastdb)}},
    })

    execute("makeblastdb", [makeblastdb, "-in", "input/S11_contigs.fasta", "-input_type", "fasta",
                             "-out", "db/S11", "-dbtype", "nucl"])
    execute("blastn", [blastn, "-query", "software/ClermonTyping/data/primers.fasta",
                        "-perc_identity", "90", "-task", "blastn", "-word_size", "6", "-outfmt", "5",
                        "-num_threads", "1", "-db", "db/S11", "-out", "results/S11_primers_blast.xml"])
    raw = execute("clermont_classifier", [sys.executable, "-B", "software/ClermonTyping/bin/clermont.py",
                   "-x", "results/S11_primers_blast.xml", "-m", "6", "-l", "3", "-s", "0",
                   "--log_level", "DEBUG", "--log_file", "logs/clermont_classifier.debug.log"],
                   stdout_file="results/S11_clermont_raw.tsv")
    fields = raw.strip().split("\t")
    if len(fields) != 4:
        raise ValueError(f"Unexpected upstream classifier output: {raw!r}")
    markers, quadruplex, supplementary = [ast.literal_eval(value) for value in fields[:3]]
    group = fields[3]
    order = ["arpA", "chuA", "yjaA", "TspE4.C2"]
    if len(quadruplex) != 4:
        raise ValueError("Malformed quadruplex output")
    result = {"sample": "S11", "phylogroup": group, "detected_amplicons": markers,
              "quadruplex_marker_order": order, "quadruplex_profile": quadruplex,
              "quadruplex_by_marker": dict(zip(order, quadruplex)),
              "supplementary_amplicons": supplementary, "ybgD_group_G_marker": "+" if "ybgD" in markers else "-",
              "method": "Unmodified ClermonTyping standalone PCR classifier", "assembly_sha256": expected,
              "scope": "In silico PCR phylogroup assignment; Mash comparison and R-generated HTML were not run",
              "all_commands_exit_zero": True}
    save_json(output / "results/S11_clermont_result.json", result)
    with (output / "results/S11_phylogroup.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "phylogroup", *order, "ybgD", "detected_amplicons", "supplementary_amplicons"])
        writer.writerow(["S11", group, *quadruplex, result["ybgD_group_G_marker"], ";".join(markers), ";".join(supplementary)])
    if output != target:
        # Source/input files are already present in the requested folder. All generated
        # database, result, diagnostic and execution records are returned verbatim.
        for file in output.rglob("*"):
            if not file.is_file():
                continue
            relative = file.relative_to(output)
            returned = target / relative
            returned.parent.mkdir(parents=True, exist_ok=True)
            if relative.parts[0] in {"input", "software"} and returned.exists():
                if digest(file) != digest(returned):
                    raise RuntimeError(f"Existing source/input differs: {relative}")
                continue
            shutil.copyfile(file, returned)
            if digest(file) != digest(returned):
                raise RuntimeError(f"Returned output checksum differs: {relative}")
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()
