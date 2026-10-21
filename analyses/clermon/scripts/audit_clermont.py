"""Extract and verify the S11 primer evidence independently from BLAST XML.

This is an audit of the unmodified upstream result, not a replacement classifier.
It checks alignment sequences against the original input, reconstructs supported
primer pairs, and records a sensitivity run using the settings printed in the
upstream help text (which differ from its actual executable defaults).
"""

from __future__ import annotations

import argparse
from run_local_clermont import ASSEMBLY_SHA256, CLASSIFIER_SHA256, PRIMER_SHA256
import ast
from collections import Counter, defaultdict
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

SOURCE_ROOT = Path(__file__).resolve().parents[1]
ROOT = SOURCE_ROOT


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def fasta(path):
    result, name, letters = {}, None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                result[name] = "".join(letters).upper()
            name, letters = line[1:].split()[0], []
        elif line.strip():
            letters.append(line.strip())
    if name is not None:
        result[name] = "".join(letters).upper()
    return result


def revcomp(sequence):
    return sequence.translate(str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN"))[::-1]


def write_tsv(path, rows, columns=None):
    if columns is None:
        columns = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    global ROOT
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--outdir", type=Path, default=SOURCE_ROOT.parents[1] / "outputs/clermon")
    args = cli.parse_args()
    ROOT = args.outdir.resolve()
    for path, expected_hash in [(ROOT / "software/ClermonTyping/bin/clermont.py", CLASSIFIER_SHA256),
                                (ROOT / "software/ClermonTyping/data/primers.fasta", PRIMER_SHA256)]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected_hash:
            raise ValueError("Staged classifier or primer differs from the pinned upstream source")
    genome = fasta(ROOT / "input/S11_contigs.fasta")
    primers = fasta(ROOT / "software/ClermonTyping/data/primers.fasta")
    result = json.loads((ROOT / "results/S11_clermont_result.json").read_text(encoding="utf-8"))
    source = ast.parse((ROOT / "software/ClermonTyping/bin/clermont.py").read_text(encoding="utf-8"))
    expected = None
    for node in ast.walk(source):
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Attribute) and t.attr == "names" for t in node.targets):
            expected = ast.literal_eval(node.value)
            break
    if expected is None:
        raise ValueError("Unable to read the unchanged upstream marker-distance dictionary")
    xml = ET.parse(ROOT / "results/S11_primers_blast.xml").getroot()
    seen_queries, hsps, accepted = set(), [], defaultdict(lambda: {"F": [], "R": []})
    for iteration in xml.findall(".//Iteration"):
        query = iteration.findtext("Iteration_query-def").split()[0]
        seen_queries.add(query)
        marker, direction = query.rsplit("_", 1)
        query_length = int(iteration.findtext("Iteration_query-len"))
        if query_length != len(primers[query]):
            raise ValueError("BLAST query length differs from the frozen primer file")
        for hit in iteration.findall("Iteration_hits/Hit"):
            contig = hit.findtext("Hit_def").split()[0]
            length = int(hit.findtext("Hit_len"))
            if length != len(genome[contig]):
                raise ValueError("BLAST database contig length differs from the frozen FASTA")
            for hsp in hit.findall("Hit_hsps/Hsp"):
                number = len(hsps) + 1
                qstart, qend = int(hsp.findtext("Hsp_query-from")), int(hsp.findtext("Hsp_query-to"))
                start, end = int(hsp.findtext("Hsp_hit-from")), int(hsp.findtext("Hsp_hit-to"))
                frame = int(hsp.findtext("Hsp_hit-frame"))
                identities = int(hsp.findtext("Hsp_identity"))
                alignment_length = int(hsp.findtext("Hsp_align-len"))
                qseq = hsp.findtext("Hsp_qseq").replace("-", "").upper()
                sseq = hsp.findtext("Hsp_hseq").replace("-", "").upper()
                observed_query = primers[query][min(qstart, qend)-1:max(qstart, qend)]
                if qstart > qend:
                    observed_query = revcomp(observed_query)
                observed_subject = genome[contig][min(start, end)-1:max(start, end)]
                if start > end:
                    observed_subject = revcomp(observed_subject)
                if qseq != observed_query or sseq != observed_subject:
                    raise ValueError(f"Aligned sequence disagrees with its input: HSP {number}")
                missing_or_mismatched = query_length - identities
                end_matches = hsp.findtext("Hsp_midline")[-3:].count("|")
                passing = missing_or_mismatched <= 6 and qend >= query_length and end_matches >= 3
                row = {"hsp_number": number, "primer": query, "marker": marker, "primer_direction": direction,
                       "contig": contig, "contig_length_bp": length, "query_length_bp": query_length,
                       "query_start": qstart, "query_end": qend, "subject_start": start, "subject_end": end,
                       "subject_strand": frame, "identical_bases": identities, "alignment_length_bp": alignment_length,
                       "identity_percent_of_alignment": f"{100*identities/alignment_length:.4f}",
                       "query_length_minus_identities": missing_or_mismatched, "three_prime_seed_matches": end_matches,
                       "passes_upstream_primer_filters": passing, "alignment_matches_input_sequences": True}
                hsps.append(row)
                if passing:
                    accepted[marker][direction].append(row)
    if seen_queries != set(primers):
        raise ValueError("Not all 32 frozen primers were searched")

    pairs = []
    for marker, strands in accepted.items():
        for forward in strands["F"]:
            for reverse in strands["R"]:
                distance = abs(forward["subject_start"] - reverse["subject_start"])
                if (forward["contig"] != reverse["contig"] or
                    forward["subject_strand"] + reverse["subject_strand"] != 0 or
                    not expected[marker]*0.9 <= distance <= expected[marker]*1.1):
                    continue
                convergent = ((forward["subject_start"] < reverse["subject_start"] and forward["subject_strand"] == 1)
                              or (forward["subject_start"] > reverse["subject_start"] and forward["subject_strand"] == -1))
                pairs.append({"marker": marker, "contig": forward["contig"], "forward_hsp_number": forward["hsp_number"],
                              "reverse_hsp_number": reverse["hsp_number"], "forward_primer_start": forward["subject_start"],
                              "reverse_primer_start": reverse["subject_start"], "forward_strand": forward["subject_strand"],
                              "reverse_strand": reverse["subject_strand"], "primer_start_distance_bp": distance,
                              "expected_start_distance_bp": expected[marker], "within_10_percent": True,
                              "primers_convergent": convergent,
                              "genomic_interval_start": min(forward["subject_start"], reverse["subject_start"]),
                              "genomic_interval_end": max(forward["subject_start"], reverse["subject_start"]),
                              "inclusive_genomic_span_bp": distance+1})
    observed_markers = {p["marker"] for p in pairs}
    if observed_markers != set(result["detected_amplicons"]):
        raise ValueError("Independent XML/sequence audit disagrees with upstream positive amplicons")
    if not all(p["primers_convergent"] for p in pairs):
        raise ValueError("At least one upstream-supported amplicon has nonconvergent primer orientation")
    write_tsv(ROOT / "results/primer_alignment_audit.tsv", hsps)
    write_tsv(ROOT / "results/valid_amplicon_pairs.tsv", pairs)
    marker_rows = []
    for marker, expected_size in expected.items():
        current = [p for p in pairs if p["marker"] == marker]
        marker_rows.append({"marker": marker, "valid_amplicon_detected": "+" if current else "-",
                            "accepted_forward_primer_hits": len(accepted[marker]["F"]),
                            "accepted_reverse_primer_hits": len(accepted[marker]["R"]),
                            "valid_primer_pairs": len(current), "expected_primer_start_distance_bp": expected_size,
                            "observed_primer_start_distances_bp": ";".join(str(p["primer_start_distance_bp"]) for p in current),
                            "positive_contigs": ";".join(sorted({p["contig"] for p in current}))})
    write_tsv(ROOT / "results/marker_summary.tsv", marker_rows)
    with (ROOT / "results/genomic_amplicon_intervals.fasta").open("w", encoding="ascii", newline="\n") as handle:
        for number, pair in enumerate(pairs, 1):
            seq = genome[pair["contig"]][pair["genomic_interval_start"]-1:pair["genomic_interval_end"]]
            if pair["forward_strand"] == -1:
                seq = revcomp(seq)
            handle.write(f">S11_{pair['marker']}_{number} contig={pair['contig']} interval={pair['genomic_interval_start']}-{pair['genomic_interval_end']} forward_primer_oriented=true\n")
            for position in range(0, len(seq), 80):
                handle.write(seq[position:position+80] + "\n")

    args = [sys.executable, "-B", "software/ClermonTyping/bin/clermont.py", "-x", "results/S11_primers_blast.xml",
            "-m", "2", "-l", "5", "-s", "0", "--log_level", "DEBUG", "--log_file", "logs/help_documented_settings.debug.log"]
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    alternative = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=120)
    (ROOT / "results/S11_help_documented_settings_raw.tsv").write_text(alternative.stdout, encoding="utf-8")
    (ROOT / "logs/help_documented_settings.stderr.txt").write_text(alternative.stderr, encoding="utf-8")
    save_json(ROOT / "provenance/validation_commands.json", [{"purpose": "Sensitivity to the help-documented mismatch=2, seed=5 settings",
              "argv": args, "cwd": str(ROOT), "started_utc": started,
              "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "exit_code": alternative.returncode}])
    if alternative.returncode != 0:
        raise RuntimeError(alternative.stderr)
    alt_fields = alternative.stdout.strip().split("\t")
    if len(alt_fields) != 4:
        raise ValueError("Malformed sensitivity-check output")
    alternative_group = alt_fields[3]
    settings_check = {"primary_effective_settings": {"mismatch": 6, "seed_length": 3, "phylogroup": result["phylogroup"]},
                      "help_documented_settings": {"mismatch": 2, "seed_length": 5, "phylogroup": alternative_group,
                                                   "detected_amplicons": ast.literal_eval(alt_fields[0])},
                      "phylogroup_agrees": alternative_group == result["phylogroup"],
                      "amplicon_sets_agree": set(ast.literal_eval(alt_fields[0])) == set(result["detected_amplicons"])}
    save_json(ROOT / "provenance/settings_sensitivity_check.json", settings_check)

    xml_parameters = {x.tag: x.text for x in xml.findall(".//Parameters/*")}
    statistics = {x.tag: x.text for x in xml.findall(".//Statistics")[0]}
    if int(statistics["Statistics_db-num"]) != len(genome) or int(statistics["Statistics_db-len"]) != sum(map(len, genome.values())):
        raise ValueError("BLAST database does not contain the complete supplied assembly")
    save_json(ROOT / "provenance/blast_xml_parameters.json", {"version": xml.findtext("BlastOutput_version"),
              "parameters": xml_parameters, "database_statistics": statistics, "query_count": len(seen_queries)})
    input_ok = hashlib.sha256((ROOT / "input/S11_contigs.fasta").read_bytes()).hexdigest() == ASSEMBLY_SHA256
    if not input_ok:
        raise ValueError("Assembly checksum differs from the frozen S11 input")
    validation = {"passed": True, "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                  "independent_parser": "Python xml.etree.ElementTree; upstream uses Bio.Blast.NCBIXML",
                  "all_hsp_sequences_match_frozen_input": True, "primer_queries_checked": len(seen_queries),
                  "hsp_alignments_checked": len(hsps), "accepted_primer_hsps": sum(len(v[s]) for v in accepted.values() for s in ["F", "R"]),
                  "supported_amplicon_pairs": len(pairs), "markers_match_upstream_result": True,
                  "all_supported_pairs_are_convergent": True, "complete_assembly_searched": True,
                  "input_contigs": len(genome), "input_bases": sum(map(len, genome.values())),
                  "phylogroup": result["phylogroup"], "documented_settings_check": settings_check,
                  "frozen_assembly_sha256_matches": input_ok}
    save_json(ROOT / "provenance/validation.json", validation)
    print(json.dumps(validation, ensure_ascii=True))


if __name__ == "__main__":
    main()
