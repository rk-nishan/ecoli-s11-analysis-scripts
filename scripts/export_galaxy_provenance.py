"""Preserve Galaxy job provenance without retrieving deleted dataset contents."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import threading
import time
from urllib.parse import quote
import xml.etree.ElementTree as ET
import zipfile

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[1] / "outputs/galaxy_history/tool_provenance"
HISTORY_ID = "11ac94870d0bb33a834c81ccf4f16ef2"
BASE = "https://usegalaxy.eu"
KEY = ""
LOCAL = threading.local()
TZ = timezone(timedelta(hours=6))


def now():
    return datetime.now(TZ).isoformat()


def redact(value):
    if isinstance(value, str):
        return value.replace(KEY, "[API key redacted]")
    if isinstance(value, dict):
        return {key: redact(item) for key, item in value.items()}
    if isinstance(value, list):
        return [redact(item) for item in value]
    return value


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(redact(value), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def session():
    if not hasattr(LOCAL, "session"):
        LOCAL.session = requests.Session()
        LOCAL.session.headers["x-api-key"] = KEY
        retry = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
        LOCAL.session.mount("https://", HTTPAdapter(max_retries=retry))
    return LOCAL.session


def fetch(endpoint, params=None, optional=False, raw=False):
    for attempt in range(4):
        try:
            response = session().get(BASE + endpoint, params=params, timeout=(20, 90), allow_redirects=False)
            if response.is_redirect:
                raise ValueError("Galaxy API redirect refused to keep credentials on the selected server")
            if optional and response.status_code in [400, 401, 403, 404, 405, 501]:
                return {"available": False, "http_status": response.status_code, "endpoint": endpoint}
            response.raise_for_status()
            if raw:
                return {"available": True, "text": redact(response.text), "language": response.headers.get("language"), "endpoint": endpoint}
            return response.json()
        except (requests.RequestException, ValueError):
            if attempt == 3:
                if optional:
                    return {"available": False, "reason": "Request failed after retries", "endpoint": endpoint}
                raise
            time.sleep(attempt + 1)


def safe(value, limit=65):
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", value).strip("_.")[:limit] or "unknown"


def decode_params(params):
    decoded = {}
    for key, value in params.items():
        if isinstance(value, str):
            try:
                value = json.loads(value)
            except (ValueError, TypeError):
                pass
        decoded[key] = value
    return decoded


def flatten(value, prefix=""):
    if isinstance(value, dict):
        if not value:
            yield prefix, "{}"
        for key, item in value.items():
            yield from flatten(item, prefix + ("|" if prefix else "") + str(key))
    elif isinstance(value, list):
        if not value:
            yield prefix, "[]"
        for index, item in enumerate(value):
            yield from flatten(item, prefix + "[%d]" % index)
    else:
        yield prefix, json.dumps(value, ensure_ascii=False)


def tsv(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def export_job(job_id, history_items):
    job = fetch("/api/jobs/" + job_id, {"full": "true"})
    if job.get("history_id") != HISTORY_ID:
        raise ValueError("Job does not belong to the requested history")
    tool_id = job.get("tool_id") or "unknown"
    parts = tool_id.split("/")
    name = parts[-2] if len(parts) > 1 else tool_id
    version = parts[-1] if len(parts) > 1 else "Galaxy_builtin"
    folder = ROOT / "jobs" / (safe(name, 30) + "__" + job_id[-12:])
    folder.mkdir(parents=True, exist_ok=True)
    save_json(folder / "job.json", job)
    save_json(folder / "parameters_raw.json", job.get("params") or {})
    params = decode_params(job.get("params") or {})
    save_json(folder / "parameters.json", params)
    parameter_rows = [{"parameter": path, "value": redact(value)} for path, value in flatten(params)]
    tsv(folder / "parameters.tsv", parameter_rows, ["parameter", "value"])
    (folder / "parameters.txt").write_text("\n".join(row["parameter"] + " = " + row["value"] for row in parameter_rows) + "\n", encoding="utf-8")
    command = job.get("command_line")
    (folder / "command.txt").write_text(redact(command or "[Command not provided by Galaxy for this job]\n"), encoding="utf-8")
    availability = {"command_line": "provided" if command else "not provided by Galaxy"}
    for key in ["stdout", "stderr", "tool_stdout", "tool_stderr", "job_stdout", "job_stderr"]:
        value = job.get(key)
        availability[key] = "not provided" if value is None else ("recorded, empty" if value == "" else "recorded, %d characters" % len(value))
        (folder / (key + ".txt")).write_text(redact(value or ""), encoding="utf-8")
    save_json(folder / "log_availability.json", availability)
    for key in ["inputs", "outputs", "output_collections", "dependencies", "job_messages"]:
        save_json(folder / (key + ".json"), job.get(key))
    display = fetch("/api/jobs/" + job_id + "/parameters_display", optional=True)
    save_json(folder / "parameters_display.json", display)
    metrics = fetch("/api/jobs/" + job_id + "/metrics", optional=True)
    save_json(folder / "runtime_metrics.json", metrics)
    timing = {key: job.get(key) for key in ["id", "history_id", "tool_id", "state", "exit_code", "create_time", "update_time", "command_version", "galaxy_version", "job_runner_name", "handler", "external_id"]}
    timing["exported_at"] = now()
    save_json(folder / "status_versions_and_timing.json", timing)
    mappings = []
    for role in ["inputs", "outputs", "output_collections"]:
        for parameter, value in (job.get(role) or {}).items():
            values = value if isinstance(value, list) else [value]
            for item in values:
                if not isinstance(item, dict):
                    continue
                source = history_items.get(item.get("id"), {})
                mappings.append({"role": role, "parameter": parameter, "id": item.get("id"), "src": item.get("src"),
                                 "history_number": source.get("hid"), "name": source.get("name"), "deleted": source.get("deleted"),
                                 "state": source.get("state"), "local_file": source.get("local_file", "")})
    tsv(folder / "input_output_mapping.tsv", mappings, ["role", "parameter", "history_number", "name", "id", "src", "deleted", "state", "local_file"])
    output_hids = sorted(source["hid"] for source in history_items.values() if source.get("creating_job") == job_id)
    row = {"tool": name, "wrapper_version": version, "runtime_version_reported": job.get("command_version") or "not separately reported",
           "job_id": job_id, "state": job.get("state"), "exit_code": job.get("exit_code"),
           "output_history_numbers": ",".join(map(str, output_hids)), "parameter_count": len(parameter_rows),
           "command_available": bool(command), "tool_stdout_characters": len(job.get("tool_stdout") or ""),
           "tool_stderr_characters": len(job.get("tool_stderr") or ""),
           "folder": folder.relative_to(ROOT).as_posix(), "tool_id": tool_id}
    print(json.dumps({"job_exported": name, "state": job.get("state"), "job": job_id, "command_saved": bool(command)}), flush=True)
    return row


def export_tool(tool_id):
    parts = tool_id.split("/")
    name = parts[-2] if len(parts) > 1 else tool_id
    version = parts[-1] if len(parts) > 1 else "Galaxy_builtin"
    folder = ROOT / "tool_definitions" / (safe(name, 35) + "__" + safe(version, 35))
    folder.mkdir(parents=True, exist_ok=True)
    encoded = quote(tool_id, safe="")
    definition = fetch("/api/tools/" + encoded, {"io_details": "true"}, optional=True)
    save_json(folder / "tool_definition_and_defaults.json", definition)
    raw = fetch("/api/tools/" + encoded + "/raw_tool_source", optional=True, raw=True)
    availability = {"tool_id": tool_id, "tool_definition_available": not (isinstance(definition, dict) and definition.get("available") is False),
                    "source_available": raw.get("available")}
    if raw.get("available"):
        language = raw.get("language") or "xml"
        (folder / ("tool_source.xml" if language == "xml" else "tool_source.txt")).write_text(raw["text"], encoding="utf-8")
        if language == "xml":
            try:
                xml = ET.fromstring(raw["text"])
                requirements = [{"name": r.text, "attributes": r.attrib} for r in xml.findall(".//requirement")]
                containers = [{"identifier": r.text, "attributes": r.attrib} for r in xml.findall(".//container")]
                save_json(folder / "declared_software_requirements.json", {"evidence": "Requirements declared by the Galaxy wrapper, not a complete runtime environment snapshot.", "requirements": requirements, "containers": containers})
            except ET.ParseError:
                availability["source_parse_note"] = "Source saved exactly as returned; XML parsing failed"
    else:
        save_json(folder / "source_unavailable.json", raw)
    citations = fetch("/api/tools/" + encoded + "/citations", optional=True)
    save_json(folder / "citations.json", citations)
    save_json(folder / "availability.json", availability)
    return availability


def main():
    global ROOT, KEY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    ROOT = args.outdir.resolve()
    KEY = os.environ.pop("GALAXY_API_KEY", "")
    if not KEY:
        raise SystemExit("Set GALAXY_API_KEY in the environment")
    if ROOT.exists():
        raise SystemExit("Select a fresh --outdir to preserve existing records")
    ROOT.mkdir(parents=True)
    history = fetch("/api/histories/" + HISTORY_ID)
    if history["name"] != "Manual Analysis S-11":
        raise ValueError("Unexpected history name")
    contents = fetch("/api/histories/%s/contents" % HISTORY_ID, {"details": "all"})
    # Metadata are preserved for complete lineage; deleted dataset sequences/contents are never fetched.
    metadata = [{key: value for key, value in item.items() if key != "peek"} for item in contents]
    save_json(ROOT / "history" / "history.json", history)
    save_json(ROOT / "history" / "history_contents_metadata.json", metadata)
    save_json(ROOT / "history" / "galaxy_server_version.json", fetch("/api/version", optional=True))
    items = {item["id"]: item for item in metadata}
    inventory = ROOT.parent / "inventory.tsv"
    if inventory.exists():
        with inventory.open(encoding="utf-8-sig", newline="") as source:
            for row in csv.DictReader(source, delimiter="\t"):
                item = items.get(row["dataset_id"])
                if item:
                    item["local_file"] = str(ROOT.parent / row["relative_path"])
    job_ids = sorted({item["creating_job"] for item in contents if item.get("history_content_type") == "dataset" and item.get("creating_job")})
    print(json.dumps({"event": "provenance_export_started", "jobs": len(job_ids)}), flush=True)
    rows = []
    errors = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(export_job, job_id, items): job_id for job_id in job_ids}
        for future in as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as error:
                errors.append({"job_id": futures[future], "type": type(error).__name__, "message": redact(str(error))})
    rows.sort(key=lambda row: (row["tool"].lower(), row["job_id"]))
    tsv(ROOT / "tool_runs.tsv", rows, ["tool", "wrapper_version", "runtime_version_reported", "job_id", "state", "exit_code", "output_history_numbers", "parameter_count", "command_available", "tool_stdout_characters", "tool_stderr_characters", "folder", "tool_id"])
    tool_ids = sorted({row["tool_id"] for row in rows})
    tool_results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(export_tool, tool_ids):
            tool_results.append(result)
    active_tools = sorted({row["tool"] for row in rows if any(item.get("creating_job") == row["job_id"] and not item.get("deleted") for item in contents)})
    summary = {"exported_at": now(), "history_name": history["name"], "history_id": HISTORY_ID,
               "scope": "Parameters, commands and available logs for every dataset-producing job in the history, including previous failed attempts. No deleted dataset contents are downloaded.",
               "jobs_expected": len(job_ids), "jobs_saved": len(rows), "tool_versions_saved": len(tool_ids),
               "tools": sorted({row["tool"] for row in rows}), "tools_with_current_outputs": active_tools,
               "successful_jobs": sum(row["state"] == "ok" for row in rows),
               "previous_failed_jobs": sum(row["state"] == "error" for row in rows),
               "jobs_with_commands": sum(row["command_available"] for row in rows),
               "tool_source_availability": tool_results, "errors": errors,
               "export_complete": len(rows) == len(job_ids) and not errors}
    save_json(ROOT / "provenance_summary.json", summary)
    text = ["# Manual Analysis S-11: tool parameters and execution records", "",
            "Saved on " + summary["exported_at"] + ".", "",
            "This package preserves Galaxy's actual job records, including previous failed attempts. Deleted dataset contents are excluded.", "",
            "- Job runs saved: **%d / %d**." % (len(rows), len(job_ids)),
            "- Tool versions represented: **%d**." % len(tool_ids),
            "- Successful jobs: **%d**; previous failed jobs: **%d**." % (summary["successful_jobs"], summary["previous_failed_jobs"]),
            "- Tools: " + ", ".join(summary["tools"]) + ".", "",
            "Open `tool_runs.tsv` for the inventory. Each `jobs/` folder contains:", "",
            "- Original and decoded parameter JSON, a readable parameter table and plain text settings.",
            "- The recorded command (`command.txt`), separate tool/job/combined stdout and stderr files, and log availability notes.",
            "- Input/output identifiers and mapping to history numbers and downloaded files.",
            "- Runtime metrics, timing/status, versions, dependencies and job messages returned by Galaxy.",
            "- A complete original job JSON response.", "",
            "`tool_definitions/` contains wrapper source, parameter definitions/defaults, declared requirements and citations where exposed by Galaxy.", "",
            "Actual job parameters are authoritative for the run. Tool definition defaults also include unselected branches and are stored separately.",
            "A missing command or empty log is recorded explicitly; no unavailable command, log, container digest or package version has been reconstructed.",
            "Declared wrapper requirements describe the tool recipe and do not prove a complete runtime package inventory.",
            "Commands can contain server-specific paths and are retained as evidence, not local scripts to execute unchanged.", "",
            "| Tool | Wrapper version | State | History outputs | Command | Folder |", "|---|---|---|---|---|---|"]
    for row in rows:
        text.append("| %s | %s | %s | %s | %s | %s |" % (row["tool"], row["wrapper_version"], row["state"], row["output_history_numbers"], "saved" if row["command_available"] else "not provided", row["folder"]))
    (ROOT / "README.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    checksums = []
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.tsv":
            with path.open("rb") as source:
                digest = hashlib.file_digest(source, "sha256").hexdigest()
            checksums.append({"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": digest})
    tsv(ROOT / "SHA256SUMS.tsv", checksums, ["path", "bytes", "sha256"])
    archive = ROOT.parent / (ROOT.name + ".zip")
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as output:
        for path in sorted(ROOT.rglob("*")):
            if path.is_file():
                output.write(path, arcname=(Path(ROOT.name) / path.relative_to(ROOT)).as_posix())
    print(json.dumps({"event": "provenance_export_finished", "summary": summary, "archive": str(archive), "files": len(checksums) + 1}), flush=True)
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
