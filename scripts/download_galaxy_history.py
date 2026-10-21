"""Download only non-deleted Galaxy history datasets; credentials stay in memory."""

from concurrent.futures import ThreadPoolExecutor, as_completed, wait, FIRST_COMPLETED
from collections import defaultdict
from datetime import datetime, timedelta, timezone
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import threading
import time
from urllib.parse import quote, urljoin, urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1] / "outputs/galaxy_history"
BASE = "https://usegalaxy.eu"
HISTORY_ID = "11ac94870d0bb33a834c81ccf4f16ef2"
API_KEY = ""
LOCAL_TIMEZONE = timezone(timedelta(hours=6))
LOCAL = threading.local()
LOCK = threading.Lock()
PROGRESS = {}
RESULTS = {}


def now():
    return datetime.now(LOCAL_TIMEZONE).isoformat()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def session():
    if not hasattr(LOCAL, "session"):
        LOCAL.session = requests.Session()
        retry = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
        LOCAL.session.mount("https://", HTTPAdapter(max_retries=retry))
    return LOCAL.session


def get(endpoint, *, stream=False, headers=None, params=None):
    url = urljoin(BASE, endpoint)
    for _ in range(6):
        if urlparse(url).scheme != "https":
            raise ValueError("Refusing a non-HTTPS download")
        request_headers = {"Accept-Encoding": "identity"} if stream else {}
        if urlparse(url).netloc == urlparse(BASE).netloc:
            request_headers["x-api-key"] = API_KEY
        request_headers.update(headers or {})
        response = session().get(url, headers=request_headers, params=params, stream=stream,
                                 allow_redirects=False, timeout=(30, 120))
        params = None
        if response.status_code in (301, 302, 303, 307, 308):
            url = urljoin(url, response.headers["Location"])
            response.close()
            continue
        response.raise_for_status()
        return response
    raise RuntimeError("Too many redirects")


def api(endpoint):
    for attempt in range(1, 6):
        try:
            with get(endpoint) as response:
                return response.json()
        except (requests.RequestException, ValueError):
            if attempt == 5:
                raise
            time.sleep(attempt * 2)


def safe_name(value, limit=150):
    value = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", str(value))
    value = re.sub(r"\s+", " ", value).strip(" .")[:limit].rstrip(" .") or "unnamed"
    if value.split(".")[0].upper() in {"CON", "PRN", "AUX", "NUL", *("COM%d" % n for n in range(1, 10)), *("LPT%d" % n for n in range(1, 10))}:
        value = "_" + value
    return value


def dataset_path(dataset):
    extension = dataset.get("extension") or dataset.get("file_ext") or "dat"
    extension = {"tabular": "tsv", "fastqsanger.gz": "fastq.gz", "fastqsanger": "fastq", "fasta.gz": "fasta.gz"}.get(extension, extension)
    name = safe_name(dataset["name"])
    if not name.lower().endswith("." + extension.lower()):
        name += "." + safe_name(extension, 25)
    return ROOT / "datasets" / ("%03d_%s" % (dataset["hid"], name))


def hash_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_sha256(dataset):
    for item in dataset.get("hashes") or []:
        if not item.get("extra_files_path") and item.get("hash_function", "").upper().replace("-", "") == "SHA256":
            return item["hash_value"].lower()
    return None


def transfer(endpoint, path, expected_size=None, expected_hash=None, progress_key=None, params=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        digest = hash_file(path)
        if expected_size is not None and path.stat().st_size != expected_size:
            raise ValueError("An existing file has an unexpected size")
        if expected_hash and digest != expected_hash:
            raise ValueError("An existing file has an unexpected checksum")
        return {"bytes": path.stat().st_size, "sha256": digest, "reused": True}
    partial = path.with_name(path.name + ".part")
    for attempt in range(1, 4):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            if expected_size is not None and offset > expected_size:
                raise ValueError("Partial download is larger than the source")
            if partial.exists() and expected_size is not None and offset == expected_size:
                pass
            else:
                range_headers = {"Range": "bytes=%d-" % offset} if offset else {}
                with get(endpoint, stream=True, headers=range_headers, params=params) as response:
                    if offset and response.status_code == 206:
                        if not response.headers.get("Content-Range", "").startswith("bytes %d-" % offset):
                            raise ValueError("Server returned an unexpected byte range")
                        mode = "ab"
                    else:
                        offset, mode = 0, "wb"
                    with partial.open(mode) as output:
                        for block in response.iter_content(chunk_size=64 * 1024):
                            if block:
                                output.write(block)
                                offset += len(block)
                                if progress_key is not None:
                                    with LOCK:
                                        PROGRESS[progress_key] = offset
                if expected_size is not None and offset != expected_size:
                    raise ValueError("Download size differs from Galaxy metadata")
            digest = hash_file(partial)
            if expected_hash and digest != expected_hash:
                partial.unlink()
                raise ValueError("Download checksum differs from Galaxy metadata")
            partial.replace(path)
            return {"bytes": path.stat().st_size, "sha256": digest, "reused": False}
        except Exception:
            if attempt == 3:
                raise
            time.sleep(attempt * 2)


def download_group(group):
    for item in group:
        if item.get("deleted") or item.get("purged"):
            raise ValueError("Deleted or purged dataset excluded from downloads")
    source = min(group, key=lambda item: item["hid"])
    target = dataset_path(source)
    download_params = {"raw": "true", "preview": "false"} if source.get("extension") == "html" else {"to_ext": "data", "preview": "false"}
    result = transfer(source.get("download_url") or "/api/datasets/%s/display" % source["id"], target,
                      int(source.get("file_size") or 0), expected_sha256(source), source["id"],
                      download_params)
    records = []
    for item in group:
        destination = dataset_path(item)
        if destination != target:
            if not destination.exists():
                shutil.copy2(target, destination)
            if destination.stat().st_size != result["bytes"] or hash_file(destination) != result["sha256"]:
                raise ValueError("Duplicate history entry failed verification")
        records.append({
            "hid": item["hid"], "name": item["name"], "dataset_id": item["id"],
            "underlying_dataset_id": item.get("dataset_id"), "state": item.get("state"),
            "visible": item.get("visible"), "deleted": False, "extension": item.get("extension"),
            "expected_bytes": int(item.get("file_size") or 0), "bytes": result["bytes"],
            "sha256": result["sha256"], "source_sha256": expected_sha256(item),
            "source_sha256_verified": bool(expected_sha256(item)),
            "relative_path": destination.relative_to(ROOT).as_posix(), "status": "complete",
            "downloaded_via_hid": source["hid"],
        })
    return records


def get_extras(dataset):
    cache_path = ROOT / "metadata" / "extra_files" / ("%03d.json" % dataset["hid"])
    entries = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else api("/api/datasets/%s/extra_files" % dataset["id"])
    write_json(cache_path, entries)
    records = []
    for item in entries:
        raw_path = item.get("path", "")
        parts = PurePosixPath(raw_path.replace("\\", "/"))
        if parts.is_absolute() or ".." in parts.parts:
            raise ValueError("Unsafe extra file path")
        target = ROOT / "extra_files" / ("%03d" % dataset["hid"])
        target = target.joinpath(*(safe_name(part) for part in parts.parts))
        if item.get("class") == "Directory":
            target.mkdir(parents=True, exist_ok=True)
        elif item.get("class") == "File":
            result = transfer("/api/datasets/%s/display" % dataset["id"], target,
                              params={"raw": "true", "preview": "false", "filename": raw_path})
            records.append({"dataset_hid": dataset["hid"], "original_path": raw_path,
                            "relative_path": target.relative_to(ROOT).as_posix(), **result})
    return records


def collect_members(collection, datasets_by_id):
    cache_path = ROOT / "metadata" / "collections" / ("%03d.json" % collection["hid"])
    detail = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else api("/api/histories/%s/contents/dataset_collections/%s" % (HISTORY_ID, collection["id"]))
    if detail.get("deleted"):
        raise ValueError("Collection became deleted during the download")
    write_json(ROOT / "metadata" / "collections" / ("%03d.json" % collection["hid"]), detail)
    members = []

    def walk(elements, parents):
        for element in elements:
            element_path = parents + [str(element.get("element_identifier", ""))]
            obj = element.get("object") or {}
            if element.get("element_type") == "hda":
                member = datasets_by_id.get(obj.get("id"))
                if not member:
                    if obj.get("deleted"):
                        continue
                    raise ValueError("Collection member is missing from the history inventory")
                members.append({"collection_hid": detail["hid"], "collection_name": detail["name"],
                                "element_path": "/".join(element_path), "dataset_hid": member["hid"],
                                "dataset_id": member["id"], "relative_path": dataset_path(member).relative_to(ROOT).as_posix()})
            else:
                if "elements" not in obj:
                    raise ValueError("Nested collection was not fully expanded by Galaxy")
                walk(obj["elements"], element_path)
    walk(detail.get("elements") or [], [])
    return members


def write_tsv(path, records, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def main():
    global ROOT, BASE, HISTORY_ID, API_KEY
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--outdir", type=Path, default=ROOT)
    cli.add_argument("--base-url", default=BASE)
    cli.add_argument("--history-id", default=HISTORY_ID)
    cli.add_argument("--inventory-only", action="store_true")
    args = cli.parse_args()
    ROOT, BASE, HISTORY_ID = args.outdir.resolve(), args.base_url.rstrip("/"), args.history_id
    if not BASE.startswith("https://"):
        raise SystemExit("An HTTPS Galaxy server is required")
    API_KEY = os.environ.pop("GALAXY_API_KEY", "")
    if not API_KEY:
        raise SystemExit("Set GALAXY_API_KEY in the environment")
    ROOT.mkdir(parents=True, exist_ok=True)
    history = api("/api/histories/" + HISTORY_ID)
    if history["name"] != "Manual Analysis S-11":
        raise ValueError("History name changed")
    all_contents = api("/api/histories/%s/contents?details=all" % HISTORY_ID)
    contents = [item for item in all_contents if not item.get("deleted")]
    datasets = sorted([item for item in contents if item.get("history_content_type") == "dataset"], key=lambda item: item["hid"])
    collections = [item for item in contents if item.get("history_content_type") == "dataset_collection"]
    excluded_count = sum(item.get("history_content_type") == "dataset" and bool(item.get("deleted")) for item in all_contents)
    if any(item.get("purged") for item in datasets):
        raise ValueError("A current dataset has been purged")
    write_json(ROOT / "metadata" / "history.json", history)
    write_json(ROOT / "metadata" / "history_contents.json", contents)
    if args.inventory_only:
        print(json.dumps({"history": history["name"], "non_deleted_datasets": len(datasets),
                          "excluded_deleted_datasets": excluded_count, "dataset_contents_downloaded": False}))
        return
    (ROOT / "datasets").mkdir(exist_ok=True)
    for item in datasets:
        write_json(ROOT / "metadata" / ("dataset_%03d.json" % item["hid"]), item)
    groups = defaultdict(list)
    for item in datasets:
        groups[item.get("dataset_id") or item["id"]].append(item)
    expected_total = sum(int(item.get("file_size") or 0) for item in datasets)
    if shutil.disk_usage(ROOT).free < expected_total + 100 * 1024 * 1024:
        raise ValueError("Not enough free space")
    print(json.dumps({"event": "started", "at": now(), "datasets": len(datasets), "collections": len(collections),
                      "excluded_deleted_datasets": excluded_count, "expected_dataset_bytes": expected_total,
                      "unique_download_bytes": sum(int(group[0].get("file_size") or 0) for group in groups.values())}), flush=True)
    errors = []
    extras = []
    members = []
    ordered_groups = sorted(groups.values(), key=lambda group: int(group[0].get("file_size") or 0))
    with ThreadPoolExecutor(max_workers=8) as download_pool, ThreadPoolExecutor(max_workers=4) as metadata_pool:
        futures = {download_pool.submit(download_group, group): ("dataset", [item["hid"] for item in group]) for group in ordered_groups}
        for item in datasets:
            futures[metadata_pool.submit(get_extras, item)] = ("extra_files", item["hid"])
        by_id = {item["id"]: item for item in datasets}
        for collection in collections:
            futures[metadata_pool.submit(collect_members, collection, by_id)] = ("collection", collection["hid"])
        last_print = 0
        while futures:
            done, _ = wait(futures, timeout=5, return_when=FIRST_COMPLETED)
            for future in done:
                kind, identifier = futures.pop(future)
                try:
                    value = future.result()
                    if kind == "dataset":
                        for item in value:
                            RESULTS[str(item["hid"])] = item
                        write_json(ROOT / "download_state.json", RESULTS)
                        print(json.dumps({"event": "datasets_complete", "hids": identifier, "completed": len(RESULTS), "total": len(datasets)}), flush=True)
                    elif kind == "extra_files":
                        extras.extend(value)
                    else:
                        members.extend(value)
                except Exception as error:
                    record = {"kind": kind, "identifier": identifier, "error_type": type(error).__name__, "message": str(error).replace(API_KEY, "[redacted]")}
                    errors.append(record)
                    print(json.dumps({"event": "error", **record}), flush=True)
            if time.monotonic() - last_print > 20:
                with LOCK:
                    transferred = sum(PROGRESS.values())
                progress = {"at": now(), "datasets_complete": len(RESULTS), "datasets_total": len(datasets),
                            "downloaded_network_bytes": transferred, "pending_operations": len(futures), "errors": len(errors)}
                write_json(ROOT / "progress.json", progress)
                print(json.dumps({"event": "progress", **progress}), flush=True)
                last_print = time.monotonic()

    records = sorted(RESULTS.values(), key=lambda item: item["hid"])
    write_tsv(ROOT / "inventory.tsv", records, ["hid", "name", "relative_path", "extension", "state", "visible", "deleted", "expected_bytes", "bytes", "sha256", "source_sha256", "source_sha256_verified", "downloaded_via_hid", "dataset_id", "underlying_dataset_id", "status"])
    write_tsv(ROOT / "collection_membership.tsv", sorted(members, key=lambda item: (item["collection_hid"], item["element_path"])),
              ["collection_hid", "collection_name", "element_path", "dataset_hid", "relative_path", "dataset_id"])
    summary = {"completed_at": now(), "history_name": history["name"], "history_id": HISTORY_ID,
               "scope": "Only non-deleted datasets, including hidden datasets and collection members. Deleted datasets excluded.",
               "expected_datasets": len(datasets), "saved_datasets": len(records),
               "saved_dataset_bytes": sum(item["bytes"] for item in records), "non_deleted_collections": len(collections),
               "collection_member_links": len(members), "extra_files": extras,
               "excluded_deleted_datasets": excluded_count, "errors": errors,
               "all_downloads_complete": not errors and len(records) == len(datasets)}
    write_json(ROOT / "download_summary.json", summary)
    text = ["Manual Analysis S-11 - downloaded Galaxy history datasets", "",
            "Scope: non-deleted datasets only; hidden datasets and collection members included.",
            "Deleted datasets excluded: %d" % excluded_count,
            "Datasets saved: %d / %d" % (len(records), len(datasets)),
            "Dataset bytes saved: %s" % format(summary["saved_dataset_bytes"], ","),
            "Collection structures saved: %d" % len(collections),
            "Extra files saved: %d" % len(extras),
            "Every downloaded file was checked against Galaxy's file size and assigned a SHA-256 checksum.",
            "Where Galaxy supplies a source SHA-256 checksum, it was checked as well.",
            "", "Files: datasets/", "Inventory and checksums: inventory.tsv", "Collection mapping: collection_membership.tsv",
            "History and dataset metadata: metadata/", "Detailed completion status: download_summary.json",
            "", "Completed at: " + summary["completed_at"], "Complete: " + str(summary["all_downloads_complete"])]
    (ROOT / "README.txt").write_text("\n".join(text) + "\n", encoding="utf-8")
    print(json.dumps({"event": "finished", **summary}), flush=True)
    if not summary["all_downloads_complete"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
