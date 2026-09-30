#!/usr/bin/env python3
"""Write SHA-256 provenance for the frozen graph and audit inputs."""
from datetime import datetime, timezone
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "validation"
META = VAL / "metadata"
GRAPH = ROOT / "backend/data/processed/expanded_graph_t900.csv"
URLS = {
    "GSE92742_": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/",
    "GSE70138_": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/",
    "hgnc_complete_set.txt": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "depmap_catalog_verification_response.html": "https://depmap.org/portal/api/download/files (returned verification HTML)",
}

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    now = datetime.now(timezone.utc).isoformat()
    records = []
    for path in [GRAPH, *sorted(META.iterdir())]:
        if not path.is_file():
            continue
        name = path.name
        url = next((v + name for k, v in URLS.items() if name.startswith(k)), "")
        if name == "depmap_catalog_verification_response.html":
            url = URLS[name]
        records.append({
            "relative_path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha(path),
            "source_url": url or "local frozen GenX graph",
            "retrieved_utc": now if path.parent == META else "frozen graph; checksum captured during audit",
        })
    target = VAL / "metadata/input_manifest.csv"
    with target.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)
    print(f"Wrote {target} with {len(records)} records")

if __name__ == "__main__":
    main()
