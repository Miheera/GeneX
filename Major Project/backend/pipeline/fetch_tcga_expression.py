"""Fetch TCGA breast-cancer expression scores for genes in the current PPI graph.

Pipeline
--------
1. Read unique Hugo gene symbols from ``data/breast_cancer_graph.csv``.
2. Resolve those symbols to Entrez IDs through the public cBioPortal API.
3. Discover and explicitly select the TCGA RNA-seq z-score molecular profile.
4. Fetch raw and z-score RNA-seq values for all mapped genes across the TCGA
   BRCA PanCancer Atlas sample list.
5. Write one per-gene summary row to ``data/processed/tcga_expression_scores.csv``.

``pct_overexpressed`` currently uses a per-gene descriptive cutoff of that
gene's median expression plus one population standard deviation.  It is a
placeholder rather than a clinically validated overexpression threshold.

Only the Python standard library is used.  Run this script manually when a
fresh API-derived output is wanted; it makes live requests to cBioPortal.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from statistics import median, pstdev
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

LOGGER = logging.getLogger(__name__)

API_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "brca_tcga_pan_can_atlas_2018"
SAMPLE_LIST_ID = "brca_tcga_pan_can_atlas_2018_all"
RAW_EXPRESSION_PROFILE_ID = "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
REQUEST_DELAY_SECONDS = 0.2
REQUEST_TIMEOUT_SECONDS = 30
OUTPUT_COLUMNS = [
    "gene_symbol",
    "entrez_id",
    "median_expression",
    "median_zscore",
    "pct_overexpressed",
    "n_samples",
]
# The graph retains the familiar breast-cancer label HER2, while cBioPortal
# identifies this gene by its current Hugo symbol, ERBB2.
GENE_SYMBOL_ALIASES = {"HER2": "ERBB2"}

_BACKEND = Path(__file__).parents[1]
DEFAULT_GRAPH_PATH = _BACKEND / "data" / "processed" / "expanded_graph_t900.csv"
DEFAULT_OUTPUT_PATH = _BACKEND / "data" / "processed" / "tcga_expression_scores.csv"


def load_graph_genes(graph_path: Path) -> list[str]:
    """Return sorted unique gene symbols used by the current graph CSV."""
    if not graph_path.exists():
        raise FileNotFoundError(f"Graph CSV not found: {graph_path}")

    genes: set[str] = set()
    with graph_path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        required = {"gene1", "gene2"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Graph CSV must contain gene1 and gene2 columns")
        for row in reader:
            genes.update(
                symbol.strip()
                for symbol in (row["gene1"], row["gene2"])
                if symbol and symbol.strip()
            )
    return sorted(genes)


def _post_json(path: str, body: Any, *, query: dict[str, str] | None = None) -> Any | None:
    """POST JSON to cBioPortal, retrying one failed request before returning None."""
    url = f"{API_BASE}{path}"
    if query:
        url = f"{url}?{urlencode(query)}"
    payload = json.dumps(body).encode("utf-8")

    for attempt in range(2):
        if attempt:
            LOGGER.warning("  Retrying cBioPortal request (%d/2): %s", attempt + 1, path)
        try:
            request = Request(
                url,
                data=payload,
                headers={
                    "Accept": "application/json",
                    "Content-Type": "application/json",
                    "User-Agent": "GeneX-TCGA-expression-pipeline/1.0",
                },
                method="POST",
            )
            with urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                if response.status != 200:
                    raise RuntimeError(f"unexpected HTTP status {response.status}")
                result = json.load(response)
            time.sleep(REQUEST_DELAY_SECONDS)
            return result
        except (HTTPError, URLError, TimeoutError, OSError, ValueError, RuntimeError) as exc:
            LOGGER.warning("  cBioPortal request failed for %s: %s", path, exc)
            if attempt == 0:
                time.sleep(REQUEST_DELAY_SECONDS)

    LOGGER.warning("  Skipping request after retry failure: %s", path)
    return None


def _get_json(path: str) -> Any | None:
    """GET JSON from cBioPortal, retrying one failed request before returning None."""
    url = f"{API_BASE}{path}"
    for attempt in range(2):
        if attempt:
            LOGGER.warning("  Retrying cBioPortal request (%d/2): %s", attempt + 1, path)
        try:
            request = Request(
                url,
                headers={"Accept": "application/json", "User-Agent": "GeneX-TCGA-expression-pipeline/1.0"},
                method="GET",
            )
            with urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                if response.status != 200:
                    raise RuntimeError(f"unexpected HTTP status {response.status}")
                result = json.load(response)
            time.sleep(REQUEST_DELAY_SECONDS)
            return result
        except (HTTPError, URLError, TimeoutError, OSError, ValueError, RuntimeError) as exc:
            LOGGER.warning("  cBioPortal request failed for %s: %s", path, exc)
            if attempt == 0:
                time.sleep(REQUEST_DELAY_SECONDS)
    LOGGER.warning("  Skipping request after retry failure: %s", path)
    return None


def discover_zscore_profiles() -> list[dict[str, str]]:
    """Return and print TCGA BRCA mRNA-expression profiles whose labels indicate z-scores."""
    response = _get_json(f"/studies/{STUDY_ID}/molecular-profiles")
    if not isinstance(response, list):
        raise ValueError("Could not retrieve molecular profiles from cBioPortal")

    candidates = []
    for profile in response:
        if not isinstance(profile, dict) or profile.get("molecularAlterationType") != "MRNA_EXPRESSION":
            continue
        profile_id = str(profile.get("molecularProfileId") or "").strip()
        name = str(profile.get("name") or "").strip()
        normalized = f"{profile_id}{name}".lower().replace("-", "").replace("_", "").replace(" ", "")
        if "zscore" not in normalized and "zscores" not in normalized:
            continue
        candidates.append({
            "molecular_profile_id": profile_id,
            "name": name,
            "description": str(profile.get("description") or "").strip(),
        })

    print("\n=== cBioPortal TCGA BRCA mRNA z-score profile candidates ===")
    if candidates:
        for candidate in candidates:
            print(f"  ID:          {candidate['molecular_profile_id']}")
            print(f"  Name:        {candidate['name']}")
            if candidate["description"]:
                print(f"  Description: {candidate['description']}")
    else:
        print("  No matching MRNA_EXPRESSION z-score profiles found.")
    print()
    return candidates


def fetch_entrez_ids(gene_symbols: list[str]) -> dict[str, int]:
    """Resolve Hugo gene symbols to Entrez IDs, warning for unmapped symbols."""
    requested_symbols = [GENE_SYMBOL_ALIASES.get(symbol.upper(), symbol) for symbol in gene_symbols]
    response = _post_json(
        "/genes/fetch",
        requested_symbols,
        query={"geneIdType": "HUGO_GENE_SYMBOL"},
    )
    if response is None:
        LOGGER.warning("Unable to resolve Entrez IDs; no expression scores will be produced")
        return {}
    if not isinstance(response, list):
        LOGGER.warning("Unexpected cBioPortal gene response format; no expression scores will be produced")
        return {}

    mappings: dict[str, int] = {}
    for item in response:
        if not isinstance(item, dict):
            continue
        symbol = str(item.get("hugoGeneSymbol") or "").strip()
        entrez_id = item.get("entrezGeneId")
        if not symbol or not isinstance(entrez_id, int):
            continue
        mappings[symbol.upper()] = entrez_id

    resolved = {
        symbol.upper(): mappings[GENE_SYMBOL_ALIASES.get(symbol.upper(), symbol.upper())]
        for symbol in gene_symbols
        if GENE_SYMBOL_ALIASES.get(symbol.upper(), symbol.upper()) in mappings
    }
    unmapped = [symbol for symbol in gene_symbols if symbol.upper() not in resolved]
    if unmapped:
        LOGGER.warning("  %d graph genes did not map to Entrez IDs: %s",
                       len(unmapped), ", ".join(unmapped))
    LOGGER.info("  Resolved %d of %d graph genes to Entrez IDs", len(resolved), len(gene_symbols))
    return resolved


def fetch_expression_values(profile_id: str, entrez_ids: list[int]) -> dict[int, list[float]]:
    """Fetch one molecular profile's values and group numeric results by Entrez ID."""
    response = _post_json(
        f"/molecular-profiles/{profile_id}/molecular-data/fetch",
        {"sampleListId": SAMPLE_LIST_ID, "entrezGeneIds": entrez_ids},
    )
    if response is None:
        LOGGER.warning("Unable to fetch TCGA expression data; no expression scores will be produced")
        return {}
    if not isinstance(response, list):
        LOGGER.warning("Unexpected cBioPortal molecular-data response format")
        return {}

    values_by_entrez: dict[int, list[float]] = defaultdict(list)
    for item in response:
        if not isinstance(item, dict):
            continue
        entrez_id = item.get("entrezGeneId")
        value = item.get("value")
        if isinstance(entrez_id, int) and isinstance(value, (int, float)):
            values_by_entrez[entrez_id].append(float(value))
    return dict(values_by_entrez)


def score_expression(
    gene_symbols: list[str],
    entrez_by_symbol: dict[str, int],
    raw_values_by_entrez: dict[int, list[float]],
    zscore_values_by_entrez: dict[int, list[float]],
) -> list[dict[str, str | int | float]]:
    """Summarize raw and z-score values for each mapped graph gene."""
    rows: list[dict[str, str | int | float]] = []
    for symbol in gene_symbols:
        entrez_id = entrez_by_symbol.get(symbol.upper())
        if entrez_id is None:
            continue
        raw_values = raw_values_by_entrez.get(entrez_id, [])
        zscore_values = zscore_values_by_entrez.get(entrez_id, [])
        if not raw_values or not zscore_values:
            LOGGER.warning("  Missing raw or z-score expression values for %s (Entrez %d)", symbol, entrez_id)
            continue
        expression_median = median(raw_values)
        zscore_median = median(zscore_values)
        overexpression_cutoff = expression_median + pstdev(raw_values)
        pct_overexpressed = 100.0 * sum(value > overexpression_cutoff for value in raw_values) / len(raw_values)
        rows.append({
            "gene_symbol": symbol,
            "entrez_id": entrez_id,
            "median_expression": round(expression_median, 6),
            "median_zscore": round(zscore_median, 6),
            "pct_overexpressed": round(pct_overexpressed, 6),
            "n_samples": len(zscore_values),
        })
    return rows


def write_scores(output_path: Path, rows: list[dict[str, str | int | float]]) -> None:
    """Write TCGA expression summaries using the project's processed-data convention."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def fetch_tcga_expression(
    graph_path: Path = DEFAULT_GRAPH_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
    *,
    zscore_profile_id: str,
) -> list[dict[str, str | int | float]]:
    """Fetch and summarize TCGA BRCA raw and z-score expression for graph genes."""
    gene_symbols = load_graph_genes(graph_path)
    LOGGER.info("Loaded %d unique graph genes from %s", len(gene_symbols), graph_path.name)

    LOGGER.info("Resolving graph genes through cBioPortal study %s ...", STUDY_ID)
    entrez_by_symbol = fetch_entrez_ids(gene_symbols)
    if not entrez_by_symbol:
        return []

    entrez_ids = sorted(set(entrez_by_symbol.values()))
    LOGGER.info("Fetching raw RNA-seq expression values for %d mapped genes ...", len(entrez_ids))
    raw_values_by_entrez = fetch_expression_values(RAW_EXPRESSION_PROFILE_ID, entrez_ids)
    LOGGER.info("Fetching z-score RNA-seq expression values from %s ...", zscore_profile_id)
    zscore_values_by_entrez = fetch_expression_values(zscore_profile_id, entrez_ids)
    rows = score_expression(gene_symbols, entrez_by_symbol, raw_values_by_entrez, zscore_values_by_entrez)
    write_scores(output_path, rows)

    print("\n=== fetch_tcga_expression summary ===")
    print(f"  Study:                 {STUDY_ID}")
    print(f"  Graph genes:           {len(gene_symbols)}")
    print(f"  Entrez-mapped genes:   {len(entrez_by_symbol)}")
    print(f"  Genes with expression: {len(rows)}")
    print(f"  Z-score profile:       {zscore_profile_id}")
    print("  Overexpression cutoff: per-gene median + 1 population std (placeholder)")
    print(f"  Output:                {output_path}")
    print()
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-path", type=Path, default=DEFAULT_GRAPH_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--discover-zscore-profile", action="store_true",
                        help="Print matching mRNA z-score profiles and exit without fetching expression data")
    parser.add_argument("--zscore-profile-id",
                        help="A reviewed mRNA-expression z-score profile ID returned by --discover-zscore-profile")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.discover_zscore_profile:
        discover_zscore_profiles()
        return
    if not args.zscore_profile_id:
        parser.error("--zscore-profile-id is required; run --discover-zscore-profile first")
    fetch_tcga_expression(
        graph_path=args.graph_path,
        output_path=args.output_path,
        zscore_profile_id=args.zscore_profile_id,
    )


if __name__ == "__main__":
    main()
