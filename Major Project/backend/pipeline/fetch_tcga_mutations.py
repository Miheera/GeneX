"""Fetch TCGA BRCA mutation frequency scores for genes in the PPI graph."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any

try:
    from .fetch_tcga_expression import (
        API_BASE,
        GENE_SYMBOL_ALIASES,
        SAMPLE_LIST_ID,
        STUDY_ID,
        _get_json,
        _post_json,
        fetch_entrez_ids,
        load_graph_genes,
    )
except ImportError:
    from fetch_tcga_expression import (
        API_BASE,
        GENE_SYMBOL_ALIASES,
        SAMPLE_LIST_ID,
        STUDY_ID,
        _get_json,
        _post_json,
        fetch_entrez_ids,
        load_graph_genes,
    )

LOGGER = logging.getLogger(__name__)
MUTATION_PROFILE_ID = f"{STUDY_ID}_mutations"
OUTPUT_COLUMNS = ["gene", "mutation_count", "mutation_frequency_pct", "sample_count"]

_BACKEND = Path(__file__).parents[1]
DEFAULT_GRAPH_PATH = _BACKEND / "data" / "processed" / "expanded_graph_t900.csv"
DEFAULT_OUTPUT_PATH = _BACKEND / "data" / "processed" / "tcga_mutation_scores.csv"


def fetch_mutations(profile_id: str, entrez_ids: list[int]) -> dict[int, set[str]]:
    """Return unique mutated sample IDs grouped by Entrez gene ID."""
    response = _post_json(
        f"/molecular-profiles/{profile_id}/mutations/fetch",
        {"sampleListId": SAMPLE_LIST_ID, "entrezGeneIds": entrez_ids},
        query={"pageSize": "10000000"},
    )
    if not isinstance(response, list):
        LOGGER.warning("Unexpected cBioPortal mutation response; no mutation scores produced")
        return {}

    samples_by_entrez: dict[int, set[str]] = {}
    for item in response:
        if not isinstance(item, dict):
            continue
        entrez_id = item.get("entrezGeneId")
        sample_id = str(item.get("sampleId") or "").strip()
        if isinstance(entrez_id, int) and sample_id:
            samples_by_entrez.setdefault(entrez_id, set()).add(sample_id)
    return samples_by_entrez


def score_mutations(
    gene_symbols: list[str],
    entrez_by_symbol: dict[str, int],
    samples_by_entrez: dict[int, set[str]],
    sample_count: int,
) -> list[dict[str, str | int | float]]:
    """Calculate per-gene unique mutated sample counts and frequencies."""
    rows = []
    for symbol in gene_symbols:
        entrez_id = entrez_by_symbol.get(symbol.upper())
        if entrez_id is None:
            continue
        mutation_count = len(samples_by_entrez.get(entrez_id, set()))
        frequency = 100.0 * mutation_count / sample_count if sample_count else 0.0
        rows.append({
            "gene": symbol,
            "mutation_count": mutation_count,
            "mutation_frequency_pct": round(frequency, 6),
            "sample_count": sample_count,
        })
    return rows


def fetch_tcga_mutations(
    graph_path: Path = DEFAULT_GRAPH_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
    *,
    mutation_profile_id: str = MUTATION_PROFILE_ID,
) -> list[dict[str, str | int | float]]:
    """Fetch and write TCGA BRCA mutation frequencies for graph genes."""
    gene_symbols = load_graph_genes(graph_path)
    entrez_by_symbol = fetch_entrez_ids(gene_symbols)
    if not entrez_by_symbol:
        return []
    sample_ids = _get_json(f"/sample-lists/{SAMPLE_LIST_ID}/sample-ids")
    if not isinstance(sample_ids, list):
        raise ValueError("Could not retrieve the cBioPortal sample list")
    samples_by_entrez = fetch_mutations(mutation_profile_id, sorted(set(entrez_by_symbol.values())))
    rows = score_mutations(gene_symbols, entrez_by_symbol, samples_by_entrez, len(sample_ids))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print("\n=== fetch_tcga_mutations summary ===")
    print(f"  Study:       {STUDY_ID}")
    print(f"  Graph genes: {len(gene_symbols)}")
    print(f"  Scored:      {len(rows)}")
    print(f"  Samples:     {len(sample_ids)}")
    print(f"  Output:      {output_path}\n")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-path", type=Path, default=DEFAULT_GRAPH_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--mutation-profile-id", default=MUTATION_PROFILE_ID)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fetch_tcga_mutations(args.graph_path, args.output_path, mutation_profile_id=args.mutation_profile_id)


if __name__ == "__main__":
    main()