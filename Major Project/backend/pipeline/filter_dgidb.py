"""Filter a local DGIdb interaction export to genes in the PPI graph.

The bundled DGIdb v5.0.11 export uses ``gene_name``, ``drug_name``,
``drug_is_approved``, and ``interaction_source_db_name``.  The parser validates
those fields rather than relying on a schema from a different DGIdb release.
"""

from __future__ import annotations

import argparse
import csv
import logging
from collections import defaultdict
from itertools import chain
from pathlib import Path

try:
    from .map_identifiers import resolve_gene_symbol
except ImportError:
    from map_identifiers import resolve_gene_symbol

LOGGER = logging.getLogger(__name__)
REQUIRED_COLUMNS = {
    "gene_name", "drug_name", "drug_is_approved", "interaction_source_db_name",
}
OUTPUT_COLUMNS = ["gene_symbol", "n_interactions", "drug_names", "has_approved_drug", "n_sources", "sources"]
UNMATCHED_COLUMNS = ["raw_gene_symbol", "resolved_gene_symbol", "n_rows", "reason"]


def _read_rows(path: Path):
    """Yield rows after skipping DGIdb metadata comment lines."""
    with path.open(newline="", encoding="utf-8-sig") as input_file:
        header = None
        for line in input_file:
            if line.startswith("#") or not line.strip():
                continue
            header = line
            break
        if header is None:
            raise ValueError(f"DGIdb input has no tabular header: {path}")
        reader = csv.DictReader(chain([header], input_file), delimiter="\t")
        available = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_COLUMNS - available)
        if missing:
            raise ValueError(
                "DGIdb input is missing required columns "
                f"{missing}; available columns are {sorted(available)}"
            )
        yield from reader


def _parse_approved(value: str) -> bool:
    """Parse the v5 export's explicit lowercase true/false approval flag."""
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no", ""}:
        return False
    raise ValueError(f"Unexpected DGIdb drug_is_approved value: {value!r}")


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def filter_dgidb(
    input_path: Path,
    graph_genes: set[str] | list[str],
    *,
    output_path: Path | None = None,
    unmatched_output_path: Path | None = None,
) -> list[dict]:
    """Return per-graph-gene DGIdb annotations and write stable CSV outputs."""
    if not input_path.exists():
        raise FileNotFoundError(f"DGIdb input does not exist: {input_path}")
    graph_genes = {str(gene) for gene in graph_genes}
    if output_path is None:
        output_path = Path(__file__).parents[1] / "data" / "processed" / "dgidb_gene_annotations.csv"
    if unmatched_output_path is None:
        unmatched_output_path = Path(__file__).parents[1] / "data" / "processed" / "dgidb_unmatched_genes.csv"

    aggregates = defaultdict(lambda: {"records": set(), "drugs": set(), "approved": False, "sources": set()})
    unmatched = defaultdict(int)
    for row in _read_rows(input_path):
        raw_gene = (row["gene_name"] or "").strip()
        if not raw_gene:
            continue
        resolved = resolve_gene_symbol(raw_gene, graph_genes)
        if resolved not in graph_genes:
            unmatched[(raw_gene, resolved, "not present in graph after synonym resolution")] += 1
            continue
        drug_name = (row["drug_name"] or "").strip()
        drug_identifier = (row.get("drug_concept_id") or "").strip() or drug_name.upper()
        interaction_types = (row.get("interaction_types") or "").strip()
        source = (row["interaction_source_db_name"] or "").strip()
        aggregate = aggregates[resolved]
        # One claim may appear in several source databases; source is therefore
        # deliberately excluded from this record identity.
        aggregate["records"].add((resolved, drug_identifier, interaction_types))
        if drug_name:
            aggregate["drugs"].add(drug_name)
        aggregate["approved"] = aggregate["approved"] or _parse_approved(row["drug_is_approved"] or "")
        if source:
            aggregate["sources"].add(source)

    annotations = []
    for gene in sorted(graph_genes):
        aggregate = aggregates[gene]
        annotations.append({
            "gene_symbol": gene,
            "n_interactions": len(aggregate["records"]),
            "drug_names": "|".join(sorted(aggregate["drugs"])),
            "has_approved_drug": bool(aggregate["approved"]),
            "n_sources": len(aggregate["sources"]),
            "sources": "|".join(sorted(aggregate["sources"])),
        })
    unmatched_rows = [
        {"raw_gene_symbol": raw, "resolved_gene_symbol": resolved, "n_rows": count, "reason": reason}
        for (raw, resolved, reason), count in sorted(unmatched.items())
    ]
    _write_csv(output_path, OUTPUT_COLUMNS, annotations)
    _write_csv(unmatched_output_path, UNMATCHED_COLUMNS, unmatched_rows)
    LOGGER.info("DGIdb annotations: graph_genes=%d matched=%d unmatched_rows=%d", len(graph_genes), len(aggregates), sum(unmatched.values()))
    return annotations


def _default_input_path() -> Path:
    backend = Path(__file__).parents[1]
    requested = backend / "data" / "raw" / "dgidb" / "interactions.tsv"
    return requested if requested.exists() else backend / "data" / "interactions.tsv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-path", type=Path, default=_default_input_path())
    parser.add_argument("--output-path", type=Path, default=Path(__file__).parents[1] / "data" / "processed" / "dgidb_gene_annotations.csv")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    try:
        import graph_loader
    except ModuleNotFoundError:
        import sys
        sys.path.insert(0, str(Path(__file__).parents[1]))
        import graph_loader
    graph = graph_loader.load_graph()
    filter_dgidb(args.input_path, set(graph), output_path=args.output_path)


if __name__ == "__main__":
    main()
