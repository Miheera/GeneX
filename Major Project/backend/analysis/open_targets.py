"""Identify high-impact graph genes with and without DGIdb interaction claims."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import networkx as nx

try:
    import perturbation
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, str(Path(__file__).parents[1]))
    import perturbation

METHODOLOGY_NOTE = (
    "A DGIdb interaction indicates that at least one aggregated source reports a drug-gene interaction claim. "
    "It does not establish clinical efficacy, suitability for this cancer type, or an effective treatment. "
    "Conversely, no DGIdb interaction in this dataset does not mean that a gene is intrinsically undruggable "
    "or untreatable; it only means no interaction claim was found in the imported DGIdb data."
)
OUTPUT_COLUMNS = [
    "group", "rank", "gene_symbol", "structural_score", "ranking_method", "n_interactions",
    "drug_names", "has_approved_drug", "n_sources", "sources", "methodology_note",
]


def _annotations(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"DGIdb annotations are not available: {path}")
    with path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        required = {"gene_symbol", "n_interactions", "drug_names", "has_approved_drug", "n_sources", "sources"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("DGIdb annotations CSV has an invalid schema")
        return {row["gene_symbol"]: row for row in reader}


def _impact_scores(graph: nx.Graph) -> dict[str, float]:
    perturbation.precompute_baseline(graph)
    return {
        str(gene): float(perturbation.run_perturbation(graph, str(gene), "knockout", 0.0)["impact_score"])
        for gene in graph
    }


def open_targets(
    graph: nx.Graph,
    annotation_path: Path,
    top_n: int = 20,
    *,
    output_path: Path | None = None,
) -> dict[str, object]:
    """Rank knockout impact and split the selected genes by DGIdb claim presence."""
    if type(top_n) is not int or top_n <= 0:
        raise ValueError("top_n must be a positive integer")
    if output_path is None:
        output_path = Path(__file__).parents[1] / "data" / "analysis_outputs" / "open_targets.csv"
    annotations = _annotations(annotation_path)
    missing = set(map(str, graph)) - set(annotations)
    if missing:
        raise ValueError(f"DGIdb annotations are missing {len(missing)} graph genes")
    scores = _impact_scores(graph)
    ranked = sorted(scores, key=lambda gene: (-scores[gene], gene))[:top_n]
    rows = []
    for rank, gene in enumerate(ranked, start=1):
        annotation = annotations[gene]
        group = "druggable" if int(annotation["n_interactions"]) > 0 else "open_targets"
        rows.append({
            "group": group, "rank": rank, "gene_symbol": gene, "structural_score": scores[gene],
            "ranking_method": "impact_score", "n_interactions": int(annotation["n_interactions"]),
            "drug_names": annotation["drug_names"], "has_approved_drug": annotation["has_approved_drug"],
            "n_sources": int(annotation["n_sources"]), "sources": annotation["sources"],
            "methodology_note": METHODOLOGY_NOTE,
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return {
        "ranking_method": "impact_score", "top_n": len(ranked),
        "druggable": [row for row in rows if row["group"] == "druggable"],
        "open_targets": [row for row in rows if row["group"] == "open_targets"],
        "methodology_note": METHODOLOGY_NOTE,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--annotations-path", type=Path, default=Path(__file__).parents[1] / "data" / "processed" / "dgidb_gene_annotations.csv")
    parser.add_argument("--output-path", type=Path, default=Path(__file__).parents[1] / "data" / "analysis_outputs" / "open_targets.csv")
    args = parser.parse_args()
    import graph_loader
    open_targets(graph_loader.load_graph(), args.annotations_path, args.top_n, output_path=args.output_path)


if __name__ == "__main__":
    main()
