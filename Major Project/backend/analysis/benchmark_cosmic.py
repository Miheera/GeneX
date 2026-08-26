"""Benchmark network rankings against filtered COSMIC breast genes.

Graph-universe construction provenance is not retained by the available
pipeline, so output rows use the required cautious ``Unknown:`` caveat unless
the caller explicitly supplies ``graph_universe_includes_cosmic``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from pathlib import Path
from typing import Literal

import networkx as nx

try:
    import graph_loader
    import perturbation
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).parents[1]))
    import graph_loader
    import perturbation

RankMethod = Literal["betweenness", "degree", "impact_score"]
OUTPUT_COLUMNS = [
    "rank_by", "top_n", "precision", "recall", "f1", "n_positives_in_graph",
    "n_positives_total", "n_positives_excluded_from_graph", "methodology_note",
]
UNKNOWN_CAVEAT = (
    "Unknown: graph-universe construction provenance does not establish whether "
    "COSMIC genes were included. Interpret this benchmark cautiously."
)
KNOWN_INCLUDED_CAVEAT = (
    "COSMIC genes were used in constructing the graph universe. This benchmark therefore does not measure "
    "de novo discovery of COSMIC drivers. It measures whether internal network structure ranks known drivers "
    "above non-drivers among genes already represented in the graph; recall may be inflated because known "
    "COSMIC genes were preferentially included."
)
KNOWN_EXCLUDED_CAVEAT = "COSMIC was not used to build the graph; this circularity does not apply based on the available pipeline configuration."
LOGGER = logging.getLogger(__name__)


def _methodology_note(includes_cosmic: bool | None) -> str:
    if includes_cosmic is None:
        return UNKNOWN_CAVEAT
    return KNOWN_INCLUDED_CAVEAT if includes_cosmic else KNOWN_EXCLUDED_CAVEAT


def _read_cosmic(path: Path) -> set[str]:
    if not path.exists():
        raise FileNotFoundError(f"COSMIC processed output does not exist: {path}")
    required = {"gene_symbol", "breast_annotated"}
    with path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("COSMIC processed output must contain gene_symbol and breast_annotated columns")
        positives = set()
        for row in reader:
            if not row.get("gene_symbol", "").strip():
                raise ValueError("COSMIC processed output contains a blank gene_symbol")
            if _parse_bool(row["breast_annotated"]):
                positives.add(row["gene_symbol"].strip())
    return positives


def _parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid boolean value in COSMIC output: {value}")


def _rankings(graph: nx.Graph) -> dict[str, dict[str, float]]:
    path_weight = "path_weight" if any("path_weight" in data for _, _, data in graph.edges(data=True)) else None
    betweenness = nx.betweenness_centrality(graph, weight=path_weight, normalized=True)
    degree = {node: float(graph.degree(node)) for node in graph}
    perturbation.precompute_baseline(graph)
    impact = {
        node: float(perturbation.run_perturbation(graph, str(node), "knockout", 0.0)["impact_score"])
        for node in graph
    }
    return {"betweenness": betweenness, "degree": degree, "impact_score": impact}


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "precision": f"{row['precision']:.6f}", "recall": f"{row['recall']:.6f}", "f1": f"{row['f1']:.6f}"})


def benchmark_cosmic(
    graph: nx.Graph,
    cosmic_path: Path,
    top_n_values: list[int] | None = None,
    *,
    output_path: Path | None = None,
    graph_universe_includes_cosmic: bool | None = None,
) -> list[dict]:
    """Benchmark all ranking methods against breast-annotated COSMIC genes."""
    top_n_values = [10, 20, 50] if top_n_values is None else top_n_values
    if not top_n_values or any(type(value) is not int or value <= 0 for value in top_n_values):
        raise ValueError("top_n_values must contain positive integers")
    if output_path is None:
        output_path = Path(__file__).parents[1] / "data" / "analysis_outputs" / "cosmic_benchmark.csv"
    positives_total = _read_cosmic(cosmic_path)
    graph_nodes = set(graph.nodes)
    positives = positives_total & graph_nodes
    LOGGER.info(
        "COSMIC breast positives total=%d in_graph=%d excluded=%d graph_nodes=%d",
        len(positives_total), len(positives), len(positives_total - positives), len(graph_nodes),
    )
    rankings = _rankings(graph)
    note = _methodology_note(graph_universe_includes_cosmic)
    rows: list[dict] = []
    for method in ("betweenness", "degree", "impact_score"):
        ordered = sorted(rankings[method], key=lambda node: (-rankings[method][node], str(node)))
        for top_n in top_n_values:
            selected = set(ordered[:top_n])
            hits = len(selected & positives)
            precision = hits / len(selected) if selected else 0.0
            recall = hits / len(positives) if positives else 0.0
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
            rows.append({
                "rank_by": method, "top_n": top_n, "precision": precision, "recall": recall, "f1": f1,
                "n_positives_in_graph": len(positives), "n_positives_total": len(positives_total),
                "n_positives_excluded_from_graph": len(positives_total - positives), "methodology_note": note,
            })
    _write_csv(output_path, rows)
    return rows


def _load_graph(path: Path | None) -> nx.Graph:
    if path is None:
        return graph_loader.load_graph()
    previous = os.environ.get("GRAPH_CSV_PATH")
    os.environ["GRAPH_CSV_PATH"] = str(path)
    try:
        return graph_loader.load_graph()
    finally:
        if previous is None:
            os.environ.pop("GRAPH_CSV_PATH", None)
        else:
            os.environ["GRAPH_CSV_PATH"] = previous


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cosmic-path", type=Path, default=Path(__file__).parents[1] / "data" / "processed" / "cosmic_breast_genes.csv")
    parser.add_argument("--top-n", type=lambda value: [int(item) for item in value.split(",")], default=[10, 20, 50])
    parser.add_argument("--graph-path", type=Path)
    parser.add_argument("--output-path", type=Path, default=Path(__file__).parents[1] / "data" / "analysis_outputs" / "cosmic_benchmark.csv")
    parser.add_argument("--cosmic-used-in-graph", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    benchmark_cosmic(_load_graph(args.graph_path), args.cosmic_path, args.top_n,
                     output_path=args.output_path,
                     graph_universe_includes_cosmic=True if args.cosmic_used_in_graph else None)


if __name__ == "__main__":
    main()
