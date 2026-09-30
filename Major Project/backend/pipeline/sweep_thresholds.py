"""Run a one-pass STRING confidence-threshold sensitivity sweep."""

from __future__ import annotations

import argparse
import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping

import networkx as nx

try:
    from ..path_cost import require_positive_betweenness_costs
except ImportError:
    from path_cost import require_positive_betweenness_costs

try:
    from .graph_builder import build_graph, write_processed_graph
except ImportError:
    from graph_builder import build_graph, write_processed_graph

LOGGER = logging.getLogger(__name__)
DEFAULT_THRESHOLDS = [400, 500, 600, 700, 800, 900]
DEFAULT_SHIPPING_THRESHOLD = 700
MANIFEST_COLUMNS = [
    "threshold", "status", "error_message", "node_count", "edge_count",
    "mean_degree", "density", "n_connected_components", "is_fully_connected",
    "processed_graph_csv",
]
REPORT_COLUMNS = [
    "threshold", "node_count", "edge_count", "mean_degree", "density",
    "n_connected_components", "is_fully_connected", "top10_betweenness_genes",
    "jaccard_vs_shipping_threshold", "status", "error_message",
]


def validate_thresholds(thresholds: Iterable[int]) -> list[int]:
    """Validate positive, unique, ascending integer thresholds."""
    values = list(thresholds)
    if any(type(value) is not int or value <= 0 for value in values):
        raise ValueError("thresholds must be positive integers")
    if len(set(values)) != len(values):
        raise ValueError("thresholds must be unique")
    if values != sorted(values):
        raise ValueError("thresholds must be sorted in ascending order")
    if not values:
        raise ValueError("at least one threshold is required")
    return values


def load_gene_universe(path: Path) -> set[str]:
    """Load gene symbols from a processed graph CSV."""
    genes: set[str] = set()
    with path.open(newline="", encoding="utf-8") as input_file:
        for row in csv.DictReader(input_file):
            genes.update(value.strip() for value in (row["gene1"], row["gene2"]) if value.strip())
    return genes


def load_string_mapping(path: Path, gene_universe: set[str]) -> dict[str, str]:
    """Map STRING protein IDs to symbols for the requested gene universe."""
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        for row in reader:
            protein_id = row.get("#string_protein_id", "").strip()
            symbol = row.get("preferred_name", "").strip()
            if protein_id and symbol in gene_universe:
                mapping[protein_id] = symbol
    return mapping


def bucket_string_edges(
    links_path: Path,
    mapping: Mapping[str, str],
    thresholds: list[int],
) -> dict[int, list[tuple[str, str, int]]]:
    """Read the raw STRING link file once and bucket qualifying mapped edges."""
    buckets: dict[int, list[tuple[str, str, int]]] = defaultdict(list)
    with links_path.open(encoding="utf-8") as input_file:
        header = input_file.readline().split()
        try:
            protein1_index = header.index("protein1")
            protein2_index = header.index("protein2")
            score_index = header.index("combined_score")
        except ValueError as exc:
            raise ValueError("STRING links file lacks the expected header") from exc

        for line in input_file:
            fields = line.split()
            if len(fields) <= score_index:
                continue
            gene1 = mapping.get(fields[protein1_index])
            gene2 = mapping.get(fields[protein2_index])
            if not gene1 or not gene2 or gene1 == gene2:
                continue
            score = int(fields[score_index])
            for threshold in thresholds:
                if score >= threshold:
                    buckets[threshold].append((gene1, gene2, score))
    return buckets


def _graph_metrics(graph: nx.Graph) -> dict[str, float | int | bool]:
    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    components = nx.number_connected_components(graph) if node_count else 0
    return {
        "node_count": node_count,
        "edge_count": edge_count,
        "mean_degree": 2.0 * edge_count / node_count if node_count else 0.0,
        "density": nx.density(graph) if node_count >= 2 else 0.0,
        "n_connected_components": components,
        "is_fully_connected": node_count > 0 and components == 1,
    }


def _top_betweenness(graph: nx.Graph) -> str:
    if not graph:
        return ""
    require_positive_betweenness_costs(graph)
    centrality = nx.betweenness_centrality(graph, weight="path_weight", normalized=True)
    ordered = sorted(centrality, key=lambda node: (-centrality[node], str(node)))
    return ",".join(str(node) for node in ordered[:10])


def _jaccard(first: set[tuple[str, str]], second: set[tuple[str, str]]) -> float | None:
    union = first | second
    if not union:
        return None
    return len(first & second) / len(union)


def _edge_set(graph: nx.Graph) -> set[tuple[str, str]]:
    return {tuple(sorted((str(first), str(second)))) for first, second in graph.edges}


def _write_rows(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run_sweep(
    *,
    thresholds: Iterable[int] = DEFAULT_THRESHOLDS,
    shipping_threshold: int = DEFAULT_SHIPPING_THRESHOLD,
    gene_universe_path: Path,
    mapping_path: Path,
    links_path: Path,
    output_root: Path,
) -> tuple[list[dict], list[dict]]:
    """Run the complete sweep and write graphs, manifest, and reports."""
    requested = validate_thresholds(thresholds)
    if type(shipping_threshold) is not int or shipping_threshold <= 0:
        raise ValueError("shipping_threshold must be a positive integer")
    setup_error = ""
    try:
        universe = load_gene_universe(gene_universe_path)
        mapping = load_string_mapping(mapping_path, universe)
        mapping_error = "identifier mapping produced no symbols from the gene universe" if universe and not mapping else ""
        buckets = bucket_string_edges(links_path, mapping, requested)
        setup_error = mapping_error
    except Exception as exc:
        LOGGER.exception("Threshold sweep input preparation failed")
        buckets = {}
        setup_error = str(exc)
    graph_dir = output_root / "sweep"
    graphs: dict[int, nx.Graph] = {}
    manifest: list[dict] = []

    for threshold in requested:
        csv_path = graph_dir / f"breast_cancer_graph_t{threshold}.csv"
        try:
            if setup_error and not buckets:
                raise RuntimeError(setup_error)
            graph = build_graph(buckets.get(threshold, []))
            graphs[threshold] = graph
            metrics = _graph_metrics(graph)
            status = "empty" if not graph else "degraded" if metrics["n_connected_components"] > 1 else "success"
            error = mapping_error or ("no mapped edges qualified at this threshold" if not graph else "" if status == "success" else "graph is disconnected")
            write_processed_graph(graph, csv_path)
        except Exception as exc:
            LOGGER.exception("Threshold %s failed", threshold)
            metrics = _graph_metrics(nx.Graph())
            status, error = "failed", str(exc)
            csv_path = graph_dir / f"breast_cancer_graph_t{threshold}.csv"
            write_processed_graph(nx.Graph(), csv_path)
        manifest.append({"threshold": threshold, "status": status, "error_message": error, **metrics,
                         "processed_graph_csv": str(csv_path)})

    _write_rows(graph_dir / "threshold_sweep_manifest.csv", MANIFEST_COLUMNS, manifest)
    successful = {threshold for threshold, graph in graphs.items() if graph}
    shipping_edges = _edge_set(graphs[shipping_threshold]) if shipping_threshold in successful else set()
    report: list[dict] = []
    for row in manifest:
        threshold = row["threshold"]
        graph = graphs.get(threshold, nx.Graph())
        report.append({**{key: row[key] for key in MANIFEST_COLUMNS if key not in {"processed_graph_csv"}},
                       "top10_betweenness_genes": _top_betweenness(graph),
                       "jaccard_vs_shipping_threshold": _jaccard(_edge_set(graph), shipping_edges) if threshold in successful and shipping_edges else None})
    analysis_output_root = output_root.parent / "analysis_outputs"
    _write_rows(analysis_output_root / "threshold_sensitivity.csv", REPORT_COLUMNS, report)

    comparisons: list[dict] = []
    for first, second in zip(requested, requested[1:]):
        comparisons.append({"comparison_type": "consecutive", "threshold_pair": f"{first}_vs_{second}",
                            "jaccard_overlap": _jaccard(_edge_set(graphs.get(first, nx.Graph())), _edge_set(graphs.get(second, nx.Graph())))})
    if shipping_threshold in successful:
        for threshold in requested:
            if threshold != shipping_threshold and threshold in successful:
                comparisons.append({"comparison_type": "vs_shipping", "threshold_pair": f"{threshold}_vs_{shipping_threshold}",
                                    "jaccard_overlap": _jaccard(_edge_set(graphs[threshold]), shipping_edges)})
    _write_rows(analysis_output_root / "threshold_top10_jaccard.csv",
                ["comparison_type", "threshold_pair", "jaccard_overlap"], comparisons)
    return manifest, report


def _parse_thresholds(value: str) -> list[int]:
    try:
        return [int(item.strip()) for item in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("thresholds must be comma-separated integers") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--thresholds", type=_parse_thresholds, default=DEFAULT_THRESHOLDS)
    parser.add_argument("--shipping-threshold", type=int, default=DEFAULT_SHIPPING_THRESHOLD)
    parser.add_argument("--gene-universe", type=Path, default=Path(__file__).parents[1] / "data" / "breast_cancer_graph.csv")
    parser.add_argument("--mapping", type=Path, default=Path(__file__).parents[3] / "major_project_dataset" / "9606.protein.info.v12.0.txt" / "9606.protein.info.v12.0.txt")
    parser.add_argument("--links", type=Path, default=Path(__file__).parents[3] / "major_project_dataset" / "9606.protein.links.v12.0.txt" / "9606.protein.links.v12.0.txt")
    parser.add_argument("--output-root", type=Path, default=Path(__file__).parents[1] / "data" / "processed")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_sweep(thresholds=args.thresholds, shipping_threshold=args.shipping_threshold,
              gene_universe_path=args.gene_universe, mapping_path=args.mapping,
              links_path=args.links, output_root=args.output_root)


if __name__ == "__main__":
    main()
