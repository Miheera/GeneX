"""Build and serialize graphs using the backend's processed CSV convention."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import networkx as nx

try:
    from .filter_string import StringEdge
except ImportError:
    from filter_string import StringEdge

CSV_COLUMNS = ["gene1", "gene2", "confidence_score", "interaction_type"]


def build_graph(edges: Iterable[StringEdge]) -> nx.Graph:
    """Build an undirected weighted graph from mapped STRING edges."""
    graph = nx.Graph()
    for gene1, gene2, score in edges:
        if gene1 and gene2 and gene1 != gene2:
            graph.add_edge(
                gene1,
                gene2,
                weight=score / 1000.0,
                path_weight=1000.0 / score if score > 0 else float("inf"),
            )
    return graph


def write_processed_graph(graph: nx.Graph, output_path: Path) -> None:
    """Write a graph in the normal processed graph CSV format."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for gene1, gene2, attributes in graph.edges(data=True):
            writer.writerow(
                {
                    "gene1": gene1,
                    "gene2": gene2,
                    "confidence_score": int(round(attributes["weight"] * 1000)),
                    "interaction_type": "string",
                }
            )
