import os
import csv
import logging
from pathlib import Path
from typing import Dict
import networkx as nx
from path_cost import confidence_to_path_cost, get_path_cost_model

logger = logging.getLogger(__name__)

DISEASE_SCORES: Dict[str, float] = {
    "TP53": 0.98, "BRCA1": 0.95, "BRCA2": 0.91, "ERBB2": 0.88,
    "PIK3CA": 0.85, "ESR1": 0.82, "MYC": 0.90, "EGFR": 0.87,
}

G: nx.Graph = nx.Graph()

_DATA_DIR = Path(__file__).parent / "data"
DEFAULT_GRAPH_CSV_PATH = _DATA_DIR / "processed" / "expanded_graph_t900.csv"
DEFAULT_NODES_CSV_PATH = _DATA_DIR / "processed" / "expanded_nodes_t900.csv"


def load_graph() -> nx.Graph:
    graph_path_override = os.environ.get("GRAPH_CSV_PATH")
    csv_path = graph_path_override or str(DEFAULT_GRAPH_CSV_PATH)
    nodes_path_override = os.environ.get("GRAPH_NODES_CSV_PATH")
    nodes_path = (
        nodes_path_override
        if nodes_path_override is not None
        else str(DEFAULT_NODES_CSV_PATH) if graph_path_override is None else None
    )
    logger.info("Loading graph from %s", csv_path)
    
    graph = nx.Graph()
    
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            g1 = row["gene1"].strip()
            g2 = row["gene2"].strip()
            if not g1 or not g2: continue
            
            weight = int(row["confidence_score"]) / 1000.0
            path_weight = confidence_to_path_cost(weight)
            
            graph.add_edge(g1, g2, weight=weight, path_weight=path_weight)

    node_file = Path(nodes_path) if nodes_path else None
    if node_file and node_file.exists():
        with node_file.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or "gene_symbol" not in reader.fieldnames:
                raise ValueError("Graph nodes CSV must contain a gene_symbol column")
            graph.add_nodes_from(
                row["gene_symbol"].strip()
                for row in reader
                if row.get("gene_symbol", "").strip()
            )
        logger.info("Loaded explicit graph nodes from %s", node_file)
    elif node_file:
        logger.warning("Graph nodes CSV not found; deriving nodes from edges: %s", node_file)
            
    for node in graph.nodes:
        graph.nodes[node]["id"] = node
        graph.nodes[node]["label"] = node
        graph.nodes[node]["degree"] = graph.degree(node)
        graph.nodes[node]["disease_score"] = DISEASE_SCORES.get(node, 0.5)

    graph.graph["path_cost_model"] = get_path_cost_model()
    return graph

def init_graph() -> None:
    global G
    G = load_graph()
