import os
import csv
import logging
from pathlib import Path
from typing import Dict
import networkx as nx

logger = logging.getLogger(__name__)

DISEASE_SCORES: Dict[str, float] = {
    "TP53": 0.98, "BRCA1": 0.95, "BRCA2": 0.91, "HER2": 0.88,
    "PIK3CA": 0.85, "ESR1": 0.82, "MYC": 0.90, "EGFR": 0.87,
}

G: nx.Graph = nx.Graph()

def load_graph() -> nx.Graph:
    csv_path = os.environ.get("GRAPH_CSV_PATH", str(Path(__file__).parent / "data" / "breast_cancer_graph.csv"))
    logger.info("Loading graph from %s", csv_path)
    
    graph = nx.Graph()
    
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            g1 = row["gene1"].strip()
            g2 = row["gene2"].strip()
            if not g1 or not g2: continue
            
            weight = int(row["confidence_score"]) / 1000.0
            path_weight = 1.0 / weight if weight > 0 else float('inf')
            
            graph.add_edge(g1, g2, weight=weight, path_weight=path_weight)
            
    for node in graph.nodes:
        graph.nodes[node]["id"] = node
        graph.nodes[node]["label"] = node
        graph.nodes[node]["degree"] = graph.degree(node)
        graph.nodes[node]["disease_score"] = DISEASE_SCORES.get(node, 0.5)

    return graph

def init_graph() -> None:
    global G
    G = load_graph()
