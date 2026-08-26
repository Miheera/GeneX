import copy
import math
import logging
from typing import Dict, List, Tuple
import networkx as nx

logger = logging.getLogger(__name__)

_baseline_efficiency: float = 0.0
_baseline_bc: Dict[str, float] = {}

def global_efficiency(graph: nx.Graph) -> float:
    n = len(graph.nodes)
    if n < 2: return 0.0
    total = 0.0
    for source in graph.nodes:
        lengths = nx.single_source_dijkstra_path_length(
            graph, source, weight='path_weight'
        )
        for target, dist in lengths.items():
            if target != source and dist > 0:
                total += 1.0 / dist
    return total / (n * (n - 1))

def precompute_baseline(G: nx.Graph) -> None:
    global _baseline_efficiency, _baseline_bc
    _baseline_efficiency = global_efficiency(G)
    _baseline_bc = nx.betweenness_centrality(G, weight='path_weight', normalized=True)

def run_perturbation(G: nx.Graph, gene_id: str, perturbation_type: str, strength: float) -> dict:
    G_prime = copy.deepcopy(G)
    
    if perturbation_type == "knockout":
        G_prime.remove_node(gene_id)
    elif perturbation_type == "inhibition":
        factor = 1.0 - strength
        for neighbor in list(G_prime.neighbors(gene_id)):
            edge = G_prime[gene_id][neighbor]
            new_weight = edge['weight'] * factor
            edge['weight'] = new_weight if new_weight > 0 else 1e-9
            edge['path_weight'] = 1.0 / new_weight if new_weight > 0 else float('inf')
    elif perturbation_type == "overexpression":
        factor = 1.0 + strength
        for neighbor in list(G_prime.neighbors(gene_id)):
            edge = G_prime[gene_id][neighbor]
            new_weight = min(edge['weight'] * factor, 1.0)
            edge['weight'] = new_weight
            edge['path_weight'] = 1.0 / new_weight

    E_baseline = _baseline_efficiency
    E_prime = global_efficiency(G_prime)
    
    efficiency_delta_pct = ((E_prime - E_baseline) / E_baseline * 100) if E_baseline > 0 else 0.0
    impact_score = round(abs(efficiency_delta_pct), 2)
    
    bc_prime = nx.betweenness_centrality(G_prime, weight='path_weight', normalized=True)
    
    delta_bc = {}
    for v in _baseline_bc:
        if v == gene_id: continue
        delta_bc[v] = bc_prime.get(v, 0.0) - _baseline_bc.get(v, 0.0)
        
    sorted_items = sorted(delta_bc.items(), key=lambda x: abs(x[1]), reverse=True)
    top_impacted = [
        {"gene": gene, "delta_bc": delta, "direction": "gain" if delta > 0 else "loss"}
        for gene, delta in sorted_items[:5]
    ]
    
    comps = list(nx.connected_components(G_prime))
    network_fragmented = len(comps) > 1
    
    if network_fragmented and comps:
        largest = max(comps, key=len)
        isolated_genes = sorted([n for n in G_prime.nodes if n not in largest])
    else:
        isolated_genes = []

    return {
        "efficiency_baseline": E_baseline,
        "efficiency_perturbed": E_prime,
        "efficiency_delta_pct": efficiency_delta_pct,
        "impact_score": impact_score,
        "perturbation_type": perturbation_type,
        "gene_id": gene_id,
        "top_impacted": top_impacted,
        "network_fragmented": network_fragmented,
        "isolated_genes": isolated_genes,
        "bc_baseline": _baseline_bc,
        "bc_perturbed": bc_prime
    }
