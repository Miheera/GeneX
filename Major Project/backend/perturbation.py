import copy
import logging
import random
from collections.abc import Collection
from typing import Dict
import networkx as nx
from path_cost import require_positive_betweenness_costs

logger = logging.getLogger(__name__)

_baseline_efficiency: float = 0.0
_baseline_bc: Dict[str, float] = {}
_baseline_sources: list[str] = []

APPROXIMATION_SAMPLES = 100
APPROXIMATION_SEED = 42
APPROXIMATION_METHOD = f"approximate (k={APPROXIMATION_SAMPLES} samples, seed={APPROXIMATION_SEED})"

def _sample_sources(graph: nx.Graph, sample_size: int = APPROXIMATION_SAMPLES) -> list[str]:
    """Select reproducible shortest-path source nodes for approximation."""
    nodes = sorted(map(str, graph.nodes))
    if not nodes:
        return []
    return random.Random(APPROXIMATION_SEED).sample(nodes, min(sample_size, len(nodes)))


def global_efficiency(
    graph: nx.Graph,
    *,
    sources: Collection[str] | None = None,
    sample_size: int = APPROXIMATION_SAMPLES,
) -> float:
    """Estimate weighted global efficiency from reproducibly sampled source nodes.

    Sampling every node (the outcome for graphs with <= 100 nodes) produces the
    exact directed-pair global-efficiency value used by the prior implementation.
    """
    n = graph.number_of_nodes()
    if n < 2:
        return 0.0
    selected_sources = [source for source in (sources or _sample_sources(graph, sample_size)) if source in graph]
    if not selected_sources:
        return 0.0
    total = 0.0
    for source in selected_sources:
        lengths = nx.single_source_dijkstra_path_length(
            graph, source, weight='path_weight'
        )
        for target, dist in lengths.items():
            if target == source:
                continue
            if dist == 0:
                raise ValueError(
                    "Global efficiency is undefined for a reachable zero-cost pair "
                    "under the selected path-cost model"
                )
            total += 1.0 / dist
    return total / (len(selected_sources) * (n - 1))


def approximate_betweenness(graph: nx.Graph) -> Dict[str, float]:
    """Return reproducible sampled-source weighted betweenness centrality."""
    if graph.number_of_nodes() <= 1:
        return {str(node): 0.0 for node in graph}
    require_positive_betweenness_costs(graph)
    return nx.betweenness_centrality(
        graph,
        k=min(APPROXIMATION_SAMPLES, graph.number_of_nodes()),
        weight="path_weight",
        normalized=True,
        seed=APPROXIMATION_SEED,
    )

def precompute_baseline(G: nx.Graph) -> None:
    global _baseline_efficiency, _baseline_bc, _baseline_sources
    _baseline_sources = _sample_sources(G)
    _baseline_efficiency = global_efficiency(G, sources=_baseline_sources)
    _baseline_bc = approximate_betweenness(G)


def knockout_impact_score(G: nx.Graph, gene_id: str) -> float:
    """Return approximate efficiency-loss impact without unnecessary centrality work."""
    if gene_id not in G:
        raise ValueError(f"Gene not found: {gene_id}")
    if not _baseline_sources:
        precompute_baseline(G)
    perturbed = G.copy()
    perturbed.remove_node(gene_id)
    efficiency = global_efficiency(perturbed, sources=_baseline_sources)
    delta_pct = ((efficiency - _baseline_efficiency) / _baseline_efficiency * 100) if _baseline_efficiency else 0.0
    return round(abs(delta_pct), 2)

def run_perturbation(G: nx.Graph, gene_id: str, perturbation_type: str, strength: float) -> dict:
    if perturbation_type != "knockout":
        raise ValueError("run_perturbation now handles only virtual topological knockout")
    if not 0 <= strength <= 1:
        raise ValueError("strength must be in [0, 1]")
    G_prime = copy.deepcopy(G)
    if gene_id not in G_prime:
        raise ValueError(f"Gene not found: {gene_id}")
    G_prime.remove_node(gene_id)

    E_baseline = _baseline_efficiency
    E_prime = global_efficiency(G_prime, sources=_baseline_sources)
    
    efficiency_delta_pct = ((E_prime - E_baseline) / E_baseline * 100) if E_baseline > 0 else 0.0
    impact_score = round(abs(efficiency_delta_pct), 2)
    
    bc_prime = approximate_betweenness(G_prime)
    
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
    largest_component_size = max((len(component) for component in comps), default=0)
    lcc_fraction = largest_component_size / G_prime.number_of_nodes() if G_prime.number_of_nodes() else 0.0
    
    if network_fragmented and comps:
        largest = max(comps, key=len)
        isolated_genes = sorted([n for n in G_prime.nodes if n not in largest])
    else:
        isolated_genes = []

    return {
        "experiment_type": "virtual_topological_knockout",
        "efficiency_baseline": E_baseline,
        "efficiency_perturbed": E_prime,
        "efficiency_delta_pct": efficiency_delta_pct,
        "impact_score": impact_score,
        "perturbation_type": perturbation_type,
        "gene_id": gene_id,
        "top_impacted": top_impacted,
        "network_fragmented": network_fragmented,
        "connected_components": len(comps),
        "lcc_size": largest_component_size,
        "lcc_fraction": lcc_fraction,
        "isolated_genes": isolated_genes,
        "bc_baseline": _baseline_bc,
        "bc_perturbed": bc_prime,
        "approximation_method": APPROXIMATION_METHOD,
    }
