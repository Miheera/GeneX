"""Attack-tolerance analysis for the breast-cancer PPI network."""

from __future__ import annotations

import csv
import random
from pathlib import Path
from statistics import mean, pstdev
from typing import Literal

import networkx as nx

try:
    from .. import perturbation
except ImportError:
    import perturbation

Strategy = Literal["degree", "betweenness", "random"]
Result = dict[str, float | str | None]
CSV_COLUMNS = [
    "fraction_removed",
    "efficiency",
    "efficiency_std",
    "lcc_fraction",
    "lcc_fraction_std",
    "approximation_method",
]


def _metrics(graph: nx.Graph, original_nodes: int, removed: int) -> Result:
    remaining = graph.number_of_nodes()
    if remaining <= 1:
        efficiency = 0.0
    else:
        efficiency = perturbation.global_efficiency(graph)

    largest_component = max((len(component) for component in nx.connected_components(graph)), default=0)
    lcc_fraction = largest_component / remaining if remaining else 0.0
    return {
        "fraction_removed": removed / original_nodes,
        "efficiency": efficiency,
        "efficiency_std": None,
        "lcc_fraction": lcc_fraction,
        "lcc_fraction_std": None,
        "approximation_method": perturbation.APPROXIMATION_METHOD,
    }


def _highest_centrality_node(graph: nx.Graph, strategy: Literal["degree", "betweenness"]):
    centrality = nx.degree_centrality(graph) if strategy == "degree" else perturbation.approximate_betweenness(graph)
    return min(centrality, key=lambda node: (-centrality[node], str(node)), default=None)


def _targeted_results(graph: nx.Graph, strategy: Literal["degree", "betweenness"]) -> list[Result]:
    original_nodes = graph.number_of_nodes()
    current = graph.copy()
    results = [_metrics(current, original_nodes, 0)]

    for removed in range(1, original_nodes):
        node = _highest_centrality_node(current, strategy)
        if node is None:
            break
        current.remove_node(node)
        results.append(_metrics(current, original_nodes, removed))
    return results


def _random_trial(graph: nx.Graph, order: list, original_nodes: int) -> list[Result]:
    current = graph.copy()
    results = [_metrics(current, original_nodes, 0)]
    for removed, node in enumerate(order, start=1):
        if current.number_of_nodes() <= 1:
            break
        current.remove_node(node)
        results.append(_metrics(current, original_nodes, removed))
    return results


def _write_csv(strategy: Strategy, results: list[Result]) -> None:
    output_path = Path(__file__).resolve().parents[1] / "data" / "analysis_outputs" / f"attack_tolerance_{strategy}.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(results)


def analyze_attack_tolerance(
    graph: nx.Graph,
    strategy: Strategy,
    n_random_trials: int = 50,
    *,
    seed: int | None = None,
) -> list[Result]:
    """Measure weighted efficiency and LCC retention under node removal."""
    if strategy not in {"degree", "betweenness", "random"}:
        raise ValueError(f"Unsupported attack strategy: {strategy}")
    if graph.number_of_nodes() == 0:
        results: list[Result] = []
    elif strategy != "random":
        results = _targeted_results(graph, strategy)
    else:
        if n_random_trials < 1:
            raise ValueError("n_random_trials must be at least 1")
        original_nodes = graph.number_of_nodes()
        rng = random.Random(seed)
        trials = []
        for _ in range(n_random_trials):
            order = list(graph.nodes)
            rng.shuffle(order)
            trials.append(_random_trial(graph, order, original_nodes))

        results = []
        for step in range(original_nodes):
            efficiencies = [trial[step]["efficiency"] for trial in trials]
            lcc_fractions = [trial[step]["lcc_fraction"] for trial in trials]
            results.append({
                "fraction_removed": step / original_nodes,
                "efficiency": mean(efficiencies),
                "efficiency_std": pstdev(efficiencies),
                "lcc_fraction": mean(lcc_fractions),
                "lcc_fraction_std": pstdev(lcc_fractions),
                "approximation_method": perturbation.APPROXIMATION_METHOD,
            })

    _write_csv(strategy, results)
    return results
