"""Confidence-weighted Random Walk with Restart on GenX's fixed graph.

Scores describe target-centered topological proximity on the STRING association
network. They are not calibrated biological response predictions.
"""

from __future__ import annotations

import math
import os
from typing import Literal

import networkx as nx

RWR_ALPHA = float(os.environ.get("RWR_ALPHA", "0.85"))
RWR_TOL = float(os.environ.get("RWR_TOL", "1e-10"))
RWR_MAX_ITER = int(os.environ.get("RWR_MAX_ITER", "1000"))
PROPAGATION_MODEL = "confidence_weighted_rwr"
PerturbationType = Literal["inhibition", "overexpression"]


def _validate_rwr_parameters(alpha: float, tol: float, max_iter: int) -> None:
    if not math.isfinite(alpha) or not 0 <= alpha < 1:
        raise ValueError("alpha must be finite and in [0, 1)")
    if not math.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    if type(max_iter) is not int or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")


def rwr(
    graph: nx.Graph,
    target_node: str,
    *,
    alpha: float = RWR_ALPHA,
    tol: float = RWR_TOL,
    max_iter: int = RWR_MAX_ITER,
) -> tuple[dict[str, float], int]:
    """Return a unit-mass personalized PageRank profile and iteration count.

    The undirected graph is interpreted as a weighted adjacency matrix using
    only the original ``weight`` confidence attribute. Dangling mass returns
    to the one-hot target restart distribution, keeping the transition process
    stochastic. Zero-weight edges are ignored.
    """
    _validate_rwr_parameters(alpha, tol, max_iter)
    if target_node not in graph:
        raise ValueError(f"Target gene not found: {target_node}")

    nodes = list(graph.nodes)
    weighted_neighbors: dict[str, list[tuple[str, float]]] = {}
    weighted_degree: dict[str, float] = {}
    for node in nodes:
        neighbors = []
        for neighbor, attributes in graph[node].items():
            weight = float(attributes.get("weight", 1.0))
            if not math.isfinite(weight) or weight < 0:
                raise ValueError(f"Invalid STRING confidence on edge {node!r}-{neighbor!r}: {weight}")
            if weight > 0:
                neighbors.append((neighbor, weight))
        weighted_neighbors[node] = neighbors
        weighted_degree[node] = sum(weight for _, weight in neighbors)

    scores = {node: float(node == target_node) for node in nodes}
    for iteration in range(1, max_iter + 1):
        dangling_mass = sum(scores[node] for node in nodes if weighted_degree[node] == 0)
        next_scores = {node: 0.0 for node in nodes}
        for node in nodes:
            degree = weighted_degree[node]
            if degree == 0:
                continue
            mass = alpha * scores[node] / degree
            for neighbor, weight in weighted_neighbors[node]:
                next_scores[neighbor] += mass * weight

        # Restart plus dangling redistribution both return mass to the target.
        next_scores[target_node] += (1.0 - alpha) + alpha * dangling_mass
        difference = sum(abs(next_scores[node] - scores[node]) for node in nodes)
        scores = next_scores
        if difference < tol:
            total = sum(scores.values())
            if total <= 0 or not math.isfinite(total):
                raise RuntimeError("RWR produced an invalid score mass")
            if not math.isclose(total, 1.0, rel_tol=1e-12, abs_tol=1e-12):
                scores = {node: value / total for node, value in scores.items()}
            return scores, iteration

    raise RuntimeError(f"RWR did not converge within max_iter={max_iter} (tol={tol})")


def compute_perturbation_propagation(
    graph: nx.Graph,
    target_node: str,
    perturbation_type: PerturbationType,
    strength: float,
    *,
    alpha: float = RWR_ALPHA,
    tol: float = RWR_TOL,
    max_iter: int = RWR_MAX_ITER,
    top_n: int = 10,
) -> dict:
    """Scale a unit RWR profile into a signed inhibition/overexpression signal."""
    if perturbation_type not in {"inhibition", "overexpression"}:
        raise ValueError(f"Unsupported propagation perturbation type: {perturbation_type}")
    if not math.isfinite(strength) or not 0 <= strength <= 1:
        raise ValueError("strength must be finite and in [0, 1]")
    if type(top_n) is not int or top_n < 1:
        raise ValueError("top_n must be a positive integer")

    unit_scores, iterations = rwr(graph, target_node, alpha=alpha, tol=tol, max_iter=max_iter)
    sign = -1.0 if perturbation_type == "inhibition" else 1.0
    response = {node: sign * strength * score for node, score in unit_scores.items()}
    ranked_nodes = sorted(
        (node for node in graph if node != target_node and abs(response[node]) > 0),
        key=lambda node: (-abs(response[node]), str(node)),
    )
    top_impacted = [
        {
            "gene": str(node),
            "rank": rank,
            "propagation_score": response[node],
            "absolute_score": abs(response[node]),
        }
        for rank, node in enumerate(ranked_nodes[:top_n], start=1)
    ] if strength > 0 else []

    return {
        "target_gene": str(target_node),
        "perturbation_type": perturbation_type,
        "perturbation_strength": float(strength),
        "propagation_model": PROPAGATION_MODEL,
        "alpha": alpha,
        "propagation_scores": {str(node): score for node, score in response.items()},
        "top_impacted_genes": top_impacted,
        "iterations": iterations,
        "converged": True,
        "ranking_strength_invariant": True,
        "interpretation": (
            "Target-centered topological proximity on a fixed confidence-weighted "
            "STRING association network; not an experimentally calibrated biological response."
        ),
    }
