import asyncio

import networkx as nx
import pytest
from pydantic import ValidationError

import graph_loader
import main
import perturbation
from models import PerturbRequest
from pipeline.graph_builder import build_graph
from propagation import compute_perturbation_propagation, rwr


def test_one_node_profile_is_unit_mass():
    graph = nx.Graph()
    graph.add_node("T")
    scores, iterations = rwr(graph, "T")
    assert scores == {"T": 1.0}
    assert iterations == 1


def test_two_node_profile_matches_closed_form():
    graph = nx.Graph()
    graph.add_edge("T", "N", weight=0.7)
    alpha = 0.85
    scores, _ = rwr(graph, "T", alpha=alpha)
    assert scores["T"] == pytest.approx(1 / (1 + alpha))
    assert scores["N"] == pytest.approx(alpha / (1 + alpha))
    assert sum(scores.values()) == pytest.approx(1)


def test_weighted_graph_uses_row_normalized_original_confidence():
    graph = nx.Graph()
    graph.add_edge("T", "A", weight=0.8)
    graph.add_edge("T", "B", weight=0.2)
    scores, _ = rwr(graph, "T", alpha=0.5)
    # First iterate from T sends 0.4 and 0.1 to A/B; subsequent returns
    # divide each leaf's mass equally over T's outgoing confidence weights.
    assert scores["A"] > scores["B"]
    assert scores["T"] == pytest.approx(2 / 3)
    assert scores["A"] == pytest.approx(4 / 15)
    assert scores["B"] == pytest.approx(1 / 15)
    assert sum(scores.values()) == pytest.approx(1)


def test_larger_alpha_increases_mass_away_from_restart_target():
    graph = nx.path_graph(["0", "1", "2"])
    nx.set_edge_attributes(graph, 1.0, "weight")
    low, _ = rwr(graph, "0", alpha=0.5)
    high, _ = rwr(graph, "0", alpha=0.9)
    assert 1 - low["0"] < 1 - high["0"]


def test_target_has_restart_contribution_and_dangling_mass_returns_to_target():
    graph = nx.Graph()
    graph.add_edge("T", "N", weight=1.0)
    graph.add_node("isolated")
    scores, _ = rwr(graph, "T")
    assert scores["T"] >= 1 - 0.85
    assert sum(scores.values()) == pytest.approx(1)


def test_isolated_target_remains_one_hot():
    graph = nx.Graph()
    graph.add_nodes_from(["T", "other"])
    scores, _ = rwr(graph, "T")
    assert scores == {"T": 1.0, "other": 0.0}


def test_zero_weight_edges_are_ignored_without_nan():
    graph = nx.Graph()
    graph.add_edge("T", "N", weight=0.0)
    scores, _ = rwr(graph, "T")
    assert scores["T"] == 1.0
    assert scores["N"] == 0.0
    assert all(value == value for value in scores.values())


def test_rwr_is_independent_of_selected_path_cost_model(monkeypatch):
    edges = [("T", "A", 900), ("A", "B", 500), ("T", "B", 700)]
    monkeypatch.setenv("PATH_COST_MODEL", "inverse")
    inverse_graph = build_graph(edges)
    inverse_scores, _ = rwr(inverse_graph, "T")
    monkeypatch.setenv("PATH_COST_MODEL", "log")
    log_graph = build_graph(edges)
    log_scores, _ = rwr(log_graph, "T")
    assert inverse_scores == pytest.approx(log_scores, abs=0)
    assert [d["weight"] for _, _, d in inverse_graph.edges(data=True)] == [
        d["weight"] for _, _, d in log_graph.edges(data=True)
    ]
    assert any(
        inverse_graph[u][v]["path_weight"] != log_graph[u][v]["path_weight"]
        for u, v in inverse_graph.edges
    )


def test_convergence_is_deterministic_and_failure_is_explicit():
    graph = nx.path_graph(4)
    nx.set_edge_attributes(graph, 1.0, "weight")
    first = rwr(graph, 0)
    second = rwr(graph, 0)
    assert first == second
    with pytest.raises(RuntimeError, match="did not converge"):
        rwr(graph, 0, max_iter=1)


@pytest.mark.parametrize(("kind", "sign"), [("inhibition", -1), ("overexpression", 1)])
def test_signed_profiles_scale_unit_profile_and_preserve_graph(kind, sign):
    graph = nx.Graph()
    graph.add_edge("T", "A", weight=0.8, path_weight=123)
    graph.add_edge("A", "B", weight=0.4, path_weight=456)
    edge_snapshot = [(u, v, dict(data)) for u, v, data in graph.edges(data=True)]
    unit, _ = rwr(graph, "T")
    result = compute_perturbation_propagation(graph, "T", kind, 0.4)
    assert result["propagation_scores"] == pytest.approx({n: sign * 0.4 * v for n, v in unit.items()})
    assert list(graph.edges(data=True)) == edge_snapshot
    assert result["top_impacted_genes"]
    assert result["top_impacted_genes"][0]["gene"] != "T"


def test_zero_strength_returns_zero_scores_and_no_ranked_nodes():
    graph = nx.path_graph(2)
    nx.set_edge_attributes(graph, 0.5, "weight")
    result = compute_perturbation_propagation(graph, 0, "inhibition", 0)
    assert all(value == 0 for value in result["propagation_scores"].values())
    assert result["top_impacted_genes"] == []


def test_absolute_ranking_does_not_depend_on_positive_strength():
    graph = nx.path_graph(4)
    nx.set_edge_attributes(graph, 1.0, "weight")
    a = compute_perturbation_propagation(graph, 0, "overexpression", 0.2)
    b = compute_perturbation_propagation(graph, 0, "overexpression", 0.8)
    assert [row["gene"] for row in a["top_impacted_genes"]] == [row["gene"] for row in b["top_impacted_genes"]]
    for node in graph:
        assert abs(b["propagation_scores"][str(node)]) == pytest.approx(4 * abs(a["propagation_scores"][str(node)]))


@pytest.mark.parametrize("strength", [-0.01, 1.01])
def test_invalid_strength_is_rejected(strength):
    graph = nx.path_graph(2)
    nx.set_edge_attributes(graph, 1.0, "weight")
    with pytest.raises(ValueError, match="strength must"):
        compute_perturbation_propagation(graph, 0, "inhibition", strength)


def test_unknown_propagation_target_is_rejected():
    with pytest.raises(ValueError, match="Target gene not found"):
        rwr(nx.Graph(), "missing")


def test_api_returns_distinct_propagation_and_knockout_schemas(monkeypatch):
    graph = nx.path_graph(["0", "1", "2"])
    nx.set_edge_attributes(graph, 0.8, "weight")
    nx.set_edge_attributes(graph, 1.25, "path_weight")
    monkeypatch.setattr(graph_loader, "G", graph)
    perturbation.precompute_baseline(graph)
    edge_snapshot = [(u, v, dict(data)) for u, v, data in graph.edges(data=True)]

    partial = asyncio.run(main.perturb_gene(PerturbRequest(
        gene_id="0", perturbation_type="inhibition", strength=0.5
    )))
    assert partial.propagation_model == "confidence_weighted_rwr"
    assert hasattr(partial, "propagation_scores")
    assert not hasattr(partial, "efficiency_baseline")
    assert [(u, v, dict(data)) for u, v, data in graph.edges(data=True)] == edge_snapshot

    knockout = asyncio.run(main.perturb_gene(PerturbRequest(
        gene_id="1", perturbation_type="knockout", strength=0.5
    )))
    assert knockout.experiment_type == "virtual_topological_knockout"
    assert hasattr(knockout, "efficiency_delta_pct")
    assert knockout.network_fragmented
    assert knockout.lcc_size == 1
    assert knockout.connected_components == 2


def test_api_unknown_target_and_request_strength_validation(monkeypatch):
    graph = nx.Graph()
    graph.add_node("T")
    monkeypatch.setattr(graph_loader, "G", graph)
    with pytest.raises(Exception) as error:
        asyncio.run(main.perturb_gene(PerturbRequest(
            gene_id="missing", perturbation_type="overexpression", strength=0.5
        )))
    assert getattr(error.value, "status_code", None) == 404
    with pytest.raises(ValidationError):
        PerturbRequest(gene_id="T", perturbation_type="inhibition", strength=2)
