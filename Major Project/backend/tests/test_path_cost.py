import math

import networkx as nx
import pytest

import perturbation
from path_cost import confidence_to_path_cost
from pipeline.graph_builder import build_graph


@pytest.mark.parametrize("w, expected", [(1.0, 1.0), (0.9, 1 / 0.9), (0.5, 2.0), (0.1, 10.0), (0.001, 1000.0)])
def test_inverse_cost_values(w, expected):
    assert confidence_to_path_cost(w, "inverse") == pytest.approx(expected)


@pytest.mark.parametrize("w, expected", [(1.0, 0.0), (0.9, -math.log(0.9)), (0.5, math.log(2)), (0.1, -math.log(0.1)), (0.001, -math.log(0.001))])
def test_log_cost_values(w, expected):
    assert confidence_to_path_cost(w, "log") == pytest.approx(expected)


@pytest.mark.parametrize("model", ["inverse", "log"])
def test_zero_confidence_maps_to_infinite_cost(model):
    assert math.isinf(confidence_to_path_cost(0, model))


@pytest.mark.parametrize("w", [-0.1, float("nan"), float("inf")])
def test_invalid_confidence_rejected(w):
    with pytest.raises(ValueError):
        confidence_to_path_cost(w, "log")


def test_log_rejects_confidence_above_one():
    with pytest.raises(ValueError, match="requires normalized confidence"):
        confidence_to_path_cost(1.001, "log")


def test_model_selection_changes_builder_cost_without_changing_confidence(monkeypatch):
    edge = [("A", "B", 500)]
    monkeypatch.setenv("PATH_COST_MODEL", "inverse")
    inverse = build_graph(edge)
    monkeypatch.setenv("PATH_COST_MODEL", "log")
    logarithmic = build_graph(edge)
    assert inverse["A"]["B"]["weight"] == logarithmic["A"]["B"]["weight"] == 0.5
    assert inverse["A"]["B"]["path_weight"] == 2
    assert logarithmic["A"]["B"]["path_weight"] == pytest.approx(math.log(2))


@pytest.mark.parametrize("model, expected", [("inverse", 2.0), ("log", math.log(2))])
def test_dijkstra_and_disconnected_efficiency(monkeypatch, model, expected):
    monkeypatch.setenv("PATH_COST_MODEL", model)
    graph = build_graph([("A", "B", 500)])
    graph.add_node("C")
    lengths = nx.single_source_dijkstra_path_length(graph, "A", weight="path_weight")
    assert lengths["B"] == pytest.approx(expected)
    assert "C" not in lengths
    assert perturbation.global_efficiency(graph, sources=["A"]) == pytest.approx(1 / (2 * expected))


def test_log_zero_cost_efficiency_and_betweenness_are_explicit_errors(monkeypatch):
    monkeypatch.setenv("PATH_COST_MODEL", "log")
    graph = build_graph([("A", "B", 1000)])
    assert graph["A"]["B"]["path_weight"] == 0
    with pytest.raises(ValueError, match="zero-cost pair"):
        perturbation.global_efficiency(graph, sources=["A"])
    with pytest.raises(ValueError, match="zero path costs"):
        perturbation.approximate_betweenness(graph)


@pytest.mark.parametrize("kind", ["inhibition", "overexpression"])
def test_structural_perturbation_api_no_longer_accepts_partial_models(kind):
    graph = build_graph([("A", "B", 500)])
    with pytest.raises(ValueError, match="only virtual topological knockout"):
        perturbation.run_perturbation(graph, "A", kind, 0.5)


def test_knockout_changes_topology_only_by_deleting_node(monkeypatch):
    monkeypatch.setenv("PATH_COST_MODEL", "log")
    graph = build_graph([("A", "B", 500), ("B", "C", 600)])
    perturbation.precompute_baseline(graph)
    result = perturbation.run_perturbation(graph, "B", "knockout", 0)
    assert "B" not in result["bc_perturbed"]
    assert result["network_fragmented"] is True
    assert result["efficiency_perturbed"] == 0
