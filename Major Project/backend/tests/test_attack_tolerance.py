import csv
from pathlib import Path

import networkx as nx

from analysis import attack_tolerance


def make_graph() -> nx.Graph:
    graph = nx.Graph()
    graph.add_edges_from(("hub", leaf) for leaf in ("a", "b", "c", "d", "e"))
    graph.add_edges_from(
        [("u", "v"), ("v", "w"), ("w", "u"), ("u", "x"), ("x", "w"), ("v", "x")]
    )
    for source, target in graph.edges:
        graph[source][target]["weight"] = 1.0
        graph[source][target]["path_weight"] = 1.0
    return graph


def test_degree_attack_recomputes_centrality_and_removes_hub_first(monkeypatch):
    graph = make_graph()
    calls = []
    original = nx.degree_centrality

    def recording_centrality(current):
        calls.append(current.number_of_nodes())
        return original(current)

    monkeypatch.setattr(attack_tolerance.nx, "degree_centrality", recording_centrality)
    results = attack_tolerance.analyze_attack_tolerance(graph, "degree")

    assert len(results) == graph.number_of_nodes()
    assert calls == list(range(graph.number_of_nodes(), 1, -1))
    assert results[1]["lcc_fraction"] == 4 / 9
    assert results[1]["efficiency_std"] is None


def test_tie_breaking_is_lexicographically_deterministic():
    graph = nx.Graph()
    graph.add_edge("a", "b", weight=1.0, path_weight=1.0)
    graph.add_edge("y", "z", weight=1.0, path_weight=1.0)

    assert attack_tolerance._highest_centrality_node(graph, "degree") == "a"


def test_random_attack_has_statistics_and_targeted_attack_is_stronger():
    graph = make_graph()
    targeted = attack_tolerance.analyze_attack_tolerance(graph, "degree")
    random_results = attack_tolerance.analyze_attack_tolerance(
        graph, "random", n_random_trials=50, seed=7
    )

    assert targeted[1]["efficiency"] < random_results[1]["efficiency"] * 0.75
    assert all(result["efficiency_std"] is not None for result in random_results)
    assert all(result["lcc_fraction_std"] is not None for result in random_results)
    assert all(result["efficiency_std"] is None for result in targeted)
    assert all(result["lcc_fraction_std"] is None for result in targeted)
    assert [result["fraction_removed"] for result in targeted] == [
        index / graph.number_of_nodes() for index in range(graph.number_of_nodes())
    ]


def test_analysis_does_not_mutate_input_and_writes_expected_csv():
    graph = make_graph()
    before = nx.node_link_data(graph)

    results = attack_tolerance.analyze_attack_tolerance(graph, "degree")

    assert nx.node_link_data(graph) == before
    output_path = Path(attack_tolerance.__file__).resolve().parents[1] / "data" / "analysis_outputs" / "attack_tolerance_degree.csv"
    assert output_path.exists()
    with output_path.open(newline="", encoding="utf-8") as output_file:
        reader = csv.DictReader(output_file)
        assert reader.fieldnames == attack_tolerance.CSV_COLUMNS
        assert len(list(reader)) == len(results)
