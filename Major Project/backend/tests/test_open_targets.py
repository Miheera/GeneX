import csv

import networkx as nx

from analysis.open_targets import METHODOLOGY_NOTE, OUTPUT_COLUMNS, open_targets


def _graph():
    graph = nx.Graph()
    for leaf in ("A", "B", "C"):
        graph.add_edge("H", leaf, weight=1.0, path_weight=1.0)
    return graph


def _annotations(path):
    path.write_text(
        "gene_symbol,n_interactions,drug_names,has_approved_drug,n_sources,sources\n"
        "A,0,,false,0,\nB,1,Drug B,true,1,Source B\nC,0,,false,0,\nH,2,Drug H,true,2,Source A|Source B\n",
        encoding="utf-8",
    )


def test_open_targets_splits_high_impact_genes_and_writes_note(tmp_path):
    graph = _graph()
    annotations = tmp_path / "annotations.csv"
    output = tmp_path / "open_targets.csv"
    _annotations(annotations)
    before = nx.node_link_data(graph)

    result = open_targets(graph, annotations, top_n=10, output_path=output)

    assert result["top_n"] == 4
    assert result["ranking_method"] == "impact_score"
    assert result["druggable"][0]["gene_symbol"] == "H"
    assert {row["gene_symbol"] for row in result["open_targets"]} == {"A", "C"}
    assert result["methodology_note"] == METHODOLOGY_NOTE
    assert nx.node_link_data(graph) == before
    with output.open(newline="", encoding="utf-8") as output_file:
        reader = csv.DictReader(output_file)
        assert reader.fieldnames == OUTPUT_COLUMNS
        rows = list(reader)
    assert len(rows) == 4
    assert all(row["methodology_note"] == METHODOLOGY_NOTE for row in rows)
