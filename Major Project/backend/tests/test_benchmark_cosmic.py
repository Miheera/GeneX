import csv
from pathlib import Path

import networkx as nx
import pytest

from analysis import benchmark_cosmic
from pipeline.filter_cosmic import filter_cosmic


def make_graph() -> nx.Graph:
    graph = nx.Graph()
    for leaf in ("A", "B", "C", "D"):
        graph.add_edge("H", leaf, weight=1.0, path_weight=1.0)
    return graph


def write_cosmic(path: Path) -> None:
    path.write_text(
        "gene_symbol,cosmic_tier,breast_annotated,matched_term,is_tier1\n"
        "H,1,True,breast,True\n"
        "B,2,True,brca,False\n"
        "X,1,True,breast,True\n"
        "A,2,False,,False\n",
        encoding="utf-8",
    )


def test_filter_cosmic_preserves_tier1_and_deduplicates(tmp_path):
    source = tmp_path / "cosmic.csv"
    source.write_text(
        "Gene Symbol,Tier,Tumour Types(Somatic),Tumour Types(Germline),Other\n"
        "H,2,Breast carcinoma,,ignored\n"
        "H,1,,BRCA,ignored\n"
        "X,1,,Lung,breast in unrelated field\n"
        "Y,2,,,ignored\n",
        encoding="utf-8",
    )
    output = tmp_path / "processed.csv"

    rows = filter_cosmic(source, output)

    by_gene = {row["gene_symbol"]: row for row in rows}
    assert by_gene["H"]["cosmic_tier"] == 1
    assert by_gene["H"]["matched_term"] == "brca|breast"
    assert by_gene["X"]["breast_annotated"] is False
    assert by_gene["Y"]["is_tier1"] is False
    assert output.exists()


def test_benchmark_rankings_metrics_caveat_and_csv(tmp_path):
    graph = make_graph()
    cosmic_path = tmp_path / "cosmic_breast_genes.csv"
    output_path = tmp_path / "analysis_outputs" / "cosmic_benchmark.csv"
    write_cosmic(cosmic_path)
    before = nx.node_link_data(graph)

    rows = benchmark_cosmic.benchmark_cosmic(
        graph, cosmic_path, [1, 3, 10], output_path=output_path
    )

    assert len(rows) == 9
    degree_top1 = next(row for row in rows if row["rank_by"] == "degree" and row["top_n"] == 1)
    assert degree_top1["precision"] == 1.0
    assert degree_top1["recall"] == 0.5
    assert degree_top1["f1"] == pytest.approx(2 / 3)
    degree_top3 = next(row for row in rows if row["rank_by"] == "degree" and row["top_n"] == 3)
    assert degree_top3["precision"] == pytest.approx(2 / 3)
    assert degree_top3["recall"] == 1.0
    assert degree_top3["f1"] == pytest.approx(0.8)
    assert degree_top3["n_positives_total"] == 3
    assert degree_top3["n_positives_in_graph"] == 2
    assert degree_top3["n_positives_excluded_from_graph"] == 1
    assert degree_top3["methodology_note"].startswith("Unknown:")

    for method in ("betweenness", "degree", "impact_score"):
        top1 = next(row for row in rows if row["rank_by"] == method and row["top_n"] == 1)
        assert top1["precision"] == 1.0
    ranking = benchmark_cosmic._rankings(graph)
    assert min(ranking["betweenness"], key=lambda gene: (-ranking["betweenness"][gene], str(gene))) == "H"
    assert min(ranking["impact_score"], key=lambda gene: (-ranking["impact_score"][gene], str(gene))) == "H"

    top10 = next(row for row in rows if row["rank_by"] == "degree" and row["top_n"] == 10)
    assert top10["precision"] == 0.4
    assert nx.node_link_data(graph) == before

    with output_path.open(newline="", encoding="utf-8") as output_file:
        reader = csv.DictReader(output_file)
        assert reader.fieldnames == benchmark_cosmic.OUTPUT_COLUMNS
        assert len(list(reader)) == 9


def test_missing_or_malformed_cosmic_fails_clearly(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        benchmark_cosmic.benchmark_cosmic(make_graph(), tmp_path / "missing.csv", [1])

    malformed = tmp_path / "malformed.csv"
    malformed.write_text("gene_symbol\nH\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must contain"):
        benchmark_cosmic.benchmark_cosmic(make_graph(), malformed, [1])
