import csv
from pathlib import Path

import pytest

from pipeline.sweep_thresholds import run_sweep


def write_fixture_files(tmp_path: Path) -> tuple[Path, Path, Path]:
    universe = tmp_path / "universe.csv"
    universe.write_text(
        "gene1,gene2,confidence_score,interaction_type\n"
        "A,B,900,string\nB,C,800,string\nC,D,700,string\nD,E,600,string\n",
        encoding="utf-8",
    )
    mapping = tmp_path / "protein.info.tsv"
    mapping.write_text(
        "#string_protein_id\tpreferred_name\n"
        "9606.PA\tA\n9606.PB\tB\n9606.PC\tC\n9606.PD\tD\n9606.PE\tE\n",
        encoding="utf-8",
    )
    links = tmp_path / "protein.links.tsv"
    links.write_text(
        "protein1 protein2 combined_score\n"
        "9606.PA 9606.PB 900\n"
        "9606.PB 9606.PC 800\n"
        "9606.PC 9606.PD 700\n"
        "9606.PD 9606.PE 600\n"
        "9606.PA 9606.PE 100\n",
        encoding="utf-8",
    )
    return universe, mapping, links


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as input_file:
        return list(csv.DictReader(input_file))


def test_threshold_sweep_outputs_metrics_and_one_pass_results(tmp_path, monkeypatch):
    universe, mapping, links = write_fixture_files(tmp_path)
    output_root = tmp_path / "data" / "processed"

    manifest, report = run_sweep(
        thresholds=[600, 700, 800],
        shipping_threshold=700,
        gene_universe_path=universe,
        mapping_path=mapping,
        links_path=links,
        output_root=output_root,
    )

    assert [row["edge_count"] for row in manifest] == [4, 3, 2]
    assert [row["node_count"] for row in manifest] == [5, 4, 3]
    assert manifest[0]["mean_degree"] == 1.6
    assert manifest[0]["density"] == 0.4
    assert manifest[0]["n_connected_components"] == 1
    assert all((output_root / "sweep" / f"breast_cancer_graph_t{threshold}.csv").exists()
               for threshold in (600, 700, 800))
    assert len(read_csv(output_root / "sweep" / "threshold_sweep_manifest.csv")) == 3
    assert len(report) == 3
    assert report[0]["top10_betweenness_genes"] == "C,B,D,A,E"

    comparisons = read_csv(output_root.parent / "analysis_outputs" / "threshold_top10_jaccard.csv")
    by_pair = {row["threshold_pair"]: row["jaccard_overlap"] for row in comparisons}
    assert by_pair["600_vs_700"] == str(3 / 4)
    assert float(by_pair["600_vs_700"]) == report[0]["jaccard_vs_shipping_threshold"]
    assert float(by_pair["700_vs_800"]) == pytest.approx(2 / 3)


def test_empty_threshold_is_recorded(tmp_path):
    universe, mapping, links = write_fixture_files(tmp_path)
    manifest, report = run_sweep(
        thresholds=[950, 1000],
        shipping_threshold=950,
        gene_universe_path=universe,
        mapping_path=mapping,
        links_path=links,
        output_root=tmp_path / "processed",
    )

    assert len(manifest) == 2
    assert all(row["status"] == "empty" for row in manifest)
    assert all(row["node_count"] == 0 and row["edge_count"] == 0 for row in manifest)
    assert all(row["top10_betweenness_genes"] == "" for row in report)


def test_input_failure_is_recorded_for_every_threshold(tmp_path):
    missing = tmp_path / "missing.tsv"
    manifest, report = run_sweep(
        thresholds=[400, 700],
        shipping_threshold=700,
        gene_universe_path=missing,
        mapping_path=missing,
        links_path=missing,
        output_root=tmp_path / "processed",
    )

    assert len(manifest) == 2
    assert all(row["status"] == "failed" for row in manifest)
    assert all(row["error_message"] for row in manifest)
    assert all(row["status"] == "failed" for row in report)
