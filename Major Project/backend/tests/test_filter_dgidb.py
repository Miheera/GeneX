import csv

from pipeline.filter_dgidb import OUTPUT_COLUMNS, UNMATCHED_COLUMNS, filter_dgidb


HEADER = (
    "gene_claim_name\tgene_concept_id\tgene_name\tdrug_claim_name\tdrug_concept_id\t"
    "drug_name\tdrug_is_approved\tdrug_is_immunotherapy\tdrug_is_antineoplastic\t"
    "interaction_source_db_name\tinteraction_source_db_version\tinteraction_types\n"
)


def test_filter_dgidb_aggregates_schema_and_resolves_erbb2_alias(tmp_path):
    source = tmp_path / "interactions.tsv"
    source.write_text(
        "# Data version: test\n" + HEADER +
        "x\tx\tERBB2\tx\trxcui:1\tTrastuzumab\ttrue\tfalse\ttrue\tOncoKB\tv\tantibody\n"
        "x\tx\tERBB2\tx\trxcui:1\tTrastuzumab\ttrue\tfalse\ttrue\tTTD\tv\tantibody\n"
        "x\tx\tTP53\tx\tdrugbank:2\tNutlin\tfalse\tfalse\ttrue\tDrugBank\tv\tinhibitor\n"
        "x\tx\tNOT_IN_GRAPH\tx\tdrugbank:3\tUnknown\tfalse\tfalse\tfalse\tSource\tv\t\n",
        encoding="utf-8",
    )
    output = tmp_path / "annotations.csv"
    unmatched = tmp_path / "unmatched.csv"

    rows = filter_dgidb(source, {"HER2", "TP53", "EMPTY"}, output_path=output, unmatched_output_path=unmatched)
    by_gene = {row["gene_symbol"]: row for row in rows}

    assert by_gene["HER2"]["n_interactions"] == 1
    assert by_gene["HER2"]["drug_names"] == "Trastuzumab"
    assert by_gene["HER2"]["has_approved_drug"] is True
    assert by_gene["HER2"]["n_sources"] == 2
    assert by_gene["HER2"]["sources"] == "OncoKB|TTD"
    assert by_gene["TP53"]["n_interactions"] == 1
    assert by_gene["TP53"]["has_approved_drug"] is False
    assert by_gene["EMPTY"]["n_interactions"] == 0

    with output.open(newline="", encoding="utf-8") as output_file:
        assert csv.DictReader(output_file).fieldnames == OUTPUT_COLUMNS
    with unmatched.open(newline="", encoding="utf-8") as unmatched_file:
        reader = csv.DictReader(unmatched_file)
        assert reader.fieldnames == UNMATCHED_COLUMNS
        unmatched_rows = list(reader)
    assert unmatched_rows[0]["raw_gene_symbol"] == "NOT_IN_GRAPH"
    assert unmatched_rows[0]["reason"] == "not present in graph after synonym resolution"
