#!/usr/bin/env python3
"""Availability-only audit of the user-provided DepMap Public 26Q1 files.

This script counts identifiers and finite/missing matrix cells only. It never
exports CRISPR effect values and never reads any GenX performance statistic.
"""
from __future__ import annotations

import csv
import hashlib
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "Major Project"
MODEL = ROOT / "Model.csv"
EFFECT = ROOT / "CRISPRGeneEffect.csv"
NODES = PROJECT / "backend" / "data" / "processed" / "expanded_nodes_t900.csv"
HGNC = PROJECT / "validation" / "metadata" / "hgnc_complete_set.txt"
OUT = ROOT / "validation" / "audit"

BREAST_LINEAGE = "breast"
BREAST_PRIMARY_DISEASE = "invasive breast carcinoma"
GENE_ID_RE = re.compile(r"\((\d+)\)\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_symbols(value: str | None) -> list[str]:
    return [x.strip() for x in (value or "").split("|") if x.strip()]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path in (MODEL, EFFECT, NODES, HGNC):
        if not path.is_file():
            raise SystemExit(f"Required input missing: {path}")

    with NODES.open(encoding="utf-8-sig", newline="") as handle:
        graph_nodes = [row["gene_symbol"].strip() for row in csv.DictReader(handle)]
    if len(graph_nodes) != 881 or len(set(graph_nodes)) != 881:
        raise SystemExit(f"Frozen node list expected 881 unique symbols; observed {len(graph_nodes)} rows / {len(set(graph_nodes))} unique")

    # HGNC exact symbol first, then exact unique alias/previous symbol.
    approved: dict[str, list[dict]] = defaultdict(list)
    aliases: dict[str, list[tuple[dict, str]]] = defaultdict(list)
    with HGNC.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("symbol"):
                approved[row["symbol"].strip().upper()].append(row)
            for field in ("alias_symbol", "prev_symbol"):
                for alias in split_symbols(row.get(field)):
                    aliases[alias.upper()].append((row, field))

    with EFFECT.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    gene_columns: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for column_index, name in enumerate(header[1:], start=1):
        match = GENE_ID_RE.search(name)
        if match:
            gene_columns[match.group(1)].append((column_index, name, name[:match.start()].strip()))
    column_duplicate_ids = {gene_id for gene_id, cols in gene_columns.items() if len(cols) > 1}

    mapping: list[dict] = []
    mapped_node_to_col: dict[str, int] = {}
    mapped_node_to_entrez: dict[str, str] = {}
    for symbol in graph_nodes:
        key = symbol.upper()
        method = ""
        hgnc_record = None
        approved_matches = approved.get(key, [])
        if len(approved_matches) == 1:
            hgnc_record, method = approved_matches[0], "exact_HGNC_approved_symbol"
        elif len(approved_matches) > 1:
            method = "ambiguous_HGNC_approved_symbol"
        else:
            alias_matches = {}
            for record, alias_kind in aliases.get(key, []):
                alias_matches[record.get("hgnc_id", "")] = (record, alias_kind)
            if len(alias_matches) == 1:
                _, (hgnc_record, alias_kind) = next(iter(alias_matches.items()))
                method = alias_kind
                method = f"exact_HGNC_{method}"
            elif len(alias_matches) > 1:
                method = "ambiguous_HGNC_alias_or_previous_symbol"
            else:
                method = "no_exact_HGNC_match"

        entrez = (hgnc_record or {}).get("entrez_id", "").strip()
        columns = gene_columns.get(entrez, []) if entrez else []
        if method.startswith("ambiguous_"):
            status, reason = "ambiguous", "GenX symbol resolves to multiple HGNC records"
        elif not hgnc_record:
            status, reason = "unmapped", "No exact approved symbol or unique exact alias/previous-symbol match in HGNC"
        elif not entrez:
            status, reason = "unmapped", "Unique HGNC record has no Entrez identifier"
        elif len(columns) == 0:
            status, reason = "absent_from_matrix", "HGNC Entrez identifier is not a CRISPRGeneEffect column"
        elif len(columns) > 1:
            status, reason = "ambiguous_matrix_column", "Entrez identifier occurs in multiple matrix columns"
        else:
            status, reason = "mapped", "Unique HGNC Entrez identifier matched to one matrix column"
            mapped_node_to_col[symbol] = columns[0][0]
            mapped_node_to_entrez[symbol] = entrez
        mapping.append({
            "genx_gene": symbol,
            "hgnc_symbol": (hgnc_record or {}).get("symbol", ""),
            "hgnc_id": (hgnc_record or {}).get("hgnc_id", ""),
            "entrez_id": entrez,
            "depmap_identifier": columns[0][1] if len(columns) == 1 else "",
            "mapping_status": status,
            "mapping_reason": reason,
            "mapping_rule": method,
            "duplicate_flag": False,
        })

    entrez_to_nodes: dict[str, list[str]] = defaultdict(list)
    for node, entrez in mapped_node_to_entrez.items():
        entrez_to_nodes[entrez].append(node)
    for record in mapping:
        record["duplicate_flag"] = len(entrez_to_nodes.get(record["entrez_id"], [])) > 1 if record["entrez_id"] else False
        if record["duplicate_flag"] and record["mapping_status"] == "mapped":
            record["mapping_status"] = "duplicate_genx_mapping"
            record["mapping_reason"] = "Multiple GenX nodes resolve to the same DepMap Entrez identifier; retained and flagged"
            mapped_node_to_col.pop(record["genx_gene"], None)
            mapped_node_to_entrez.pop(record["genx_gene"], None)
    write_csv(OUT / "depmap_genx_gene_mapping.csv", list(mapping[0]), mapping)

    with MODEL.open(encoding="utf-8-sig", newline="") as handle:
        model_rows = list(csv.DictReader(handle))
    model_counts = Counter((row.get("ModelID") or "").strip() for row in model_rows)
    matrix_ids = set()
    duplicate_matrix_ids = set()
    cohort_records: list[dict] = []
    cohort_ids: list[str] = []
    included_disease_counts = Counter()
    lineage_counts = Counter()
    for row in model_rows:
        model_id = (row.get("ModelID") or "").strip()
        lineage = (row.get("OncotreeLineage") or "").strip()
        disease = (row.get("OncotreePrimaryDisease") or "").strip()
        is_breast_lineage = lineage.casefold() == BREAST_LINEAGE
        include = is_breast_lineage and disease.casefold() == BREAST_PRIMARY_DISEASE
        if include:
            cohort_ids.append(model_id)
            included_disease_counts[disease] += 1
        lineage_counts[lineage or "<blank>"] += 1
        if not model_id:
            reason = "Missing ModelID"
        elif model_counts[model_id] > 1:
            reason = "Duplicate ModelID in Model.csv; cohort eligibility retained by annotation but identifier ambiguous"
        elif include:
            reason = "Included: OncotreeLineage=Breast and OncotreePrimaryDisease=Invasive Breast Carcinoma"
        elif is_breast_lineage:
            reason = "Excluded: Breast lineage but primary disease is not Invasive Breast Carcinoma (e.g. in situ, NOS, or non-cancerous)"
        elif disease.casefold() == BREAST_PRIMARY_DISEASE:
            reason = "Excluded: primary disease matches breast carcinoma but lineage is not Breast"
        else:
            reason = "Excluded: does not meet the frozen breast-cancer metadata rule"
        cohort_records.append({
            "model_id": model_id,
            "cell_line_name": row.get("CellLineName", ""),
            "depmap_model_type": row.get("DepmapModelType", ""),
            "oncotree_lineage": lineage,
            "oncotree_primary_disease": disease,
            "oncotree_subtype": row.get("OncotreeSubtype", ""),
            "oncotree_code": row.get("OncotreeCode", ""),
            "tissue_origin": row.get("TissueOrigin", ""),
            "valid_model_id": bool(model_id),
            "model_id_duplicate": model_counts[model_id] > 1 if model_id else False,
            "breast_cancer_inclusion": include,
            "exclusion_reason": "" if include else reason,
        })
    write_csv(OUT / "depmap_breast_model_metadata.csv", list(cohort_records[0]), cohort_records)

    # Numeric completeness: inspect the mapped GenX columns only; never emit scores.
    cohort_id_set = set(cohort_ids)
    expected_mapped = len(mapped_node_to_col)
    model_stats: dict[str, dict] = {
        mid: {"seen_in_matrix": False, "duplicate_matrix_row": False, "numeric": 0, "missing": 0, "invalid": 0}
        for mid in cohort_ids
    }
    breast_genes_with_numeric: set[str] = set()
    matrix_row_count = 0
    with EFFECT.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            matrix_row_count += 1
            if not row:
                continue
            model_id = row[0].strip()
            if model_id not in cohort_id_set:
                continue
            stat = model_stats[model_id]
            if stat["seen_in_matrix"]:
                stat["duplicate_matrix_row"] = True
            stat["seen_in_matrix"] = True
            for symbol, column_index in mapped_node_to_col.items():
                value = row[column_index].strip() if column_index < len(row) else ""
                if not value or value.casefold() in {"na", "nan", "null", "none"}:
                    stat["missing"] += 1
                    continue
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    stat["invalid"] += 1
                    continue
                if not math.isfinite(number):
                    stat["invalid"] += 1
                    continue
                stat["numeric"] += 1
                breast_genes_with_numeric.add(symbol)

    details = []
    for model_id in cohort_ids:
        stat = model_stats[model_id]
        absent_mapped = expected_mapped - stat["numeric"] - stat["missing"] - stat["invalid"] if stat["seen_in_matrix"] else 0
        details.append({
            "model_id": model_id,
            "model_present_in_effect_matrix": stat["seen_in_matrix"],
            "duplicate_matrix_row": stat["duplicate_matrix_row"],
            "mapped_genx_gene_columns": expected_mapped,
            "numeric_dependency_cells": stat["numeric"],
            "missing_cells": stat["missing"] + (expected_mapped if not stat["seen_in_matrix"] else 0),
            "invalid_nonfinite_or_nonnumeric_cells": stat["invalid"],
            "unobserved_due_to_short_row_cells": max(absent_mapped, 0),
            "eligible_model_with_any_numeric_genx_dependency": stat["numeric"] > 0,
        })
    write_csv(OUT / "depmap_numeric_completeness.csv", list(details[0]), details)

    breast_models = len(cohort_ids)
    valid_model_ids = sum(bool((row.get("ModelID") or "").strip()) for row in model_rows)
    unique_matrix_models = matrix_row_count - len(duplicate_matrix_ids)
    breast_matrix_rows = sum(model_stats[x]["seen_in_matrix"] for x in cohort_ids)
    numeric_cells = sum(model_stats[x]["numeric"] for x in cohort_ids)
    missing_cells = sum(model_stats[x]["missing"] + (expected_mapped if not model_stats[x]["seen_in_matrix"] else 0) for x in cohort_ids)
    invalid_cells = sum(model_stats[x]["invalid"] for x in cohort_ids)
    duplicate_rows = sum(model_stats[x]["duplicate_matrix_row"] for x in cohort_ids)
    # Matrix IDs are unique-checked in a dedicated pass without inspecting values.
    with EFFECT.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        seen = set()
        for row in reader:
            if row:
                if row[0] in seen:
                    duplicate_matrix_ids.add(row[0])
                seen.add(row[0])
    unique_matrix_models = matrix_row_count - len(duplicate_matrix_ids)
    n_mapped_unique = sum(record["mapping_status"] == "mapped" for record in mapping)
    n_duplicate_genx = sum(record["mapping_status"] == "duplicate_genx_mapping" for record in mapping)
    n_any_numeric = len(breast_genes_with_numeric)
    all_genx_cells = len(graph_nodes) * breast_models
    mapped_cells_expected = n_mapped_unique * breast_models
    metrics = [
        ("GenX nodes", len(graph_nodes)),
        ("GenX genes with unique HGNC/Entrez mapping", sum(record["mapping_status"] in {"mapped", "absent_from_matrix"} for record in mapping)),
        ("GenX genes present as unique columns in DepMap matrix", n_mapped_unique),
        ("GenX genes absent from DepMap matrix", sum(record["mapping_status"] == "absent_from_matrix" for record in mapping)),
        ("GenX genes without unambiguous HGNC/Entrez mapping", sum(record["mapping_status"] == "unmapped" or record["mapping_status"].startswith("ambiguous") for record in mapping)),
        ("GenX genes with numeric dependency data in breast cohort", n_any_numeric),
        ("Total DepMap models in Model.csv", len(model_rows)),
        ("Model.csv rows with valid model IDs", valid_model_ids),
        ("Breast-cancer models (frozen metadata rule)", breast_models),
        ("Breast models represented in CRISPRGeneEffect.csv", breast_matrix_rows),
        ("Breast models with at least one numeric mapped GenX dependency value", sum(model_stats[x]["numeric"] > 0 for x in cohort_ids)),
        ("GenX x breast matrix cells (881 x cohort size)", all_genx_cells),
        ("Mapped GenX x breast cells expected", mapped_cells_expected),
        ("Numeric mapped GenX x breast cells", numeric_cells),
        ("Missing mapped GenX x breast cells (including models without a matrix row)", missing_cells),
        ("Missing cells across full 881 x breast cohort grid", all_genx_cells - numeric_cells),
        ("Mapped columns expected in breast models present in matrix", expected_mapped * breast_matrix_rows),
        ("Missing mapped cells among breast models present in matrix", expected_mapped * breast_matrix_rows - numeric_cells),
        ("Invalid/nonfinite/nonnumeric mapped GenX x breast cells", invalid_cells),
        ("Matrix rows (models)", matrix_row_count),
        ("Unique matrix model IDs", unique_matrix_models),
        ("Duplicate matrix model IDs", len(duplicate_matrix_ids)),
        ("Duplicate matrix gene Entrez identifiers", len(column_duplicate_ids)),
        ("Duplicate GenX to DepMap gene mappings", n_duplicate_genx),
        ("Potential DepMap validation targets (mapped and >=1 numeric breast model)", n_any_numeric),
        ("Breast-lineage models excluded by exact primary-disease rule", sum(1 for row in cohort_records if row["oncotree_lineage"].casefold() == "breast" and not row["breast_cancer_inclusion"])),
    ]
    write_csv(OUT / "depmap_numeric_completeness_summary.csv", ["quantity", "result"], [{"quantity": k, "result": v} for k, v in metrics])

    # File inventory: local release files have no internal release label; provenance is
    # the user-provided official 26Q1 download plus format/count cross-checks.
    inventory = []
    for path, role in ((MODEL, "model metadata"), (EFFECT, "Chronos CRISPR gene effect matrix"), (HGNC, "identifier crosswalk"), (NODES, "frozen GenX node list")):
        inventory.append({"filename": path.name, "path": str(path), "role": role,
                          "bytes": path.stat().st_size, "sha256": sha256(path),
                          "release": "DepMap Public 26Q1 (user-supplied download provenance)" if path in (MODEL, EFFECT) else "not applicable"})
    write_csv(OUT / "depmap_input_manifest.csv", list(inventory[0]), inventory)

    print("DepMap audit complete; no score values exported and no performance statistics calculated.")
    print(f"Model rows={len(model_rows)}; CRISPR rows={matrix_row_count}; unique CRISPR IDs={unique_matrix_models}")
    print(f"Breast cohort={breast_models}; mapped GenX genes={n_mapped_unique}; numeric genes={n_any_numeric}; numeric cells={numeric_cells}/{mapped_cells_expected}")
    print(f"Potential targets={n_any_numeric}; duplicate matrix model IDs={len(duplicate_matrix_ids)}; duplicate Entrez columns={len(column_duplicate_ids)}")


if __name__ == "__main__":
    main()
