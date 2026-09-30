#!/usr/bin/env python3
"""Frozen GenX virtual-knockout vs DepMap dependency validation.

The `freeze` phase creates the cohort, target eligibility file, degree strata,
and analysis manifest using metadata/data-availability only. The `analyze`
phase verifies those frozen inputs, calculates the unchanged GenX knockout
efficiency metric and prespecified statistics, and writes validation outputs.
No GenX source, graph, or raw DepMap file is modified.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "Major Project"
BACKEND = PROJECT / "backend"
MODEL = ROOT / "Model.csv"
EFFECT = ROOT / "CRISPRGeneEffect.csv"
NODES = BACKEND / "data" / "processed" / "expanded_nodes_t900.csv"
GRAPH = BACKEND / "data" / "processed" / "expanded_graph_t900.csv"
GENX_MAPPING = ROOT / "validation" / "audit" / "depmap_genx_gene_mapping.csv"
COHORT_PATH = ROOT / "validation" / "frozen_depmap_validation_cohort.csv"
TARGETS_PATH = ROOT / "validation" / "frozen_depmap_target_population.csv"
MANIFEST_PATH = ROOT / "validation" / "depmap_validation_analysis_manifest.json"
RESULTS_DIR = ROOT / "validation" / "results"
FIGURES_DIR = ROOT / "validation" / "figures"
TARGET_TABLE = ROOT / "validation" / "depmap_validation_target_table.csv"

MIN_MODELS = 10
SENSITIVITY_MIN_MODELS = 20
BOOTSTRAP_B = 10_000
PERMUTATION_B = 10_000
DEGREE_PERMUTATION_B = 10_000
BOOTSTRAP_SEED = 20260929
PERMUTATION_SEED = 20260930
DEGREE_PERMUTATION_SEED = 20260931
DEGREE_BINS = 10
BREAST_LINEAGE = "breast"
BREAST_DISEASE = "invasive breast carcinoma"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value: str) -> float | None:
    text = (value or "").strip()
    if not text or text.casefold() in {"na", "nan", "null", "none"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def degree_bins(items: list[tuple[str, int]], requested: int = DEGREE_BINS) -> tuple[dict[str, int], dict]:
    """Partition whole equal-degree groups into nearly equal contiguous bins.

    Dynamic programming minimizes the sum of squared deviations from equal bin
    size, with degrees sorted ascending and ties kept together. The lower split
    index wins exact objective ties. If fewer than 10 distinct degrees exist,
    use that many nonempty bins.
    """
    by_degree: dict[int, list[str]] = defaultdict(list)
    for gene, degree in items:
        by_degree[int(degree)].append(gene)
    degrees = sorted(by_degree)
    k = min(requested, len(degrees))
    if not k:
        return {}, {"requested_bins": requested, "actual_bins": 0, "degree_groups": []}
    groups = [sorted(by_degree[d]) for d in degrees]
    counts = [len(group) for group in groups]
    prefix = [0]
    for count in counts:
        prefix.append(prefix[-1] + count)
    total = prefix[-1]
    target_size = total / k
    # dp[b][j] = minimum cost to partition first j degree groups into b bins.
    dp = [[math.inf] * (len(groups) + 1) for _ in range(k + 1)]
    prev = [[-1] * (len(groups) + 1) for _ in range(k + 1)]
    dp[0][0] = 0.0
    for b in range(1, k + 1):
        for j in range(b, len(groups) + 1):
            for i in range(b - 1, j):
                size = prefix[j] - prefix[i]
                candidate = dp[b - 1][i] + (size - target_size) ** 2
                if candidate < dp[b][j] - 1e-12:
                    dp[b][j] = candidate
                    prev[b][j] = i
    intervals = []
    j = len(groups)
    for b in range(k, 0, -1):
        i = prev[b][j]
        if i < 0:
            raise RuntimeError("Could not construct deterministic degree strata")
        intervals.append((i, j))
        j = i
    intervals.reverse()
    result = {}
    details = []
    for label, (i, j) in enumerate(intervals, start=1):
        members = [gene for group in groups[i:j] for gene in group]
        for gene in members:
            result[gene] = label
        details.append({
            "bin": label,
            "min_degree": degrees[i],
            "max_degree": degrees[j - 1],
            "target_count": len(members),
            "genes": members,
        })
    return result, {"requested_bins": requested, "actual_bins": k, "tie_rule": "whole equal-degree groups; dynamic programming minimizes sum of squared bin-size deviations; lower split index resolves exact ties", "degree_groups": details}


def verify_frozen_inputs() -> dict:
    required = (MODEL, EFFECT, NODES, GRAPH, GENX_MAPPING)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required input(s) missing: " + ", ".join(missing))
    hashes = {str(path.relative_to(ROOT)): sha256(path) for path in required}
    expected = "D9E4DC10A4A91154662D1A3C4EF047ACB7145658EBE76D6344B07C15BC6BBF91"
    if hashes[str(GRAPH.relative_to(ROOT))].upper() != expected:
        raise RuntimeError("Frozen graph SHA-256 mismatch; stopping before validation")
    return hashes


def freeze() -> None:
    hashes = verify_frozen_inputs()
    with MODEL.open(encoding="utf-8-sig", newline="") as stream:
        model_rows = list(csv.DictReader(stream))
    with EFFECT.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        matrix_header = next(reader)
        matrix_ids = set()
        for row in reader:
            if row:
                if row[0] in matrix_ids:
                    raise RuntimeError(f"Duplicate CRISPRGeneEffect model ID: {row[0]}")
                matrix_ids.add(row[0])

    cohort = []
    validation_ids = set()
    metadata_eligible = 0
    for row in model_rows:
        model_id = (row.get("ModelID") or "").strip()
        lineage = (row.get("OncotreeLineage") or "").strip()
        disease = (row.get("OncotreePrimaryDisease") or "").strip()
        metadata_ok = lineage.casefold() == BREAST_LINEAGE and disease.casefold() == BREAST_DISEASE
        has_matrix = model_id in matrix_ids
        included = bool(metadata_ok and has_matrix)
        if metadata_ok:
            metadata_eligible += 1
        if included:
            validation_ids.add(model_id)
        if not model_id:
            reason = "Excluded: missing ModelID"
        elif metadata_ok and has_matrix:
            reason = "Included: frozen Breast/Invasive Breast Carcinoma metadata rule and CRISPRGeneEffect model row present"
        elif metadata_ok:
            reason = "Excluded from validation measurements: meets frozen 83-model biological cohort rule, but no CRISPRGeneEffect row (33-model availability gap)"
        elif lineage.casefold() == BREAST_LINEAGE:
            reason = "Excluded: Breast lineage, but primary disease is not Invasive Breast Carcinoma"
        elif disease.casefold() == BREAST_DISEASE:
            reason = "Excluded: primary disease matches, but lineage is not Breast"
        else:
            reason = "Excluded: does not meet frozen breast-cancer metadata rule"
        cohort.append({
            "model_id": model_id,
            "cell_line_name": row.get("CellLineName", ""),
            "disease": row.get("OncotreePrimaryDisease", ""),
            "tissue": lineage,
            "tissue_origin": row.get("TissueOrigin", ""),
            "lineage": lineage,
            "primary_disease": disease,
            "metadata_cohort_eligible": metadata_ok,
            "has_crispr_gene_effect_row": has_matrix,
            "validation_inclusion_status": "included" if included else "excluded",
            "inclusion_exclusion_reason": reason,
        })
    if metadata_eligible != 83 or len(validation_ids) != 50:
        raise RuntimeError(f"Frozen cohort counts changed: metadata={metadata_eligible}, represented={len(validation_ids)}")
    write_csv(COHORT_PATH, list(cohort[0]), cohort)

    mapping = [r for r in csv_rows(GENX_MAPPING) if r["mapping_status"] == "mapped" and r["duplicate_flag"].casefold() == "false"]
    if len(mapping) != 858:
        raise RuntimeError(f"Expected 858 uniquely mapped matrix genes, observed {len(mapping)}")
    depmap_id_to_col = {}
    for i, name in enumerate(matrix_header[1:], start=1):
        if name:
            depmap_id_to_col[name] = i
    missing_cols = [r["depmap_identifier"] for r in mapping if r["depmap_identifier"] not in depmap_id_to_col]
    if missing_cols:
        raise RuntimeError(f"Frozen gene mapping no longer matches matrix columns: {missing_cols[:3]}")

    finite_counts = {r["genx_gene"]: 0 for r in mapping}
    cohort_matrix_rows = 0
    with EFFECT.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        next(reader)
        for row in reader:
            if not row or row[0] not in validation_ids:
                continue
            cohort_matrix_rows += 1
            for record in mapping:
                value = finite_float(row[depmap_id_to_col[record["depmap_identifier"]]])
                if value is not None:
                    finite_counts[record["genx_gene"]] += 1
    if cohort_matrix_rows != 50:
        raise RuntimeError(f"Expected 50 represented breast models, observed {cohort_matrix_rows}")

    # Baseline degree is graph-only. Freeze degree bins before computing any
    # target-level dependency aggregation or GenX-vs-DepMap statistic.
    sys.path.insert(0, str(BACKEND))
    os.environ["PATH_COST_MODEL"] = "inverse"
    import graph_loader
    graph = graph_loader.load_graph()
    if graph.number_of_nodes() != 881 or graph.number_of_edges() != 9681:
        raise RuntimeError(f"Frozen graph dimensions changed: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    eligible_genes = [r["genx_gene"] for r in mapping if finite_counts[r["genx_gene"]] >= MIN_MODELS]
    degree_data = [(gene, int(graph.degree(gene))) for gene in eligible_genes]
    bins, bin_manifest = degree_bins(degree_data)

    records = []
    for r in sorted(mapping, key=lambda x: x["genx_gene"]):
        gene = r["genx_gene"]
        count = finite_counts[gene]
        eligible = count >= MIN_MODELS
        records.append({
            "genx_gene": gene,
            "depmap_gene_identifier": r["depmap_identifier"],
            "finite_breast_model_dependency_count": count,
            "baseline_graph_degree": int(graph.degree(gene)),
            "degree_stratum": bins.get(gene, ""),
            "primary_analysis_eligible": eligible,
            "exclusion_reason": "" if eligible else f"Only {count} finite values; frozen minimum is {MIN_MODELS}",
        })
    eligible_n = sum(finite_counts[g] >= MIN_MODELS for g in finite_counts)
    write_csv(TARGETS_PATH, list(records[0]), records)
    if eligible_n < 3:
        raise RuntimeError(f"Only {eligible_n} targets pass the frozen threshold; stop before correlation analysis")

    package_versions = {
        "python": platform.python_version(),
        "networkx": nx.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }
    manifest = {
        "status": "FROZEN_BEFORE_ASSOCIATION_ANALYSIS",
        "freeze_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "execution_date_local": datetime.now().astimezone().date().isoformat(),
        "analysis_script": str(Path(__file__).resolve()),
        "inputs_sha256": hashes,
        "graph": {"path": str(GRAPH), "sha256": hashes[str(GRAPH.relative_to(ROOT))], "nodes": 881, "edges": 9681, "path_cost_model": "inverse (existing default; explicitly fixed for run)"},
        "depmap": {"release": "Public 26Q1, user-provided official downloads", "model_file": str(MODEL), "gene_effect_file": str(EFFECT), "matrix_model_rows": len(matrix_ids), "breast_metadata_models": metadata_eligible, "breast_models_with_matrix_rows": len(validation_ids)},
        "cohort_definition": {"OncotreeLineage": "exact case-insensitive Breast", "OncotreePrimaryDisease": "exact case-insensitive Invasive Breast Carcinoma", "matrix_row_required_for_validation": True, "eligible_model_ids": sorted(validation_ids), "cohort_file": str(COHORT_PATH)},
        "target_mapping": "Exact HGNC approved symbol; otherwise exact unique alias_symbol/prev_symbol; map via HGNC Entrez ID to exact unique DepMap SYMBOL (Entrez ID) matrix column; no fuzzy mapping.",
        "primary_target_eligibility": {"minimum_finite_breast_model_values": MIN_MODELS, "eligible_target_count": eligible_n, "target_file": str(TARGETS_PATH), "all_unique_mapped_genes_retained_in_file": len(mapping)},
        "primary_outcome": "D_t = -median(finite raw Chronos gene-effect scores across the 50 represented breast models); missing values are not imputed.",
        "primary_predictor": "I_t = abs(unrounded signed efficiency_delta_pct) from the frozen virtual node-deletion efficiency calculation; also export app impact_score rounded to 2 decimals for reproducibility.",
        "knockout_method": {"implementation": "backend/graph_loader.py load_graph + backend/perturbation.py precompute_baseline/global_efficiency; for each target copy graph, remove target, reuse cached 100 baseline efficiency source IDs; signed ΔE% = 100*(Eprime-E0)/E0.", "source_sample_count": 100, "source_seed": 42, "global_efficiency_formula": "sum reciprocal finite weighted shortest distances / (selected_sources*(n-1)); disconnected pairs contribute zero", "strength": "not applicable to virtual node deletion"},
        "primary_statistic": {"name": "Spearman rank correlation", "two_sided": True, "predictor": "I_t", "outcome": "D_t", "p_value": "asymptotic scipy.stats.spearmanr p-value; empirical target-label permutation p-value also reported"},
        "bootstrap": {"unit": "GenX target", "replicates": BOOTSTRAP_B, "seed": BOOTSTRAP_SEED, "interval": "95% percentile CI"},
        "permutation_null": {"procedure": "permute dependency labels across eligible targets, hold structural impact fixed", "replicates": PERMUTATION_B, "seed": PERMUTATION_SEED, "empirical_two_sided_p": "(1 + count(|rho_null| >= |rho_observed|))/(B+1)"},
        "degree_baseline": "Spearman correlation of frozen baseline unweighted node degree and D_t.",
        "degree_matched_null": {"procedure": "permute dependency labels within fixed degree strata; compute Spearman(I_t,D_perm)", "replicates": DEGREE_PERMUTATION_B, "seed": DEGREE_PERMUTATION_SEED, "strata": bin_manifest},
        "partial_spearman_sensitivity": "Average midranks; regress ranked I and ranked D separately on intercept plus ranked baseline degree using numpy.linalg.lstsq; Pearson correlation of residuals; two-sided t reference with df=N-3.",
        "mean_median_sensitivity": "Same ≥10-model target population; D_mean = -mean(finite Chronos); Spearman with I_t.",
        "model_count_sensitivity": "Same mapping and cohort; target threshold ≥20 finite values; report N and Spearman(I_t,D_t).",
        "degree_baseline_and_strata_fixed_without_dependency_values": True,
        "posthoc_optimization": "prohibited",
        "no_biological_validity_claim": True,
        "versions": package_versions,
    }
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"FROZEN: metadata cohort={metadata_eligible}; matrix cohort={len(validation_ids)}; mapped targets={len(mapping)}; primary-eligible targets={eligible_n}; threshold={MIN_MODELS}.")
    print(f"Degree strata fixed before association analysis: {bin_manifest['actual_bins']} bins; no dependency aggregations or GenX scores calculated.")
    print(f"Manifest: {MANIFEST_PATH}")


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    result = stats.spearmanr(x, y, alternative="two-sided", nan_policy="raise")
    return float(result.statistic), float(result.pvalue)


def empirical_two_sided(observed: float, null: np.ndarray) -> float:
    valid = null[np.isfinite(null)]
    return float((1 + np.count_nonzero(np.abs(valid) >= abs(observed))) / (len(valid) + 1))


def rank_partial(x: np.ndarray, y: np.ndarray, covariate: np.ndarray) -> tuple[float, float]:
    rx, ry, rz = (stats.rankdata(v, method="average") for v in (x, y, covariate))
    design = np.column_stack([np.ones(len(rz)), rz])
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    rho = float(np.corrcoef(ex, ey)[0, 1])
    df = len(x) - 3
    if df <= 0 or abs(rho) >= 1:
        p = 0.0 if abs(rho) >= 1 else float("nan")
    else:
        t_value = rho * math.sqrt(df / max(1e-300, 1.0 - rho * rho))
        p = float(2 * stats.t.sf(abs(t_value), df))
    return rho, p


def load_dependency_values(primary_records: list[dict[str, str]], cohort_ids: set[str]) -> dict[str, list[float]]:
    with EFFECT.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        col_by_name = {name: i for i, name in enumerate(header)}
        columns = {r["genx_gene"]: col_by_name[r["depmap_gene_identifier"]] for r in primary_records}
        values = {gene: [] for gene in columns}
        seen = set()
        for row in reader:
            if not row or row[0] not in cohort_ids:
                continue
            if row[0] in seen:
                raise RuntimeError(f"Duplicate matrix model during analysis: {row[0]}")
            seen.add(row[0])
            for gene, col in columns.items():
                value = finite_float(row[col])
                if value is not None:
                    values[gene].append(value)
        if seen != cohort_ids:
            raise RuntimeError(f"Matrix cohort mismatch: observed {len(seen)} expected {len(cohort_ids)}")
    return values


def run_genx_scores(targets: list[str]) -> dict[str, dict[str, float]]:
    os.environ["PATH_COST_MODEL"] = "inverse"
    sys.path.insert(0, str(BACKEND))
    import graph_loader
    import perturbation

    graph = graph_loader.load_graph()
    if graph.graph.get("path_cost_model") != "inverse":
        raise RuntimeError(f"Expected frozen default inverse costs; observed {graph.graph.get('path_cost_model')}")
    if graph.number_of_nodes() != 881 or graph.number_of_edges() != 9681:
        raise RuntimeError("Frozen graph dimensions changed")
    perturbation.precompute_baseline(graph)
    e0 = perturbation._baseline_efficiency
    sources = list(perturbation._baseline_sources)
    if not e0 > 0 or not sources:
        raise RuntimeError("Invalid frozen baseline efficiency or source list")
    scores = {}
    for index, gene in enumerate(targets, start=1):
        if gene not in graph:
            raise RuntimeError(f"Eligible target is not a frozen graph node: {gene}")
        perturbed = graph.copy()
        perturbed.remove_node(gene)
        e_prime = perturbation.global_efficiency(perturbed, sources=sources)
        delta_pct = 100.0 * (e_prime - e0) / e0 if e0 > 0 else 0.0
        scores[gene] = {
            "efficiency_baseline": e0,
            "efficiency_perturbed": e_prime,
            "efficiency_delta_pct": delta_pct,
            "structural_impact": abs(delta_pct),
            "frozen_impact_score_rounded_2dp": round(abs(delta_pct), 2),
            "baseline_degree": float(graph.degree(gene)),
        }
        if index % 25 == 0 or index == len(targets):
            print(f"GenX knockout efficiencies: {index}/{len(targets)}", flush=True)
    return scores


def plot_association(x: np.ndarray, y: np.ndarray, path: Path, title: str, xlabel: str, ylabel: str, rho: float, ci: tuple[float, float] | None = None, p: float | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
    ax.scatter(x, y, s=22, alpha=0.7, linewidths=0, color="#2369a1")
    slope, intercept, _, _ = stats.theilslopes(y, x, alpha=0.95)
    line_x = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    ax.plot(line_x, intercept + slope * line_x, color="#c14b35", linewidth=2, label="Theil–Sen monotone linear trend")
    subtitle = f"N={len(x)}   Spearman ρ={rho:.3f}"
    if ci is not None:
        subtitle += f"   95% bootstrap CI [{ci[0]:.3f}, {ci[1]:.3f}]"
    if p is not None:
        subtitle += f"   permutation p={p:.4g}"
    ax.set_title(title + "\n" + subtitle, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def analyze() -> None:
    if not MANIFEST_PATH.is_file() or not TARGETS_PATH.is_file() or not COHORT_PATH.is_file():
        raise RuntimeError("Run freeze phase first; no validation statistics have been calculated")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if manifest.get("status") != "FROZEN_BEFORE_ASSOCIATION_ANALYSIS":
        raise RuntimeError("Analysis manifest is not in frozen pre-analysis state")
    current_hashes = verify_frozen_inputs()
    if current_hashes != manifest["inputs_sha256"]:
        raise RuntimeError("Frozen input checksum changed after cohort/target freeze")

    cohort = csv_rows(COHORT_PATH)
    cohort_ids = {r["model_id"] for r in cohort if r["validation_inclusion_status"] == "included"}
    if len(cohort_ids) != 50:
        raise RuntimeError(f"Frozen analysis cohort expected 50 models, got {len(cohort_ids)}")
    all_targets = csv_rows(TARGETS_PATH)
    primary_records = [r for r in all_targets if r["primary_analysis_eligible"].casefold() == "true"]
    if len(primary_records) != manifest["primary_target_eligibility"]["eligible_target_count"]:
        raise RuntimeError("Frozen eligible target count changed")

    # Aggregate experimental measures and export the full target table before
    # running the primary correlation, bootstrap, or permutation tests.
    dep_values = load_dependency_values(primary_records, cohort_ids)
    for record in primary_records:
        vals = np.asarray(dep_values[record["genx_gene"]], dtype=float)
        n_finite = len(vals)
        expected_count = int(record["finite_breast_model_dependency_count"])
        if n_finite != expected_count:
            raise RuntimeError(f"Frozen finite count mismatch for {record['genx_gene']}: {n_finite} vs {expected_count}")
        record["chronos_median"] = float(np.median(vals))
        record["dependency_oriented_negative_median_chronos"] = -float(np.median(vals))
        record["chronos_mean"] = float(np.mean(vals))
        record["dependency_oriented_negative_mean_chronos"] = -float(np.mean(vals))
        record["chronos_sd_sample"] = float(np.std(vals, ddof=1)) if n_finite > 1 else 0.0
        record["chronos_min"] = float(np.min(vals))
        record["chronos_max"] = float(np.max(vals))
    genx_scores = run_genx_scores([r["genx_gene"] for r in primary_records])
    target_rows = []
    for r in primary_records:
        g = genx_scores[r["genx_gene"]]
        target_rows.append({
            "genx_gene": r["genx_gene"],
            "depmap_gene_identifier": r["depmap_gene_identifier"],
            "genx_structural_knockout_delta_efficiency_pct": g["efficiency_delta_pct"],
            "genx_absolute_knockout_impact_unrounded": g["structural_impact"],
            "genx_impact_score_frozen_round_2dp": g["frozen_impact_score_rounded_2dp"],
            "efficiency_baseline": g["efficiency_baseline"],
            "efficiency_after_virtual_knockout": g["efficiency_perturbed"],
            "baseline_graph_degree": int(g["baseline_degree"]),
            "degree_stratum": r["degree_stratum"],
            "finite_breast_model_dependency_count": len(dep_values[r["genx_gene"]]),
            "median_chronos": r["chronos_median"],
            "dependency_score_negative_median_chronos": r["dependency_oriented_negative_median_chronos"],
            "mean_chronos": r["chronos_mean"],
            "dependency_score_negative_mean_chronos": r["dependency_oriented_negative_mean_chronos"],
            "chronos_sd_sample": r["chronos_sd_sample"],
            "minimum_chronos": r["chronos_min"],
            "maximum_chronos": r["chronos_max"],
        })
    table_fields = list(target_rows[0])
    write_csv(TARGET_TABLE, table_fields, target_rows)

    # Primary result: Spearman over target-level observations.
    x = np.asarray([r["genx_absolute_knockout_impact_unrounded"] for r in target_rows], dtype=float)
    y = np.asarray([r["dependency_score_negative_median_chronos"] for r in target_rows], dtype=float)
    degree = np.asarray([r["baseline_graph_degree"] for r in target_rows], dtype=float)
    rho, asymptotic_p = spearman(x, y)

    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    boot = np.empty(BOOTSTRAP_B, dtype=float)
    for i in range(BOOTSTRAP_B):
        indexes = bootstrap_rng.integers(0, len(x), size=len(x))
        boot[i], _ = spearman(x[indexes], y[indexes])
    valid_boot = boot[np.isfinite(boot)]
    ci_low, ci_high = np.quantile(valid_boot, [0.025, 0.975])

    perm_rng = np.random.default_rng(PERMUTATION_SEED)
    null = np.empty(PERMUTATION_B, dtype=float)
    for i in range(PERMUTATION_B):
        null[i], _ = spearman(x, perm_rng.permutation(y))
    permutation_p = empirical_two_sided(rho, null)

    degree_rho, degree_p = spearman(degree, y)
    partial_rho, partial_p = rank_partial(x, y, degree)
    mean_y = np.asarray([r["dependency_score_negative_mean_chronos"] for r in target_rows], dtype=float)
    mean_rho, mean_p = spearman(x, mean_y)

    # Strata were frozen in the target-population file before response values
    # were aggregated. Permute y only within each fixed degree bin.
    bins_by_label: dict[int, list[int]] = defaultdict(list)
    for i, r in enumerate(target_rows):
        bins_by_label[int(r["degree_stratum"])].append(i)
    degree_rng = np.random.default_rng(DEGREE_PERMUTATION_SEED)
    degree_null = np.empty(DEGREE_PERMUTATION_B, dtype=float)
    for b in range(DEGREE_PERMUTATION_B):
        permuted = y.copy()
        for indexes in bins_by_label.values():
            permuted[indexes] = degree_rng.permutation(y[indexes])
        degree_null[b], _ = spearman(x, permuted)
    degree_perm_p = empirical_two_sided(rho, degree_null)

    stricter = [r for r in all_targets if int(r["finite_breast_model_dependency_count"]) >= SENSITIVITY_MIN_MODELS]
    stricter_rho = stricter_p = float("nan")
    if len(stricter) >= 3:
        stricter_dep = load_dependency_values(stricter, cohort_ids)
        sx = np.asarray([genx_scores[r["genx_gene"]]["structural_impact"] for r in stricter], dtype=float)
        sy = np.asarray([-np.median(stricter_dep[r["genx_gene"]]) for r in stricter], dtype=float)
        stricter_rho, stricter_p = spearman(sx, sy)

    # Save all numerical outputs and null draws for a reproducible record.
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    primary_result = [{
        "N": len(x), "spearman_rho": rho, "asymptotic_p": asymptotic_p,
        "bootstrap_ci_lower": ci_low, "bootstrap_ci_upper": ci_high,
        "bootstrap_B": BOOTSTRAP_B, "bootstrap_valid_replicates": len(valid_boot),
        "bootstrap_seed": BOOTSTRAP_SEED, "permutation_p": permutation_p,
        "permutation_B": PERMUTATION_B, "permutation_seed": PERMUTATION_SEED,
        "degree_matched_permutation_p": degree_perm_p,
        "degree_matched_permutation_B": DEGREE_PERMUTATION_B,
        "degree_matched_permutation_seed": DEGREE_PERMUTATION_SEED,
    }]
    write_csv(RESULTS_DIR / "depmap_primary_result.csv", list(primary_result[0]), primary_result)
    sensitivities = [
        {"analysis": "degree_baseline", "N": len(x), "spearman_rho": degree_rho, "p_value": degree_p, "details": "Spearman(baseline unweighted degree, -median Chronos)"},
        {"analysis": "partial_spearman_controlling_degree", "N": len(x), "spearman_rho": partial_rho, "p_value": partial_p, "details": "Pearson correlation of residualized average midranks; control=average-ranked degree; scipy/numpy least squares"},
        {"analysis": "mean_chronos_sensitivity", "N": len(x), "spearman_rho": mean_rho, "p_value": mean_p, "details": "Same primary target set; outcome=-mean Chronos"},
        {"analysis": "minimum_20_finite_models", "N": len(stricter), "spearman_rho": stricter_rho, "p_value": stricter_p, "details": "Same mappings/cohort; outcome=-median Chronos; threshold >=20"},
        {"analysis": "degree_matched_permutation_null", "N": len(x), "spearman_rho": rho, "p_value": degree_perm_p, "details": f"Within {len(bins_by_label)} frozen degree strata; B={DEGREE_PERMUTATION_B}; seed={DEGREE_PERMUTATION_SEED}"},
    ]
    write_csv(RESULTS_DIR / "depmap_sensitivity_results.csv", list(sensitivities[0]), sensitivities)
    write_csv(RESULTS_DIR / "depmap_permutation_null.csv", ["permutation", "spearman_rho"], [{"permutation": i + 1, "spearman_rho": float(v)} for i, v in enumerate(null)])
    write_csv(RESULTS_DIR / "depmap_degree_matched_null.csv", ["permutation", "spearman_rho"], [{"permutation": i + 1, "spearman_rho": float(v)} for i, v in enumerate(degree_null)])
    write_csv(RESULTS_DIR / "depmap_bootstrap_rhos.csv", ["replicate", "spearman_rho"], [{"replicate": i + 1, "spearman_rho": (float(v) if np.isfinite(v) else "")} for i, v in enumerate(boot)])

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plot_association(x, y, FIGURES_DIR / "depmap_primary_association.png", "GenX virtual knockout impact and breast-cancer dependency", "GenX virtual knockout structural impact |ΔE| (%)", "Experimental dependency (−median Chronos)", rho, (float(ci_low), float(ci_high)), permutation_p)
    plot_association(degree, y, FIGURES_DIR / "depmap_degree_baseline.png", "Baseline degree and breast-cancer dependency", "Frozen STRING graph baseline degree", "Experimental dependency (−median Chronos)", degree_rho, p=degree_p)
    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ax.hist(null, bins=45, color="#7ca6c8", edgecolor="white")
    ax.axvline(rho, color="#bd3c32", linewidth=2, label=f"Observed ρ = {rho:.3f}")
    ax.axvline(-rho, color="#bd3c32", linewidth=1, linestyle="--", alpha=0.65, label="Two-sided reference")
    ax.set_xlabel("Spearman ρ under target-label permutation")
    ax.set_ylabel("Permutation count")
    ax.set_title(f"Target-level permutation null (B={PERMUTATION_B}; empirical p={permutation_p:.4g})")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(FIGURES_DIR / "depmap_permutation_null.png", dpi=300)
    plt.close(fig)

    # Readable results narrative, strictly evidence-calibrated.
    summary = f'''# GenX virtual knockout vs DepMap dependency: first validation

## Frozen population and data coverage

The metadata rule selected 83 breast-cancer models (`OncotreeLineage = Breast`; `OncotreePrimaryDisease = Invasive Breast Carcinoma`). Of those, 50 have rows in `CRISPRGeneEffect.csv` and contribute measurements; the 50 are not the complete set of metadata-eligible DepMap breast models. The other 33 lack matrix rows. Among the 50 represented models, 619 mapped GenX-gene/model cells are missing. Targets were eligible at ≥{MIN_MODELS} finite measurements, fixed before score calculation; **primary N = {len(x)}**. No missing values were imputed.

## Primary association

Predictor: absolute signed GenX virtual-knockout efficiency change, `|ΔE|`, computed with the existing implementation and inverse path-cost default. Outcome: `−median(Chronos)` across each target's finite values in the represented breast models. Spearman ρ = **{rho:.6f}**; 95% target-bootstrap percentile CI **[{ci_low:.6f}, {ci_high:.6f}]** (B={BOOTSTRAP_B}, seed={BOOTSTRAP_SEED}); asymptotic two-sided p = **{asymptotic_p:.6g}**; target-label permutation empirical two-sided p = **{permutation_p:.6g}** (B={PERMUTATION_B}, seed={PERMUTATION_SEED}).

## Prespecified sensitivities

- Degree baseline: Spearman( baseline degree, `−median(Chronos)` ) = {degree_rho:.6f}, p={degree_p:.6g}.
- Degree-adjusted partial Spearman (rank residualization): ρpartial={partial_rho:.6f}, p={partial_p:.6g}.
- Mean-based outcome on the same targets: Spearman ρ={mean_rho:.6f}, p={mean_p:.6g}.
- ≥20 finite-model sensitivity: N={len(stricter)}, Spearman ρ={stricter_rho:.6f}, p={stricter_p:.6g}.
- Within-degree-stratum permutation: empirical two-sided p={degree_perm_p:.6g} (B={DEGREE_PERMUTATION_B}, seed={DEGREE_PERMUTATION_SEED}; {len(bins_by_label)} fixed strata).

## Interpretation boundary

The estimates describe an association between the frozen GenX structural-impact metric and DepMap Chronos knockout growth/fitness effects in this target set and these 50 represented breast models. Positive association with uncertainty/permutation evidence would be supporting evidence for this specific relationship; a small or uncertain estimate indicates weak evidence; a near-zero estimate compatible with its null indicates no detectable relationship under this design. No arbitrary pass threshold is used. This analysis does not establish causal STRING signaling, equivalence of graph fragmentation to essentiality, perfect essentiality prediction, mechanistic viability prediction, patient-response prediction, or clinical utility. Degree and degree-aware results are comparisons/sensitivities, not winners.

## Reproducibility

Analysis manifest frozen before association calculation: `{MANIFEST_PATH}`. Full target table: `{TARGET_TABLE}`. Inputs and graph were SHA-256 verified against the frozen manifest. Raw DepMap files and GenX source/graph/specification were not modified. No L1000 data were used.
'''
    (RESULTS_DIR / "depmap_validation_summary.md").write_text(summary, encoding="utf-8")

    # Record the post-analysis state without altering the frozen protocol fields.
    manifest["analysis_completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["result_files"] = {
        "target_table": str(TARGET_TABLE),
        "primary_result": str(RESULTS_DIR / "depmap_primary_result.csv"),
        "sensitivity_results": str(RESULTS_DIR / "depmap_sensitivity_results.csv"),
        "bootstrap_draws": str(RESULTS_DIR / "depmap_bootstrap_rhos.csv"),
        "permutation_draws": str(RESULTS_DIR / "depmap_permutation_null.csv"),
        "degree_matched_draws": str(RESULTS_DIR / "depmap_degree_matched_null.csv"),
        "figures": [str(FIGURES_DIR / name) for name in ("depmap_primary_association.png", "depmap_degree_baseline.png", "depmap_permutation_null.png")],
    }
    manifest["observed_results_added_after_freeze"] = {"N": len(x), "spearman_rho": rho, "bootstrap_ci": [float(ci_low), float(ci_high)], "permutation_p": permutation_p}
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"PRIMARY N={len(x)}; Spearman rho={rho:.6f}; bootstrap CI=({ci_low:.6f}, {ci_high:.6f}); permutation p={permutation_p:.6g}")
    print(f"Degree baseline rho={degree_rho:.6f}; degree-matched null p={degree_perm_p:.6g}; mean sensitivity rho={mean_rho:.6f}; >=20 N={len(stricter)} rho={stricter_rho:.6f}")
    print(f"Outputs: {RESULTS_DIR}; figures: {FIGURES_DIR}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in {"freeze", "analyze"}:
        raise SystemExit("Usage: python run_depmap_knockout_validation.py freeze|analyze")
    freeze() if sys.argv[1] == "freeze" else analyze()
