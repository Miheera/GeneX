"""Compare network knockout impact scores with TCGA tumor expression activity."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import networkx as nx
from scipy.stats import spearmanr

try:
    import perturbation
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, str(Path(__file__).parents[1]))
    import perturbation

METHODOLOGY_NOTE = (
    "Structural importance is the absolute percentage change in approximate weighted global efficiency after a gene knockout "
    f"({perturbation.APPROXIMATION_METHOD}). "
    "The primary tumor-activity measure is median RNA-seq z-score across the TCGA BRCA PanCancer Atlas sample list; "
    "raw median expression is retained for comparison. Spearman rho measures monotonic association only; it does not "
    "establish that expression causes network importance."
)
OUTPUT_COLUMNS = [
    "gene_symbol", "entrez_id", "structural_score", "median_expression",
    "median_zscore", "pct_overexpressed", "n_samples", "zscore_spearman_rho",
    "zscore_spearman_pvalue", "raw_spearman_rho", "raw_spearman_pvalue",
    "mutation_frequency_pct", "prognostic_significance",
    "mutation_spearman_rho", "mutation_spearman_pvalue",
    "prognostic_spearman_rho", "prognostic_spearman_pvalue",
]

_BACKEND = Path(__file__).parents[1]
DEFAULT_EXPRESSION_PATH = _BACKEND / "data" / "processed" / "tcga_expression_scores.csv"
DEFAULT_MUTATION_PATH = _BACKEND / "data" / "processed" / "tcga_mutation_scores.csv"
DEFAULT_SURVIVAL_PATH = _BACKEND / "data" / "analysis_outputs" / "survival_analysis.csv"
DEFAULT_OUTPUT_PATH = _BACKEND / "data" / "analysis_outputs" / "structural_expression_correlation.csv"


def _read_expression_scores(path: Path) -> dict[str, dict[str, str]]:
    """Read and validate the processed TCGA expression-score file."""
    if not path.exists():
        raise FileNotFoundError(f"TCGA expression scores are not available: {path}")
    required = {
        "gene_symbol", "entrez_id", "median_expression", "median_zscore",
        "pct_overexpressed", "n_samples",
    }
    with path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("TCGA expression scores CSV has an invalid schema")
        scores: dict[str, dict[str, str]] = {}
        for row in reader:
            symbol = row["gene_symbol"].strip()
            if not symbol:
                continue
            try:
                raw_expression = float(row["median_expression"])
                zscore_expression = float(row["median_zscore"])
                float(row["pct_overexpressed"])
                int(row["n_samples"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid expression score for {symbol}") from exc
            if not math.isfinite(raw_expression) or not math.isfinite(zscore_expression):
                raise ValueError(f"Non-finite expression score for {symbol}")
            scores[symbol] = row
    return scores


def _structural_scores(graph: nx.Graph) -> dict[str, float]:
    """Calculate approximate knockout efficiency impact without centrality overhead."""
    perturbation.precompute_baseline(graph)
    return {
        str(gene): float(perturbation.knockout_impact_score(graph, str(gene)))
        for gene in graph
    }


def _read_optional_scores(path: Path, required: set[str]) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid schema in {path}")
        return {row["gene"].strip(): row for row in reader if row.get("gene", "").strip()}


def _correlation(values_a: list[float], values_b: list[float]) -> tuple[float, float]:
    result = spearmanr(values_a, values_b)
    if not (math.isfinite(result.statistic) and math.isfinite(result.pvalue)):
        raise ValueError("Spearman correlation is undefined for the available values")
    return round(float(result.statistic), 6), round(float(result.pvalue), 6)


def analyze_structural_expression(
    graph: nx.Graph,
    expression_path: Path = DEFAULT_EXPRESSION_PATH,
    *,
    mutation_path: Path = DEFAULT_MUTATION_PATH,
    survival_path: Path = DEFAULT_SURVIVAL_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
) -> dict[str, object]:
    """Join TCGA expression measures to network impact and calculate Spearman statistics."""
    expression_scores = _read_expression_scores(expression_path)
    mutation_scores = _read_optional_scores(mutation_path, {"gene", "mutation_frequency_pct"})
    survival_scores = _read_optional_scores(survival_path, {"gene", "logrank_p_value"})
    structural_scores = _structural_scores(graph)
    shared_genes = sorted(set(structural_scores) & set(expression_scores))
    if len(shared_genes) < 2:
        raise ValueError("At least two genes shared by the graph and expression scores are required")

    rows: list[dict[str, str | int | float]] = []
    for gene in shared_genes:
        expression = expression_scores[gene]
        rows.append({
            "gene_symbol": gene,
            "entrez_id": int(expression["entrez_id"]),
            "structural_score": round(structural_scores[gene], 6),
            "median_expression": float(expression["median_expression"]),
            "median_zscore": float(expression["median_zscore"]),
            "pct_overexpressed": float(expression["pct_overexpressed"]),
            "n_samples": int(expression["n_samples"]),
        })

    structural_values = [float(row["structural_score"]) for row in rows]
    zscore_result = spearmanr(structural_values, [float(row["median_zscore"]) for row in rows])
    raw_result = spearmanr(structural_values, [float(row["median_expression"]) for row in rows])
    if not all(math.isfinite(value) for value in (
        zscore_result.statistic, zscore_result.pvalue, raw_result.statistic, raw_result.pvalue,
    )):
        raise ValueError("Spearman correlation is undefined for the available expression values")

    zscore_rho = round(float(zscore_result.statistic), 6)
    zscore_pvalue = round(float(zscore_result.pvalue), 6)
    raw_rho = round(float(raw_result.statistic), 6)
    raw_pvalue = round(float(raw_result.pvalue), 6)
    mutation_genes = sorted(set(structural_scores) & set(mutation_scores))
    mutation_rho = mutation_pvalue = None
    if len(mutation_genes) >= 2:
        mutation_rho, mutation_pvalue = _correlation(
            [structural_scores[gene] for gene in mutation_genes],
            [float(mutation_scores[gene]["mutation_frequency_pct"]) for gene in mutation_genes],
        )
    prognostic_genes = []
    for gene in sorted(set(structural_scores) & set(survival_scores)):
        pvalue = float(survival_scores[gene]["logrank_p_value"])
        if pvalue > 0 and math.isfinite(pvalue):
            prognostic_genes.append((gene, -math.log10(pvalue)))
    prognostic_rho = prognostic_pvalue = None
    if len(prognostic_genes) >= 2:
        prognostic_rho, prognostic_pvalue = _correlation(
            [structural_scores[gene] for gene, _ in prognostic_genes],
            [score for _, score in prognostic_genes],
        )
    for row in rows:
        gene = str(row["gene_symbol"])
        row.update({
            "zscore_spearman_rho": zscore_rho,
            "zscore_spearman_pvalue": zscore_pvalue,
            "raw_spearman_rho": raw_rho,
            "raw_spearman_pvalue": raw_pvalue,
            "mutation_frequency_pct": float(mutation_scores[gene]["mutation_frequency_pct"]) if gene in mutation_scores else "",
            "prognostic_significance": next((score for candidate, score in prognostic_genes if candidate == gene), ""),
            "mutation_spearman_rho": mutation_rho if mutation_rho is not None else "",
            "mutation_spearman_pvalue": mutation_pvalue if mutation_pvalue is not None else "",
            "prognostic_spearman_rho": prognostic_rho if prognostic_rho is not None else "",
            "prognostic_spearman_pvalue": prognostic_pvalue if prognostic_pvalue is not None else "",
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    return {
        "structural_metric": "knockout_impact_score",
        "approximation_method": perturbation.APPROXIMATION_METHOD,
        "expression_metric": "median_zscore",
        "spearman_rho": zscore_rho,
        "spearman_pvalue": zscore_pvalue,
        "raw_expression_correlation": {
            "expression_metric": "median_expression",
            "spearman_rho": raw_rho,
            "spearman_pvalue": raw_pvalue,
        },
        "mutation_correlation": {"spearman_rho": mutation_rho, "spearman_pvalue": mutation_pvalue, "n_genes": len(mutation_genes)},
        "prognostic_correlation": {"spearman_rho": prognostic_rho, "spearman_pvalue": prognostic_pvalue, "n_genes": len(prognostic_genes)},
        "n_genes": len(rows),
        "methodology_note": METHODOLOGY_NOTE,
        "points": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression-path", type=Path, default=DEFAULT_EXPRESSION_PATH)
    parser.add_argument("--mutation-path", type=Path, default=DEFAULT_MUTATION_PATH)
    parser.add_argument("--survival-path", type=Path, default=DEFAULT_SURVIVAL_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    args = parser.parse_args()
    import graph_loader

    result = analyze_structural_expression(
        graph_loader.load_graph(), args.expression_path, mutation_path=args.mutation_path,
        survival_path=args.survival_path, output_path=args.output_path
    )
    print("\n=== structural-expression correlation summary ===")
    print(f"  Shared genes: {result['n_genes']}")
    print(f"  Z-score rho:  {result['spearman_rho']:.6f} (p={result['spearman_pvalue']:.6f})")
    raw = result["raw_expression_correlation"]
    print(f"  Raw rho:      {raw['spearman_rho']:.6f} (p={raw['spearman_pvalue']:.6f})")
    mutation = result["mutation_correlation"]
    prognostic = result["prognostic_correlation"]
    print(f"  Mutation rho: {mutation['spearman_rho']} (p={mutation['spearman_pvalue']})")
    print(f"  Survival rho: {prognostic['spearman_rho']} (p={prognostic['spearman_pvalue']})")
    print(f"  Output:       {args.output_path}")
    print()


if __name__ == "__main__":
    main()
