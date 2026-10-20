#!/usr/bin/env python3
"""Create the paired publication figure for the frozen GenX/DepMap comparison.

Run ``python validate_depmap.py`` first to refresh the target-level points and
summary. DepMap Public 26Q1 provides Chronos scores in this analysis; the plotted
outcome is ``-median Chronos`` so higher values mean greater dependency.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "validation" / "results"
DEFAULT_POINTS = RESULTS / "depmap_tissue_rwr_points.csv"
DEFAULT_SUMMARY = RESULTS / "depmap_tissue_rwr_comparison.csv"
DEFAULT_OUTPUT = RESULTS / "depmap_correlation_figure"

X_COLUMNS = (
    "unweighted_topological_efficiency_drop_pct",
    "tissue_weighted_rwr_l1_flux",
)
Y_COLUMN = "dependency_score_negative_median_chronos"
PREDICTORS = (
    "unweighted_topological_efficiency_drop",
    "tissue_weighted_rwr_l1_flux",
)


def load_inputs(points_path: Path, summary_path: Path) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    points = pd.read_csv(points_path)
    required = {*X_COLUMNS, Y_COLUMN}
    missing = required.difference(points.columns)
    if missing:
        raise ValueError(f"Points CSV is missing required columns: {sorted(missing)}")
    points = points[list(X_COLUMNS) + [Y_COLUMN]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(points) < 3:
        raise ValueError("At least three complete target rows are required to plot regressions")

    with summary_path.open(newline="", encoding="utf-8-sig") as source:
        rows = list(csv.DictReader(source))
    by_predictor = {row["predictor"]: row for row in rows}
    missing_stats = set(PREDICTORS).difference(by_predictor)
    if missing_stats:
        raise ValueError(f"Comparison CSV is missing predictors: {sorted(missing_stats)}")
    return points, by_predictor


def plot(points_path: Path = DEFAULT_POINTS, summary_path: Path = DEFAULT_SUMMARY,
         output_base: Path = DEFAULT_OUTPUT) -> tuple[Path, Path]:
    points, summary = load_inputs(points_path, summary_path)
    sns.set_theme(style="whitegrid", context="paper", font="DejaVu Sans")
    plt.rcParams.update({
        "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": "#d9e0e8", "grid.linewidth": 0.65,
    })
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.7), sharey=True,
                             constrained_layout=True)
    specs = [
        (X_COLUMNS[0], PREDICTORS[0], "A  Unweighted topology",
         r"Unweighted topological efficiency drop ($\Delta E\%$)", "#426b8a"),
        (X_COLUMNS[1], PREDICTORS[1], "B  Tissue-weighted RWR",
         r"Geometric tissue-weighted RWR flux magnitude ($\|\Delta r\|_1$)", "#b44b43"),
    ]

    for ax, (x_col, predictor, title, xlabel, color) in zip(axes, specs):
        row = summary[predictor]
        sns.regplot(data=points, x=x_col, y=Y_COLUMN, ax=ax, ci=95,
                    scatter_kws={"s": 17, "alpha": 0.43, "color": color,
                                 "edgecolor": "none", "rasterized": True},
                    line_kws={"color": color, "linewidth": 2},
                    color=color, truncate=False)
        # Structural efficiency changes are zero for most targets and have a
        # small number of large positive values; a symlog x-axis exposes the
        # dense near-zero observations without discarding those targets.
        if predictor == PREDICTORS[0]:
            ax.set_xscale("symlog", linthresh=0.05)
        rho = float(row["spearman_rho"])
        lo = float(row["bootstrap_ci_lower"])
        hi = float(row["bootstrap_ci_upper"])
        if predictor == PREDICTORS[0]:
            stats_text = (f"N = {len(points)}\nSpearman $\\rho$ = {rho:+.4f}\n"
                          f"$p$ = {float(row['two_sided_p_value']):.4f} (not significant)\n"
                          f"95% bootstrap CI [{lo:+.4f}, {hi:+.4f}]")
        else:
            stats_text = (f"N = {len(points)}\nSpearman $\\rho$ = {rho:+.4f}\n"
                          f"FDR $q$ = {float(row['fdr_q_value']):.2e}\n"
                          f"95% bootstrap CI [{lo:+.4f}, {hi:+.4f}]")
        ax.text(0.04, 0.96, stats_text, transform=ax.transAxes, ha="left", va="top",
                fontsize=9, bbox={"boxstyle": "round,pad=0.45", "facecolor": "white",
                                  "edgecolor": "#d5dce5", "alpha": 0.96})
        ax.set_title(title, loc="left", weight="bold", pad=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("" if ax is axes[1] else "DepMap dependency (−median Chronos)")
        ax.set_axisbelow(True)
        ax.grid(True, alpha=0.8)

    fig.suptitle("GenX network perturbation metrics and breast cancer dependency",
                 fontsize=14, weight="bold")
    fig.text(0.5, -0.025,
             "Shaded bands show 95% OLS confidence intervals for the mean; reported intervals are\n"
             "1,000-resample bootstrap confidence intervals for Spearman correlation. "
             "Panel A uses a symlog x-axis; FDR uses Benjamini–Hochberg across three predictors.",
             ha="center", va="top", fontsize=8.5, color="#526173")

    output_base.parent.mkdir(parents=True, exist_ok=True)
    png_path, svg_path = output_base.with_suffix(".png"), output_base.with_suffix(".svg")
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved 300 DPI PNG: {png_path}")
    print(f"Saved vector SVG: {svg_path}")
    print(f"Rows plotted: {len(points)}")
    return png_path, svg_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plot(args.points, args.summary, args.output_base)


if __name__ == "__main__":
    main()
