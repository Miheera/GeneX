# GenX Workbench: validation and interpretation

> **GenX is an Interactive Exploratory Network Analysis & Computational Hypothesis Generation Workbench.** Its outputs prioritize candidate target combinations for downstream wet-lab validation (for example, dual-CRISPR or Perturb-seq); they are not clinical diagnostic predictors.

This document records the current frozen DepMap analysis, the parameter sensitivity results, the scope of the shipped calculations, and the verification status. GenX is a breast-cancer-focused research prototype. Its scores are network-derived hypotheses, not experimentally calibrated responses.

## Empirical validation summary

The frozen comparison uses DepMap Public 26Q1 Chronos gene-effect data. For each eligible target, the phenotype is the negative median Chronos value across the 50 breast cancer cell lines represented in the gene-effect matrix; larger values therefore indicate stronger dependency. The gene-level analysis includes 849 targets with at least 10 finite model values. The metadata query identified 83 breast models, but only 50 had matrix rows. No missing values were imputed.

| Predictor | N targets | Spearman ρ | Two-sided p | 95% bootstrap CI for ρ | BH FDR q | Interpretation |
|---|---:|---:|---:|---:|---:|---|
| Unweighted topological efficiency drop | 849 | +0.0146 | 0.6719 | [−0.0494, +0.0859] | 0.8767 | No detectable association in this comparison |
| Geometric tissue-weighted RWR L₁ flux magnitude | 849 | +0.1706 | 5.73 × 10⁻⁷ | [+0.1017, +0.2360] | 1.72 × 10⁻⁶ | Positive retrospective association; not a causal or prospective prediction |
| Frozen inverse-STRING efficiency drop | 849 | +0.0053 | 0.8767 | [−0.0616, +0.0714] | 0.8767 | No detectable association in this comparison |

The multiple-testing family contains these three predictor associations and uses Benjamini–Hochberg correction. The published figure is generated from the frozen per-target points and this summary table by `python plot_depmap_correlation.py`; it writes a 300-DPI PNG and a vector SVG under `validation/results/`. The bootstrap intervals above are from the current checked-in 1,000-resample output. They differ slightly from earlier draft intervals, so the checked-in results are the source of truth.

The positive RWR association is an exploratory retrospective association within this frozen target set and cohort. It does not establish that STRING edges represent causal signaling, that RWR flux is a biological response, or that the score predicts an individual cell-line or patient outcome. The structural efficiency comparator is near zero in this evaluation.

## Parameter sensitivity regimes

`validation/results/sensitivity_metrics.csv` contains 36 RWR rank-stability settings across four STRING thresholds (800, 850, 900, 950), three restart values (0.70, 0.85, 0.90), and three expression transforms (min-max, log2-TPM, rank), compared with the reference configuration. **15/36** RWR settings have rank R² greater than 0.85; the remaining settings show that rankings are not uniformly stable. The efficiency sampling sweep contains 16 settings; 3/16 exceed rank R² of 0.85. These are rank-agreement statistics against the reference configuration, not predictive-model R² or evidence that the reference is biologically optimal.

GenX designates geometric expression weighting, RWR α = 0.85, and STRING combined score c ≥ 900 as its primary **operational regime** because they define the checked-in graph and default analysis parameters and provide a high-confidence, reproducible baseline for exploration. This designation is a prespecified software default, not a claim that sensitivity testing selected a universal optimum. Expression weighting uses the available TCGA tumor profile; an equivalent CCLE breast-expression comparison is not reported because no CCLE expression profile is included in the repository. Sensitivity output should be regenerated with `python sensitivity_analysis.py` after changing the code or inputs.

## Tool positioning

| Platform | Established focus | GenX relationship and boundary |
|---|---|---|
| Cytoscape | Broad network analysis ecosystem and configurable visualization | GenX offers a focused web workbench that links one breast-cancer graph to differential RWR and virtual knockout calculations. It complements rather than replaces Cytoscape's broader ecosystem. No subsecond performance claim is made here without a dedicated benchmark. |
| Open Targets Platform | Disease-target evidence aggregation and prioritization | GenX shows local network topology and imported annotation context. It does not query or reproduce Open Targets' evidence platform or its disease-association score. |
| GEARS / PDGrapher | Learned perturbation-response prediction models | GenX does not implement these models, a GNN, or an equivalent prediction engine. Its graph is an exploratory interface for its own network calculations; it is not an explainability layer for their outputs. |

## Verification status

The current backend integration suite has **59/59 passing pytest cases** covering graph ingestion, RWR solvers, perturbation handlers, and safe DGIdb annotation joins. Run from the repository root:

```powershell
$env:PYTHONPATH = (Resolve-Path backend).Path
python -m pytest -x -q --confcutdir=. --basetemp=validation\results\.pytest-temp backend\tests test_backend.py
```

These tests exercise API handlers and simulation behavior in-process; they should not be read as a deployment, concurrency, or production-load certification. The validation and sensitivity scripts are separate analyses with their own frozen inputs and generated outputs.

## 11 methodological boundaries

1. **Association network:** STRING combined scores represent gene association evidence; they are not necessarily physical interactions, direction, activation, or causal effects.
2. **Disease scope:** The shipped graph and cohort analyses are breast-cancer-focused. They do not establish performance in other cancers or tissues.
3. **Expression context:** Tissue weighting uses an expression profile to contextualize edges. It does not model regulation, protein abundance, or causal tissue-specific signaling.
4. **RWR meaning:** RWR scores describe diffusion over a chosen graph operator. They are not biological time courses, kinetic rates, or dose-response curves.
5. **Baseline choice:** Differential scores compare target-seeded diffusion with a uniform-seed baseline. This baseline is a computational reference, not an experimentally measured healthy or untreated state.
6. **Inhibition and overexpression:** These modes apply a defined sign convention to a fixed-graph differential score. They do not physically modify gene expression or simulate downstream molecular mechanisms.
7. **Virtual knockout:** Node deletion and resulting efficiency or betweenness changes measure graph structure. They are not direct estimates of cell death or gene essentiality.
8. **Combinatorial synergy:** A non-additive network score is a mathematical comparison of single- and dual-seed shifts. It is not evidence of synthetic lethality without experimental validation.
9. **DepMap association:** The reported association is retrospective and aggregated across 50 represented breast models and 849 eligible targets. It is not an independent prospective validation, causal test, or patient-level prediction.
10. **Annotations:** DGIdb, COSMIC, and other imported records provide evidence or interaction claims. Their presence does not establish drug approval, target engagement, efficacy, or suitability for a patient.
11. **Intended use:** GenX is for interactive exploratory network analysis and computational hypothesis generation to prioritize downstream experiments such as dual-CRISPR or Perturb-seq. It is not a deterministic ODE kinetic simulator, clinical diagnostic predictor, or treatment recommendation system.

## Reproduction files

- Points: `validation/results/depmap_tissue_rwr_points.csv`
- Association statistics: `validation/results/depmap_tissue_rwr_comparison.csv`
- Sensitivity report: `validation/results/sensitivity_metrics.csv`
- Figure generator: `plot_depmap_correlation.py`
- Validation implementation and provenance: `validate_depmap.py`, `validation/depmap_validation_analysis_manifest.json`
