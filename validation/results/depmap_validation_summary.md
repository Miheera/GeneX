# GenX virtual knockout vs DepMap dependency: first validation

## Frozen population and data coverage

The metadata rule selected 83 breast-cancer models (`OncotreeLineage = Breast`; `OncotreePrimaryDisease = Invasive Breast Carcinoma`). Of those, 50 have rows in `CRISPRGeneEffect.csv` and contribute measurements; the 50 are not the complete set of metadata-eligible DepMap breast models. The other 33 lack matrix rows. Among the 50 represented models, 619 mapped GenX-gene/model cells are missing. Targets were eligible at ≥10 finite measurements, fixed before score calculation; **primary N = 849**. No missing values were imputed.

## Primary association

Predictor: absolute signed GenX virtual-knockout efficiency change, `|ΔE|`, computed with the existing implementation and inverse path-cost default. Outcome: `−median(Chronos)` across each target's finite values in the represented breast models. Spearman ρ = **0.005332**; 95% target-bootstrap percentile CI **[-0.064021, 0.074643]** (B=10000, seed=20260929); asymptotic two-sided p = **0.876725**; target-label permutation empirical two-sided p = **0.872913** (B=10000, seed=20260930).

## Prespecified sensitivities

- Degree baseline: Spearman( baseline degree, `−median(Chronos)` ) = 0.163044, p=1.79084e-06.
- Degree-adjusted partial Spearman (rank residualization): ρpartial=-0.038939, p=0.257344.
- Mean-based outcome on the same targets: Spearman ρ=0.014098, p=0.681659.
- ≥20 finite-model sensitivity: N=849, Spearman ρ=0.005332, p=0.876725. This matches the primary N because none of the 858 mapped targets had 10–19 finite measurements (nine had fewer than 10).
- Within-degree-stratum null: observed ρ=0.005332; permutation distribution mean=0.052085, median=0.052014, central 95% interval [-0.011128, 0.114764] (B=10000, seed=20260931; 10 fixed strata). The recorded absolute-value exceedance fraction is 0.964204; because this conditional null is shifted positive, that fraction is descriptive and is not a conventional conditional test of added value.

## Interpretation boundary

The primary estimate is near zero (ρ=0.0053), with a target-bootstrap interval spanning small negative and positive associations and a high target-label permutation p-value. Under this prespecified design, these results provide no detectable association between GenX structural impact and the aggregated DepMap dependency outcome; they do not establish that the true relationship is exactly zero. No arbitrary pass threshold is used. This analysis does not establish causal STRING signaling, equivalence of graph fragmentation to essentiality, perfect essentiality prediction, mechanistic viability prediction, patient-response prediction, or clinical utility. The degree baseline is positively associated with dependency (ρ=0.1630), whereas the partial rank association controlling degree is small and negative (ρpartial=−0.0389; p=0.2573). The within-degree-stratum permutation distribution is shifted positive (mean 0.0521; median 0.0520; central 95% interval [−0.0111, 0.1148]) and the observed ρ=0.0053 lies within it, below its median. The reported absolute-value exceedance fraction (0.9642) compares null |ρ| to |observed ρ| around zero; because the conditional null is not centered at zero, do not interpret this as a conventional conditional p-value for GenX information beyond degree. The full 10,000 draws are retained in `depmap_degree_matched_null.csv` and summarized in `depmap_degree_matched_null_summary.csv`. Degree and degree-aware results are comparisons/sensitivities, not winners. Targets are the analysis units as prespecified, but genes can be biologically dependent; the target bootstrap/permutation does not make them independent experimental replicates.

## Reproducibility

Analysis manifest frozen before association calculation: `D:\GenX\validation\depmap_validation_analysis_manifest.json`. Post-run checksums and environment are in `D:\GenX\validation\depmap_validation_run_record.json`. Full target table: `D:\GenX\validation\depmap_validation_target_table.csv`. Inputs and graph were SHA-256 verified against the frozen manifest. Raw DepMap files and GenX source/graph/specification were not modified. No L1000 data were used.
