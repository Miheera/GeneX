# GenX frozen-model validation data-availability audit

## Scope

This is a numeric data-availability and QC audit only. It does not calculate any GenX-vs-experiment statistic, compare targets, rank genes, or modify the frozen GenX implementation, graph, or methodology. The L1000 metadata audit is preserved separately under `D:/GenX/Major Project/validation`; the 21.33 GB official GSE92742 Level 5 matrix was not downloaded. Do not treat L1000 metadata counts as numeric completeness.

## Inputs

| File | Use | Release / provenance |
|---|---|---|
| `D:/GenX/Model.csv` | Model/cancer metadata and frozen breast cohort | User-provided official DepMap Public 26Q1 download |
| `D:/GenX/CRISPRGeneEffect.csv` | Model-by-gene Chronos effect numeric completeness only | User-provided official DepMap Public 26Q1 download |
| `D:/GenX/Major Project/backend/data/processed/expanded_nodes_t900.csv` | 881 frozen GenX node symbols | Frozen GenX source input |
| `D:/GenX/Major Project/validation/metadata/hgnc_complete_set.txt` | Exact symbol/alias to HGNC Entrez crosswalk | Existing HGNC crosswalk used by prior metadata audit |

The DepMap CSV headers do not embed a release identifier. Release provenance is the user-provided official 26Q1 download context. Local matrix dimensions (1,208 unique model rows) are consistent with the 26Q1 file context documented on the [official DepMap Community Forum](https://forum.depmap.org/t/how-should-crisprinferredmodelefficacy-csv-identifiers-be-mapped-in-depmap-public-26q1/4653); the [26Q1 release announcement](https://forum.depmap.org/t/announcing-the-26q1-release/4606) documents that release's CRISPR pipeline update. File byte sizes and SHA-256 checksums are in `depmap_input_manifest.csv`. If future provenance requires independent cryptographic confirmation from DepMap, obtain the official release manifest or download receipt; do not silently replace these files.

## Frozen cohort rule

Include a model exactly when both conditions hold in `Model.csv`:

1. `OncotreeLineage` equals `Breast` (case-insensitive exact match); and
2. `OncotreePrimaryDisease` equals `Invasive Breast Carcinoma` (case-insensitive exact match).

This is a deterministic primary-disease cancer cohort and does not rely on dependency values. All models are kept in `depmap_breast_model_metadata.csv`, with inclusion flag and reason. The official [DepMap metadata description](https://forum.depmap.org/t/depmap-metadata/3694) defines `OncotreeLineage`, `OncotreePrimaryDisease`, and related `Model.csv` fields.

## Identifier mapping

For each graph symbol, match exact HGNC approved symbol first. If none, use an exact HGNC `alias_symbol` or `prev_symbol` only when it resolves to one HGNC record. Take its HGNC Entrez ID and match that integer exactly to the matrix column identifier formatted as `SYMBOL (EntrezID)`. No fuzzy mapping is used. Missing HGNC mappings, absent matrix columns, duplicate Entrez columns, and multiple GenX nodes mapping to one Entrez ID are separately flagged in `depmap_genx_gene_mapping.csv`.

## Numeric completeness

Only the uniquely mapped GenX columns are scanned. Every cell in eligible breast model rows is classified as finite numeric, missing (blank/NA/NaN/null), or invalid (non-numeric or non-finite). Dependency values are never written to audit outputs. Missing whole model rows count as missing cells for the mapped matrix grid. The full 881-node by cohort grid is separately shown with unmatched GenX genes included in the missing denominator. A target is potentially eligible for future testing if it maps uniquely to a matrix gene column and at least one included breast model with a matrix row has a finite numeric value for that gene.

Chronos sign documentation: the [official Chronos README](https://github.com/broadinstitute/chronos/blob/master/README.md) says 0 indicates no growth change, negative indicates growth/viability loss, and positive indicates growth gain; DepMap's conventional normalization centers nonessential genes around 0 and common essential genes around -1. This audit records but does not transform any score. Future positive-is-stronger orientation, if selected in a frozen protocol, would negate the raw gene-effect value.

## Run

From PowerShell at `D:\GenX`:

```powershell
python validation\audit\run_depmap_audit.py
```

The script uses only the four input files above. It scans mapped cells to classify numeric availability; it exports counts and flags only, never effect values. It makes no statistical or validation-performance calculations.

## Outputs

- `depmap_breast_model_metadata.csv`: one row for every Model.csv model; annotations, inclusion flag, exclusion reason.
- `depmap_genx_gene_mapping.csv`: one row per frozen GenX node; HGNC/Entrez mapping, exact matrix match, status, duplicate flag.
- `depmap_numeric_completeness.csv`: one row per eligible breast model; presence/duplicate status and numeric/missing/invalid mapped-gene cell counts.
- `depmap_numeric_completeness_summary.csv`: requested cohort, mapping, matrix completeness, and future eligibility counts.
- `depmap_input_manifest.csv`: source paths, file sizes, release provenance note, SHA-256 checksums.
- `run_depmap_audit.py`: reproducible audit script.

The report at `D:/GenX/Major Project/validation/validation_data_availability_report.md` also preserves prior L1000 metadata counts and explicitly states that L1000 numeric Level 5 availability is not established and biological validation is deferred.

This availability audit did not calculate biological validation statistics. A subsequent first DepMap validation is recorded separately at `D:/GenX/validation/results/depmap_validation_summary.md`; its frozen protocol and run record are at `D:/GenX/validation/depmap_validation_analysis_manifest.json` and `D:/GenX/validation/depmap_validation_run_record.json`.
