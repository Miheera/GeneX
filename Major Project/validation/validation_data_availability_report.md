# GenX v2 validation data-availability report

**Audit date:** 2026-09-29 (local).  
**Scope:** availability, coverage and provenance only. No GenX performance analysis was run.

**Follow-up:** The separate first DepMap knockout validation was run after the availability audit was frozen. Its prespecified population, analysis outputs, and interpretation are documented in [the validation summary](</D:/GenX/validation/results/depmap_validation_summary.md>). The availability counts below remain the data-availability record; the L1000 numeric matrix remains unavailable locally and was not used.

## Summary

The planned MCF7 L1000 experiment has substantial metadata-level target coverage. In the pinned GSE92742 metadata, 463 distinct frozen GenX genes have at least one MCF7 genetic-perturbation profile: 455 have an MCF7 shRNA profile, 246 have cDNA overexpression, and 238 occur in both classes (using individual `trt_sh` target coverage for the latter overlap). Each profile has 161 directly measured landmark genes that overlap GenX. All mapped target genes are in the `100+` response-coverage bucket. The profile-level table preserves all records, rather than collapsing multiple profiles.

This establishes availability in the release metadata, not numeric completeness of any signature cell. The official Level 5 matrix exists at GEO but was **not downloaded locally** because its compressed size is approximately **21.33 GB** and would require substantial storage and processing solely for this availability audit. No other assay level, release, or third-party copy was substituted. L1000 numeric availability is therefore **not established**, and L1000 biological validation is deferred. Same-cell/time control signatures and profile metrics are present, but exact control IDs matched to each profile cannot be reconstructed from `sig_info` alone.

The DepMap Public 26Q1 files supplied by the user were audited locally. `Model.csv` contains 2,154 unique model IDs; `CRISPRGeneEffect.csv` contains 1,208 unique model rows. Under the metadata-only cohort rule (`OncotreeLineage = Breast` and `OncotreePrimaryDisease = Invasive Breast Carcinoma`), 83 models qualify and 50 have a matrix row. Of 881 GenX nodes, 879 have a unique HGNC/Entrez crosswalk; 858 are represented by unique gene columns in the CRISPR matrix, and all 858 have at least one finite numeric value among the 50 breast models with dependency rows. These are availability counts only. The downloaded CSVs do not contain a release label internally; 26Q1 provenance is the user's official-download provenance, corroborated by the 26Q1 release records and the 1,208-model matrix row count. No older or unofficial matrix was substituted.

## Frozen graph provenance

The graph referenced is `backend/data/processed/expanded_graph_t900.csv`: 881 nodes and 9,681 edges. SHA-256: `D9E4DC10A4A91154662D1A3C4EF047ACB7145658EBE76D6344B07C15BC6BBF91`. The graph was not copied or modified.

## L1000 release and files

Primary source is [GEO GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742), the expanded LINCS L1000 release associated with the original L1000/CMap paper. Official supplemental directory: [GSE92742 supplemental files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/). Profile layer selected for the planned validation: `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz`, a replicate-collapsed Level 5 consensus signature matrix. The official source reports a compressed content length of 21,328,033,748 bytes (about 21.33 GB). It was not downloaded, per the storage/processing scope for this audit. Metadata availability is established; numeric feature availability is not established. Do not interpret metadata/profile counts as proof of numeric completeness.

Downloaded metadata files are listed and SHA-256 hashed in `metadata/input_manifest.csv`: `sig_info`, `sig_metrics`, `gene_info`, `gene_info_delta_landmark`, `pert_info`, `cell_info`, the release README and SHA512 manifest. LINCS signature metadata and quality metrics were used for availability counts; no signature response values were loaded.

The alternative [GEO GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) was inspected separately. Its MCF7 profile metadata contain 13,476 profile rows: 12,737 `trt_cp` and 739 `ctl_vehicle`, with no MCF7 genetic knockdown or overexpression profiles. It is therefore not a source for the planned MCF7 genetic perturbation test, and was not combined with GSE92742.

## L1000 coverage results

Counts below are from actual GSE92742 metadata and the frozen graph. Target counts refer to unique GenX genes; profile counts refer to distinct Level 5 `sig_id` records.

| Category | Count |
|---|---:|
| Frozen GenX nodes | 881 |
| Nodes mapped uniquely to HGNC | 879 |
| HGNC-unmapped nodes | 2 |
| L1000 directly measured landmark genes | 978 |
| Directly measured landmarks mapping to GenX | 161 |
| Non-landmark Level 5 feature IDs mapping to GenX (inferred features; not primary) | 586 |
| All-cell-line GenX genes in genetic perturbagen catalog (GSE92742 `pert_info`) | 538 |
| All-cell-line GenX knockdown targets (`trt_sh`/`trt_sh.cgs`) | 513 |
| All-cell-line GenX cDNA overexpression targets (`trt_oe`) | 318 |
| All-cell-line targets represented in both classes | 293 |
| MCF7 genetic signature profile rows, all targets/types | 29,639 |
| MCF7 profile rows mapping to GenX targets | 4,740 |
| MCF7 GenX individual shRNA targets (`trt_sh`) | 455 |
| MCF7 GenX consensus shRNA targets (`trt_sh.cgs`) | 453 |
| MCF7 GenX `trt_sh.css` targets | 0 |
| MCF7 GenX cDNA overexpression targets (`trt_oe`) | 246 |
| MCF7 target genes with individual shRNA and cDNA overexpression (`trt_sh` ∩ `trt_oe`) | 238 |
| MCF7 target genes with consensus shRNA and cDNA overexpression (`trt_sh.cgs` ∩ `trt_oe`) | 237 |
| MCF7 targets with same-cell/time control signatures of an intervention-appropriate type | 455 knockdown; 246 overexpression |
| Mapped MCF7 target profiles with quality metadata | 4,740 / 4,740 |
| Mapped MCF7 target profiles with `distil_nsample` ≥2 | 4,718 / 4,740 |
| Mapped MCF7 target profiles with non-sentinel `distil_cc_q75` | 4,718 / 4,740 |

**Landmark response-gene coverage per MCF7 target:** the Level 5 feature universe contains 161 directly measured landmark genes mapped to the frozen graph. Coverage bins across the 463 target genes are: 0–9: **0**; 10–19: **0**; 20–49: **0**; 50–99: **0**; 100+: **463**. This is feature-universe coverage from metadata, not a check that every matrix value is finite. The profile table reports the count on every profile row.

All 29,639 genetic MCF7 profile rows have a same-cell/time control signature in the metadata; all also have a same-cell/time control *type* consistent with the intervention category. This is contextual control availability, not a profile-to-control ID link. In the profile quality fields, 29,639/29,639 have metric rows; `distil_cc_q75` is non-sentinel for 29,479/29,639. `distil_nsample` ranges from 1 to 20, with 160 profiles at 1, 7,955 at 2, and 19,248 at 3. Every GenX target gene has at least one mapped profile with at least two replicates and non-sentinel `distil_cc_q75`; no quality or replicate cutoff was applied to eligibility.

**Identifier mapping:** GenX graph identifiers are gene symbols. Mapping uses the downloaded HGNC complete set: exact approved symbol first, otherwise exact unique alias/previous symbol. The two unmapped graph symbols are `LOC102723407` and `RAD51D-2`. Four mapped symbols resolve through a unique previous symbol (`CENPJ`, `H4-16`, `NCL`, `SLC9A3R1`). Landmark features are linked by the L1000 Entrez ID `pr_gene_id` to HGNC Entrez ID, not by fuzzy symbol matching. The complete 881-row mapping and landmark overlap tables are provided.

### Exclusions / unresolved eligibility

- MCF7 genetic profiles whose target does not map to one frozen GenX node remain in `l1000_target_coverage.csv` with an explicit mapping status and exclusion reason; they are not counted in GenX target totals.
- The two unmapped GenX graph symbols cannot be tested as targets through the frozen symbol-to-HGNC crosswalk without a documented resolution.
- The 161 directly measured response-node overlap is constant across targets, so the provisional 0–9 through 100+ buckets do not distinguish targets in this release.
- The exact matched control profile IDs are not available in the downloaded signature metadata. Level 5 is an already processed consensus signature; the audit cannot reconstruct each original profile/control pairing.
- Numeric matrix value completeness and whether every signature has finite data for all 161 nodes remain unchecked because the Level 5 matrix was not downloaded.
- No minimum response-gene threshold was applied. The report provides the requested distribution for deciding one before performance is examined.

## DepMap Public 26Q1 local-file audit

### Release and input verification

The audit used only the two files the user placed at the repository root: `D:/GenX/Model.csv` and `D:/GenX/CRISPRGeneEffect.csv`. Their sizes and SHA-256 hashes are in `D:/GenX/validation/audit/depmap_input_manifest.csv`. The supplied filenames and CSV headers themselves contain no release tag, so they cannot independently prove the quarter. Their provenance is the user's record that these are official Public 26Q1 downloads; this is consistent with the official [26Q1 release announcement](https://forum.depmap.org/t/announcing-the-26q1-release/4606) and the 1,208 model IDs in the CRISPR matrix, a count specifically referenced for 26Q1 in this [DepMap forum exchange](https://forum.depmap.org/t/how-should-crisprinferredmodelefficacy-csv-identifiers-be-mapped-in-depmap-public-26q1/4653). No older release or mirror was used. The 26Q1 announcement documents new CRISPR screens and a Chronos library-correction update, while noting that the released file format is unchanged.

- `Model.csv`: 2,154 rows; 2,154 nonblank unique `ModelID`s.
- `CRISPRGeneEffect.csv`: 18,532 columns including the first model-ID column (18,531 gene columns); 1,208 rows and 1,208 unique model IDs; no duplicate model IDs.
- Gene column identifiers use `symbol (Entrez ID)`. HGNC's exact approved symbol is matched first; otherwise the exact alias/previous symbol must resolve uniquely. HGNC Entrez ID is then matched exactly to the parenthesized DepMap column identifier. No fuzzy matching is used.
- The files are readable and matrix rows have the expected width. No raw file was modified. Input hashes and exact byte sizes are listed in the manifest.

### Frozen breast-model cohort

The deterministic inclusion rule uses metadata only: include a model iff `OncotreeLineage` is exactly `Breast` and `OncotreePrimaryDisease` is exactly `Invasive Breast Carcinoma` (case-insensitive exact text match). The DepMap documentation defines these as lineage and primary disease fields in `Model.csv` ([DepMap metadata documentation](https://forum.depmap.org/t/depmap-metadata/3694)). This includes 83 of 2,154 models. There are 96 models in the broader Breast lineage; the other 13 are excluded by the primary-disease rule (in situ, NOS, or non-cancerous categories). The cohort was set without consulting dependency values. All 2,154 models, including cohort exclusions with reason, are retained in `depmap_breast_model_metadata.csv`.

### GenX gene coverage and numeric completeness

| Quantity | Result |
|---|---:|
| GenX nodes | 881 |
| GenX genes with unique HGNC/Entrez mapping | 879 |
| GenX genes present as unique columns in DepMap matrix | 858 |
| GenX genes absent from DepMap matrix | 21 |
| GenX genes without an unambiguous HGNC/Entrez mapping | 2 |
| GenX genes with at least one numeric dependency value in the eligible breast cohort | 858 |
| Total DepMap models in `Model.csv` | 2,154 |
| Breast-cancer models under the frozen rule | 83 |
| Breast models represented in `CRISPRGeneEffect.csv` | 50 |
| GenX × breast cells (881 × 83) | 73,123 |
| Numeric cells among mapped GenX columns and breast models | 42,281 |
| Missing cells across the full 881 × 83 grid | 30,842 |

The 42,281 numeric cells are 42,281/71,214 expected cells for the 858 genes represented in the matrix across all 83 metadata-eligible models (59.37%). Of the 83 breast models, 33 have no row in `CRISPRGeneEffect.csv`; among the 50 represented models, 619 mapped GenX cells are missing. There were zero invalid/nonfinite/nonnumeric values and zero duplicate matrix model IDs or duplicate Entrez gene columns in the audited cells. The 21 genes without a matrix column and the two unmapped graph symbols remain unavailable; they are not silently filled or estimated. Exact per-model counts are in `depmap_numeric_completeness.csv`; per-gene resolution is in `depmap_genx_gene_mapping.csv`.

### Chronos score interpretation (not transformed)

`CRISPRGeneEffect.csv` contains Chronos gene-effect estimates from CRISPR knockout screens. Chronos describes these as relative change in growth rate after successful knockout: zero indicates no change, negative indicates loss of viability/growth, and positive indicates a growth gain ([Chronos official README](https://github.com/broadinstitute/chronos/blob/master/README.md)). DepMap's established scaling places nonessential genes around 0 and common essential genes around -1; lower/more negative values therefore indicate stronger dependency. This is a growth/fitness phenotype, not a direct measurement of GenX structural fragmentation. The future planned orientation, if the protocol adopts positive-is-stronger, would be `-raw_gene_effect`; **no transformation has been applied in this availability audit**. The 26Q1 release announcement says the revised library correction affects `CRISPRGeneEffect.csv` but not its file format.

### Future knockout-validation population feasibility

Using only the frozen mapping/cohort/data-availability conditions, **858 GenX targets** meet the potential eligibility rule: a unique HGNC/Entrez mapping, a unique gene-effect matrix column, and at least one finite numeric value among the 50 breast-cancer models that have matrix rows. This is population feasibility, not evidence that `|ΔE|` predicts dependency. The available CRISPR measure is a knockout growth/fitness response; whether structural impact is expected to associate with it is a separate validation hypothesis. The full 83-model metadata cohort is preserved, while 33 models cannot contribute a matrix value because they are absent from this CRISPR matrix.

## Provenance / leakage

STRING v12 integrates experiments, curated databases, co-expression and text-mining evidence; STRING confidence is a combined association score. [STRING score documentation](https://string-db.org/help/scores/) The graph file contains only combined edge scores, so a source-publication-to-edge audit cannot be performed from the shipped graph alone.

| Dataset | Release/publication timing | Leakage determination |
|---|---|---|
| LINCS GSE92742 | GEO supplemental files dated 2017; L1000 paper 2017 | Potential overlap with literature/evidence channels; **not determinable from available metadata**. |
| DepMap Public 26Q1 | Release announcement 2026; underlying screens/publications can predate release | Potential historical evidence overlap; **not determinable from available metadata**. Local matrix metadata does not expose screen-publication provenance by gene. |
| Project Score | Paper issue 2021; underlying screens/publications may be earlier | Potential historical evidence overlap; **not determinable from available metadata**. No current data matrix inspected. |

No claim of “no leakage” is made. Full provenance rows are in `audit/validation_provenance_audit.csv`.

## Answers to the requested questions

### 1. How many GenX genes can be tested with MCF7 L1000 genetic perturbation?

Metadata-level population: **455 knockdown targets** (`trt_sh`; **453** also have gene-level consensus `trt_sh.cgs`), **246 overexpression targets**, and **238 targets with both individual shRNA and overexpression** (**237** for `trt_sh.cgs` plus overexpression). Each has 161 directly measured GenX landmark response genes in the Level 5 feature universe. All have at least one profile with ≥2 replicates and non-sentinel `distil_cc_q75`; exact per-profile control IDs and numeric matrix completeness remain unverified.

### 2. How many GenX genes can potentially be tested against breast-cancer DepMap dependency data?

**858 GenX targets** have an unambiguous HGNC/Entrez mapping to a unique `CRISPRGeneEffect.csv` column and at least one finite numeric dependency value among the 50 eligible breast models with matrix rows. The cohort itself contains 83 models under the predefined metadata rule; 33 do not have CRISPRGeneEffect rows. This establishes data availability only, not a relationship between structural knockout impact and dependency.

### 3. Is the planned minimum viable experiment feasible?

- **RWR/L1000 availability:** Metadata availability is established (455 individual shRNA, 246 overexpression, and 238 shared targets; 161 GenX landmark features per candidate target in the feature universe). Numeric Level 5 availability is **not established** because the official approximately 21.33 GB matrix was not downloaded. L1000 biological validation is deferred.
- **Virtual knockout/DepMap availability:** Yes, the local Public 26Q1 files provide a usable availability population: 83 breast-cancer models by metadata; 50 have dependency rows; 858 GenX targets have numeric values in at least one such model. This makes a future dependency comparison data-feasible, but does not establish that Chronos dependency is a valid surrogate for GenX structural impact or that an association exists.

### 4. Current status

No target selection based on outcomes, dependency ranking, or GenX performance has occurred. The DepMap branch is ready for a separately preregistered/frozen validation protocol; L1000 remains a future option pending numeric matrix access and a prespecified profile-selection plan.

## Work not performed / final checks

- No correlation, AUROC, AUPRC, NDCG, p-value, bootstrap or GenX performance statistic was calculated.
- No target was excluded based on a GenX outcome; no outcomes were inspected.
- At the time of this availability audit, no validation statistics had been calculated. The subsequent DepMap validation is reported separately in `D:/GenX/validation/results/depmap_validation_summary.md`. No frozen GenX implementation files, graph, or methodology specification were edited. Availability-audit outputs are in `D:/GenX/validation/audit/`; the later validation artifacts are in `D:/GenX/validation/`.
- Existing GenX tests were requested and are reported in the completion message. No tests were added.
- DepMap reproduction steps, cohort and mapping rules, file hashes, and outputs are documented in `D:/GenX/validation/audit/README.md` and `depmap_input_manifest.csv`. Existing LINCS metadata reproduction details remain in `Major Project/validation/README.md`, `scripts/`, and `metadata/input_manifest.csv`.
