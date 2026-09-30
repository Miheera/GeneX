# GenX validation data-availability audit

## Purpose

This directory audits the availability and identifier coverage of independent perturbation data for the frozen GenX v2 model. It does not evaluate GenX performance and does not alter GenX.

## Frozen inputs

- Graph: `../backend/data/processed/expanded_graph_t900.csv` (881 nodes, 9,681 edges).
- SHA-256: `D9E4DC10A4A91154662D1A3C4EF047ACB7145658EBE76D6344B07C15BC6BBF91`.
- Frozen model and algorithm files were not edited by this audit. The repository already had unrelated modified/untracked files before this task; this audit added files only under `validation/`.

## Datasets inspected

- **Primary LINCS release:** GEO GSE92742. Retrieved only small official metadata/annotation files: `sig_info`, `sig_metrics`, `gene_info`, `gene_info_delta_landmark`, `pert_info`, `cell_info`, release README and GEO SHA512 manifest. The Level 5 signature matrix (`GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz`, listed as 20 GB) was not downloaded.
- **Earlier LINCS release:** GEO GSE70138. Retrieved only `sig_info`, `sig_metrics`, `gene_info`, `pert_info`, `cell_info` metadata. It was audited separately and not pooled with GSE92742.
- **Gene identifiers:** HGNC complete-set TSV from the official HGNC public download. Exact approved symbol was preferred, then unique alias/previous symbol; no fuzzy matching.
- **DepMap:** current release observed on the official portal was Public 26Q1. The official bulk-file catalog endpoint returned a human-verification page, recorded as `metadata/depmap_catalog_verification_response.html`. No DepMap matrix was obtained and no DepMap coverage counts are reported.
- **Project Score:** official Sanger documentation was reviewed; direct API calls failed in this environment. No Project Score matrix was substituted for DepMap.

Access date: 2026-09-29 local time. SHA-256 values and exact input paths/URLs are in `metadata/input_manifest.csv`.

## Eligibility and counting rules

- The frozen graph symbols are mapped to HGNC by exact approved symbol, otherwise exact unique alias/previous symbol. Ambiguous/unmapped symbols are not guessed.
- L1000 landmark identity is taken from GSE92742 `gene_info` (`pr_is_lm=1`) and matched to the graph through Entrez Gene ID. The directly measured landmark response universe is therefore the same for each Level 5 profile; profile-specific numeric value completeness was not checked because the large matrix was not downloaded.
- GSE92742 MCF7 profile rows are retained individually. `trt_sh`, `trt_sh.cgs`, `trt_sh.css`, and `trt_oe` remain distinct. No profiles, doses, times or shRNA reagents are silently collapsed.
- `distil_nsample` and quality fields are reported as supplied. No quality threshold, minimum response-gene threshold, or outcome-based target selection was applied.
- Same-cell/time control signature presence is reported, but `sig_info` does not identify an exact control-signature ID matched to every perturbation profile. Level 5 profiles are processed signatures; exact pairing is left as a limitation rather than inferred.
- DepMap eligibility is not evaluated until official pinned-release files are available. Blank values in its coverage CSV mean unavailable/not assessed, not zero.

## Generated files

- `audit/l1000_genx_landmark_overlap.csv`
- `audit/l1000_target_coverage.csv`
- `audit/depmap_genx_coverage.csv`
- `audit/gene_mapping_audit.csv`
- `audit/validation_provenance_audit.csv`
- `validation_data_availability_report.md`
- `metadata/input_manifest.csv`
- `scripts/audit_l1000.py`
- `scripts/write_manifest.py`

## Reproduce

From the repository root, after the listed metadata files are present:

```powershell
python "Major Project\validation\scripts\audit_l1000.py"
python "Major Project\validation\scripts\write_manifest.py"
```

The L1000 script uses only Python's standard library. It reads the frozen graph by reference; it does not copy or rewrite it.

## Deliberately not done

- No expression response values were loaded or correlated with GenX.
- No dependency values were loaded, summarized or compared.
- No p-values, AUC metrics, ranking-performance statistics or validation claims were calculated.
- No GenX graph, edge weight, path cost, knockout implementation, RWR equation, API or dashboard was changed.
- No target was selected or excluded based on GenX performance.

## Next step

To complete the knockout coverage audit, obtain the official DepMap Public 26Q1 `CRISPRGeneEffect.csv` and `Model.csv` via the official current-release downloads page after its human-verification step. To verify numerical per-profile L1000 coverage, obtain the official Level 5 signature matrix or an official profile-filtered export. Freeze the response-gene minimum, profile selection across times/doses, treatment of raw shRNA versus `trt_sh.cgs`, and control/quality rules before inspecting any GenX-versus-response result.
