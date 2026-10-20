# GenX

> GenX is an **Interactive Exploratory Network Analysis & Computational Hypothesis Generation Workbench**. Its outputs prioritize candidate target combinations for downstream wet-lab validation, rather than serving as clinical diagnostic predictors.

GenX is a research prototype for exploring a human gene-association network and examining how network structure changes under virtual perturbations. Its current shipped network is breast-cancer-focused. The application is not a clinical tool and its scores are not experimentally calibrated biological response predictions.

## Run GenX

Requirements: Python 3.10 or newer.

```powershell
python -m pip install -r "backend\requirements.txt"
.\run_genx.ps1
```

The launcher starts the backend and opens <http://localhost:8000>. The interactive API documentation is at <http://localhost:8000/docs>.

## Project layout

- `frontend/` contains the browser entry point, D3 application, and styles. FastAPI serves only this directory as the web UI.
- `backend/` contains the FastAPI app and core graph modules, plus focused `analysis/`, `pipeline/`, `tests/`, and `data/` directories.
- `validation/` contains DepMap validation; `validation/lincs/` contains the separate LINCS data-availability audit.
- `research/` contains research papers; `docs/` contains project documentation and reports.

## Accounts

GenX has JWT-based sign-in and account creation. Accounts are associated with an institution name and an institution-issued email address; common personal email providers are rejected. Passwords are stored as salted PBKDF2 hashes, and the JWT is held in an HttpOnly session cookie. Local account data and the development signing key are stored in the ignored `backend/data/private/` directory. Set `GENX_AUTH_DATA_DIR` to use a different private directory.

This prototype does not yet verify email ownership or institutional affiliation, send password resets, or provide multi-factor authentication. A production deployment needs email verification, HTTPS (`GENX_COOKIE_SECURE=true`), a managed `GENX_JWT_SECRET` of at least 32 bytes, rate limiting, account recovery, and an institution/domain verification policy.

## What is implemented

- Interactive D3 view of the processed human STRING association network.
- **Virtual topological knockout:** removes one gene from a graph copy and reports structural metrics such as weighted efficiency, betweenness changes, and connected components.
- **Fixed-graph propagation:** inhibition and overexpression use differential Random Walk with Restart over the expression-weighted graph, relative to a uniform-seed baseline. Scores are model-derived signal shifts, not measured biological responses.
- Attack-tolerance, STRING-threshold sensitivity, COSMIC comparison, structural-expression comparison, and DGIdb interaction-claim annotations.

## Positioning

| Tool | Primary strength | GenX relationship |
|---|---|---|
| Cytoscape | Broad desktop ecosystem for network analysis and highly configurable visual styling | GenX is a focused web workbench that keeps an interactive breast-cancer graph connected to its fixed-network differential RWR and virtual knockout calculations; it is not a replacement for Cytoscape's analysis ecosystem. |
| Open Targets Platform | Large-scale disease-target evidence aggregation and target prioritization | GenX presents local network topology, perturbation scenarios, and imported annotation evidence; it does not reproduce Open Targets' evidence platform or disease-association scoring. |
| GEARS / PDGrapher | Learned models for predicting perturbation-associated expression or phenotype changes | GenX does not implement a GNN or claim equivalent prediction. Its graph view can serve as an explainable exploratory interface around network-derived hypotheses and, later, externally validated model outputs. |

GenX's niche is the integrated, web-native exploratory interface for fixed-network RWR, structural knockout analysis, annotation context, and subnetwork inspection. Its outputs remain computational hypotheses for downstream experiments.

The current graph contains 881 nodes and 9,681 edges from STRING v12 human data at a combined-score threshold of 900. Method details and interpretation limits are in [the v2 methodology specification](docs/GENX_V2_METHODOLOGY_SPEC.md).

For a deterministic raw-source edge rebuild, set `GRAPH_REBUILD_FROM_RAW=1` before startup. GenX maps STRING protein IDs through the v12 protein-info table, restricts both endpoints to the frozen 881-node `expanded_nodes_t900.csv` universe, keeps the maximum score for duplicate gene pairs, applies the 900 confidence threshold recorded in `seed_genes.json`, and reapplies expression weighting. The checked-in primary seed list is included as graph provenance; the historical 881-node universe remains explicitly frozen in the nodes artifact rather than being re-expanded from those 13 genes. The bundled source files rebuild to 881 nodes and 9,681 edges. Default startup uses the processed CSV graph to avoid parsing the full STRING links archive on every launch.

## Current scope and limitations

- The cancer selector retains future-scope options, but the current network and analysis outputs are breast-cancer-specific. Only the breast network is served by the API.
- The dashboard does not provide a validated disease-association score; that value is shown as unavailable.
- The DGIdb view shows interaction claims alongside a structural-impact ranking. It does not use the Open Targets Platform and does not predict drug response or efficacy.
- Knockout and propagation analyses are network calculations, not experimental perturbations. In the current DepMap evaluation, the knockout structural-impact score did not show a detectable association with dependency (N=849, Spearman rho=0.0053). See [the validation summary](validation/results/depmap_validation_summary.md).
- If the graph API is unavailable or an unsupported cancer is selected, the dashboard displays an error rather than substituting synthetic network data.

## Data and validation

Primary data artifacts and generated analysis outputs are under `backend/data/`. Pipeline scripts are under `backend/pipeline/`; source notes and provenance are in `backend/data/DATA_SOURCES.md` and `validation/lincs/validation_data_availability_report.md`.

The reproducible DepMap evaluation script and frozen inputs/outputs are under `validation/`. Its current cohort includes 50 breast models with dependency rows out of 83 metadata-eligible models. See the validation summary for the design, missingness, and interpretation.

The frozen cohort maps exact approved HGNC symbols (or exact unique aliases/previous symbols) through Entrez IDs to exact unique DepMap `SYMBOL (Entrez ID)` columns. Of 858 uniquely mapped genes, 849 meet the requirement for at least 10 finite values among the 50 represented breast cancer models. Exact model IDs and inclusion decisions are in `validation/depmap_validation_analysis_manifest.json` and `validation/frozen_depmap_validation_cohort.csv`; gene-level mapping decisions are in `validation/audit/depmap_genx_gene_mapping.csv`.

Run the exploratory comparison between sampled-source unweighted efficiency loss, tissue-weighted RWR flux, and the frozen DepMap dependency outcome from the repository root:

```powershell
python validate_depmap.py
```

This uses the frozen target-level DepMap table and the current geometric TCGA-expression-weighted graph. It reports Spearman rho, two-sided p-values, deterministic 95% bootstrap intervals (1,000 resamples), and Benjamini-Hochberg q-values across tested predictor associations. An optional `--ccle-expression path.csv` compares a CCLE breast-expression profile with columns `gene_symbol` and TPM/expression. No CCLE profile is currently included, so that comparison only runs when supplied. Run `python sensitivity_analysis.py` for STRING threshold, RWR alpha, source sample-size, and expression-transform ranking stability; reported rank R² is not predictive model R².

### Current DepMap results

| Predictor | N | Spearman ρ | p | 95% bootstrap CI | BH q |
|---|---:|---:|---:|---:|---:|
| Unweighted efficiency drop | 849 | +0.0146 | 0.6719 | [−0.0494, +0.0859] | 0.8767 |
| Geometric tissue-weighted RWR L₁ flux | 849 | +0.1706 | 5.73 × 10⁻⁷ | [+0.1017, +0.2360] | 1.72 × 10⁻⁶ |

The outcomes use negative median Chronos across 50 breast cancer cell lines represented in the matrix. These are retrospective cohort associations, not causal or patient-level predictions. The sensitivity report has 15/36 RWR configurations above rank R² = 0.85 against the reference; therefore the default geometric weighting, α = 0.85, and STRING c ≥ 900 setting is an operational reference regime, not a universally optimal setting. Use `python plot_depmap_correlation.py` to write the publication PNG and SVG under `validation/results/`. Full results, the tool positioning matrix, test status, and the 11-point methodological boundary are in [DOCUMENTATION.md](DOCUMENTATION.md).

### 11-point methodological boundary

1. STRING combined scores describe gene association evidence, not necessarily physical, directed, activating, or causal interactions.
2. The shipped graph and analyses are breast-cancer-focused; they do not establish performance in other diseases or tissues.
3. Expression-weighted edges contextualize the graph with a selected profile; they do not model regulation, protein abundance, or causal signaling.
4. RWR describes diffusion over a graph operator, not biological time, kinetic rates, or dose-response.
5. Differential RWR uses a uniform-seed computational baseline, not an experimentally measured untreated state.
6. Inhibition and overexpression apply score sign conventions on a fixed graph; they do not physically perturb expression.
7. Virtual knockout removes a node and measures graph structure; efficiency and betweenness changes are not cell death or essentiality measurements.
8. Combinatorial synergy is non-additivity of calculated network scores, not demonstrated synthetic lethality.
9. DepMap metrics are retrospective cohort associations across 50 represented breast models and 849 eligible targets, not prospective, causal, or patient-level predictions.
10. DGIdb, COSMIC, and other imported annotations are evidence records; they do not establish approval, target engagement, efficacy, or patient suitability.
11. GenX prioritizes hypotheses for downstream wet-lab work such as dual-CRISPR or Perturb-seq. It is not a deterministic ODE kinetic simulator, clinical diagnostic predictor, or treatment recommendation system.

Before redistributing GenX or using it in a commercial product, review the applicable terms and attribution requirements for each incorporated dataset and third-party library.
