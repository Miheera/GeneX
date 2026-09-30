# GenX codebase forensic methodology report

Scope: executable code and tracked data in `Major Project/backend` as present on 2026-09-26. This report describes behavior only. It does not assess the method or compare it with publications. Line numbers below refer to the files as inspected in this checkout.

## 1. Runtime graph and data provenance

### Runtime graph actually loaded

`backend/graph_loader.py:18-19,22-70` sets the default graph input to `data/processed/expanded_graph_t900.csv` and explicit node input to `data/processed/expanded_nodes_t900.csv`. `GRAPH_CSV_PATH` can replace the edge CSV; `GRAPH_NODES_CSV_PATH` can replace the node CSV. If the edge path is overridden and no node-path override is supplied, the loader does not load the default node CSV. It reads `gene1`, `gene2`, and `confidence_score`; blank endpoint rows are skipped. Each edge row adds an undirected edge to `nx.Graph`. `interaction_type` is not read. If the nodes CSV exists, `gene_symbol` values are added, including nodes without edges; otherwise node set is derived from edges. For each node the loader adds attributes `id`, `label`, `degree` (unweighted degree), and `disease_score` (a hard-coded lookup in `graph_loader.py:10-15`, default 0.5).

The checked-in default artifacts have **881 node rows and 9,681 edge rows**, 9,681 distinct unordered endpoint pairs, all 881 nodes incident to an edge; confidence scores range 900–999. These are file row counts, not a count obtained by starting the API. At runtime, NetworkX also collapses repeated endpoint pairs if they are present in a custom input CSV. `main.py:48-52` initializes this graph at application startup and calls `precompute_baseline`.

Nodes in the processed artifacts are gene symbols. The source-generation code maps STRING protein IDs to `preferred_name` values using `9606.protein.info.v12.0.txt` columns `#string_protein_id` and `preferred_name` (`pipeline/expand_via_string.py:78-101`). It excludes labels beginning `ENSP`, unmapped protein IDs, and self-symbol edges. It uppercases / alias-resolves symbols only in `pipeline/map_identifiers.py`; the STRING expansion function itself does not invoke that resolver. The loader uses symbols exactly as present in CSV other than trimming whitespace.

### STRING construction lineage

`pipeline/expand_via_string.py:33-38,209-282` declares source paths `data/processed/seed_genes.csv`, gzipped STRING v12.0 human protein info, and gzipped STRING v12.0 protein links. Its default combined-score threshold is 900. It first maps seed symbols to STRING protein IDs, scans the links for one-hop partners where `combined_score >= threshold`, then creates the node set `seed symbols ∪ mapped partners`. A second link scan emits edges only when both protein IDs are in that final set and the score meets the same threshold. This is an induced subgraph over the selected seed-plus-one-hop symbols (including qualifying partner-to-partner links), not an all-protein graph. It sorts endpoint symbols to deduplicate unordered pairs and omits same-symbol edges. Output edge fields are `gene1`, `gene2`, raw integer `confidence_score`, `interaction_type="string"`.

The relevant raw STRING fields actually read by expansion are `protein1`, `protein2`, and `combined_score` from the links file and `#string_protein_id`, `preferred_name` from protein info. No other STRING score column is read in this routine. The generic `filter_string_links` helper (`pipeline/filter_string.py:10-15`) likewise filters tuple element 3 against a caller-supplied threshold; it does not itself read a file. The default runtime artifact name `expanded_graph_t900.csv` corresponds to a 900 threshold by filename and its observed scores are >=900; the loader itself does not apply a threshold.

The intended upstream seed producer is `pipeline/filter_disgenet.py:24-28,83-130`: by default it reads DisGeNET XLSX exports for CUI `C0006142` and `C0678222`, selects the top 30 rows per query by descending GDA `score` (tie by symbol), unions the symbols, keeps the higher score for overlaps, and writes `seed_genes.csv`. The seed file itself is not listed among tracked files, so the exact seed set used for the checked-in expansion is not directly recorded in that input artifact. A separate root `major_project_dataset` directory has three placeholder 1-byte legacy paths and is not the default expansion path (`data/DATA_SOURCES.md`).

The separate `pipeline/sweep_thresholds.py:20-21,48-99,142-215` can construct graphs at thresholds 400, 500, 600, 700, 800, 900 (default). It reads the gene universe from `data/breast_cancer_graph.csv`, maps STRING IDs to symbols found in that universe, and reads `combined_score`. It uses `graph_builder.build_graph`; that builder filters empty/self edges and builds an undirected `nx.Graph` (`graph_builder.py:19-30`). The sweep's default shipping comparison threshold is 700; this is distinct from the default API-loaded `_t900` graph. Sweep summary outputs include node/edge counts, mean degree `2m/n`, simple-graph density, connected component count, connectivity flag, weighted betweenness top ten, and edge-set Jaccard comparisons (`sweep_thresholds.py:101-130,171-215`).

## 2. Edge weights and path distances

At graph build time, `pipeline/graph_builder.py:19-30` stores for a score `c`:

\[
w = c/1000,\qquad d = 1000/c\quad(c>0),
\]

where `weight=w`, `path_weight=d`. For `c<=0`, `path_weight=∞`. The serialization function writes `confidence_score = round(1000*weight)` (`graph_builder.py:33-46`).

The runtime loader independently computes `weight=int(confidence_score)/1000` and `path_weight=1/weight` for positive weight, otherwise infinity (`graph_loader.py:35-45`). Thus for positive integer confidence it is algebraically the same reciprocal mapping. There is no graph-wide normalization operation; division by 1000 rescales the score into a 0-to-1-like numeric weight, while the shortest-path cost is its reciprocal. API graph output exposes `weight` rounded to four decimal places (`main.py:130-137`); metrics use the unrounded stored `path_weight`.

Code comments/docstrings explain that the graph is built as weighted and the runtime metric is weighted global efficiency; the expansion comments describe filtering on combined score and the graph-builder docs say “weighted graph.” No explicit code comment supplies a biological or statistical rationale for taking a reciprocal. The direct implemented effect is that higher score yields smaller path distance.

## 3. Network metrics

### Weighted global efficiency

`perturbation.global_efficiency` (`backend/perturbation.py:26-52`) takes graph `G`, chooses supplied `sources` or `_sample_sources`, and for each source runs `nx.single_source_dijkstra_path_length(G, source, weight='path_weight')`. For graph size `n>=2`, it sums `1/d_G(i,j)` for reachable targets `j != i` with positive distance, then divides by `k(n-1)`, where `k` is the number of selected sources:

\[
E(G)=\frac{1}{k(n-1)}\sum_{i\in S}\sum_{j\ne i,\;d_G(i,j)<\infty}\frac{1}{d_G(i,j)}.
\]

For `n<2`, no selected sources, result is 0. This is weighted by `path_weight`; disconnected pairs contribute nothing because Dijkstra returns only reachable nodes. The denominator retains all `n-1` possible other nodes, including unreachable ones. When all nodes are selected, the code comment says the calculation is exact over directed ordered pairs; on the undirected graph, each pair is visited from both ends. For sampling, sources are drawn from the current graph but denominator uses its current node count.

The baseline value is stored in `_baseline_efficiency` during app startup; `global_efficiency` is also used after perturbation, for knockout-impact scores, and in attack-tolerance curves. No average shortest path length routine or metric was found in executable backend code.

### Shortest paths

The only shortest-path calculation for perturbation efficiency is single-source Dijkstra (`perturbation.py:45-47`) using `path_weight`, not `weight`. Thus paths minimize sums of reciprocal confidence weights. The efficiency calculation uses the reciprocal of those shortest distances. No path lists are returned to API consumers.

### Weighted betweenness

`perturbation.approximate_betweenness` (`perturbation.py:54-64`) calls `nx.betweenness_centrality(graph, k=min(100,n), weight='path_weight', normalized=True, seed=42)`. NetworkX computes node betweenness based on weighted shortest paths; output is normalized node centrality. Source sampling is controlled by the supplied seed and `k`. For `n<=1`, it returns zeros. The code does not calculate a custom betweenness formula itself. It saves baseline values, recalculates values for a perturbed graph, forms signed change `bc_prime[v]-bc_baseline[v]` for surviving baseline keys (excluding the perturbed gene), sorts by absolute change, and returns at most five (`perturbation.py:111-122`). `direction` is “gain” only if signed delta > 0, else “loss” (including exact zero).

Other betweenness call sites: threshold sweep calculates exact normalized weighted betweenness with `path_weight`, no `k` (`sweep_thresholds.py:115-120`); COSMIC benchmarking calculates normalized weighted betweenness, apparently without approximation (`analysis/benchmark_cosmic.py:79-88`). These results are used respectively for top-10 threshold report and rank-method comparison.

### Degree, connectivity, fragmentation, related metrics

The runtime graph's displayed/stored `degree` is unweighted neighbor count from `graph.degree(node)` (`graph_loader.py:61-66`), returned via graph endpoint (`main.py:118-127`) and used for frontend display. `nx.degree_centrality` is used to select the next node in degree-targeted attack tolerance (`analysis/attack_tolerance.py:49-52`); the node with highest value is removed, ties resolved by lexicographically smallest stringified node. For a simple graph, the library degree centrality is degree divided by `n-1`; code delegates the calculation to NetworkX.

Perturbation fragmentation is `len(nx.connected_components(G_prime)) > 1` (`perturbation.py:124-131`). `isolated_genes` is the sorted set of all surviving nodes outside the largest connected component, not merely degree-zero nodes. There is no separate fragment deletion: components remain in the graph. In the attack-tolerance analysis, `lcc_fraction = largest_component_size / remaining_node_count`, or 0 if none (`analysis/attack_tolerance.py:30-45`). This differs from perturbation's `isolated_genes` field.

Threshold sweep additionally computes mean degree `2m/n`, density delegated to `nx.density`, number of connected components and `is_fully_connected` (`pipeline/sweep_thresholds.py:101-112`). No average path length calculation was found.

## 4. Sampling and approximation

Constants in `perturbation.py:16-24` are `APPROXIMATION_SAMPLES=100`, seed 42. `_sample_sources` converts node IDs to strings, sorts them, then uses a fresh `random.Random(42).sample(nodes, min(100,n))`; selection is reproducible for a given node set and Python random implementation. Baseline efficiency uses these sampled source nodes, saved in `_baseline_sources`. Perturbed efficiency passes the *same original baseline source list* into `global_efficiency`; sources absent after knockout are discarded by membership filtering, but denominator is `len(selected_sources)*(n_perturbed-1)`. Consequently, after removing the sampled gene, fewer sources may be summed, while the denominator uses the reduced graph size and filtered source count.

Betweenness separately samples `k=min(100,n)` sources through NetworkX's `betweenness_centrality(..., seed=42)`, with weighted distance and normalization. It is not passed `_baseline_sources`; only seed and sample count are shared. Thus the efficiency and betweenness approximations are separate. Efficiency sampling is a source average; betweenness uses NetworkX's sampled-node approximation. If `n<=100`, all nodes are sampled by both pathways. The efficiency docstring/comment says that this matches the earlier exact directed-pair calculation for graphs of <=100 nodes. Sampling appears in the code/docstrings as approximation to reduce path calculations; this phrasing is descriptive, while exact rationale beyond runtime reduction is not separately encoded.

## 5. Perturbation models

Input `PerturbRequest` (`models.py:19-34`) allows `knockout`, `inhibition`, `overexpression`; default strength is 0.5 and validated inclusive range [0,1]. The frontend slider sends its value divided by 100 (`app.js:339-343`). `run_perturbation` deep-copies the graph first (`perturbation.py:85-88`).

**Knockout.** `G_prime.remove_node(gene_id)` (`perturbation.py:88-90`) removes the node and all incident edges (NetworkX graph operation). Remaining disconnected components are retained; efficiency skips unreachable pairs and fragmentation reports component membership.

**Inhibition.** For each edge incident to the selected node, code sets `factor=1-strength`, then `new_weight = weight*factor`; if `new_weight<=0`, stored `weight` becomes `1e-9`. `path_weight` is then `1/new_weight` (`perturbation.py:91-97`). For positive result this is exactly `w'=w(1-s)` and `d'=1/w'`. At `s=1`, displayed weight is floored to `1e-9`, not zero, so distance becomes `1e9`. There is no clipping of `s` in this function; API validation constrains it to [0,1], while direct callers can pass other values. Only incident edges change.

**Overexpression.** For each incident edge, code sets `factor=1+strength`; `new_weight=min(weight*factor,1.0)` and `path_weight=1/new_weight` (`perturbation.py:98-103`). This implements exactly `w'=min(w(1+s),1)` and reciprocal distance, provided weight is positive. It caps at 1; no lower bound is applied. API request validation constrains `s` to [0,1], but direct function calls are not clipped. Only incident edges change.

For invalid direct `perturbation_type`, no branch changes the deep-copied graph, but the function continues to calculate metrics and returns the unmodified graph result; API request typing rejects values outside its literal set.

## 6. Perturbation outputs and exact equations

At startup `precompute_baseline(G)` saves sampled sources, baseline efficiency `E0`, and baseline betweenness (`perturbation.py:66-71`). `run_perturbation` computes `E' = global_efficiency(G_prime, sources=_baseline_sources)`. It then computes:

\[
\Delta E\% = 100(E'-E_0)/E_0 \quad(E_0>0),\qquad 0\quad(E_0\le0),
\]

and `impact_score = round(abs(ΔE%), 2)` (`perturbation.py:105-109`). So both expressions asked about are implemented, with a zero-baseline branch and two-decimal rounding only for impact score; `efficiency_delta_pct` is not rounded in this function. The separate helper `knockout_impact_score` uses the same signed percentage formula with baseline truthiness, then returns absolute value rounded to two decimals (`perturbation.py:73-83`).

Betweenness change is absolute difference in normalized centrality values, signed: `delta_bc[v]=bc_prime[v]-bc_baseline[v]` for surviving nodes. It is not percentage change. The top-five list is ranked by `abs(delta_bc)` descending; tie ordering is not explicitly specified. The output carries the full baseline and perturbed centrality dictionaries as well as top-five list.

Fragmentation output uses connected components of perturbed graph; flag true exactly when component count exceeds one. If fragmented, `isolated_genes` contains nodes outside the largest component, sorted; if not fragmented, it is empty. Ties between largest components are resolved by iteration order through Python `max` over NetworkX component iteration.

`knockout_impact_score` is called by structural-expression analysis and used to produce per-gene `structural_score`; that analysis correlates it with expression/mutation/survival fields but does not rank results (`analysis/expression_correlation.py:74-80,108-180`). `impact_score` is the ranking metric for open-target output (`analysis/open_targets.py:41-84`) and one method in COSMIC benchmarking (`analysis/benchmark_cosmic.py:79-88`).

## 7. Ranking and generated interpretations

The perturbation endpoint does not rank genes globally. It ranks *other nodes' betweenness changes* within a single perturbation by absolute change, top five; direction is generated from the sign as above. Knockout impact helper itself returns an absolute score, no direction.

Open-target analysis computes knockout impact score for every graph node and sorts descending score, then gene symbol ascending for ties; truncates at caller-supplied `top_n` (API default 20, must be positive). It assigns `rank` starting at 1 and groups rows according to DGIdb annotation availability, not score threshold (`analysis/open_targets.py:41-84`; route `main.py:235-251`). COSMIC ranking computes three methods—betweenness, degree, impact score—and sorts each descending with stringified node tie-break (`analysis/benchmark_cosmic.py:79-136`). Degree is raw degree; impact is absolute knockout efficiency percent change; betweenness is weighted normalized centrality. No thresholds are applied to these rankings apart from chosen top-N. Threshold filtering applies only to graph generation.

Attack tolerance repeatedly removes a node chosen by highest recomputed degree centrality or weighted approximate betweenness; random mode shuffles node order per trial (50 trials in API, `main.py:186-197`; default seed argument is `None`). It outputs mean efficiency and LCC retention by removed fraction; random mode outputs population standard deviation, targeted modes have `None` standard deviation (`analysis/attack_tolerance.py:30-128`).

Labels generated in the perturbation response are `gain` / `loss`, `network_fragmented`, `isolated_genes`. Node `disease_score` values are hard-coded (`graph_loader.py:10-15`), with default 0.5; code does not derive that field from the graph or calculate it as a network metric. Frontend copy and mock results exist, but mock values are not backend calculations; the frontend requests the backend endpoint and uses a mock on request failure (`app.js:339-363,906-926`).

## 8. End-to-end calculation pipeline

`seed gene symbols + STRING v12 protein-info mapping + STRING v12 protein-links (protein1, protein2, combined_score) → one-hop seed expansion at combined_score >=900 → induced edge CSV and node CSV → API loads undirected NetworkX graph → confidence_score/1000 edge weight and reciprocal path_weight → deterministic sampled-source baseline efficiency and sampled weighted normalized betweenness → copy graph and remove target node or rescale incident edge weights → recalculate efficiency and betweenness → signed percentage efficiency difference, absolute impact score, signed betweenness differences, and component summary → return JSON; separate analyses sort per-gene knockout impact / centrality scores into rankings and CSV/API outputs.`

For the checked-in runtime graph, the beginning of this pipeline is represented by processed artifacts; the actual historical seed set is not available in the tracked paths inspected, so generation cannot be replayed exactly from this checkout without that input.

## 9. Evidence index

| Calculation / operation | File and function | Main lines | Variables / expression |
|---|---|---:|---|
| DisGeNET seed extraction | `backend/pipeline/filter_disgenet.py`, `_read_xlsx`, `_top_n`, `filter_disgenet` | 24–28, 35–91, 94–130 | top 30 per CUI by `score`; union symbols; `seed_genes.csv` |
| STRING source mapping and expansion | `backend/pipeline/expand_via_string.py`, `load_full_string_mapping`, `find_one_hop_partners`, `extract_induced_subgraph`, `expand_via_string` | 33–38, 78–194, 209–282 | `protein1`, `protein2`, `combined_score`, `threshold`, `seed_pids`, `all_syms`, `edges` |
| Undirected graph, builder weights | `backend/pipeline/graph_builder.py`, `build_graph` | 19–30 | `weight=score/1000`, `path_weight=1000/score` |
| Runtime graph CSV loading | `backend/graph_loader.py`, `load_graph` | 18–68 | `gene1`, `gene2`, `confidence_score`, `weight`, `path_weight` |
| API baseline setup and output | `backend/main.py`, `lifespan`, `get_graph`, `perturb_gene`, `health_check` | 48–52, 100–137, 152–181, 276–285 | graph globals and response fields |
| Source selection and efficiency | `backend/perturbation.py`, `_sample_sources`, `global_efficiency` | 16–52 | `APPROXIMATION_SAMPLES`, `APPROXIMATION_SEED`, `selected_sources`, `dist`, sum reciprocal distance |
| Betweenness | `backend/perturbation.py`, `approximate_betweenness` | 54–64 | `k`, `weight='path_weight'`, `normalized=True`, `seed=42` |
| Baseline cache | `backend/perturbation.py`, `precompute_baseline` | 66–71 | `_baseline_sources`, `_baseline_efficiency`, `_baseline_bc` |
| Knockout helper / impact | `backend/perturbation.py`, `knockout_impact_score` | 73–83 | remove node; signed percent; absolute rounded score |
| Perturbation, output deltas, fragmentation | `backend/perturbation.py`, `run_perturbation` | 85–146 | `strength`, `new_weight`, `E_prime`, `efficiency_delta_pct`, `impact_score`, `delta_bc`, components |
| Request constraints | `backend/models.py`, `PerturbRequest` | 19–34 | literal modes; strength `[0,1]` |
| Attack tolerance | `backend/analysis/attack_tolerance.py`, `_metrics`, `_highest_centrality_node`, `analyze_attack_tolerance` | 30–128 | removal order, efficiency, LCC fraction, trial mean / population SD |
| Threshold-sweep graph metrics | `backend/pipeline/sweep_thresholds.py`, `_graph_metrics`, `_top_betweenness` | 101–120 | `2m/n`, density, components, weighted betweenness |
| Open-target ranking | `backend/analysis/open_targets.py`, `_impact_scores`, `open_targets` | 41–84 | descending impact, symbol tie-break, top-N |
| COSMIC ranking calculation | `backend/analysis/benchmark_cosmic.py`, `_rankings`, benchmark function | 79–136 | weighted BC, raw degree, knockout impact; per-method sorting |
| Structural expression comparisons | `backend/analysis/expression_correlation.py`, `_structural_scores`, `analyze_structural_expression` | 74–80, 108–180 | shared genes, Spearman rho/p-value |

## 10. Uncertainties and limits of what code establishes

- The exact seed gene CSV used to generate the checked-in graph is absent from the inspected tracked file inventory. Therefore the precise initial seed set and historical derivation of the particular `_t900` artifact cannot be confirmed from this checkout alone.
- Runtime loader behavior is established, but the lineages of already-generated CSV artifacts are not encoded as metadata in each CSV. The `_t900` filename, observed score range, and builder convention are consistent with threshold 900, but a runtime loader does not enforce that provenance.
- STRING `preferred_name` is used as the node label by code; the code alone does not establish whether every returned symbol represents a distinct gene or a particular protein-to-gene biological mapping beyond the file mapping.
- NetworkX is delegated shortest-path and centrality operations. GenX specifies arguments but does not implement NetworkX's internal betweenness algorithm in its own source. In addition, version-specific internal behavior cannot be reconstructed from GenX source alone.
- `APPROXIMATION_METHOD` is a label and does not capture the different source-selection mechanics for efficiency versus NetworkX betweenness.
- The baseline caches are module globals initialized at API startup. Calling perturbation helpers in another context before `precompute_baseline`, or reusing them with a different graph after a baseline exists, is not made graph-specific by the code.
- No average shortest path length calculation was located in the backend source inspected. The report does not infer one from efficiency.
- Node disease scores are hard-coded input attributes. Their derivation is not present in the runtime graph construction code.
- Other data analyses (external source retrieval, survival, mutation, expression) are present, but they are not required to define the central graph-perturbation equations. Correlation calculations are noted above; their upstream dataset-generation formulas are outside this report's graph-metric trace.
