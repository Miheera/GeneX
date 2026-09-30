# GenX v2 Methodology Specification

**Status:** frozen description of the repository implementation inspected on 2026-09-29. This is an implementation specification, not a claim of biological validation.

**Authority note:** this specification follows the current executable code. Earlier forensic and path-cost notes still contain v1 descriptions of edge-rescaled inhibition/overexpression and should not be used to describe the v2 partial-perturbation route.

## 1. Executive definition

GenX v2 uses one fixed, undirected STRING association graph for two distinct experiments with different outputs and interpretations:

1. **Virtual Topological Knockout** removes a selected node from a graph copy and measures structural consequences using shortest-path and connectivity metrics.
2. **Fixed-Graph Confidence-Weighted RWR** retains the selected node and every edge, keeps STRING confidence fixed, and propagates a target-centered score over the weighted graph.

The knockout asks how graph structure changes on node deletion. The propagation experiment reports target-centered topological proximity on an unchanged graph. Neither experiment, by itself, measures biological viability, gene essentiality, causal effects, or experimentally observed perturbation response.

## 2. Network construction and loaded graph

### Source and scope

The expansion pipeline is configured for **STRING v12.0 human data**. Its input names use the STRING organism prefix 9606: 9606.protein.info.v12.0.txt.gz and 9606.protein.links.v12.0.txt.gz, under backend/data/raw/string/. The relevant link columns are protein1, protein2, and combined_score; protein IDs are mapped to preferred gene symbols using the protein-info file.

The default expansion threshold is **900**, applied inclusively (combined_score >= 900). The pipeline starts with seed_genes.csv, finds qualifying one-hop partners, then extracts qualifying links whose two mapped endpoints are in the seed-plus-partner set. Thus it produces the induced subgraph on that selected symbol set, rather than the complete human STRING network. The default runtime graph is backend/data/processed/expanded_graph_t900.csv, with its explicit node list in expanded_nodes_t900.csv.

The inspected shipped graph contains:

| Property | Shipped graph |
|---|---:|
| Nodes | 881 |
| Edges | 9,681 |
| confidence_score minimum / median / maximum | 900 / 958 / 999 |
| Self-loops | 0 |
| Duplicate unordered endpoint pairs | 0 |

The processed edge CSV stores gene1, gene2, integer confidence_score, and interaction_type. The expansion pipeline drops same-symbol edges and deduplicates unordered symbol pairs. Runtime loading constructs a NetworkX simple undirected Graph; duplicate rows in an arbitrary override CSV would overwrite edge attributes in that graph, and the loader itself does not reject self-loops. The no-loop/no-duplicate counts above describe the shipped file, not every possible override.

### Meaning of STRING confidence

GenX stores the score as an association-confidence quantity. STRING combined scores represent confidence in an association between proteins, not automatically a calibrated probability of physical interaction, causal influence, activity, or signal transmission. STRING combines functional associations and physical interactions; the resulting GenX graph is not a signed causal regulatory graph.

For each loaded edge, GenX computes:

    weight = confidence_score / 1000.0

The graph loader keeps this normalized STRING association-confidence value in weight. The graph builder serializes it back to integer confidence_score when writing processed CSV. interaction_type is written by the expansion pipeline but is not used by the runtime loader for propagation or path calculations.

## 3. Confidence weight and shortest-path cost are distinct

weight is the normalized STRING association-confidence value. path_weight is a derived cost used by shortest-path structural calculations. The two attributes are not interchangeable:

    STRING combined score
           │
           └──► weight = confidence_score / 1000
                     ├──► row-normalized transitions ──► RWR
                     └──► path-cost transform ─────────► shortest-path metrics

RWR reads weight directly and normalizes the incident weights at each node. RWR does not read path_weight, nor does it apply 1/w or −log(w) to its transition probabilities. Path-cost selection does not define the RWR transition matrix.

## 4. Path-cost models for structural analyses

PATH_COST_MODEL selects the derived path_weight through the centralized confidence_to_path_cost function. The default is inverse; valid selections are inverse and log.

### Inverse (default)

    d_ij = 1 / w_ij

For positive normalized weights, this is a positive cost. A zero weight maps to infinite cost.

### Log (experimental)

    d_ij = −log(w_ij)

For 0 < w ≤ 1, the cost is finite and nonnegative. w=1 gives zero cost; w=0 maps to infinite cost. Negative or non-finite confidence is rejected, and the log transform rejects w>1.

Both transformations are nonnegative on the shipped graph’s domain 0.9 ≤ w ≤ 0.999. They affect Dijkstra distances, weighted efficiency, and weighted betweenness, not RWR. Under log, GenX explicitly raises if global efficiency encounters a reachable zero-cost pair, because its reciprocal-distance term would be undefined. Weighted betweenness also rejects zero-cost edges. The shipped t900 graph has no score-1000 edge, so it has no zero-cost log edges. RWR remains independent of this setting because it consumes weight.

## 5. Virtual Topological Knockout

For target node t, run_perturbation makes a deep copy and removes the target and its incident edges:

    G′ = G \ {t}

The canonical graph is not mutated. Knockout is the only perturbation type handled by this structural routine. The API still accepts a strength value in [0,1], but the deletion operation does not use its magnitude.

### Baseline and perturbed efficiency

At backend startup, the full graph’s baseline efficiency and weighted betweenness are cached. Efficiency uses up to 100 reproducibly selected sources: node IDs are sorted and sampled with Python seed 42. Knockout reuses the cached source list, dropping sources absent from G′.

For each selected source i, Dijkstra computes distances using path_weight. Reachable targets j ≠ i contribute 1/d(i,j); unreachable pairs contribute nothing. A reachable zero-distance pair raises an error. For n ≥ 2 and k surviving selected sources, the implementation calculates:

    E(G) = [Σ_i Σ_(j ≠ i, reachable) 1/d(i,j)] / [k(n−1)]

If fewer than two nodes remain, or no sources remain, efficiency is zero. The denominator includes all n−1 possible targets, including unreachable ones.

The signed percentage change and structural impact are:

    ΔE% = 100 × (E′ − E0) / E0       when E0 > 0; otherwise 0
    impact_score = round(abs(ΔE%), 2)

efficiency_delta_pct is returned without rounding by the routine; impact_score is rounded to two decimals. These are structural/topological consequences of node removal, not biological dose-response estimates.

### Weighted betweenness and change ranking

Baseline and perturbed node betweenness are computed by NetworkX betweenness_centrality using weight="path_weight", normalized output, k=min(100,n), and seed 42 (with zero values returned for graphs of at most one node). For each surviving node other than the target:

    delta_bc(v) = bc_perturbed(v) − bc_baseline(v)

The response includes the five nodes with greatest absolute delta_bc. Positive changes are labeled gain; zero and negative changes are labeled loss by the current implementation.

### Connectivity and returned structural fields

Connected components are computed on G′ without deleting any resulting components. The response reports:

    connected_components = number of components in G′
    lcc_size = size of the largest component (0 if empty)
    lcc_fraction = lcc_size / number of nodes in G′ (0 if empty)
    network_fragmented = connected_components > 1

isolated_genes is the sorted set of surviving nodes outside the largest component when fragmented; it is not limited to degree-zero nodes. If largest components tie, NetworkX component iteration order determines which one max selects. The response also returns bc_baseline, bc_perturbed, approximation_method, target ID, perturbation type, and experiment_type="virtual_topological_knockout".

Knockout measures graph structure under node removal. It is not a direct measurement of biological viability or gene essentiality.

## 6. Fixed-Graph Confidence-Weighted RWR

Inhibition and overexpression call compute_perturbation_propagation, not the structural knockout routine. They read the original graph and do not copy-and-mutate it: nodes, topology, weight, and path_weight remain unchanged.

### Transition probabilities and iteration

For an undirected edge with normalized association confidence w_ij, the implementation forms outgoing weighted neighbor lists and the row-normalized transition probabilities:

    P_ij = w_ij / Σ_k w_ik

The target t receives a one-hot restart vector q_t, with q_t(t)=1 and zero elsewhere. The iteration is:

    r_new = (1−alpha)q_t + alpha Pᵀr

P is row-normalized as a source-to-destination transition matrix. r is represented as a column vector of node mass, so incoming mass at node j is the sum of r_i P_ij; this is why the update uses Pᵀ.

The implementation iterates through sparse neighbor lists rather than creating a dense matrix. It ignores zero-weight edges. Nodes with no positive-weight neighbors are dangling: their current mass is reassigned to the target restart vector. This preserves total mass and gives an isolated target a valid one-hot profile rather than NaNs.

Defaults are read from environment variables when the propagation module is imported:

| Parameter | Environment variable | Default |
|---|---|---:|
| Restart/damping parameter alpha | RWR_ALPHA | 0.85 |
| L1 convergence tolerance | RWR_TOL | 1e-10 |
| Maximum iterations | RWR_MAX_ITER | 1000 |

The stopping condition is ||r_new−r||₁ < tol. If it is not met within max_iter, the routine raises an explicit error. After convergence, the unit profile is normalized only if its total differs from 1 beyond the configured floating-point closeness check. The API returns the iteration count and converged=true; it does not return the residual or tolerance/max-iteration settings.

For a graph of |V| nodes and |E| edges, neighbor-list construction is O(|V|+|E|) and each iteration is O(|V|+|E|) with O(|V|+|E|) stored adjacency. Results are deterministic for a fixed graph ordering and parameter set, subject to floating-point arithmetic.

## 7. Inhibition and overexpression signal

The unit RWR profile is scaled after convergence. The API validates strength in [0,1].

### Inhibition

    S_i = −s r_i

The negative sign is a modeling convention for representing an inhibitory perturbation signal; it is not a claim that the underlying STRING association network contains signed inhibitory interactions. STRING confidence and graph topology are unchanged.

### Overexpression

    S_i = +s r_i

The positive sign represents the direction assigned to an overexpression perturbation signal; it does not mean STRING edges encode activation. STRING confidence and graph topology are unchanged, and the target is not deleted.

The score is a confidence-weighted RWR score for target-centered topological proximity on the fixed STRING association graph. It is not gene-expression change, protein-activity change, causal influence, probability of biological response, drug sensitivity, cell viability, or an experimentally calibrated dose-response.

### Strength invariance

For this linear single-seed model:

    S_i(s) = s S_i(1)       for both perturbation signs
    S_i(1) = +r_i            for overexpression
    S_i(1) = −r_i            for inhibition

Thus strength scales score magnitude but does not change ranking by absolute score for positive strengths. At strength zero, every signed score is zero and the returned ranked list is empty. This is a mathematical property of the implementation, not evidence for a linear biological dose-response relationship. Matched experimental perturbation-response data would be needed for calibration.

## 8. Propagation ranking

The full signed vector is returned as propagation_scores, including the target and zero-score nodes. top_impacted_genes excludes the target and excludes nodes whose absolute score is zero. It returns at most ten nodes, sorted by descending abs(propagation_score); ties are ordered by the string form of the node ID. Each item contains gene, one-based rank, signed propagation_score, and absolute_score. The response also includes ranking_strength_invariant=true and an interpretation string explicitly disclaiming experimentally calibrated response.

## 9. API and dashboard

### Request

POST /api/perturb accepts the same request shape for the two experiment types:

    {
      "gene_id": "TP53",
      "perturbation_type": "inhibition",
      "strength": 0.5
    }

perturbation_type is one of knockout, inhibition, or overexpression; strength defaults to 0.5 and is constrained to [0,1]. An unknown target returns HTTP 404. Request validation errors use FastAPI/Pydantic validation behavior. Computation errors are returned as HTTP 500 with a computation-error detail.

### Knockout response schema

The knockout response contains:

    {
      "experiment_type": "virtual_topological_knockout",
      "efficiency_baseline": 0.0,
      "efficiency_perturbed": 0.0,
      "efficiency_delta_pct": 0.0,
      "impact_score": 0.0,
      "perturbation_type": "knockout",
      "gene_id": "TP53",
      "top_impacted": [
        {"gene": "GENE", "delta_bc": 0.0, "direction": "gain"}
      ],
      "network_fragmented": false,
      "connected_components": 1,
      "lcc_size": 0,
      "lcc_fraction": 0.0,
      "isolated_genes": [],
      "bc_baseline": {"GENE": 0.0},
      "bc_perturbed": {"GENE": 0.0},
      "approximation_method": "approximate (k=100 samples, seed=42)"
    }

The values above illustrate field types and names; they are not a data result. Actual betweenness maps contain graph-node keys. top_impacted.direction is gain or loss.

### Inhibition and overexpression response schema

Both partial perturbation types have this response shape; the type and score signs differ:

    {
      "target_gene": "TP53",
      "perturbation_type": "inhibition",
      "perturbation_strength": 0.5,
      "propagation_model": "confidence_weighted_rwr",
      "alpha": 0.85,
      "propagation_scores": {"TP53": -0.1, "GENE": -0.02},
      "top_impacted_genes": [
        {"gene": "GENE", "rank": 1, "propagation_score": -0.02, "absolute_score": 0.02}
      ],
      "iterations": 97,
      "converged": true,
      "ranking_strength_invariant": true,
      "interpretation": "Target-centered topological proximity on a fixed confidence-weighted STRING association network; not an experimentally calibrated biological response."
    }

The example values are schematic; the response uses the field names above. Partial responses do not include efficiency, ΔE%, structural impact, betweenness change, fragmentation, or LCC fields.

### Dashboard

For knockout, the dashboard shows efficiency before/after, ΔE%, structural impact, and a note with component count and LCC size/fraction; the chart shows the largest betweenness changes. For inhibition/overexpression, it shows target, perturbation, magnitude, “Confidence-weighted RWR,” a fixed-graph/topological-proximity caveat, and the ranked propagated nodes with signed scores. It does not label partial perturbation results as efficiency change, structural impact, fragmentation, or betweenness change.

If the perturbation request fails, the dashboard clears prior perturbation results and displays an error. Separately, graph loading still falls back to a synthetic demo graph when the graph API is unavailable; this is not an RWR result and should not be confused with the shipped network.

## 10. Other analyses and interpretation boundaries

The following remain structural analyses, not validation of the RWR profile:

* Attack tolerance repeatedly removes nodes and reports global efficiency and LCC fraction. Degree and weighted betweenness strategies choose removal order; random removal aggregates trials.
* Open Targets ranking uses the knockout impact_score and joins DGIdb interaction annotations.
* COSMIC benchmarking compares degree, weighted betweenness, and knockout impact rankings with the filtered COSMIC gene set.
* Structural-expression analysis relates knockout structural scores to expression and available mutation/prognostic data.
* Threshold sensitivity reports graph size/connectivity and weighted betweenness at different confidence thresholds.

These outputs retain their structural interpretation. They are not automatically RWR rankings or biological validation of the propagation experiment.

## 11. Methodological evidence status and limitations

Virtual node deletion has precedent as a topological or virtual-knockout analysis. Random Walk with Restart/Personalized PageRank has precedent as a ranking/proximity method on biological interaction networks. GenX uses RWR here as a topology-conditioned ranking method. That precedent does not validate GenX’s complete pipeline or establish biological perturbation prediction on unsigned STRING.

Known limitations:

1. STRING is an association network rather than a signed causal regulatory network and combines different association/evidence types.
2. The loaded graph is undirected; direction and activation/inhibition are not encoded in its edges.
3. RWR measures topology-conditioned proximity, not mechanistic perturbation dynamics.
4. Inhibition and overexpression signs are conventions applied to the target-centered score; STRING edges do not encode those signs.
5. Strength scales scores linearly but does not change their absolute ranking in this formulation.
6. No experimental calibration maps user-selected strength to biological effect.
7. Structural knockout and fixed-graph propagation answer different questions and have different outputs.
8. Path-cost selection affects shortest-path structural analyses but not RWR.
9. Whether propagated rankings correspond to actual inhibition/overexpression response, whether the graph captures the relevant mechanism, and whether score magnitude corresponds to dose remain unvalidated.

Independent biological validation is a separate research stage. Implementation and mathematical checks do not establish predictive validity.

## 12. Frozen model summary

    NETWORK
    G = (V, E, W)
    w_ij = STRING confidence_score / 1000

    STRUCTURAL KNOCKOUT
    G′ = G \ {t}
    structural metrics = f(G′), using path_weight for shortest paths

    RWR
    P_ij = w_ij / Σ_k w_ik           for positive weighted degree
    r_new = (1−alpha)q + alpha Pᵀr
    unit profile is normalized to total mass 1 after convergence
    dangling-node mass is reassigned to q

    INHIBITION
    S_i = −s r_i

    OVEREXPRESSION
    S_i = +s r_i

    RANKING
    exclude target; rank nonzero propagated nodes by |S_i|, descending

**This document represents the GenX v2 methodology as implemented at the time of generation. Future biological-validation experiments should evaluate this frozen specification rather than silently modifying the model.**
