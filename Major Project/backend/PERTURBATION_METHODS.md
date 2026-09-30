# Perturbation methods

## Static network

STRING v12 combined confidence is retained as the static edge attribute
`weight = confidence_score / 1000`. It is an association-confidence score and
is not treated as a calibrated transmission probability. Shortest-path
`path_weight` remains a separate attribute, selected by `PATH_COST_MODEL`.

## Virtual topological knockout

Knockout removes the selected node from a graph copy and recomputes structural
metrics: weighted global efficiency, weighted betweenness, efficiency change,
and connected-component/LCC size and fraction information. This is a **Virtual Topological
Knockout**. Its structural ΔE% is not a biological dose-response prediction.
Existing `inverse` and experimental `log` path-cost options apply only to these
shortest-path structural metrics.

## Fixed-graph inhibition and overexpression

Partial inhibition and overexpression leave every node, edge, STRING
confidence, and path cost unchanged. They compute personalized Random Walk
with Restart using confidence-weighted transition probabilities derived from
`weight`, not `path_weight`. The default is `alpha=0.85`, `tol=1e-10`,
`max_iter=1000`; set `RWR_ALPHA`, `RWR_TOL`, and `RWR_MAX_ITER` before starting
the backend to configure them.

The unit target-centered profile is converted into a signed signal:

* inhibition: `response = -strength * r_unit`
* overexpression: `response = +strength * r_unit`

Scores represent target-centered topological proximity/propagation on a fixed
confidence-weighted STRING association network. They are not experimentally
calibrated biological responses, expression predictions, activity changes,
causal effects, or proof of biological information flow. In this linear
single-seed model, `score(s) = s * score(1)`, so strength scales magnitude but
does not change the absolute ranking. Matched perturbation-response data would
be required to calibrate strength to biological response.

The UI keeps structural knockout results separate from partial-perturbation
propagation results. COSMIC, Open Targets, expression, and attack-tolerance
analyses remain structural analyses; they are not validation of RWR scores.
