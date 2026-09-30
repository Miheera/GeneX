# Experimental path-cost models

GenX retains STRING association confidence as `confidence_score` and normalized
`weight = confidence_score / 1000`. `path_weight` is derived separately for
shortest-path calculations and is selected with `PATH_COST_MODEL`:

* `inverse` (default): `d = 1 / w`, preserving the existing implementation.
* `log`: `d = -log(w)`, an experimental nonnegative transformation. The
  identity `-log(Π wᵢ) = Σ -log(wᵢ)` explains the additive path-cost form when
  inputs are multiplicative reliability-like quantities.

STRING combined scores are association-confidence scores, not automatically
calibrated probabilities. STRING includes functional associations as well as
physical interactions. The log transformation therefore does not establish
biological information flow or validate a biological interpretation of paths.

Select a model before starting the backend, for example in PowerShell:

```powershell
$env:PATH_COST_MODEL = "inverse"  # or "log"
python main.py
```

Restart the backend after changing the setting so the graph and cached baseline
metrics are rebuilt consistently. The dashboard continues to use `weight` for
edge display; it does not display `path_weight`.

Both models require nonnegative edge costs for Dijkstra. In `log`, `w = 1`
produces a zero-cost edge. A reachable zero-cost pair makes GenX's existing
`1 / shortest_distance` global-efficiency term undefined, so efficiency raises
a clear error instead of applying an epsilon or silently assigning a value.
NetworkX weighted betweenness requires strictly positive edge costs; GenX
likewise rejects zero-cost edges for that metric. Overexpression can saturate
`w` at 1 and trigger both guards. No canonical downstream scientific reports
are regenerated when switching models.
