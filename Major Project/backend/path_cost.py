"""Central conversion from normalized STRING confidence to path cost.

``log`` is an experimental cost transform, not a claim that STRING scores are
calibrated probabilities or that paths model biological information flow.
"""

from __future__ import annotations

import math
import os

PATH_COST_MODEL_ENV = "PATH_COST_MODEL"
DEFAULT_PATH_COST_MODEL = "inverse"
VALID_PATH_COST_MODELS = frozenset({"inverse", "log"})


def get_path_cost_model() -> str:
    model = os.environ.get(PATH_COST_MODEL_ENV, DEFAULT_PATH_COST_MODEL).strip().lower()
    if model not in VALID_PATH_COST_MODELS:
        raise ValueError(f"{PATH_COST_MODEL_ENV} must be one of {sorted(VALID_PATH_COST_MODELS)}; got {model!r}")
    return model


def confidence_to_path_cost(w: float, model: str | None = None) -> float:
    """Convert normalized STRING confidence ``w`` to a nonnegative cost.

    Zero confidence maps to infinite cost (no finite-cost route), matching the
    existing loader's handling. Negative confidence is invalid. The log model
    requires the normalized domain ``0 <= w <= 1``; ``w == 1`` yields zero.
    """
    selected = get_path_cost_model() if model is None else model.strip().lower()
    if selected not in VALID_PATH_COST_MODELS:
        raise ValueError(f"Unknown path-cost model: {selected!r}")
    value = float(w)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Normalized confidence must be finite and nonnegative; got {w!r}")
    if value == 0:
        return math.inf
    if selected == "inverse":
        return 1.0 / value
    if value > 1:
        raise ValueError(f"The log path-cost model requires normalized confidence <= 1; got {value}")
    return -math.log(value)


def require_positive_betweenness_costs(graph) -> None:
    """NetworkX weighted betweenness requires strictly positive edge weights."""
    zero_edges = [(u, v) for u, v, data in graph.edges(data=True) if data.get("path_weight", 1.0) == 0]
    if zero_edges:
        raise ValueError(
            "Weighted betweenness is undefined/unsupported for zero path costs; "
            f"found {len(zero_edges)} zero-cost edge(s), e.g. {zero_edges[0]!r}"
        )
