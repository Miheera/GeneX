"""Filtering helpers for STRING protein links."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TypeAlias

StringEdge: TypeAlias = tuple[str, str, int]


def filter_string_links(
    edges: Iterable[StringEdge], combined_score_threshold: int
) -> list[StringEdge]:
    """Return STRING links whose combined score meets the explicit threshold."""
    return [edge for edge in edges if edge[2] >= combined_score_threshold]


# Compatibility alias for callers that use the shorter pipeline name.
filter_string = filter_string_links
