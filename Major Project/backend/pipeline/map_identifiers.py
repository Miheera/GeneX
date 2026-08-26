"""Shared gene-symbol normalization used by pipeline data sources."""

from __future__ import annotations

from collections.abc import Collection

# Keep aliases in one shared place.  Both names resolve to the HGNC symbol; the
# resolver then selects the representation used by a particular graph.
GENE_SYMBOL_ALIASES = {"HER2": "ERBB2", "ERBB2": "ERBB2"}


def resolve_gene_symbol(symbol: str, available_genes: Collection[str] | None = None) -> str:
    """Return a graph-compatible symbol, resolving known aliases deterministically."""
    normalized = symbol.strip().upper()
    if not normalized:
        return ""
    available = {str(gene).upper(): str(gene) for gene in available_genes or ()}
    if normalized in available:
        return available[normalized]

    canonical = GENE_SYMBOL_ALIASES.get(normalized, normalized)
    if canonical in available:
        return available[canonical]
    for alias, alias_canonical in sorted(GENE_SYMBOL_ALIASES.items()):
        if alias_canonical == canonical and alias in available:
            return available[alias]
    return canonical
