"""One-hop STRING expansion from a seed gene set.

Pipeline
--------
1. Read data/processed/seed_genes.csv -> seed gene symbols.
2. Load STRING v12.0 protein-info flat file (gzipped) to map
   gene symbols -> STRING protein IDs (forward mapping) and
   protein IDs -> gene symbols (reverse mapping, for all proteins).
3. Scan STRING v12.0 protein-links file (gzipped) once:
   - Collect all one-hop partners of any seed protein with combined_score >= threshold.
4. Build final node set: seed genes U one-hop partners (mapped to gene symbols).
5. Scan STRING links file a second time to extract the induced subgraph:
   all edges where BOTH endpoints are in the final node set AND combined_score >= threshold.
   (This naturally includes partner-partner edges.)
6. Write:
   - data/processed/expanded_nodes.csv  (gene_symbol, is_seed)
   - data/processed/expanded_edges.csv  (gene1, gene2, confidence_score, interaction_type)
     interaction_type is always "string" to match sweep-graph convention.

Only stdlib + gzip (no live API calls).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)

DEFAULT_THRESHOLD = 900
DEFAULT_SEED_PATH = Path(__file__).parents[1] / "data" / "processed" / "seed_genes.csv"
DEFAULT_INFO_PATH = Path(__file__).parents[1] / "data" / "raw" / "string" / "9606.protein.info.v12.0.txt.gz"
DEFAULT_LINKS_PATH = Path(__file__).parents[1] / "data" / "raw" / "string" / "9606.protein.links.v12.0.txt.gz"
DEFAULT_NODES_OUT = Path(__file__).parents[1] / "data" / "processed" / "expanded_nodes.csv"
DEFAULT_EDGES_OUT = Path(__file__).parents[1] / "data" / "processed" / "expanded_edges.csv"

NODES_COLUMNS = ["gene_symbol", "is_seed"]
EDGES_COLUMNS = ["gene1", "gene2", "confidence_score", "interaction_type"]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _open(path: Path):
    """Open a plain or gzipped text file for reading."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _is_ensp(symbol: str) -> bool:
    """Return True if symbol looks like a raw STRING/Ensembl protein ID rather than a gene name.

    STRING's preferred_name falls back to the ENSP accession when no HGNC
    symbol is available.  These IDs are not useful as gene-symbol node labels.
    """
    return symbol.startswith("ENSP")


def load_seed_genes(path: Path) -> set[str]:
    """Read gene symbols from seed_genes.csv."""
    if not path.exists():
        raise FileNotFoundError(f"Seed gene file not found: {path}")
    genes: set[str] = set()
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sym = row.get("gene_symbol", "").strip()
            if sym:
                genes.add(sym)
    return genes


def load_full_string_mapping(info_path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Return (symbol->protein_id, protein_id->symbol) for ALL human proteins.

    Columns in 9606.protein.info.v12.0.txt:
        #string_protein_id  preferred_name  protein_size  annotation
    """
    sym_to_id: dict[str, str] = {}
    id_to_sym: dict[str, str] = {}
    ensp_only = 0
    with _open(info_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = row.get("#string_protein_id", "").strip()
            sym = row.get("preferred_name", "").strip()
            if not pid or not sym:
                continue
            if _is_ensp(sym):
                # Keep these temporarily so the later graph filter can count what it removes.
                ensp_only += 1
            sym_to_id[sym] = pid
            id_to_sym[pid] = sym
    LOGGER.info("Loaded %d STRING protein mappings (%d ENSP-only entries)",
                len(sym_to_id), ensp_only)
    return sym_to_id, id_to_sym


def find_one_hop_partners(
    links_path: Path,
    seed_pids: set[str],
    id_to_sym: dict[str, str],
    threshold: int,
) -> tuple[set[str], set[str], set[str]]:
    """Return mapped partners plus ENSP-only and unmapped partner protein IDs."""
    partner_syms: set[str] = set()
    ensp_partner_pids: set[str] = set()
    unmapped_partner_pids: set[str] = set()
    with _open(links_path) as fh:
        header = fh.readline().split()
        try:
            p1_idx = header.index("protein1")
            p2_idx = header.index("protein2")
            sc_idx = header.index("combined_score")
        except ValueError as exc:
            raise ValueError("STRING links file lacks expected header columns") from exc
        for line in fh:
            fields = line.split()
            if len(fields) <= sc_idx:
                continue
            score = int(fields[sc_idx])
            if score < threshold:
                continue
            p1 = fields[p1_idx]
            p2 = fields[p2_idx]
            if p1 in seed_pids:
                sym = id_to_sym.get(p2)
                if not sym:
                    unmapped_partner_pids.add(p2)
                elif _is_ensp(sym):
                    ensp_partner_pids.add(p2)
                else:
                    partner_syms.add(sym)
            if p2 in seed_pids:
                sym = id_to_sym.get(p1)
                if not sym:
                    unmapped_partner_pids.add(p1)
                elif _is_ensp(sym):
                    ensp_partner_pids.add(p1)
                else:
                    partner_syms.add(sym)
    return partner_syms, ensp_partner_pids, unmapped_partner_pids


def extract_induced_subgraph(
    links_path: Path,
    node_pids: set[str],
    id_to_sym: dict[str, str],
    threshold: int,
) -> tuple[list[tuple[str, str, int]], set[tuple[str, str]], set[tuple[str, str]]]:
    """Extract valid induced edges and retain dropped ENSP/unmapped edge IDs for logging."""
    edges: list[tuple[str, str, int]] = []
    seen: set[tuple[str, str]] = set()
    ensp_edges_dropped: set[tuple[str, str]] = set()
    unmapped_edges_dropped: set[tuple[str, str]] = set()
    with _open(links_path) as fh:
        header = fh.readline().split()
        try:
            p1_idx = header.index("protein1")
            p2_idx = header.index("protein2")
            sc_idx = header.index("combined_score")
        except ValueError as exc:
            raise ValueError("STRING links file lacks expected header columns") from exc
        for line in fh:
            fields = line.split()
            if len(fields) <= sc_idx:
                continue
            score = int(fields[sc_idx])
            if score < threshold:
                continue
            p1 = fields[p1_idx]
            p2 = fields[p2_idx]
            if p1 not in node_pids or p2 not in node_pids:
                continue
            sym1 = id_to_sym.get(p1)
            sym2 = id_to_sym.get(p2)
            if not sym1 or not sym2:
                unmapped_edges_dropped.add(tuple(sorted((p1, p2))))
                continue
            if _is_ensp(sym1) or _is_ensp(sym2):
                ensp_edges_dropped.add(tuple(sorted((p1, p2))))
                continue
            if sym1 == sym2:
                continue
            key = tuple(sorted((sym1, sym2)))
            if key not in seen:
                seen.add(key)
                edges.append((key[0], key[1], score))
    return edges, ensp_edges_dropped, unmapped_edges_dropped


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------

def expand_via_string(
    seed_path: Path = DEFAULT_SEED_PATH,
    info_path: Path = DEFAULT_INFO_PATH,
    links_path: Path = DEFAULT_LINKS_PATH,
    threshold: int = DEFAULT_THRESHOLD,
    nodes_output: Path = DEFAULT_NODES_OUT,
    edges_output: Path = DEFAULT_EDGES_OUT,
) -> tuple[list[dict], list[dict]]:
    """Run the full seed-and-expand pipeline."""

    # 1. Load seeds
    seed_syms = load_seed_genes(seed_path)
    ensp_seed_syms = {sym for sym in seed_syms if _is_ensp(sym)}
    seed_syms -= ensp_seed_syms
    LOGGER.info("Loaded %d seed genes", len(seed_syms))

    # 2. Load full STRING mapping
    LOGGER.info("Loading STRING protein info from %s ...", info_path.name)
    sym_to_id, id_to_sym = load_full_string_mapping(info_path)

    # Map seeds -> protein IDs (warn on misses)
    seed_pids: set[str] = set()
    unmapped: set[str] = set()
    for sym in seed_syms:
        pid = sym_to_id.get(sym)
        if pid:
            seed_pids.add(pid)
        else:
            unmapped.add(sym)
    if unmapped:
        LOGGER.warning("  %d seed genes not found in STRING info: %s",
                       len(unmapped), ", ".join(sorted(unmapped)))

    # 3. First pass: find one-hop partners
    LOGGER.info("Pass 1/2: scanning STRING links for one-hop partners (threshold=%d) ...", threshold)
    partner_syms, ensp_partner_pids, unmapped_partner_pids = find_one_hop_partners(
        links_path, seed_pids, id_to_sym, threshold
    )
    # Remove seeds themselves from partner set to keep labels clean
    partner_syms -= seed_syms
    LOGGER.info("  Found %d one-hop partner genes", len(partner_syms))

    # 4. Build final node set
    all_syms = seed_syms | partner_syms
    all_pids = {sym_to_id[sym] for sym in all_syms if sym in sym_to_id}
    # Include rejected partners while extracting so every candidate edge that has
    # an ENSP-only or absent STRING mapping is counted and logged before it is
    # excluded from output CSVs.
    all_pids |= ensp_partner_pids | unmapped_partner_pids
    LOGGER.info("  Final node set: %d genes (%d seeds + %d partners)",
                len(all_syms), len(seed_syms), len(partner_syms))

    # 5. Second pass: extract induced subgraph
    LOGGER.info("Pass 2/2: extracting induced subgraph ...")
    edges, ensp_edges_dropped, unmapped_edges_dropped = extract_induced_subgraph(
        links_path, all_pids, id_to_sym, threshold
    )
    LOGGER.info("  Extracted %d edges", len(edges))

    # 6. Write outputs
    node_rows = [
        {"gene_symbol": sym, "is_seed": sym in seed_syms}
        for sym in sorted(all_syms)
    ]
    edge_rows = [
        {"gene1": g1, "gene2": g2, "confidence_score": sc, "interaction_type": "string"}
        for g1, g2, sc in sorted(edges)
    ]
    _write_csv(nodes_output, NODES_COLUMNS, node_rows)
    _write_csv(edges_output, EDGES_COLUMNS, edge_rows)

    # Summary
    degrees: dict[str, int] = {}
    for g1, g2, _ in edges:
        degrees[g1] = degrees.get(g1, 0) + 1
        degrees[g2] = degrees.get(g2, 0) + 1
    avg_degree = sum(degrees.values()) / len(degrees) if degrees else 0.0

    scores = [sc for _, _, sc in edges]
    n_ge900 = sum(1 for s in scores if s >= 900)
    n_800   = sum(1 for s in scores if 800 <= s < 900)
    n_700   = sum(1 for s in scores if 700 <= s < 800)

    print("\n=== expand_via_string summary ===")
    print(f"  Threshold:           {threshold} (STRING combined_score)")
    print(f"  Seed genes:          {len(seed_syms)}")
    if unmapped:
        print(f"  Seeds not in STRING: {len(unmapped)}  [{', '.join(sorted(unmapped))}]")
    if ensp_seed_syms:
        print(f"  ENSP-ID seed nodes dropped: {len(ensp_seed_syms)}")
    print(f"  One-hop partners:    {len(partner_syms)}")
    print(f"  Total nodes:         {len(all_syms)}")
    print(f"  Total edges:         {len(edges)}")
    print(f"  Average degree:      {avg_degree:.2f}")
    print(f"  Score distribution:")
    print(f"    >= 900 (very high): {n_ge900}  ({100*n_ge900/len(edges):.1f}%)" if edges else "    (no edges)")
    print(f"    800–899 (high):     {n_800}  ({100*n_800/len(edges):.1f}%)" if edges else "")
    print(f"    700–799 (medium):   {n_700}  ({100*n_700/len(edges):.1f}%)" if edges else "")
    print(f"  ENSP-ID nodes dropped: {len(ensp_seed_syms) + len(ensp_partner_pids)}")
    print(f"  ENSP-ID edges dropped: {len(ensp_edges_dropped)}")
    print(f"  Unmapped-ID nodes dropped: {len(unmapped_partner_pids)}")
    print(f"  Unmapped-ID edges dropped: {len(unmapped_edges_dropped)}")
    if unmapped_edges_dropped:
        examples = ", ".join("--".join(edge) for edge in sorted(unmapped_edges_dropped)[:5])
        print(f"  Unmapped-ID edge examples: {examples}")
    print(f"  Nodes output:        {nodes_output}")
    print(f"  Edges output:        {edges_output}")
    print()

    return node_rows, edge_rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-path", type=Path, default=DEFAULT_SEED_PATH)
    parser.add_argument("--info-path", type=Path, default=DEFAULT_INFO_PATH)
    parser.add_argument("--links-path", type=Path, default=DEFAULT_LINKS_PATH)
    parser.add_argument("--threshold", type=int, default=DEFAULT_THRESHOLD,
                        help="STRING combined_score threshold (default: 900)")
    parser.add_argument("--nodes-output", type=Path, default=DEFAULT_NODES_OUT)
    parser.add_argument("--edges-output", type=Path, default=DEFAULT_EDGES_OUT)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    expand_via_string(
        seed_path=args.seed_path,
        info_path=args.info_path,
        links_path=args.links_path,
        threshold=args.threshold,
        nodes_output=args.nodes_output,
        edges_output=args.edges_output,
    )


if __name__ == "__main__":
    main()
