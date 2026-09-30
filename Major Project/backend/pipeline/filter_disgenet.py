"""Extract and union top-N seed genes from DisGeNET XLSX exports.

Design notes
------------
* Input format: DisGeNET web-UI XLSX exports.  The first row is the column
  header; the score column is named ``score`` (GDA score, 0-1).
* For genes that appear in both CUI queries, the row with the higher ``score``
  is kept and both source CUIs are recorded in ``source_cui``.
* Output: ``data/processed/seed_genes.csv``
  Columns: gene_symbol, score, source_cui, in_both
"""

from __future__ import annotations

import argparse
import csv
import logging
from dataclasses import dataclass
from pathlib import Path

import openpyxl

LOGGER = logging.getLogger(__name__)

DEFAULT_TOP_N = 30
OUTPUT_COLUMNS = ["gene_symbol", "score", "source_cui", "in_both"]

_BACKEND = Path(__file__).parents[1]
_DISGENET_RAW = _BACKEND / "data" / "raw" / "disgenet"
DEFAULT_INPUTS: list[tuple[str, Path]] = [
    ("C0006142", _DISGENET_RAW / "search_C0006142" / "search_result_C0006142.xlsx"),
    ("C0678222", _DISGENET_RAW / "search_C0678222" / "search_result_C0678222.xlsx"),
]
DEFAULT_OUTPUT = _BACKEND / "data" / "processed" / "seed_genes.csv"


@dataclass
class GDARow:
    gene_symbol: str
    score: float
    cui: str


def _read_xlsx(path: Path, cui: str) -> list[GDARow]:
    """Return every GDA row from a DisGeNET XLSX export."""
    if not path.exists():
        raise FileNotFoundError(f"DisGeNET XLSX not found: {path}")
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows_iter = ws.iter_rows(values_only=True)
    header: list[str] = []
    for raw_row in rows_iter:
        header = [str(cell).strip() if cell is not None else "" for cell in raw_row]
        if any(header):
            break
    if not header:
        raise ValueError(f"XLSX has no header row: {path}")
    # The XLSX header uses display labels, not API field names.
    # Try API name first, then the shortened display label as fallback.
    def _find_col(candidates: list[str]) -> int:
        for name in candidates:
            try:
                return header.index(name)
            except ValueError:
                continue
        raise ValueError(
            f"XLSX missing gene-symbol column (tried {candidates}); "
            f"found columns: {header}"
        )

    symbol_idx = _find_col(["geneSymbol", "Gene"])
    score_idx  = _find_col(["score", "Score"])

    records: list[GDARow] = []
    for raw_row in rows_iter:
        gene = raw_row[symbol_idx]
        score = raw_row[score_idx]
        if gene is None or score is None:
            continue
        gene_str = str(gene).strip()
        try:
            score_float = float(score)
        except (TypeError, ValueError):
            LOGGER.debug("Skipping non-numeric score for gene %s: %r", gene_str, score)
            continue
        if gene_str:
            records.append(GDARow(gene_symbol=gene_str, score=score_float, cui=cui))
    wb.close()
    LOGGER.info("  Read %d GDA rows from %s (CUI %s)", len(records), path.name, cui)
    return records


def _top_n(records: list[GDARow], n: int) -> list[GDARow]:
    return sorted(records, key=lambda r: (-r.score, r.gene_symbol))[:n]


def filter_disgenet(
    inputs: list[tuple[str, Path]],
    top_n: int = DEFAULT_TOP_N,
    output_path: Path = DEFAULT_OUTPUT,
) -> list[dict]:
    """Union top-N seeds from each CUI query; keep higher score for duplicates."""
    per_query: dict[str, list[GDARow]] = {}
    for cui, path in inputs:
        records = _read_xlsx(path, cui)
        top = _top_n(records, top_n)
        per_query[cui] = top
        LOGGER.info("  Top %d from CUI %s: %s",
                    top_n, cui,
                    ", ".join(r.gene_symbol for r in top[:5]) + (" ..." if len(top) > 5 else ""))
    best: dict[str, dict] = {}
    for cui, top_records in per_query.items():
        for row in top_records:
            sym = row.gene_symbol
            if sym not in best:
                best[sym] = {"score": row.score, "cuis": {cui}}
            else:
                best[sym]["cuis"].add(cui)
                if row.score > best[sym]["score"]:
                    best[sym]["score"] = row.score
    output_rows = []
    for sym in sorted(best, key=lambda s: (-best[s]["score"], s)):
        entry = best[sym]
        output_rows.append({
            "gene_symbol": sym,
            "score": round(entry["score"], 6),
            "source_cui": "|".join(sorted(entry["cuis"])),
            "in_both": len(entry["cuis"]) > 1,
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(output_rows)
    cuis = list(per_query.keys())
    symbols_per = {cui: {r.gene_symbol for r in rows} for cui, rows in per_query.items()}
    all_sets = list(symbols_per.values())
    overlap = set.intersection(*all_sets) if len(all_sets) >= 2 else set()
    print("\n=== filter_disgenet summary ===")
    for cui in cuis:
        print(f"  CUI {cui}: {len(per_query[cui])} genes selected (top {top_n} by GDA score)")
    if overlap:
        print(f"  Overlap (in both queries): {len(overlap)} genes  [{', '.join(sorted(overlap)[:10])}{'...' if len(overlap) > 10 else ''}]")
    else:
        print("  Overlap (in both queries): 0 genes")
    print(f"  Final unique seed genes:  {len(output_rows)}")
    print(f"  Output: {output_path}")
    print()
    return output_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--input", nargs=2, action="append", metavar=("CUI", "XLSX_PATH"), dest="inputs")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    inputs: list[tuple[str, Path]] = (
        [(cui, Path(p)) for cui, p in args.inputs] if args.inputs else DEFAULT_INPUTS
    )
    filter_disgenet(inputs, top_n=args.top_n, output_path=args.output_path)


if __name__ == "__main__":
    main()
