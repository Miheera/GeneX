"""Filter a manually downloaded COSMIC Cancer Gene Census export.

The current repository does not retain graph-universe provenance, so downstream
benchmarks must report that provenance as unknown unless explicitly supplied.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)
BREAST_MATCH_TERMS = ("breast", "brca", "mammary")
OUTPUT_COLUMNS = ["gene_symbol", "cosmic_tier", "breast_annotated", "matched_term", "is_tier1"]


def _tier_number(value: str) -> int | None:
    digits = "".join(character for character in value if character.isdigit())
    return int(digits) if digits else None


def filter_cosmic(
    input_path: Path,
    output_path: Path | None = None,
    *,
    breast_match_terms: tuple[str, ...] = BREAST_MATCH_TERMS,
) -> list[dict[str, str | bool | int | None]]:
    """Create a deduplicated breast-annotated COSMIC gene table."""
    if not input_path.exists():
        raise FileNotFoundError(f"COSMIC input file does not exist: {input_path}")
    if output_path is None:
        output_path = Path(__file__).parents[1] / "data" / "processed" / "cosmic_breast_genes.csv"
    if not breast_match_terms or any(not term.strip() for term in breast_match_terms):
        raise ValueError("breast_match_terms must contain non-empty terms")

    required = {"Gene Symbol", "Tier", "Tumour Types(Somatic)", "Tumour Types(Germline)"}
    rows_seen = 0
    records: dict[str, dict[str, str | bool | int | None]] = {}
    with input_path.open(newline="", encoding="utf-8-sig") as input_file:
        reader = csv.DictReader(input_file)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"COSMIC input is missing required columns: {', '.join(missing)}")
        for row in reader:
            rows_seen += 1
            gene = (row.get("Gene Symbol") or "").strip()
            if not gene:
                continue
            tier = _tier_number(row.get("Tier") or "")
            annotations = " ".join(
                (row.get(column) or "") for column in ("Tumour Types(Somatic)", "Tumour Types(Germline)")
            ).lower()
            matched = sorted({term.lower() for term in breast_match_terms if term.lower() in annotations})
            current = records.get(gene)
            if current is None or (
                tier is not None and (current["cosmic_tier"] is None or tier < current["cosmic_tier"])
            ):
                previous_terms = set()
                previous_breast = False
                if current is not None:
                    previous_terms = set(filter(None, str(current["matched_term"]).split("|")))
                    previous_breast = bool(current["breast_annotated"])
                records[gene] = {
                    "gene_symbol": gene,
                    "cosmic_tier": tier,
                    "breast_annotated": bool(matched) or previous_breast,
                    "matched_term": "|".join(sorted(previous_terms | set(matched))),
                    "is_tier1": tier == 1,
                }
            else:
                current["breast_annotated"] = bool(current["breast_annotated"] or matched)
                terms = set(filter(None, str(current["matched_term"]).split("|"))) | set(matched)
                current["matched_term"] = "|".join(sorted(terms))
                current["is_tier1"] = bool(current["is_tier1"] or tier == 1)

    result = [records[gene] for gene in sorted(records)]
    breast_count = sum(bool(row["breast_annotated"]) for row in result)
    tier1_count = sum(bool(row["is_tier1"]) for row in result)
    LOGGER.info(
        "COSMIC rows=%d breast_annotated=%d tier1=%d unique_genes=%d",
        rows_seen, breast_count, tier1_count, len(result),
    )
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_path", type=Path)
    parser.add_argument(
        "--output-path", type=Path,
        default=Path(__file__).parents[1] / "data" / "processed" / "cosmic_breast_genes.csv",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    filter_cosmic(args.input_path, args.output_path)


if __name__ == "__main__":
    main()
