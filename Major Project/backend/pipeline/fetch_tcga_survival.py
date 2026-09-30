"""Join TCGA BRCA expression z-scores to patient OS and run log-rank tests."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path
from statistics import mean, median, pstdev

try:
    from .fetch_tcga_expression import (
        SAMPLE_LIST_ID,
        STUDY_ID,
        _post_json,
        fetch_entrez_ids,
        load_graph_genes,
    )
except ImportError:
    from fetch_tcga_expression import (
        SAMPLE_LIST_ID,
        STUDY_ID,
        _post_json,
        fetch_entrez_ids,
        load_graph_genes,
    )

LOGGER = logging.getLogger(__name__)
OS_STATUS = "OS_STATUS"
OS_MONTHS = "OS_MONTHS"
DEFAULT_ZSCORE_PROFILE_ID = f"{STUDY_ID}_rna_seq_v2_mrna_median_all_sample_Zscores"
OUTPUT_COLUMNS = ["gene", "split_method", "logrank_p_value", "hazard_direction", "n_high", "n_low"]

_BACKEND = Path(__file__).parents[1]
DEFAULT_GRAPH_PATH = _BACKEND / "data" / "processed" / "expanded_graph_t900.csv"
DEFAULT_OUTPUT_PATH = _BACKEND / "data" / "analysis_outputs" / "survival_analysis.csv"
MOLECULAR_BATCH_SIZE = 5


def _clinical_records(patient_ids: list[str]) -> dict[str, dict[str, str]]:
    response = _post_json(
        f"/studies/{STUDY_ID}/clinical-data/fetch",
        {"ids": patient_ids, "attributeIds": [OS_STATUS, OS_MONTHS]},
        query={"clinicalDataType": "PATIENT"},
    )
    records: dict[str, dict[str, str]] = {}
    if not isinstance(response, list):
        return records
    for item in response:
        if not isinstance(item, dict):
            continue
        patient_id = str(item.get("patientId") or "").strip()
        attribute = str(item.get("clinicalAttributeId") or "").strip()
        value = str(item.get("value") or "").strip()
        if patient_id and attribute in {OS_STATUS, OS_MONTHS} and value:
            records.setdefault(patient_id, {})[attribute] = value
    return records


def _event(value: str) -> bool:
    normalized = value.upper()
    return "DECEASED" in normalized or normalized.startswith("1")


def _split(values: list[float]) -> tuple[str, list[bool]]:
    LOGGER.info(
        "  z-score distribution: n=%d min=%.3f median=%.3f mean=%.3f sd=%.3f max=%.3f",
        len(values), min(values), median(values), mean(values), pstdev(values), max(values),
    )
    negative_extreme = sum(value <= -2 for value in values)
    positive_extreme = sum(value >= 2 for value in values)
    if negative_extreme >= 10 and positive_extreme >= 10:
        return "cBioPortal z-score threshold +/-2", [value >= 2 for value in values]
    cutoff = median(values)
    return "median z-score split", [value > cutoff for value in values]


def analyze_survival(
    gene_symbols: list[str],
    entrez_by_symbol: dict[str, int],
    values_by_entrez: dict[int, list[tuple[str, str, float]]],
    clinical: dict[str, dict[str, str]],
) -> list[dict[str, str | int | float]]:
    """Run per-gene log-rank tests on complete patient-level observations."""
    from lifelines.statistics import logrank_test

    rows = []
    for symbol in gene_symbols:
        entrez_id = entrez_by_symbol.get(symbol.upper())
        observations = []
        for sample_id, patient_id, value in values_by_entrez.get(entrez_id, []):
            patient_data = clinical.get(patient_id, {})
            try:
                duration = float(patient_data[OS_MONTHS])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(duration) and duration >= 0 and OS_STATUS in patient_data:
                observations.append((patient_id, value, duration, _event(patient_data[OS_STATUS])))
        if len(observations) < 4:
            LOGGER.warning("  Too few complete survival observations for %s", symbol)
            continue
        # A patient can have multiple samples; use the first stable observation per patient.
        unique = {patient_id: (value, duration, event) for patient_id, value, duration, event in observations}
        patient_ids = list(unique)
        values = [unique[patient_id][0] for patient_id in patient_ids]
        split_method, high_flags = _split(values)
        high = [unique[patient_id] for patient_id, is_high in zip(patient_ids, high_flags) if is_high]
        low = [unique[patient_id] for patient_id, is_high in zip(patient_ids, high_flags) if not is_high]
        if len(high) < 2 or len(low) < 2:
            LOGGER.warning("  Split produced too few groups for %s", symbol)
            continue
        result = logrank_test(
            [row[1] for row in high], [row[1] for row in low],
            event_observed_A=[row[2] for row in high], event_observed_B=[row[2] for row in low],
        )
        high_event_rate = sum(row[2] for row in high) / len(high)
        low_event_rate = sum(row[2] for row in low) / len(low)
        if high_event_rate != low_event_rate:
            direction = "high-expression-worse" if high_event_rate > low_event_rate else "high-expression-better"
        else:
            direction = "high-expression-worse" if median(row[1] for row in high) < median(row[1] for row in low) else "high-expression-better"
        rows.append({
            "gene": symbol,
            "split_method": split_method,
            "logrank_p_value": round(float(result.p_value), 12),
            "hazard_direction": direction,
            "n_high": len(high),
            "n_low": len(low),
        })
    return rows


def fetch_tcga_survival(
    graph_path: Path = DEFAULT_GRAPH_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
    *,
    zscore_profile_id: str = DEFAULT_ZSCORE_PROFILE_ID,
) -> list[dict[str, str | int | float]]:
    print(f"Loading graph genes from: {graph_path}")
    gene_symbols = load_graph_genes(graph_path)
    print(f"Loaded {len(gene_symbols)} unique genes from {graph_path.name}")
    entrez_by_symbol = fetch_entrez_ids(gene_symbols)
    if not entrez_by_symbol:
        return []
    values_by_entrez: dict[int, list[tuple[str, str, float]]] = {}
    patient_ids = set()
    entrez_ids = sorted(set(entrez_by_symbol.values()))
    total_batches = (len(entrez_ids) + MOLECULAR_BATCH_SIZE - 1) // MOLECULAR_BATCH_SIZE
    for batch_number, start in enumerate(range(0, len(entrez_ids), MOLECULAR_BATCH_SIZE), 1):
        batch_ids = entrez_ids[start:start + MOLECULAR_BATCH_SIZE]
        LOGGER.info("Fetching expression batch %d/%d (%d genes)", batch_number, total_batches, len(batch_ids))
        molecular_data = _post_json(
            f"/molecular-profiles/{zscore_profile_id}/molecular-data/fetch",
            {"sampleListId": SAMPLE_LIST_ID, "entrezGeneIds": batch_ids},
        )
        if not isinstance(molecular_data, list):
            LOGGER.warning("Skipping expression batch %d/%d after request failure", batch_number, total_batches)
            continue
        for item in molecular_data:
            if not isinstance(item, dict):
                continue
            entrez_id = item.get("entrezGeneId")
            sample_id = str(item.get("sampleId") or "").strip()
            patient_id = str(item.get("patientId") or "").strip()
            value = item.get("value")
            if isinstance(entrez_id, int) and sample_id and patient_id and isinstance(value, (int, float)):
                values_by_entrez.setdefault(entrez_id, []).append((sample_id, patient_id, float(value)))
                patient_ids.add(patient_id)
    LOGGER.info("Fetching clinical survival data once for %d patients", len(patient_ids))
    clinical = _clinical_records(sorted(patient_ids))
    LOGGER.info("Received clinical records for %d patients", len(clinical))
    rows = analyze_survival(gene_symbols, entrez_by_symbol, values_by_entrez, clinical)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n=== fetch_tcga_survival summary ===\n  Genes scored: {len(rows)}\n  Output:       {output_path}\n")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-path", type=Path, default=DEFAULT_GRAPH_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--zscore-profile-id", default=DEFAULT_ZSCORE_PROFILE_ID)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fetch_tcga_survival(args.graph_path, args.output_path, zscore_profile_id=args.zscore_profile_id)


if __name__ == "__main__":
    main()