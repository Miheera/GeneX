#!/usr/bin/env python3
"""Metadata-only LINCS L1000 coverage audit for frozen GenX v2.

This deliberately does not load expression signatures or compare them with GenX.
It inventories GSE92742 Level 5 profile metadata and landmark coverage only.
"""
from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "validation"
META = VAL / "metadata"
OUT = VAL / "audit"
MAP = VAL / "mappings"
GRAPH = ROOT / "backend/data/processed/expanded_nodes_t900.csv"

def rows(path: Path, compressed: bool = False, delimiter: str = "\t"):
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as f:
        yield from csv.DictReader(f, delimiter=delimiter)

def write_csv(path: Path, fields, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)

def split_names(value):
    return [x.strip() for x in (value or "").split("|") if x.strip()]

def main():
    graph_nodes = [r["gene_symbol"] for r in rows(GRAPH, delimiter=",")]
    hgnc_rows = list(rows(META / "hgnc_complete_set.txt"))
    approved = defaultdict(list)
    aliases = defaultdict(list)
    entrez = defaultdict(list)
    for h in hgnc_rows:
        approved[h["symbol"].upper()].append(h)
        for field in ("alias_symbol", "prev_symbol"):
            for name in split_names(h.get(field)):
                aliases[name.upper()].append((h, field))
        if h.get("entrez_id"):
            entrez[h["entrez_id"]].append(h)

    def map_hgnc(name):
        key = (name or "").strip().upper()
        exact = approved.get(key, [])
        if len(exact) == 1:
            return exact[0], "approved_symbol", "unique"
        if len(exact) > 1:
            return None, "approved_symbol", "ambiguous"
        found = {x[0]["hgnc_id"]: x for x in aliases.get(key, [])}
        if len(found) == 1:
            h, source = next(iter(found.values()))
            return h, source, "unique_alias"
        if len(found) > 1:
            return None, "alias_or_previous_symbol", "ambiguous"
        return None, "none", "unmapped"

    graph_map = {}
    mapping_audit = []
    for node in graph_nodes:
        h, method, status = map_hgnc(node)
        graph_map[node] = (h, method, status)
        mapping_audit.append({
            "genx_id": node, "genx_symbol": node, "external_database": "HGNC",
            "external_id": h["hgnc_id"] if h else "",
            "mapping_method": method,
            "mapping_status": status,
            "ambiguity_flag": status == "ambiguous",
            "notes": (f"Entrez={h.get('entrez_id','')}; Ensembl={h.get('ensembl_gene_id','')}" if h else "No fuzzy matching; exact approved symbol then exact alias/previous symbol."),
        })

    graph_by_entrez = defaultdict(list)
    for node, (h, _, status) in graph_map.items():
        if h and status in ("unique", "unique_alias") and h.get("entrez_id"):
            graph_by_entrez[h["entrez_id"]].append(node)

    lincs_genes = list(rows(META / "GSE92742_Broad_LINCS_gene_info.txt.gz", compressed=True))
    landmark = [g for g in lincs_genes if g.get("pr_is_lm") == "1"]
    lm_id_to_node = {}
    l1000_overlap = []
    overlap_by_node = defaultdict(list)
    for g in landmark:
        hit = graph_by_entrez.get(g["pr_gene_id"], [])
        if len(hit) == 1:
            node = hit[0]
            lm_id_to_node[g["pr_gene_id"]] = node
            overlap_by_node[node].append(g)
            h = graph_map[node][0]
            l1000_overlap.append({
                "genx_node_id": node,
                "gene_symbol": node,
                "stable_gene_id": h.get("hgnc_id", "") if h else "",
                "entrez_gene_id": g["pr_gene_id"],
                "ensembl_gene_id": h.get("ensembl_gene_id", "") if h else "",
                "l1000_gene_id": g["pr_gene_id"],
                "mapping_status": "unique_entrez_exact",
                "landmark_status": "directly_measured_landmark",
            })
    # Level 5 includes the 12k feature universe, where non-landmarks are
    # inferred rather than directly measured. Keep this count separate.
    feature_node_ids = set()
    for g in lincs_genes:
        if g.get("pr_is_lm") == "1":
            continue
        hits = graph_by_entrez.get(g["pr_gene_id"], [])
        if len(hits) == 1:
            feature_node_ids.add(hits[0])
    inferred_overlap_count = len(feature_node_ids)

    # Include all graph nodes so non-overlap and mapping exclusions are explicit.
    overlap_nodes = {r["genx_node_id"] for r in l1000_overlap}
    for node in graph_nodes:
        if node in overlap_nodes:
            continue
        h, _, status = graph_map[node]
        l1000_overlap.append({
            "genx_node_id": node, "gene_symbol": node,
            "stable_gene_id": h.get("hgnc_id", "") if h else "",
            "entrez_gene_id": h.get("entrez_id", "") if h else "",
            "ensembl_gene_id": h.get("ensembl_gene_id", "") if h else "",
            "l1000_gene_id": "",
            "mapping_status": status if not h else "mapped_hgnc_but_not_landmark",
            "landmark_status": "not_in_direct_landmark_overlap",
        })
    write_csv(OUT / "l1000_genx_landmark_overlap.csv", [
        "genx_node_id", "gene_symbol", "stable_gene_id", "entrez_gene_id",
        "ensembl_gene_id", "l1000_gene_id", "mapping_status", "landmark_status"
    ], l1000_overlap)

    # A signature row remains a distinct profile. trt_sh.cgs/css are retained
    # separately from individual shRNA reagent profiles; none are collapsed.
    genetic_types = {"trt_sh", "trt_sh.cgs", "trt_sh.css", "trt_oe"}
    control_types = {"ctl_vehicle", "ctl_vector", "ctl_untrt", "ctl_vector.cns", "ctl_vehicle.cns"}
    all_mcf7 = []
    control_by_context = defaultdict(int)
    control_by_type_context = defaultdict(int)
    sig_meta = META / "GSE92742_Broad_LINCS_sig_info.txt.gz"
    for r in rows(sig_meta, compressed=True):
        if r.get("cell_id") != "MCF7":
            continue
        ctx = (r.get("cell_id", ""), r.get("pert_time", ""), r.get("pert_time_unit", ""))
        if r.get("pert_type") in control_types:
            control_by_context[ctx] += 1
            control_by_type_context[(r.get("pert_type"),) + ctx] += 1
    for r in rows(sig_meta, compressed=True):
        if r.get("cell_id") == "MCF7" and r.get("pert_type") in genetic_types:
            all_mcf7.append(r)

    target_names = {r.get("pert_iname", "") for r in all_mcf7}
    selected_sig_ids = {x["sig_id"] for x in all_mcf7}
    metric_by_sig = {}
    for r in rows(META / "GSE92742_Broad_LINCS_sig_metrics.txt.gz", compressed=True):
        if r.get("sig_id") and r["sig_id"] in selected_sig_ids:
            metric_by_sig[r["sig_id"]] = r

    profile_rows = []
    target_type_profiles = defaultdict(set)
    target_coverage = {}
    for r in all_mcf7:
        h, map_method, map_status = map_hgnc(r.get("pert_iname", ""))
        node = ""
        if h and map_status in ("unique", "unique_alias"):
            possible = [n for n, (gh, _, gs) in graph_map.items()
                        if gh and gs in ("unique", "unique_alias") and gh["hgnc_id"] == h["hgnc_id"]]
            if len(possible) == 1:
                node = possible[0]
        metric = metric_by_sig.get(r["sig_id"], {})
        n_landmark_overlap = len(overlap_nodes)  # constant directly measured feature universe count
        if node:
            target_type_profiles[(node, r["pert_type"])].add(r["sig_id"])
            target_coverage[node] = n_landmark_overlap
        ctx = (r.get("cell_id", ""), r.get("pert_time", ""), r.get("pert_time_unit", ""))
        required_controls = ("ctl_vector", "ctl_vector.cns") if r.get("pert_type") == "trt_oe" else ("ctl_untrt", "ctl_vehicle", "ctl_vehicle.cns")
        required_control_count = sum(control_by_type_context.get((ct,) + ctx, 0) for ct in required_controls)
        profile_rows.append({
            "genx_node_id": node,
            "gene_symbol": node or r.get("pert_iname", ""),
            "stable_gene_id": h.get("hgnc_id", "") if h else "",
            "entrez_gene_id": h.get("entrez_id", "") if h else "",
            "l1000_target_id": r.get("pert_id", ""),
            "pert_iname": r.get("pert_iname", ""),
            "target_mapping_method": map_method,
            "target_mapping_status": map_status if node else ("not_in_genx" if h else map_status),
            "perturbation_type": r.get("pert_type", ""),
            "cell_line": r.get("cell_id", ""),
            "time": r.get("pert_time", ""),
            "time_unit": r.get("pert_time_unit", ""),
            "dose": r.get("pert_dose", ""),
            "dose_unit": r.get("pert_dose_unit", ""),
            "profile_id": r.get("sig_id", ""),
            "replicate_count": metric.get("distil_nsample", ""),
            "control_signature_count_same_cell_time": control_by_context.get(ctx, 0),
            "same_cell_time_control_available": control_by_context.get(ctx, 0) > 0,
            "appropriate_control_type_same_cell_time_count": required_control_count,
            "appropriate_control_context_available": required_control_count > 0,
            "control_pair_identifiable": "not_determinable_from_sig_info; Level 5 is processed signature",
            "profile_quality_distil_cc_q75": metric.get("distil_cc_q75", ""),
            "profile_quality_distil_ss": metric.get("distil_ss", ""),
            "profile_quality_tas": metric.get("tas", ""),
            "profile_quality_is_exemplar": metric.get("is_exemplar", ""),
            "profile_quality_metadata_present": bool(metric),
            "directly_measured_overlap_count": n_landmark_overlap,
            "inferred_overlap_count": inferred_overlap_count,
            "level5_signature_metadata_present": True,
            "eligible": "pending exact matched-control rule and predeclared response-gene threshold" if node else False,
            "exclusion_reason": "" if node else ("perturbagen does not map uniquely to a GenX node" if map_status == "ambiguous" else "perturbagen not in frozen GenX graph"),
        })
    write_csv(OUT / "l1000_target_coverage.csv", list(profile_rows[0].keys()) if profile_rows else [], profile_rows)

    # Append L1000 external identifier records to the mapping audit.
    for row in l1000_overlap:
        if row["l1000_gene_id"]:
            mapping_audit.append({
                "genx_id": row["genx_node_id"], "genx_symbol": row["gene_symbol"],
                "external_database": "LINCS_L1000_GSE92742",
                "external_id": row["l1000_gene_id"], "mapping_method": "HGNC Entrez ID exact match",
                "mapping_status": "unique_direct_landmark", "ambiguity_flag": False,
                "notes": "pr_is_lm=1 in GSE92742 gene_info; direct landmark, not inferred feature.",
            })
    write_csv(OUT / "gene_mapping_audit.csv", [
        "genx_id", "genx_symbol", "external_database", "external_id",
        "mapping_method", "mapping_status", "ambiguity_flag", "notes"
    ], mapping_audit)

    # Broad DepMap is explicitly marked unavailable until its official pinned
    # release files can be obtained; no coverage or effect values are guessed.
    depmap_rows = []
    for node in graph_nodes:
        h, _, status = graph_map[node]
        depmap_rows.append({
            "genx_node_id": node, "gene_symbol": node,
            "stable_gene_id": h.get("hgnc_id", "") if h else "",
            "entrez_gene_id": h.get("entrez_id", "") if h else "",
            "depmap_gene_id": "", "breast_model_count": "",
            "nonmissing_effect_count": "", "median_gene_effect": "",
            "eligible": "not_assessed",
            "exclusion_reason": "DepMap 26Q1 files not retrieved; official portal endpoint returned human-verification page.",
            "audit_status": "blocked_pending_official_files",
        })
    write_csv(OUT / "depmap_genx_coverage.csv", [
        "genx_node_id", "gene_symbol", "stable_gene_id", "entrez_gene_id",
        "depmap_gene_id", "breast_model_count", "nonmissing_effect_count",
        "median_gene_effect", "eligible", "exclusion_reason", "audit_status"
    ], depmap_rows)

    provenance_rows = [
        {"dataset": "LINCS L1000", "release": "GSE92742; GEO supplemental files dated 2017-12-13; Level 5 release n=473647 signatures",
         "publication": "Subramanian et al., A Next Generation Connectivity Map", "publication_date": "2017-11-30",
         "STRING_version": "STRING v12.0 graph frozen by GenX",
         "potential_overlap": "possible through STRING text-mining/experimental evidence; exact edge provenance not assessed",
         "evidence": "GSE92742 data and signatures predate the current STRING release; STRING integrates evidence channels including experiments and literature text mining.",
         "determination": "Not determinable from available metadata.",
         "limitation": "GenX graph contains combined scores, not per-edge provenance; no source-paper-to-edge audit was possible."},
        {"dataset": "DepMap Public", "release": "26Q1 (current release observed 2026-09-29; files not downloaded)",
         "publication": "DepMap release notes; release data are published separately from papers", "publication_date": "2026-04-01 release announcement",
         "STRING_version": "STRING v12.0 graph frozen by GenX",
         "potential_overlap": "possible at the level of historical publications/evidence; release itself postdates STRING v12.0",
         "evidence": "Official portal release notes announce 26Q1; STRING combines experiments, curated databases and text mining.",
         "determination": "Not determinable from available metadata.",
         "limitation": "Official download endpoint returned a human-verification page; no 26Q1 matrix or per-edge provenance was inspected."},
        {"dataset": "Project Score / Cell Model Passports", "release": "Current Sanger data version not retrieved; historic Project Score v2 described in 2021 paper",
         "publication": "Dwane et al., Project Score database", "publication_date": "2021-01-08",
         "STRING_version": "STRING v12.0 graph frozen by GenX",
         "potential_overlap": "possible through historical experimental evidence or literature; not established",
         "evidence": "Official documentation describes CRISPR knockout fitness data and a public API; direct API requests failed in this environment.",
         "determination": "Not determinable from available metadata.",
         "limitation": "No Sanger target/model matrix or provenance table was retrieved; it is not used as a substitute for DepMap."},
    ]
    write_csv(OUT / "validation_provenance_audit.csv", list(provenance_rows[0].keys()), provenance_rows)

    # Deterministic availability summary (no expression values or performance stats).
    nodes_by_type = defaultdict(set)
    profiles_by_type = defaultdict(int)
    for r in profile_rows:
        if r["genx_node_id"]:
            nodes_by_type[r["perturbation_type"]].add(r["genx_node_id"])
            profiles_by_type[r["perturbation_type"]] += 1
    knock_types = {"trt_sh", "trt_sh.cgs", "trt_sh.css"}
    knock_nodes = set().union(*(nodes_by_type[x] for x in knock_types))
    oe_nodes = nodes_by_type["trt_oe"]
    union_nodes = knock_nodes | oe_nodes
    print(f"GenX nodes: {len(graph_nodes)}")
    print(f"HGNC unique-approved-or-alias mapped graph nodes: {sum(1 for n in graph_nodes if graph_map[n][0] and graph_map[n][2] in ('unique','unique_alias'))}")
    print(f"L1000 landmarks: {len(landmark)}; directly measured unique GenX overlap: {len(overlap_nodes)}")
    print(f"MCF7 Level 5 genetic signature profiles: {len(profile_rows)}")
    print("MCF7 GenX perturbation target unique genes by type:", {k: len(v) for k,v in sorted(nodes_by_type.items())})
    print("MCF7 GenX knockdown unique targets (including cgs/css):", len(knock_nodes))
    print("MCF7 GenX cDNA overexpression unique targets:", len(oe_nodes))
    print("MCF7 GenX targets in either modality:", len(union_nodes))
    print("Targets both knockdown and overexpression:", len(knock_nodes & oe_nodes))
    print("Landmark coverage buckets (target genes with profiles; same feature universe per profile):")
    print({"0-9": sum(1 for n in union_nodes if target_coverage.get(n,0) < 10),
           "10-19": sum(1 for n in union_nodes if 10 <= target_coverage.get(n,0) < 20),
           "20-49": sum(1 for n in union_nodes if 20 <= target_coverage.get(n,0) < 50),
           "50-99": sum(1 for n in union_nodes if 50 <= target_coverage.get(n,0) < 100),
           "100+": sum(1 for n in union_nodes if target_coverage.get(n,0) >= 100)})
    print("Profiles with matched-control identity: not determinable from downloaded sig_info; metadata records only Level 5 processed signatures and same-cell/time control-signature counts.")
    print("Profiles with quality metadata:", len(metric_by_sig), "/", len(profile_rows))
    for col in ("profile_quality_distil_cc_q75", "profile_quality_distil_ss", "profile_quality_tas", "profile_quality_is_exemplar"):
        valid = sum(r.get(col, "") not in ("", "-666", "-666.0") for r in profile_rows)
        print(f"Profiles with non-sentinel {col}: {valid} / {len(profile_rows)}")
    valid_reps = [int(r["replicate_count"]) for r in profile_rows if str(r.get("replicate_count", "")).isdigit()]
    print("Replicate-count values:", dict(sorted(__import__('collections').Counter(valid_reps).items())))
    print("Profiles with a control signature at same cell and time:", sum(bool(r["same_cell_time_control_available"]) for r in profile_rows), "/", len(profile_rows), "(context-level availability, not explicit per-profile pairing)")
    print("Profiles with a same-time control type appropriate to intervention:", sum(bool(r["appropriate_control_context_available"]) for r in profile_rows), "/", len(profile_rows), "(context-level only; not an explicit per-profile pairing)")
    print("Profiles with replicate_count >=2 / >=3:", sum(str(r["replicate_count"]).isdigit() and int(r["replicate_count"]) >= 2 for r in profile_rows), "/", sum(str(r["replicate_count"]).isdigit() and int(r["replicate_count"]) >= 3 for r in profile_rows))
    gene_type_profiles = defaultdict(list)
    for r in profile_rows:
        if r["genx_node_id"]:
            gene_type_profiles[(r["genx_node_id"], r["perturbation_type"])].append(r)
    print("GenX targets by perturbation type with any profile having appropriate same-time control context:", {
        typ: len({g for (g,t), rs in gene_type_profiles.items() if t == typ and any(x["appropriate_control_context_available"] for x in rs)})
        for typ in sorted(nodes_by_type)
    })
    print("GenX targets by perturbation type with at least one profile having >=2 replicates:", {
        typ: len({g for (g,t), rs in gene_type_profiles.items() if t == typ and any(str(x["replicate_count"]).isdigit() and int(x["replicate_count"]) >= 2 for x in rs)})
        for typ in sorted(nodes_by_type)
    })

    # Full compendium target inventory, separate from the MCF7 profile counts.
    global_targets = defaultdict(set)
    for r in rows(META / "GSE92742_Broad_LINCS_pert_info.txt.gz", compressed=True):
        if r.get("pert_type") not in genetic_types:
            continue
        h, _, st = map_hgnc(r.get("pert_iname", ""))
        if not h or st not in ("unique", "unique_alias"):
            continue
        matched = [n for n, (gh, _, gs) in graph_map.items()
                   if gh and gs in ("unique", "unique_alias") and gh["hgnc_id"] == h["hgnc_id"]]
        if len(matched) == 1:
            global_targets[r["pert_type"]].add(matched[0])
    global_sh = set().union(*(global_targets[x] for x in knock_types))
    print("All-cell-line GenX perturbation targets by type:", {k: len(v) for k,v in sorted(global_targets.items())})
    print("All-cell-line unique GenX genetic perturbation targets:", len(set().union(*global_targets.values())))
    print("All-cell-line unique GenX knockdown targets:", len(global_sh))
    print("All-cell-line unique GenX cDNA overexpression targets:", len(global_targets['trt_oe']))
    print("All-cell-line GenX targets with both knockdown and cDNA overexpression:", len(global_sh & global_targets['trt_oe']))
    print("Inferred non-landmark Level 5 feature IDs overlapping mapped GenX nodes:", inferred_overlap_count)

    # Inspect the earlier GSE70138 profile metadata independently. Do not merge
    # it into the pinned GSE92742 coverage population.
    early_sig = META / "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
    if early_sig.exists():
        early_cells = defaultdict(int)
        early_mcf7_types = defaultdict(int)
        for r in rows(early_sig, compressed=True):
            early_cells[r.get("cell_id", "")] += 1
            if r.get("cell_id") == "MCF7":
                early_mcf7_types[r.get("pert_type", "")] += 1
        print("GSE70138 MCF7 profile rows by perturbation type:", dict(sorted(early_mcf7_types.items())))
        print("GSE70138 has MCF7 genetic knockdown/overexpression profiles:", sum(early_mcf7_types[t] for t in genetic_types))

if __name__ == "__main__":
    main()
