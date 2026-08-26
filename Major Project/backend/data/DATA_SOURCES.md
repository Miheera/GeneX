# Data Sources

This file documents how to re-download the large raw data files that are **excluded from git**
(see `/.gitignore`). All other data (cosmic, disgenet, processed outputs) are tracked normally.

---

## 1. `major_project_dataset/` (~611 MB) — root-level legacy copy

This is a legacy/root-level mirror of the STRING and DisGeNET source files kept for convenience.
It is **not** used directly by the pipeline (the pipeline reads from `data/raw/` below).

| Filename | Source |
|----------|--------|
| `9606.protein.info.v12.0.txt` | STRING DB v12.0 |
| `9606.protein.links.v12.0.txt` | STRING DB v12.0 |
| `curated_gene_disease_associations.txt` | DisGeNET |

**Download instructions:**

```bash
# STRING DB v12.0 — human protein info and links
# Navigate to: https://string-db.org/cgi/download?species_text=Homo+sapiens
# Download "protein.info" and "protein.links" for organism 9606 (Homo sapiens), version 12.0.
# Place the unzipped .txt files into:  major_project_dataset/<filename>/<filename>

# DisGeNET curated associations
# Navigate to: https://www.disgenet.org/downloads
# Download "Curated gene-disease associations" (curated_gene_disease_associations.tsv or .txt).
# Place into: major_project_dataset/curated_gene_disease_associations.txt/curated_gene_disease_associations.txt
```

---

## 2. `backend/data/raw/string/` (~95 MB) — pipeline input files

These are the compressed STRING DB files consumed directly by the data pipeline scripts
(`pipeline/filter_string.py`, `pipeline/map_identifiers.py`, etc.).

| Filename | Version | Note |
|----------|---------|------|
| `9606.protein.info.v12.0.txt.gz` | STRING v12.0 | Protein metadata |
| `9606.protein.links.v12.0.txt.gz` | STRING v12.0 | Protein interaction links |
| `9606.protein.actions.v11.0.txt.gz` | **STRING v11.0** | Action-type annotations — v12.0 does **not** publish this file; v11.0 is the latest available |

**Download instructions:**

```bash
# Go to: https://string-db.org/cgi/download?species_text=Homo+sapiens
# Select organism: Homo sapiens (taxon ID 9606)
#
# For v12.0 files (protein.info and protein.links):
#   Download the .txt.gz files directly — no extraction needed; the pipeline reads gzipped input.
#
# For protein.actions (v11.0 only):
#   Switch the version selector to v11.0 on the same download page, or use:
#   https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz
#
# Place all three .gz files into:  Major Project/backend/data/raw/string/
```

---

## Tracked raw data (no action needed)

The following small raw data files **are** committed to git and require no manual download:

| Path | Size | Source |
|------|------|--------|
| `data/raw/cosmic/cancer_gene_census.csv` | ~185 KB | [COSMIC](https://cancer.sanger.ac.uk/census) |
| `data/raw/disgenet/search_C0006142.*` | ~32 KB | [DisGeNET](https://www.disgenet.org) |
| `data/raw/disgenet/search_C0678222.*` | ~32 KB | [DisGeNET](https://www.disgenet.org) |
