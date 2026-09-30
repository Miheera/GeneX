"""
main.py — FastAPI application for Cancer Gene Network Explorer.

Startup sequence:
  1. Load graph from CSV (graph_loader.init_graph)
  2. Precompute baseline metrics (perturbation.precompute_baseline)
  3. Serve requests

Routes:
  GET  /api/graph/{cancer_type}  → GraphResponse
  POST /api/perturb              → PerturbResponse
"""
from __future__ import annotations

import csv
import logging
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

import graph_loader
import perturbation as perturb_module
from analysis.attack_tolerance import analyze_attack_tolerance
from analysis.expression_correlation import analyze_structural_expression
from analysis.open_targets import open_targets
from models import (
    EdgeModel,
    GraphResponse,
    NodeModel,
    PerturbRequest,
    PerturbResponse,
    PropagationResponse,
)
from propagation import compute_perturbation_propagation

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
)
logger = logging.getLogger(__name__)


# ── Lifespan (replaces deprecated @app.on_event) ─────────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    # ── Startup ──
    logger.info("=== Cancer Gene Network Explorer — startup ===")
    graph_loader.init_graph()                          # loads CSV → G
    perturb_module.precompute_baseline(graph_loader.G) # cache baseline metrics
    logger.info("=== Startup complete — ready to serve requests ===")
    yield
    # ── Shutdown ──
    logger.info("=== Shutting down ===")


# ── App ───────────────────────────────────────────────────────────────────────
app = FastAPI(
    title="Cancer Gene Network Explorer API",
    description=(
        "FastAPI backend for the Cancer Gene Network Explorer research dashboard. "
        "Provides gene interaction network data and perturbation simulation results."
    ),
    version="1.0.0",
    lifespan=lifespan,
)

# ── CORS ──────────────────────────────────────────────────────────────────────
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",   # React / Next.js dev server
        "http://localhost:5500",   # VS Code Live Server
        "http://127.0.0.1:5500",  # Live Server alternate
        "http://localhost:5173",   # Vite dev server
        "http://127.0.0.1:3000",
        "http://localhost:8000",   # Same-origin (static files served here)
        "http://127.0.0.1:8000",
        "null",                    # file:// origin (local dev fallback)
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)




# ─────────────────────────────────────────────────────────────────────────────
# ENDPOINT 1 — GET /api/graph/{cancer_type}
# ─────────────────────────────────────────────────────────────────────────────
@app.get(
    "/api/graph/{cancer_type}",
    response_model=GraphResponse,
    summary="Fetch the full gene interaction graph",
    tags=["Graph"],
)
async def get_graph(cancer_type: str) -> GraphResponse:
    """
    Return the full graph as a node-link JSON suitable for D3 rendering.

    Currently only **breast_cancer** is supported.  Other cancer types will
    return a 404 until additional CSV files are added.
    """
    supported = {"breast_cancer"}
    if cancer_type not in supported:
        raise HTTPException(
            status_code=404,
            detail=(
                f"Cancer type '{cancer_type}' is not supported. "
                f"Supported types: {sorted(supported)}"
            ),
        )

    G = graph_loader.G
    if G.number_of_nodes() == 0:
        raise HTTPException(status_code=503, detail="Graph data is not available.")

    nodes = [
        NodeModel(
            id=node,
            label=attrs.get("label", node),
            degree=attrs.get("degree", G.degree(node)),
            disease_score=attrs.get("disease_score", 0.0),
        )
        for node, attrs in G.nodes(data=True)
    ]

    edges = [
        EdgeModel(
            source=str(u),
            target=str(v),
            weight=round(data.get("weight", 0.0), 4),
        )
        for u, v, data in G.edges(data=True)
    ]

    return GraphResponse(nodes=nodes, edges=edges)


# ─────────────────────────────────────────────────────────────────────────────
# ENDPOINT 2 — POST /api/perturb
# ─────────────────────────────────────────────────────────────────────────────
@app.post(
    "/api/perturb",
    response_model=PerturbResponse | PropagationResponse,
    summary="Run a gene perturbation simulation",
    tags=["Perturbation"],
)
async def perturb_gene(request: PerturbRequest) -> PerturbResponse | PropagationResponse:
    """
    Run a virtual topological knockout or fixed-graph confidence-weighted RWR.
    """
    G = graph_loader.G

    # ── 404 if gene not in graph ──
    if request.gene_id not in G:
        raise HTTPException(
            status_code=404,
            detail=f"Gene not found: {request.gene_id}",
        )

    # ── Run perturbation ──
    try:
        if request.perturbation_type == "knockout":
            result = perturb_module.run_perturbation(
                G=G,
                gene_id=request.gene_id,
                perturbation_type="knockout",
                strength=request.strength,
            )
            return PerturbResponse(**result)
        result = compute_perturbation_propagation(
            graph=G,
            target_node=request.gene_id,
            perturbation_type=request.perturbation_type,
            strength=request.strength,
        )
    except Exception as exc:
        logger.exception("Perturbation computation failed for gene '%s'", request.gene_id)
        raise HTTPException(
            status_code=500,
            detail=f"Computation error: {str(exc)}",
        )

    return PropagationResponse(**result)


@app.get("/api/analysis/attack-tolerance", tags=["Analysis"])
async def attack_tolerance_analysis(
    strategy: Literal["degree", "betweenness", "random"],
):
    """Run attack-tolerance analysis on a copy of the loaded graph."""
    try:
        return analyze_attack_tolerance(
            graph_loader.G,
            strategy=strategy,
            n_random_trials=50,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Attack-tolerance analysis failed")
        raise HTTPException(status_code=500, detail=f"Analysis error: {str(exc)}") from exc


@app.get("/api/analysis/threshold-sensitivity", tags=["Analysis"])
async def threshold_sensitivity():
    """Return the pre-generated threshold sensitivity report."""
    report_path = Path(__file__).parent / "data" / "analysis_outputs" / "threshold_sensitivity.csv"
    if not report_path.exists():
        raise HTTPException(status_code=404, detail="Threshold sensitivity report is not available.")
    with report_path.open(newline="", encoding="utf-8") as report_file:
        return list(csv.DictReader(report_file))


@app.get("/api/analysis/cosmic-benchmark", tags=["Analysis"])
async def cosmic_benchmark():
    """Return the pre-generated COSMIC benchmark report."""
    report_path = Path(__file__).parent / "data" / "analysis_outputs" / "cosmic_benchmark.csv"
    if not report_path.exists():
        raise HTTPException(
            status_code=404,
            detail="COSMIC benchmark is not available. Run analysis/benchmark_cosmic.py first.",
        )
    with report_path.open(newline="", encoding="utf-8") as report_file:
        rows = []
        for row in csv.DictReader(report_file):
            row["top_n"] = int(row["top_n"])
            for key in ("precision", "recall", "f1"):
                row[key] = float(row[key])
            for key in ("n_positives_in_graph", "n_positives_total", "n_positives_excluded_from_graph"):
                row[key] = int(row[key])
            rows.append(row)
        return rows


@app.get("/api/analysis/open-targets", tags=["Analysis"])
async def open_target_analysis(top_n: int = 20):
    """Return high-impact genes split by DGIdb interaction-claim coverage."""
    if top_n <= 0:
        raise HTTPException(status_code=422, detail="top_n must be a positive integer")
    annotations_path = Path(__file__).parent / "data" / "processed" / "dgidb_gene_annotations.csv"
    try:
        return open_targets(graph_loader.G, annotations_path, top_n)
    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=404,
            detail="DGIdb annotations are not available. Run pipeline/filter_dgidb.py first.",
        ) from exc
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Open-target analysis failed")
        raise HTTPException(status_code=500, detail=f"Analysis error: {str(exc)}") from exc


@app.get("/api/analysis/structural-expression", tags=["Analysis"])
async def structural_expression_analysis():
    """Compare knockout network impact with TCGA BRCA median expression z-scores."""
    expression_path = Path(__file__).parent / "data" / "processed" / "tcga_expression_scores.csv"
    try:
        return analyze_structural_expression(graph_loader.G, expression_path)
    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=404,
            detail="TCGA expression scores are not available. Run pipeline/fetch_tcga_expression.py first.",
        ) from exc
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Structural-expression analysis failed")
        raise HTTPException(status_code=500, detail=f"Analysis error: {str(exc)}") from exc


# ─────────────────────────────────────────────────────────────────────────────
# HEALTH CHECK
# ─────────────────────────────────────────────────────────────────────────────
@app.get("/health", tags=["Utility"])
async def health_check() -> dict:
    G = graph_loader.G
    return {
        "status": "ok",
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "baseline_efficiency": round(perturb_module._baseline_efficiency, 6),
        "approximation_method": perturb_module.APPROXIMATION_METHOD,
    }


# ───────────────────────────────────────────────────────────────────────────────
# STATIC FRONTEND  — MUST be last so /api/* routes match first
# ───────────────────────────────────────────────────────────────────────────────
# Serve index.html / app.js / style.css from Major Project/ at /
_frontend_dir = Path(__file__).parents[1]  # Major Project/
app.mount("/", StaticFiles(directory=str(_frontend_dir), html=True), name="frontend")


# ───────────────────────────────────────────────────────────────────────────────
# DEV ENTRY POINT
# ───────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import uvicorn

    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
