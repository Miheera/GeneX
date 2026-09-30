from __future__ import annotations
from typing import List, Literal, Dict
from pydantic import BaseModel, Field, field_validator

class NodeModel(BaseModel):
    id: str
    label: str
    degree: int
    disease_score: float

class EdgeModel(BaseModel):
    source: str
    target: str
    weight: float

class GraphResponse(BaseModel):
    nodes: List[NodeModel]
    edges: List[EdgeModel]

class PerturbRequest(BaseModel):
    gene_id: str = Field(..., description="Gene symbol to perturb")
    perturbation_type: Literal["knockout", "inhibition", "overexpression"]
    strength: float = Field(default=0.5, ge=0.0, le=1.0)

    @field_validator("gene_id")
    @classmethod
    def gene_id_must_not_be_empty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("gene_id must not be empty")
        return v.strip()

class ImpactedGene(BaseModel):
    gene: str
    delta_bc: float
    direction: Literal["gain", "loss"]

class PerturbResponse(BaseModel):
    experiment_type: str = "virtual_topological_knockout"
    efficiency_baseline: float
    efficiency_perturbed: float
    efficiency_delta_pct: float
    impact_score: float
    perturbation_type: str
    gene_id: str
    top_impacted: List[ImpactedGene]
    network_fragmented: bool
    connected_components: int
    lcc_size: int
    lcc_fraction: float
    isolated_genes: List[str]
    bc_baseline: Dict[str, float]
    bc_perturbed: Dict[str, float]
    approximation_method: str

class PropagatedGene(BaseModel):
    gene: str
    rank: int
    propagation_score: float
    absolute_score: float

class PropagationResponse(BaseModel):
    target_gene: str
    perturbation_type: Literal["inhibition", "overexpression"]
    perturbation_strength: float
    propagation_model: str
    alpha: float
    propagation_scores: Dict[str, float]
    top_impacted_genes: List[PropagatedGene]
    iterations: int
    converged: bool
    ranking_strength_invariant: bool
    interpretation: str
