"""Pydantic request/response models."""
from __future__ import annotations

from pydantic import BaseModel, Field


class Label(BaseModel):
    """An isotope label, e.g. ``{"isotope": "13C", "count": 2}``."""

    isotope: str = Field(..., examples=["13C"], description="Isotope, e.g. 13C, 15N, 2H, 18O")
    count: int = Field(1, ge=1, description="Number of labelled atoms")


class DistributionRequest(BaseModel):
    formula: str = Field(
        ...,
        min_length=1,
        examples=["C6H12O6", "C6H12O6[13C]2"],
        description="Molecular formula; labels may be embedded as [13C]2 etc.",
    )
    labels: list[Label] = Field(
        default_factory=list,
        description="Optional structured labels, merged with embedded ones",
    )
    top_n: int = Field(10, ge=1, le=1000, description="Number of peaks to return")


class Peak(BaseModel):
    mass: float = Field(..., description="Peak mass in Da")
    intensity: float = Field(..., description="Intensity, base peak = 100")


class DistributionResponse(BaseModel):
    formula: str
    monoisotopic_mass: float = Field(..., description="Monoisotopic mass in Da")
    peaks: list[Peak]
    total_peaks: int = Field(..., description="Peaks computed before top_n truncation")
    cached: bool = Field(..., description="True if served from the precomputed cache")


class ErrorDetail(BaseModel):
    error: str
    message: str
    position: int = Field(..., description="Character index in `formula`, or index in `labels`")
