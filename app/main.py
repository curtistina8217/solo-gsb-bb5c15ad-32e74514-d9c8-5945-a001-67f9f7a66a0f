"""FastAPI application: mass-spectrometry isotope-distribution precomputation."""
from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from . import calculator, db
from .parser import FormulaError, parse_formula, parse_label_string
from .schemas import (
    DistributionRequest,
    DistributionResponse,
    ErrorDetail,
    Peak,
)


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Dependencies (isotope reference data) and DB schema are initialized here,
    # so starting the service is the only step required.
    db.init_db()
    app.state.isotopes = db.load_isotopes()
    yield


app = FastAPI(
    title="MS Precompute API",
    version="1.0.0",
    description=(
        "Computes theoretical isotope distributions (mass-intensity polynomials "
        "convolved element by element) and monoisotopic masses for molecular "
        "formulas with optional isotope labels."
    ),
    lifespan=lifespan,
)


@app.exception_handler(FormulaError)
async def formula_error_handler(_request: Request, exc: FormulaError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content=ErrorDetail(
            error="invalid_formula", message=exc.message, position=exc.position
        ).model_dump(),
    )


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.get("/api/v1/elements")
def list_elements(request: Request) -> dict:
    """List supported elements and their isotopes (mass, natural abundance)."""
    return {"elements": request.app.state.isotopes}


@app.post(
    "/api/v1/isotope-distribution",
    response_model=DistributionResponse,
    responses={422: {"model": ErrorDetail}},
)
def isotope_distribution(payload: DistributionRequest, request: Request):
    """Compute the isotope distribution and monoisotopic mass of a formula."""
    isotopes: dict[str, list[dict]] = request.app.state.isotopes
    valid_elements = set(isotopes)
    isotope_masses = {
        el: {iso["mass_number"] for iso in lst} for el, lst in isotopes.items()
    }

    parsed = parse_formula(payload.formula, valid_elements, isotope_masses)

    # Merge structured labels from the request body with embedded ones.
    labels = dict(parsed.labels)
    label_positions = dict(parsed.label_positions)
    for index, label in enumerate(payload.labels):
        element, mass_number = parse_label_string(
            label.isotope, index, valid_elements, isotope_masses
        )
        key = (element, mass_number)
        labels[key] = labels.get(key, 0) + label.count
        label_positions.setdefault(key, index)
    for (element, _mass), count in labels.items():
        total = parsed.composition.get(element, 0)
        if count > total:
            raise FormulaError(
                f"label count ({count}) exceeds total {element} atom count ({total})",
                label_positions[(element, _mass)],
            )
    parsed.labels = labels

    cache_key = parsed.canonical_key()
    cached = db.get_cached(cache_key)
    if cached is not None:
        return DistributionResponse(
            formula=payload.formula,
            monoisotopic_mass=cached["monoisotopic_mass"],
            peaks=[Peak(**p) for p in cached["peaks"][: payload.top_n]],
            total_peaks=len(cached["peaks"]),
            cached=True,
        )

    # Build per-element polynomial terms; labelled atoms contribute a fixed
    # mass offset (their isotope is treated as pure, abundance 1.0).
    element_terms = []
    mass_offset = 0.0
    monoisotopic_mass = 0.0
    for element, total_count in parsed.composition.items():
        iso_list = isotopes[element]
        labelled = sum(
            count for (el, _m), count in parsed.labels.items() if el == element
        )
        natural_count = total_count - labelled
        if natural_count > 0:
            element_terms.append(
                ([(iso["mass"], iso["abundance"]) for iso in iso_list], natural_count)
            )
            monoisotopic_mass += natural_count * min(i["mass"] for i in iso_list)
        for (el, mass_number), count in parsed.labels.items():
            if el != element:
                continue
            label_mass = next(
                i["mass"] for i in iso_list if i["mass_number"] == mass_number
            )
            mass_offset += count * label_mass
            monoisotopic_mass += count * label_mass

    peaks = calculator.compute_distribution(element_terms)
    peaks = calculator.normalize(peaks, base=100.0)

    stored = [
        {"mass": round(m + mass_offset, 6), "intensity": intensity}
        for m, i in sorted(peaks, key=lambda p: p[1], reverse=True)
        if (intensity := round(i, 6)) > 0.0
    ]
    monoisotopic_mass = round(monoisotopic_mass, 6)
    db.store_calculation(cache_key, payload.formula, monoisotopic_mass, stored)

    return DistributionResponse(
        formula=payload.formula,
        monoisotopic_mass=monoisotopic_mass,
        peaks=[Peak(**p) for p in stored[: payload.top_n]],
        total_peaks=len(stored),
        cached=False,
    )
