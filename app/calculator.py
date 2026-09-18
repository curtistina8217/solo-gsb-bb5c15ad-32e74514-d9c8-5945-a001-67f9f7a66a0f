"""Isotope-distribution engine.

The isotope pattern of each element is treated as a mass-intensity polynomial
(sum of ``abundance * x**mass`` terms).  Element polynomials are exponentiated
to the atom count (binary powering) and convolved element by element.  Peaks
closer than ``MERGE_TOLERANCE`` Da are merged after every step, and the
working distribution is pruned so computation stays bounded.
"""
from __future__ import annotations

MERGE_TOLERANCE = 1e-4   # Da; peaks closer than this are merged
REL_THRESHOLD = 1e-9     # drop peaks below this fraction of the current max
MAX_PEAKS = 2000         # cap on peaks kept per intermediate distribution

Peak = tuple[float, float]  # (mass, intensity)


def merge_peaks(peaks: list[Peak], tol: float = MERGE_TOLERANCE) -> list[Peak]:
    """Merge peaks whose mass difference is below *tol* Da.

    Merged peaks keep the intensity-weighted mean mass and the summed
    intensity.  Input need not be sorted.
    """
    merged: list[list[float]] = []
    for mass, intensity in sorted(peaks):
        if merged and mass - merged[-1][0] < tol:
            m0, i0 = merged[-1]
            total = i0 + intensity
            merged[-1] = [(m0 * i0 + mass * intensity) / total, total]
        else:
            merged.append([mass, intensity])
    return [(m, i) for m, i in merged]


def convolve(a: list[Peak], b: list[Peak]) -> list[Peak]:
    """Multiply two mass-intensity polynomials (distribution convolution)."""
    out: list[Peak] = []
    for m1, i1 in a:
        for m2, i2 in b:
            out.append((m1 + m2, i1 * i2))
    return merge_peaks(out)


def prune(peaks: list[Peak]) -> list[Peak]:
    """Drop negligible peaks and cap the distribution size."""
    if not peaks:
        return peaks
    top = max(i for _, i in peaks)
    peaks = [p for p in peaks if p[1] >= top * REL_THRESHOLD]
    if len(peaks) > MAX_PEAKS:
        peaks = sorted(peaks, key=lambda p: p[1], reverse=True)[:MAX_PEAKS]
        peaks.sort()
    return peaks


def power_distribution(isotopes: list[Peak], count: int) -> list[Peak]:
    """Raise a single-atom distribution to *count* atoms (binary powering)."""
    result: list[Peak] = [(0.0, 1.0)]
    base = list(isotopes)
    k = count
    while k:
        if k & 1:
            result = prune(convolve(result, base))
        k >>= 1
        if k:
            base = prune(convolve(base, base))
    return result


def compute_distribution(
    element_terms: list[tuple[list[Peak], int]],
) -> list[Peak]:
    """Convolve per-element distributions into the molecular distribution.

    ``element_terms`` is a list of ``(isotope_distribution, atom_count)``.
    Returns merged, pruned peaks sorted by mass.
    """
    dist: list[Peak] = [(0.0, 1.0)]
    for isotopes, count in element_terms:
        if count <= 0:
            continue
        element_dist = power_distribution(isotopes, count)
        dist = prune(convolve(dist, element_dist))
    return merge_peaks(dist)


def normalize(peaks: list[Peak], base: float = 100.0) -> list[Peak]:
    """Normalize so the base (most intense) peak equals *base*."""
    if not peaks:
        return peaks
    top = max(i for _, i in peaks)
    if top <= 0:
        return peaks
    factor = base / top
    return [(m, i * factor) for m, i in peaks]
