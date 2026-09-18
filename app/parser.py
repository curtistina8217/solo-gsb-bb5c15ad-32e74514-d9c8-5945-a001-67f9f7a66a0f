"""Molecular-formula parser with isotope-label support and position-aware errors.

Grammar (whitespace is not allowed)::

    formula  := term+
    term     := element count? | label count?
    element  := [A-Z] [a-z]?
    label    := "[" mass-number element "]"        e.g. [13C], [15N], [2H]
    count    := [0-9]+                            (must be >= 1 when present)

Examples: ``C6H12O6``, ``C6H12O6[13C]2``, ``C2H5[2H]1O``.  A label states how
many of that element's atoms are isotopically labelled, so ``C6[13C]2`` means
six carbon atoms in total, two of which are carbon-13.

Every validation failure raises :class:`FormulaError` carrying the character
position (0-based) of the offending token so clients can point at the exact
location of the problem.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


class FormulaError(ValueError):
    """Validation error carrying the position of the offending token."""

    def __init__(self, message: str, position: int):
        self.message = message
        self.position = position
        super().__init__(f"{message} (at position {position})")


@dataclass
class ParsedFormula:
    # element symbol -> total atom count (labelled atoms included)
    composition: dict[str, int] = field(default_factory=dict)    # (element, mass_number) -> labelled atom count
    labels: dict[tuple[str, int], int] = field(default_factory=dict)
    # (element, mass_number) -> position of first occurrence, for error reports
    label_positions: dict[tuple[str, int], int] = field(default_factory=dict)

    def canonical_key(self) -> str:
        """Stable cache key independent of input ordering."""
        comp = "".join(f"{el}{n}" for el, n in sorted(self.composition.items()))
        lab = ",".join(
            f"{mass}{el}:{n}" for (el, mass), n in sorted(self.labels.items())
        )
        return f"{comp}|{lab}"


_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")
_LABEL_RE = re.compile(r"\[(\d+)([A-Z][a-z]?)\]")
_COUNT_RE = re.compile(r"\d+")

MAX_ATOMS = 10_000  # guard against pathological inputs


def _read_count(formula: str, i: int) -> tuple[int, int]:
    """Read an optional count starting at *i*; returns (count, next_index)."""
    m = _COUNT_RE.match(formula, i)
    if not m:
        return 1, i
    count = int(m.group(0))
    if count < 1:
        raise FormulaError("atom count must be at least 1", i)
    return count, m.end()


def parse_formula(
    formula: str,
    valid_elements: set[str],
    isotope_masses: dict[str, set[int]],
) -> ParsedFormula:
    """Parse *formula* into a :class:`ParsedFormula`.

    ``valid_elements`` is the set of supported element symbols and
    ``isotope_masses`` maps each element to the mass numbers of its known
    isotopes (used to validate labels such as ``[13C]``).
    """
    if not formula:
        raise FormulaError("formula is empty", 0)

    parsed = ParsedFormula()
    i, n = 0, len(formula)
    while i < n:
        ch = formula[i]
        if ch == "[":
            m = _LABEL_RE.match(formula, i)
            if not m:
                raise FormulaError(
                    "malformed isotope label, expected e.g. [13C]", i
                )
            mass_number = int(m.group(1))
            element = m.group(2)
            if element not in valid_elements:
                raise FormulaError(
                    f"unknown element symbol '{element}' in label", i + 1 + len(m.group(1))
                )
            if mass_number not in isotope_masses.get(element, set()):
                known = ", ".join(str(x) for x in sorted(isotope_masses.get(element, ())))
                raise FormulaError(
                    f"unknown isotope {mass_number}{element} (known: {known})", i
                )
            count, i = _read_count(formula, m.end())
            key = (element, mass_number)
            parsed.labels[key] = parsed.labels.get(key, 0) + count
            parsed.label_positions.setdefault(key, m.start())
        elif ch.isdigit():
            raise FormulaError("count without a preceding element or label", i)
        elif ch.isupper():
            m = _ELEMENT_RE.match(formula, i)
            element = m.group(0)
            if element not in valid_elements:
                raise FormulaError(f"unknown element symbol '{element}'", i)
            count, i = _read_count(formula, m.end())
            parsed.composition[element] = parsed.composition.get(element, 0) + count
        else:
            raise FormulaError(f"unexpected character '{ch}'", i)

    total = sum(parsed.composition.values())
    if total > MAX_ATOMS:
        raise FormulaError(
            f"formula too large ({total} atoms, limit {MAX_ATOMS})", 0
        )

    for (element, _mass), count in parsed.labels.items():
        if count > parsed.composition[element]:
            raise FormulaError(
                f"label count ({count}) exceeds total {element} atom count "
                f"({parsed.composition[element]})",
                parsed.label_positions[(element, _mass)],
            )
    return parsed


_LABEL_BODY_RE = re.compile(r"^(\d+)([A-Z][a-z]?)$")


def parse_label_string(
    label: str,
    index: int,
    valid_elements: set[str],
    isotope_masses: dict[str, set[int]],
) -> tuple[str, int]:
    """Validate a structured label string such as ``"13C"``.

    ``index`` is the label's position in the request's ``labels`` list and is
    reported as the error position.
    """
    m = _LABEL_BODY_RE.match(label or "")
    if not m:
        raise FormulaError(
            f"malformed isotope label '{label}', expected e.g. '13C'", index
        )
    mass_number, element = int(m.group(1)), m.group(2)
    if element not in valid_elements:
        raise FormulaError(f"unknown element symbol '{element}' in label", index)
    if mass_number not in isotope_masses.get(element, set()):
        known = ", ".join(str(x) for x in sorted(isotope_masses.get(element, ())))
        raise FormulaError(
            f"unknown isotope {mass_number}{element} (known: {known})", index
        )
    return element, mass_number
