"""Stdlib-only sanity checks for the core engine (no FastAPI required).

Run: python3 tests/test_core.py
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app import calculator, db  # noqa: E402
from app.parser import FormulaError, parse_formula, parse_label_string  # noqa: E402

# Use a throwaway database for the tests.
_tmp = tempfile.mkdtemp()
db.DATABASE_PATH = os.path.join(_tmp, "test.db")
db.init_db()
ISOTOPES = db.load_isotopes()
VALID = set(ISOTOPES)
ISO_MASSES = {el: {i["mass_number"] for i in lst} for el, lst in ISOTOPES.items()}

failures = []


def check(name, cond, extra=""):
    print(("PASS" if cond else "FAIL"), name, extra)
    if not cond:
        failures.append(name)


# --- parser ---------------------------------------------------------------
p = parse_formula("C6H12O6", VALID, ISO_MASSES)
check("parse plain formula", p.composition == {"C": 6, "H": 12, "O": 6})

p = parse_formula("C6H12O6[13C]2", VALID, ISO_MASSES)
check("parse embedded label", p.labels == {("C", 13): 2})
check("composition holds total counts", p.composition == {"C": 6, "H": 12, "O": 6})

p = parse_formula("CH3CH3", VALID, ISO_MASSES)
check("repeated elements accumulate", p.composition == {"C": 2, "H": 6})


def expect_error(name, formula, position):
    try:
        parse_formula(formula, VALID, ISO_MASSES)
        check(name, False, "no error raised")
    except FormulaError as e:
        check(name, e.position == position, f"pos={e.position} want={position} msg={e.message}")


expect_error("unknown element position", "C6Xy", 2)
expect_error("lowercase start position", "c6H12", 0)
expect_error("stray count position", "12C", 0)
expect_error("bad label position", "C6[C13]", 2)
expect_error("unknown isotope position", "C6[99C]2", 2)
expect_error("zero count position", "C0H2", 1)
expect_error("label exceeds atoms", "C2[13C]3", 2)

el, mn = parse_label_string("13C", 0, VALID, ISO_MASSES)
check("body label ok", (el, mn) == ("C", 13))
try:
    parse_label_string("99C", 3, VALID, ISO_MASSES)
    check("body label bad isotope", False)
except FormulaError as e:
    check("body label bad isotope", e.position == 3)

# --- calculator -----------------------------------------------------------
# Single 13C atom: distribution must equal the isotope table exactly.
c_isotopes = [(i["mass"], i["abundance"]) for i in ISOTOPES["C"]]
dist = calculator.power_distribution(c_isotopes, 1)
check("single atom distribution", dist == calculator.merge_peaks(c_isotopes))

# C2: 12C2 / 12C13C / 13C2 with binomial intensities.
dist = calculator.power_distribution(c_isotopes, 2)
intens = sorted((round(i, 6) for _, i in dist), reverse=True)
expect = sorted(
    [0.9893**2, 2 * 0.9893 * 0.0107, 0.0107**2], reverse=True
)
check(
    "C2 binomial intensities",
    all(abs(a - b) < 1e-6 for a, b in zip(intens, expect)),
    f"{intens} vs {expect}",
)

# Merge tolerance: peaks 5e-5 apart merge, 5e-4 apart do not.
merged = calculator.merge_peaks([(10.0, 1.0), (10.00005, 1.0), (20.0, 1.0), (20.0005, 1.0)])
check("merge below 0.0001 Da", len(merged) == 3 and abs(merged[0][0] - 10.000025) < 1e-9)

# Glucose C6H12O6: monoisotopic mass and M+1 shift.
terms = [
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["C"]], 6),
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["H"]], 12),
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["O"]], 6),
]
dist = calculator.normalize(calculator.compute_distribution(terms))
mono = 6 * 12.0 + 12 * 1.00782503223 + 6 * 15.99491461957
check("glucose monoisotopic", abs(dist[0][0] - mono) < 1e-4, f"{dist[0][0]:.6f} vs {mono:.6f}")
check("base peak normalized to 100", abs(max(i for _, i in dist) - 100.0) < 1e-9)
m_plus_1 = min(dist, key=lambda p: abs(p[0] - (mono + 1.00335483507)))
check("glucose M+1 position", abs(m_plus_1[0] - (mono + 1.00335483507)) < 1e-3)
check("glucose M+1 intensity ~6.5%", abs(m_plus_1[1] - 6.49) < 0.3, f"{m_plus_1[1]:.2f}%")

# Larger molecule: bovine insulin C254H377N65O75S6 — known mono mass 5729.6009.
terms = [
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["C"]], 254),
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["H"]], 377),
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["N"]], 65),
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["O"]], 75),
    ([(i["mass"], i["abundance"]) for i in ISOTOPES["S"]], 6),
]
dist = calculator.normalize(calculator.compute_distribution(terms))
mono_ins = dist[0][0]
check("insulin monoisotopic ~5729.6009", abs(mono_ins - 5729.6009) < 5e-4, f"{mono_ins:.4f}")
# For insulin the base peak is no longer the monoisotopic peak.
check("insulin base peak not monoisotopic", dist[0][1] < 100.0)

# --- db cache -------------------------------------------------------------
db.store_calculation("C6H12O6|", "C6H12O6", 180.063388, [{"mass": 180.063388, "intensity": 100.0}])
cached = db.get_cached("C6H12O6|")
check("cache round-trip", cached and cached["monoisotopic_mass"] == 180.063388)
check("cache miss", db.get_cached("nope") is None)

print()
if failures:
    print("FAILED:", ", ".join(failures))
    sys.exit(1)
print(f"All checks passed ({db.DATABASE_PATH})")
