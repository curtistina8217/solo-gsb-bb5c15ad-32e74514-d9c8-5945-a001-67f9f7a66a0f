"""SQLite persistence: isotope reference table and calculation cache.

The database is created and seeded from :mod:`app.isotopes` at application
startup (see the lifespan handler in :mod:`app.main`), so no manual
initialization step is required.
"""
from __future__ import annotations

import json
import os
import sqlite3
import threading

from .isotopes import ISOTOPES

DATABASE_PATH = os.environ.get("DATABASE_PATH", "msprecompute.db")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS isotopes (
    element     TEXT    NOT NULL,
    mass_number INTEGER NOT NULL,
    mass        REAL    NOT NULL,
    abundance   REAL    NOT NULL,
    PRIMARY KEY (element, mass_number)
);
CREATE TABLE IF NOT EXISTS calculations (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    cache_key         TEXT NOT NULL UNIQUE,
    formula           TEXT NOT NULL,
    monoisotopic_mass REAL NOT NULL,
    peaks_json        TEXT NOT NULL,
    created_at        TEXT NOT NULL DEFAULT (datetime('now'))
);
"""

_lock = threading.Lock()


def _connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DATABASE_PATH, timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    return conn


def init_db() -> None:
    """Create tables and seed the isotope reference data if empty."""
    directory = os.path.dirname(os.path.abspath(DATABASE_PATH))
    os.makedirs(directory, exist_ok=True)
    with _lock, _connect() as conn:
        conn.executescript(_SCHEMA)
        (count,) = conn.execute("SELECT COUNT(*) FROM isotopes").fetchone()
        if count == 0:
            conn.executemany(
                "INSERT INTO isotopes (element, mass_number, mass, abundance)"
                " VALUES (?, ?, ?, ?)",
                ISOTOPES,
            )


def load_isotopes() -> dict[str, list[dict]]:
    """Return ``{element: [{mass_number, mass, abundance}, ...]}`` from SQLite."""
    with _lock, _connect() as conn:
        rows = conn.execute(
            "SELECT element, mass_number, mass, abundance"
            " FROM isotopes ORDER BY element, mass_number"
        ).fetchall()
    table: dict[str, list[dict]] = {}
    for row in rows:
        table.setdefault(row["element"], []).append(
            {
                "mass_number": row["mass_number"],
                "mass": row["mass"],
                "abundance": row["abundance"],
            }
        )
    return table


def get_cached(cache_key: str) -> dict | None:
    """Return the cached calculation for *cache_key*, or ``None``."""
    with _lock, _connect() as conn:
        row = conn.execute(
            "SELECT formula, monoisotopic_mass, peaks_json"
            " FROM calculations WHERE cache_key = ?",
            (cache_key,),
        ).fetchone()
    if row is None:
        return None
    return {
        "formula": row["formula"],
        "monoisotopic_mass": row["monoisotopic_mass"],
        "peaks": json.loads(row["peaks_json"]),
    }


def store_calculation(
    cache_key: str, formula: str, monoisotopic_mass: float, peaks: list[dict]
) -> None:
    """Persist a computed distribution for future (precomputed) lookups."""
    with _lock, _connect() as conn:
        conn.execute(
            "INSERT OR IGNORE INTO calculations"
            " (cache_key, formula, monoisotopic_mass, peaks_json)"
            " VALUES (?, ?, ?, ?)",
            (cache_key, formula, monoisotopic_mass, json.dumps(peaks)),
        )
