"""Persist every candidate result to SQLite. Never lose a simulation.

Ported from python-scripts/src/database.py with the columns renamed from
optical/mechanical-gap-first to mechanical/optical-gap-first (mechanical is
the implemented backend in this project; optical is a future stub -- see
src/optical_comsol_2d.py) and a distinct env var / default path so this
project's results never collide with python-scripts/results/runs.sqlite.

The column list is fixed; every other record key rides in the `record` JSON
blob. That is not a limitation worth working around -- SQLite's JSON1 functions
make the blob queryable, e.g.

    SELECT id, score,
           json_extract(record, '$.mech_parity')          AS parity,
           json_extract(record, '$.mechanical_gap_oddz')  AS gap_oddz
    FROM candidates WHERE parity = 'oddz' ORDER BY score DESC;

Prefer that over adding a column. `CREATE TABLE IF NOT EXISTS` does NOT migrate
an existing runs2d.sqlite, so a new column needs an `ALTER TABLE` guarded by a
`PRAGMA table_info(candidates)` check, and every already-saved row would read
NULL for it.
"""
from __future__ import annotations

import json
import os
import sqlite3
import time
from typing import Dict, List

_DB = os.environ.get(
    "OMC2D_DB_PATH",
    os.path.join(os.path.dirname(__file__), "..", "results", "runs2d.sqlite"))

_SCHEMA = """
CREATE TABLE IF NOT EXISTS candidates (
    id TEXT PRIMARY KEY,
    ts REAL,
    u TEXT,
    params TEXT,
    mech_backend TEXT,
    optical_backend TEXT,
    status TEXT,
    mechanical_gap REAL,
    optical_gap REAL,
    mechanical_center_frequency REAL,
    optical_center_frequency REAL,
    score REAL,
    record TEXT
);
"""


def _conn(path=_DB):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    c = sqlite3.connect(path)
    c.execute(_SCHEMA)
    return c


def _jsonl_path(path):
    return os.path.splitext(path)[0] + ".jsonl"


def save_result(rec: Dict, path=_DB) -> None:
    try:
        _save_sqlite(rec, path)
    except sqlite3.OperationalError:
        # filesystem doesn't support sqlite locking (some network/FUSE mounts);
        # fall back to append-only JSONL so we never lose a simulation.
        os.makedirs(os.path.dirname(_jsonl_path(path)), exist_ok=True)
        with open(_jsonl_path(path), "a") as fh:
            fh.write(json.dumps(rec) + "\n")


def _save_sqlite(rec: Dict, path=_DB) -> None:
    c = _conn(path)
    with c:
        c.execute(
            """INSERT OR REPLACE INTO candidates
               (id, ts, u, params, mech_backend, optical_backend, status,
                mechanical_gap, optical_gap, mechanical_center_frequency,
                optical_center_frequency, score, record)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (rec["id"], time.time(), json.dumps(rec.get("u")),
             json.dumps(rec.get("params")), rec.get("mech_backend"),
             rec.get("optical_backend"), rec.get("status"),
             rec.get("mechanical_gap"), rec.get("optical_gap"),
             rec.get("mechanical_center_frequency"),
             rec.get("optical_center_frequency"), rec.get("score"),
             json.dumps(rec)))
    c.close()


def load_completed_results(path=_DB) -> List[Dict]:
    out: List[Dict] = []
    if os.path.exists(path):
        try:
            c = _conn(path)
            rows = c.execute("SELECT record FROM candidates").fetchall()
            c.close()
            out.extend(json.loads(r[0]) for r in rows)
        except sqlite3.OperationalError:
            pass
    jl = _jsonl_path(path)
    if os.path.exists(jl):
        with open(jl) as fh:
            out.extend(json.loads(line) for line in fh if line.strip())
    return out


def best_result(path=_DB):
    res = [r for r in load_completed_results(path)
           if r.get("status") == "success"]
    return max(res, key=lambda r: r["score"]) if res else None
