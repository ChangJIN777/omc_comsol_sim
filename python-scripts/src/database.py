"""Persist every candidate result to SQLite. Never lose a simulation."""
from __future__ import annotations

import json
import os
import sqlite3
import time
from typing import Dict, List

_DB = os.environ.get(
    "OMC_DB_PATH",
    os.path.join(os.path.dirname(__file__), "..", "results", "runs.sqlite"))

_SCHEMA = """
CREATE TABLE IF NOT EXISTS candidates (
    id TEXT PRIMARY KEY,
    ts REAL,
    u TEXT,
    params TEXT,
    optical_backend TEXT,
    mech_backend TEXT,
    status TEXT,
    optical_gap REAL,
    mechanical_gap REAL,
    optical_center_frequency REAL,
    mechanical_center_frequency REAL,
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
               (id, ts, u, params, optical_backend, mech_backend, status,
                optical_gap, mechanical_gap, optical_center_frequency,
                mechanical_center_frequency, score, record)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (rec["id"], time.time(), json.dumps(rec.get("u")),
             json.dumps(rec.get("params")), rec.get("optical_backend"),
             rec.get("mech_backend"), rec.get("status"),
             rec.get("optical_gap"), rec.get("mechanical_gap"),
             rec.get("optical_center_frequency"),
             rec.get("mechanical_center_frequency"), rec.get("score"),
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
