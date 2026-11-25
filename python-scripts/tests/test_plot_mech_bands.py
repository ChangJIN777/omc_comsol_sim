"""Tests of scripts/plot_mech_bands.py (read-only mechanical band analysis).

Everything runs on SYNTHETIC records in a temp SQLite file (rib sym/antisym
dicts and hole lists with fy), deleted afterwards. The real databases are
not opened. No COMSOL/MPh, no network.

Run:  python tests/test_plot_mech_bands.py
"""
import ast
import csv
import hashlib
import io
import json
import os
import socket
import sqlite3
import sys
import tempfile
import time
from contextlib import redirect_stderr, redirect_stdout

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import numpy as np

import plot_mech_bands as pmb
from bandgap import Gap, all_gaps, interval_gap, largest_gap
from database import _SCHEMA
from objective_rib import _gap_to_dict

GHZ = 1e9
F_T = 8e9
K = np.linspace(1e-3, 1.0, 10)


# ---- fixtures ---------------------------------------------------------------

def _bands(centres_ghz, amp=0.15):
    """[n_k, n_b] Hz: flat-ish bands with a little dispersion."""
    return np.array([[(c + amp * np.sin(np.pi * k + j)) * GHZ
                      for j, c in enumerate(centres_ghz)] for k in K])


def rib_record(cid, mech_gap=None, freqs=None, with_bands=True, **kw):
    sym = _bands([1.0, 3.0, 5.0, 9.0, 13.0]) if freqs is None else freqs["sym"]
    anti = _bands([2.0, 4.0, 7.6, 11.0, 14.0]) if freqs is None else freqs["antisym"]
    gaps = {"sym": interval_gap([sym], F_T), "antisym": interval_gap([anti], F_T),
            "complete": interval_gap([sym, anti], F_T)}
    gp = gaps["complete"]
    rec = dict(id=cid, u=[0.5] * 4, celltype="rib", status="success",
               mech_backend="comsol", optical_backend="none",
               params=dict(a=600e-9, w=1300e-9, t=280e-9, s=100e-9, th=250e-9,
                           r1=40e-9, r2=40e-9),
               mech_gap_mode="complete", score=0.1, score_version=4,
               score_weights=dict(gmin_o=0.2, gmin_m=0.2, lambda_optical=10.0,
                                  lambda_mech=10.0, lambda_optical_freq=0.4,
                                  lambda_light_cone=2000.0,
                                  lambda_mech_window=2000.0,
                                  mech_window_GHz=[5.0, 10.0],
                                  optical_freq_tol_THz=5.0,
                                  optical_freq_power=8.0),
               mech_gaps={n: _gap_to_dict(g) for n, g in gaps.items()},
               mechanical_gap=gp.normalized_gap if mech_gap is None else mech_gap,
               mechanical_center_frequency=gp.f_center)
    if with_bands:
        to_list = lambda b: [[None if not np.isfinite(x) else float(x)  # noqa: E731
                              for x in row] for row in b]
        rec.update(_k_m=K.tolist(), _freqs_m=dict(sym=to_list(sym),
                                                 antisym=to_list(anti)))
    rec.update(kw)
    return rec


def hole_freqs():
    f = _bands([2.0, 4.0, 6.0, 10.0, 12.0])
    fy = np.tile([0.9, 0.9, 0.2, 0.8, 0.9], (K.size, 1))
    return f, fy


def hole_record(cid, legacy=False, **kw):
    f, fy = hole_freqs()
    gb = (pmb.legacy_breathing_gap if legacy else pmb.breathing_gap)(f, fy)
    gl = largest_gap(f, f_window=(2.5e9, 20e9))
    rec = dict(id=cid, u=[0.5] * 5, status="success", mech_backend="comsol",
               optical_backend="comsol",
               params=dict(a=550e-9, w=700e-9, t=300e-9, hx=120e-9, hy=250e-9),
               G_m=gb.normalized_gap if gb.found else 0.0,
               mechanical_gap=gb.normalized_gap if gb.found else 0.0,
               mechanical_center_frequency=gb.f_center if gb.found else 0.0,
               f_m_GHz=gb.f_center * 1e-9 if gb.found else 0.0,
               G_m_largest=gl.normalized_gap, f_m_largest_GHz=gl.f_center * 1e-9,
               fy_computed=True, score=0.3,
               _k_m=K.tolist(), _freqs_m=f.tolist(), _fy_arr=fy.tolist())
    rec.update(kw)
    return rec


def make_db(path, records):
    """A database like src/database.py writes, with chosen save times."""
    con = sqlite3.connect(path)
    con.execute(_SCHEMA)
    for ts, rec in records:
        con.execute(
            "INSERT INTO candidates (id, ts, u, params, optical_backend, "
            "mech_backend, status, optical_gap, mechanical_gap, "
            "optical_center_frequency, mechanical_center_frequency, score, "
            "record) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (rec["id"], ts, json.dumps(rec.get("u")),
             json.dumps(rec.get("params")), rec.get("optical_backend"),
             rec.get("mech_backend"), rec.get("status"), None,
             rec.get("mechanical_gap"), None,
             rec.get("mechanical_center_frequency"), rec.get("score"),
             json.dumps(rec)))
    con.commit()
    con.close()


def _md5(path):
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def _run(argv):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        rc = pmb.main(argv)
    return rc, out.getvalue(), err.getvalue()


def _summary(outdir, cid):
    with open(os.path.join(outdir, cid, "summary.json")) as fh:
        return json.load(fh)


def _close(a, b):
    return abs(a - b) <= 1e-12 * max(1.0, abs(b))


# ---- tests ------------------------------------------------------------------

def test_breathing_gap_is_the_one_of_run_opt_comsol():
    """The port equals scripts/run_opt_comsol.breathing_gap (taken from its
    source with ast: importing that script would import the COMSOL client)."""
    path = os.path.join(_ROOT, "scripts", "run_opt_comsol.py")
    with open(path) as fh:
        tree = ast.parse(fh.read())
    consts, func = {}, None
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            name = getattr(node.targets[0], "id", None)
            if name in ("FY_THRESH", "MECH_F_MIN", "MECH_F_MAX"):
                consts[name] = ast.literal_eval(node.value)
        if isinstance(node, ast.FunctionDef) and node.name == "breathing_gap":
            func = node
    assert consts == dict(FY_THRESH=pmb.FY_THRESH, MECH_F_MIN=pmb.MECH_F_MIN,
                          MECH_F_MAX=pmb.MECH_F_MAX), consts
    ns = dict(np=np, all_gaps=all_gaps, Gap=Gap, **consts)
    exec(compile(ast.Module(body=[func], type_ignores=[]), path, "exec"), ns)
    original = ns["breathing_gap"]
    rng = np.random.default_rng(1)
    for _ in range(200):
        n_b = int(rng.integers(3, 9))
        f = np.sort(rng.uniform(1e9, 20e9, (7, n_b)), axis=1)
        if rng.random() < 0.3:
            f[rng.integers(0, 7), rng.integers(0, n_b)] = np.nan
        fy = rng.uniform(0, 1, f.shape)
        assert original(f, fy) == pmb.breathing_gap(f, fy)
    f, fy = hole_freqs()
    assert pmb.breathing_gap(f, fy) != pmb.legacy_breathing_gap(f, fy)


def test_rib_gaps_are_recomputed_with_interval_gap_and_compared():
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "runs_rib.sqlite3"), os.path.join(tmp, "out")
        good, bad = rib_record("cand_good"), rib_record("cand_bad", mech_gap=0.5)
        make_db(db, [(1.0, good), (2.0, bad)])
        rc, stdout, _ = _run(["--db", db, "--all", "--outdir", out])
        assert rc == 0, stdout
        s = _summary(out, "cand_good")
        sym = np.array(good["_freqs_m"]["sym"], dtype=float)
        anti = np.array(good["_freqs_m"]["antisym"], dtype=float)
        for name, sets in (("sym", [sym]), ("antisym", [anti]),
                           ("complete", [sym, anti])):
            want = interval_gap(sets, F_T)
            got = s["gaps"][name]
            assert want.found and got["found"]
            assert _close(got["normalized_gap"], want.normalized_gap)
            assert _close(got["f_lower_GHz"], want.f_lower * 1e-9)
            assert _close(got["f_upper_GHz"], want.f_upper * 1e-9)
            assert s["comparison"][name]["matches"], s["comparison"][name]
        assert s["scored_gap"]["name"] == "complete"
        assert s["matches_record"] is True
        assert s["has_bands"] and s["n_k"] == K.size
        assert s["n_bands"] == {"sym": 5, "antisym": 5}
        assert s["window_GHz"] == [5.0, 10.0]
        fc = s["scored_gap"]["f_center_GHz"]
        assert s["in_window"] == (5.0 <= fc <= 10.0)
        assert s["penalties"]["recomputed"]["pen_window_m"] >= 0.0
        # tampered record -> flagged
        sb = _summary(out, "cand_bad")
        assert sb["matches_record"] is False
        assert "normalized_gap" in sb["comparison"]["mechanical_gap"]["differs"]
        assert "DIFFERS" in stdout and "match" in stdout
        # files
        for cid in ("cand_good", "cand_bad"):
            for name in ("band_diagram.png", "bands.csv", "summary.json"):
                assert os.path.getsize(os.path.join(out, cid, name)) > 0
        with open(os.path.join(out, "cand_good", "bands.csv")) as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 2 * K.size * 5
        assert {r["sector"] for r in rows} == {"sym", "antisym"}
        with open(os.path.join(out, "mech_index.csv")) as fh:
            index = list(csv.DictReader(fh))
        assert [r["record_id"] for r in index] == ["cand_good", "cand_bad"]
        assert [r["matches_record"] for r in index] == ["1", "0"]
        assert os.path.getsize(os.path.join(out, "mech_summary.png")) > 0


def test_rib_window_penalty_uses_objective_rib():
    from objective_rib import mech_terms
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "r.sqlite3"), os.path.join(tmp, "out")
        rec = rib_record("cand_w")
        g = rec["mech_gaps"]["complete"]
        terms = mech_terms(g["normalized_gap"], g["f_center"], True,
                           rec["score_weights"])
        rec["pen_window_m"] = terms["pen_window_m"] + 1.0     # wrong on purpose
        make_db(db, [(1.0, rec)])
        rc, stdout, _ = _run(["--db", db, "cand_w", "--outdir", out, "--no-plot"])
        assert rc == 0
        s = _summary(out, "cand_w")
        assert _close(s["penalties"]["recomputed"]["pen_window_m"],
                      terms["pen_window_m"])
        assert s["penalties"]["weights_source"] == "record score_weights"
        assert any("pen_window_m differs" in n for n in s["notes"])


def test_hole_records_breathing_gap_and_legacy_definition():
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "runs.sqlite3"), os.path.join(tmp, "out")
        make_db(db, [(1.0, hole_record("opt_now")),
                     (2.0, hole_record("opt_old", legacy=True)),
                     (3.0, hole_record("opt_bad", mechanical_gap=0.77, G_m=0.77))])
        rc, stdout, _ = _run(["--db", db, "--all", "--outdir", out])
        assert rc == 0, stdout
        f, fy = hole_freqs()
        s = _summary(out, "opt_now")
        want = pmb.breathing_gap(f, fy)
        assert _close(s["gaps"]["breathing"]["normalized_gap"], want.normalized_gap)
        assert _close(s["gaps"]["largest"]["normalized_gap"],
                      largest_gap(f, f_window=(2.5e9, 20e9)).normalized_gap)
        assert _close(s["gaps"]["interval"]["normalized_gap"],
                      interval_gap([f], F_T).normalized_gap)
        assert s["matches_record"] is True
        assert s["comparison"]["G_m_largest"]["matches"]
        assert s["n_bands"] == {"sym": 5}
        old = _summary(out, "opt_old")
        assert old["matches_record"] == "legacy"
        assert "breathing_legacy" in old["gaps"]
        assert _summary(out, "opt_bad")["matches_record"] is False
        with open(os.path.join(out, "opt_now", "bands.csv")) as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == K.size * 5 and all(r["fy"] for r in rows)
        assert "legacy def" in stdout


def test_missing_arrays_and_nan_bands():
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "r.sqlite3"), os.path.join(tmp, "out")
        sym = _bands([1.0, 3.0, 5.0, 9.0, 13.0])
        anti = _bands([2.0, 4.0, 7.6, 11.0, 14.0])
        sym[3, 4] = np.nan                      # a mode not returned
        anti[:, 0] = np.nan
        make_db(db, [(1.0, rib_record("cand_none", with_bands=False)),
                     (2.0, rib_record("cand_nan", freqs=dict(sym=sym, antisym=anti))),
                     (3.0, dict(hole_record("cand_hole"), _k_m=None, _freqs_m=None,
                                _fy_arr=None))])
        rc, stdout, _ = _run(["--db", db, "--all", "--outdir", out])
        assert rc == 0, stdout
        s = _summary(out, "cand_none")
        assert not s["has_bands"] and s["gaps"] == {}
        assert any("no band arrays" in n for n in s["notes"])
        assert s["scored_gap"]["recorded_only"]            # from mech_gaps
        assert os.path.getsize(os.path.join(out, "cand_none", "band_diagram.png")) > 0
        with open(os.path.join(out, "cand_none", "bands.csv")) as fh:
            assert len(fh.read().strip().splitlines()) == 1    # header only
        s = _summary(out, "cand_nan")
        assert s["matches_record"] is True
        assert any("NaN" in n for n in s["notes"])
        with open(os.path.join(out, "cand_nan", "bands.csv")) as fh:
            rows = list(csv.DictReader(fh))
        assert sum(1 for r in rows if r["freq_GHz"] == "") == 1 + K.size
        s = _summary(out, "cand_hole")
        assert not s["has_bands"] and s["matches_record"] is None
        assert any("src/objective.py" in n for n in s["notes"])


def test_databases_are_opened_read_only():
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "r.sqlite3"), os.path.join(tmp, "out")
        make_db(db, [(1.0, rib_record("cand_a")), (2.0, hole_record("opt_b"))])
        before, files = _md5(db), sorted(os.listdir(tmp))
        rc, _, _ = _run(["--db", db, "--all", "--status", "any",
                         "--mech-backend", "any", "--outdir", out])
        assert rc == 0
        assert _md5(db) == before
        assert sorted(os.listdir(tmp)) == sorted(files + ["out"])   # no journal
        con = pmb.connect_ro(db)
        try:
            try:
                con.execute("DELETE FROM candidates")
                raise AssertionError("a read-only connection could write")
            except sqlite3.OperationalError:
                pass
        finally:
            con.close()
        assert _md5(db) == before


def test_out_json_falls_back_to_sqlite_by_id():
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "r.sqlite3"), os.path.join(tmp, "out")
        a, b = rib_record("cand_a"), rib_record("cand_only_json")
        # the same id evaluated again later: SQLite holds the last evaluation
        h_db = hole_record("opt_again", fy_computed=False)
        h_js = hole_record("opt_again", mechanical_gap=0.0123)
        make_db(db, [(1.0, a), (2.0, h_db)])
        slim = [{k: v for k, v in r.items() if not k.startswith("_")}
                for r in (a, b, h_js)]
        js = os.path.join(tmp, "opt_results.json")
        with open(js, "w") as fh:
            json.dump(slim, fh)
        before = (_md5(db), _md5(js))
        rc, stdout, stderr = _run(["--json", js, "--db", db, "--all",
                                   "--outdir", out, "--no-plot"])
        assert rc == 0, stdout + stderr
        assert "no band arrays" in stderr and "2 of 3" in stderr
        s = _summary(out, "cand_a")
        assert s["has_bands"] and s["matches_record"] is True
        assert any("band arrays from" in n for n in s["notes"])
        s = _summary(out, "cand_only_json")
        assert not s["has_bands"]
        s = _summary(out, "opt_again")       # the SQLite record, consistent
        assert s["has_bands"] and s["matches_record"] is True
        assert s["extra"]["fy_computed"] is False
        assert any("evaluated again" in n for n in s["notes"])
        assert (_md5(db), _md5(js)) == before


def test_selection_and_filters():
    day = 86400.0
    t0 = time.mktime((2026, 11, 20, 12, 0, 0, 0, 0, -1))
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "r.sqlite3"), os.path.join(tmp, "out")
        make_db(db, [
            (t0, rib_record("cand_old")),
            (t0 + 2 * day, rib_record("cand_new")),
            (t0 + 3 * day, rib_record("cand_stub", mech_backend="surrogate_stub")),
            (t0 + 4 * day, rib_record("cand_fail", status="mech_failed")),
            (t0 + 5 * day, rib_record("cand_small", mech_gap=0.01)),
        ])

        def ids(*argv):
            rc, stdout, _ = _run(["--db", db, "--outdir", out, "--no-plot"]
                                 + list(argv))
            assert rc in (0, 1), stdout
            return [line.split()[0] for line in stdout.splitlines()
                    if line.startswith("cand_")]
        assert ids("--all") == ["cand_old", "cand_new", "cand_small"]
        assert ids("--latest", "2") == ["cand_small", "cand_new"]
        assert ids("--all", "--min-gap", "0.05") == ["cand_old", "cand_new"]
        assert ids("--all", "--since", "2026-11-21") == ["cand_new", "cand_small"]
        assert ids("--all", "--status", "any") == [
            "cand_old", "cand_new", "cand_fail", "cand_small"]
        assert "cand_stub" in ids("--all", "--mech-backend", "any")
        assert ids("cand_fail") == ["cand_fail"]          # explicit id: no filter
        assert ids("cand_sm") == ["cand_small"]           # unique prefix
        rc, _, err = _run(["--db", db, "cand_", "--outdir", out])
        assert rc == 2 and "ambiguous" in err
        rc, _, err = _run(["--db", os.path.join(tmp, "nope.sqlite3"), "--all",
                           "--outdir", out])
        assert rc == 2 and "not found" in err


def test_network_is_blocked_while_the_script_runs():
    with pmb.network_blocked():
        try:
            socket.create_connection(("127.0.0.1", 9))
            raise AssertionError("connection was not blocked")
        except pmb.AnalysisError:
            pass
    assert socket.create_connection is not None


TESTS = [test_breathing_gap_is_the_one_of_run_opt_comsol,
         test_rib_gaps_are_recomputed_with_interval_gap_and_compared,
         test_rib_window_penalty_uses_objective_rib,
         test_hole_records_breathing_gap_and_legacy_definition,
         test_missing_arrays_and_nan_bands,
         test_databases_are_opened_read_only,
         test_out_json_falls_back_to_sqlite_by_id,
         test_selection_and_filters,
         test_network_is_blocked_while_the_script_runs]

if __name__ == "__main__":
    for test in TESTS:
        print(test.__name__)
        test()
    print("all %d tests passed" % len(TESTS))
