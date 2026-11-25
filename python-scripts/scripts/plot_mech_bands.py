#!/usr/bin/env python3
"""Mechanical band diagrams and gaps of saved records (read-only analysis).

Every evaluation is saved by src/database.save_result to SQLite (table
`candidates`, full JSON in the `record` column):
  hole cell  results/runs.sqlite3       (env OMC_DB_PATH overrides)
  rib cell   results/runs_rib.sqlite3   (env OMC_DB_PATH_RIB overrides)
plus a JSONL fallback next to each file. The out_json files of the loops
(results/opt_results_*.json) hold the same records WITHOUT the `_`-prefixed
band arrays; with --json the arrays are looked up in SQLite by id.

What a record must hold for a band diagram (otherwise: summary only, noted):
  rib   _k_m (k a/pi), _freqs_m = {"sym": [n_k][n_b], "antisym": ...} [Hz]
        (y-even study 'mech sym', y-odd 'mech antisym'), mech_gap_mode,
        mech_gaps, mechanical_gap, mechanical_center_frequency,
        score_version, score_weights, pen_window_m, pen_size_m
        -- written by src/objective_rib.py (run_one.py, run_loop.py,
        run_opt_comsol_rib.py, run_batch_rib.py).
  hole  _k_m, _freqs_m = [n_k][n_b] [Hz] (the run's study_mech, 'mech sym'
        by default; not recorded), _fy_arr (u_y energy fraction per mode,
        all 1.0 when fy_computed is false), mechanical_gap / G_m,
        mechanical_center_frequency, G_m_largest, f_m_largest_GHz
        -- written by scripts/run_opt_comsol.py only. Hole records of
        run_one.py / run_loop.py (src/objective.py) store no band arrays.

Gaps are recomputed with the repo's own definitions and compared with the
recorded ones (mismatches are flagged):
  rib   bandgap.interval_gap per sector and over both sectors (complete),
        target targets.yaml mechanical.target_frequency_GHz (8 GHz), rel_tol
        0.5, exactly as objective_rib.evaluate_rib_candidate; the scored gap
        is the one of mech_gap_mode. The window penalty is
        objective_rib.mech_terms with the record's score_weights (current
        targets.yaml weights when the record has none).
  hole  the breathing gap of scripts/run_opt_comsol.py (breathing_gap: the
        largest bandgap.all_gaps gap, centre in 2.5-20 GHz, with an edge band
        of mean fy >= 0.5) and its plain largest_gap in the same window.
        breathing_gap is re-implemented here (that script imports the COMSOL
        client); tests/test_plot_mech_bands.py checks it against the source.
The 5-10 GHz window (rib: score_weights.mech_window_GHz; hole: MECH_F_MIN /
MECH_F_MAX, informational) and the 8 GHz target are reported and drawn.

READ-ONLY: the databases are opened with sqlite3 URI mode=ro, no result
file is written, no COMSOL/MPh is imported and network connections are
blocked while the script runs. Output goes to --outdir only.

Examples (from python-scripts/):
  python scripts/plot_mech_bands.py cand_812956efbbba            # one id (prefix ok)
  python scripts/plot_mech_bands.py --celltype rib --latest 5
  python scripts/plot_mech_bands.py --celltype hole --all --min-gap 0.2
  python scripts/plot_mech_bands.py --celltype rib --all --since 2026-11-20
  python scripts/plot_mech_bands.py --json results/opt_results_rib_tidy3d.json --all
  python scripts/plot_mech_bands.py --db /tmp/copy.sqlite3 --all --status any
  python scripts/plot_mech_bands.py --celltype rib --all --outdir /tmp/mech

Selection: record ids (or unique prefixes), --latest N (newest save time
first) or --all. Without --celltype / --db both default databases are read.
The filters (--status, default success; --mech-backend, default comsol;
--since; --min-gap on the recorded mechanical_gap) apply to --latest / --all;
ids given explicitly are always analysed.

Outputs, per record in <outdir>/<record id>/: band_diagram.png, bands.csv
(sector, k_a_over_pi, band, freq_GHz, fy) and summary.json; with several
records also <outdir>/mech_index.csv and <outdir>/mech_summary.png (gap vs
centre and vs each geometry variable, coloured by score). Default outdir:
results/mech_analysis.
"""
import argparse
import contextlib
import csv
import datetime
import json
import math
import os
import pathlib
import socket
import sqlite3
import sys
import warnings

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_ROOT, "src"))

from bandgap import Gap, all_gaps, interval_gap, largest_gap   # noqa: E402

DEFAULT_OUTDIR = os.path.join(_ROOT, "results", "mech_analysis")
TARGETS_YAML = os.path.join(_ROOT, "configs", "targets.yaml")
DEFAULT_TARGET_GHZ = 8.0
DEFAULT_WINDOW_GHZ = (5.0, 10.0)
REL_TOL = 0.5                 # interval_gap default, as objective_rib uses it
# scripts/run_opt_comsol.py (hole loop); kept equal by the tests
MECH_F_MIN, MECH_F_MAX = 5e9, 10e9
FY_THRESH = 0.5
SECTORS = ("sym", "antisym")
SECTOR_LABEL = {"sym": "y-even (mech sym)", "antisym": "y-odd (mech antisym)"}
RTOL = 1e-6                   # recorded vs recomputed

# colours as scripts/plot_tidy3d_bands.py; sectors also differ by marker
# shape and fill, so nothing is told apart by colour alone
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
SURFACE, CONE = "#fcfcfb", "#e9e8e4"
SEQ_RAMP = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
SECTOR_COLOR = {"sym": "#257de0", "antisym": "#f47230"}
SECTOR_MARKER = {"sym": "o", "antisym": "^"}
GAP = "#1baf7a"


class AnalysisError(RuntimeError):
    """A problem to show the user without a traceback."""


@contextlib.contextmanager
def network_blocked():
    """No connection can be opened while the block runs (local analysis)."""
    def refuse(*args, **kw):
        raise AnalysisError("network access attempted: this script works on "
                            "local files only")
    saved = (socket.socket.connect, socket.socket.connect_ex,
             socket.create_connection)
    socket.socket.connect = refuse
    socket.socket.connect_ex = refuse
    socket.create_connection = refuse
    try:
        yield
    finally:
        (socket.socket.connect, socket.socket.connect_ex,
         socket.create_connection) = saved


# ---- configuration ----------------------------------------------------------

def load_targets():
    """configs/targets.yaml, or None when it cannot be read."""
    try:
        import yaml
        with open(TARGETS_YAML) as fh:
            return yaml.safe_load(fh)
    except Exception:
        return None


def target_hz(targets):
    try:
        return float(targets["mechanical"]["target_frequency_GHz"]) * 1e9
    except Exception:
        return DEFAULT_TARGET_GHZ * 1e9


def default_db(celltype):
    """The database of a cell type, as src/database.py resolves it."""
    if celltype == "hole":
        return os.environ.get("OMC_DB_PATH",
                              os.path.join(_ROOT, "results", "runs.sqlite3"))
    return os.environ.get("OMC_DB_PATH_RIB",
                          os.path.join(_ROOT, "results", "runs_rib.sqlite3"))


# ---- reading (read-only) ----------------------------------------------------

def connect_ro(path):
    """sqlite3 connection that cannot write (URI mode=ro)."""
    if not os.path.isfile(path):
        raise AnalysisError(f"database not found: {path}")
    uri = pathlib.Path(path).resolve().as_uri() + "?mode=ro"
    return sqlite3.connect(uri, uri=True)


def read_db(path):
    """Every record of a database: list of dict(rec, ts, source).

    The JSONL fallback next to it (database._jsonl_path) is read too; its
    records have no save time and are used only for ids not in SQLite.
    """
    out, seen = [], set()
    con = connect_ro(path)
    try:
        rows = con.execute("SELECT id, ts, record FROM candidates").fetchall()
    except sqlite3.DatabaseError as err:
        raise AnalysisError(f"cannot read {path}: {err}")
    finally:
        con.close()
    for cid, ts, raw in rows:
        try:
            rec = json.loads(raw)
        except (TypeError, ValueError):
            continue
        out.append(dict(rec=rec, ts=ts, source=path))
        seen.add(cid)
    jl = os.path.splitext(path)[0] + ".jsonl"
    if os.path.isfile(jl):
        with open(jl) as fh:
            for line in fh:
                if not line.strip():
                    continue
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue
                if rec.get("id") not in seen:
                    out.append(dict(rec=rec, ts=None, source=jl))
                    seen.add(rec.get("id"))
    return out


def db_lookup(paths):
    """id -> entry over several databases (the first one wins)."""
    table = {}
    for p in paths:
        if os.path.isfile(p):
            for e in read_db(p):
                table.setdefault(e["rec"].get("id"), e)
    return table


MECH_KEYS = ("status", "mech_gap_mode", "fy_computed", "mechanical_gap",
             "mechanical_center_frequency")


def _record_differences(a, b):
    """Mechanical fields that differ between two records of the same id."""
    out = []
    for key in MECH_KEYS:
        x, y = a.get(key), b.get(key)
        if isinstance(x, (int, float)) and isinstance(y, (int, float)) \
                and not isinstance(x, bool) and not isinstance(y, bool):
            if not math.isclose(x, y, rel_tol=RTOL, abs_tol=1e-12):
                out.append("%s %.6g vs %.6g" % (key, x, y))
        elif x != y:
            out.append("%s %s vs %s" % (key, x, y))
    return out


def read_json(path, db_paths):
    """Records of an out_json; band arrays (stripped there) from SQLite."""
    if not os.path.isfile(path):
        raise AnalysisError(f"out_json not found: {path}")
    with open(path) as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        data = data.get("results") or list(data.values())
    table = db_lookup(db_paths)
    out, n_missing = [], 0
    for rec in data:
        if not isinstance(rec, dict):
            continue
        e = dict(rec=rec, ts=None, source=path, note=None)
        hit = table.get(rec.get("id"))
        if rec.get("_freqs_m") is None:
            if hit is not None and hit["rec"].get("_freqs_m") is not None:
                diff = _record_differences(rec, hit["rec"])
                if diff:
                    # the id was evaluated again (hole ids depend on the
                    # geometry only) and SQLite kept the last evaluation:
                    # its arrays do not belong to this entry, so analyse
                    # the SQLite record, whose values match its arrays
                    e.update(rec=hit["rec"], note=(
                        "out_json entry differs from the %s record of this id "
                        "(%s): the id was evaluated again and SQLite keeps the "
                        "last evaluation only; the SQLite record is analysed"
                        % (os.path.basename(hit["source"]), "; ".join(diff))))
                else:
                    e.update(rec=dict(rec, **{k: v for k, v in hit["rec"].items()
                                              if k.startswith("_")}),
                             note="band arrays from %s (out_json has none)"
                                  % os.path.basename(hit["source"]))
            else:
                n_missing += 1
        if hit is not None:
            e["ts"] = hit["ts"]
        out.append(e)
    print("note: %s holds no band arrays (the loops strip `_` keys); %d of %d "
          "records found them in %s" % (
              os.path.basename(path), len(out) - n_missing, len(out),
              ", ".join(os.path.basename(p) for p in db_paths
                        if os.path.isfile(p)) or "no database"),
          file=sys.stderr)
    return out


def parse_since(text):
    try:
        return datetime.datetime.fromisoformat(text).timestamp()
    except ValueError:
        raise AnalysisError(f"--since: not a date (YYYY-MM-DD[ HH:MM]): {text!r}")


def keep(entry, args, since_ts):
    rec = entry["rec"]
    if args.status != "any" and rec.get("status") != args.status:
        return False
    if args.mech_backend != "any" and rec.get("mech_backend") != args.mech_backend:
        return False
    if since_ts is not None and (entry["ts"] is None or entry["ts"] < since_ts):
        return False
    if args.min_gap is not None:
        g = rec.get("mechanical_gap")
        if g is None or not g >= args.min_gap:
            return False
    return True


def select(entries, args):
    """The entries to analyse, in output order."""
    if args.ids:
        out = []
        for want in args.ids:
            hits = [e for e in entries if e["rec"].get("id") == want]
            if not hits:
                hits = [e for e in entries
                        if str(e["rec"].get("id", "")).startswith(want)]
            if not hits:
                raise AnalysisError(f"no record with id {want!r}")
            if len({e["rec"]["id"] for e in hits}) > 1:
                raise AnalysisError(
                    f"id prefix {want!r} is ambiguous: "
                    + ", ".join(sorted({e["rec"]["id"] for e in hits})[:8]))
            out.append(hits[0])
        return out
    since_ts = parse_since(args.since) if args.since else None
    kept = [e for e in entries if keep(e, args, since_ts)]
    if args.latest is not None:
        if args.json:      # out_json order is evaluation order
            kept = kept[-args.latest:][::-1]
        else:
            kept = sorted(kept, key=lambda e: e["ts"] or 0.0,
                          reverse=True)[:args.latest]
    elif not args.json:
        kept = sorted(kept, key=lambda e: e["ts"] or 0.0)
    return kept


# ---- analysis ---------------------------------------------------------------

def celltype_of(rec):
    if rec.get("celltype"):
        return rec["celltype"]
    p = rec.get("params") or {}
    return "rib" if "th" in p or "s" in p else "hole"


def _arr(x):
    """[n_k, n_b] float array (None -> NaN), or None."""
    if x is None:
        return None
    try:
        a = np.array([[np.nan if v is None else v for v in row] for row in x],
                     dtype=float)
    except (TypeError, ValueError):
        return None
    return a if a.ndim == 2 and a.size else None


def bands_of(rec):
    """(k, {sector: [n_k, n_b] Hz}, fy or None, note); k None if absent."""
    k = rec.get("_k_m")
    fm = rec.get("_freqs_m")
    if k is None or fm is None:
        return None, {}, None, "no band arrays in the record (summary only)"
    k = np.asarray(k, dtype=float)
    if isinstance(fm, dict):
        bands = {s: _arr(fm.get(s)) for s in SECTORS if fm.get(s) is not None}
        bands = {s: b for s, b in bands.items() if b is not None}
        fy = None
    else:
        b = _arr(fm)
        bands = {"sym": b} if b is not None else {}
        fy = _arr(rec.get("_fy_arr"))
        if fy is not None and (b is None or fy.shape != b.shape):
            fy = None
    for s, b in list(bands.items()):
        if b.shape[0] != k.size:
            return None, {}, None, ("band array %s has %d k points, _k_m %d "
                                    "(summary only)" % (s, b.shape[0], k.size))
    if not bands:
        return None, {}, None, "band arrays empty or unreadable (summary only)"
    return k, bands, fy, None


def legacy_breathing_gap(freqs_all, fy_arr, fy_thresh=FY_THRESH,
                         f_min=MECH_F_MIN, f_max=MECH_F_MAX):
    """The breathing gap as run_opt_comsol.py computed it before its
    breathing_gap was rewritten (2026-07-01): modes of fy < fy_thresh set to
    NaN point by point, then largest_gap(nan_safe=True), so BOTH edges had
    to be breathing. Used only to explain mismatches of old records."""
    freqs = np.where(np.asarray(fy_arr, dtype=float) >= fy_thresh,
                     np.asarray(freqs_all, dtype=float), np.nan)
    with warnings.catch_warnings():      # all-NaN bands are expected here
        warnings.simplefilter("ignore", RuntimeWarning)
        return largest_gap(freqs, f_window=(f_min * 0.5, f_max * 2.0),
                           nan_safe=True)


def breathing_gap(freqs_all, fy_arr, fy_thresh=FY_THRESH,
                  f_min=MECH_F_MIN, f_max=MECH_F_MAX):
    """scripts/run_opt_comsol.breathing_gap (same code; that module imports
    the COMSOL client, so it is not imported here)."""
    freqs_all = np.asarray(freqs_all, dtype=float)
    fy_arr = np.asarray(fy_arr, dtype=float)
    n_bands = freqs_all.shape[1]
    fy_mean = np.nanmean(fy_arr, axis=0)
    window = (f_min * 0.5, f_max * 2.0)
    for gp in all_gaps(freqs_all, f_window=window, nan_safe=True):
        i = gp.lower_band
        if not (0 <= i < n_bands - 1):
            continue
        if max(fy_mean[i], fy_mean[i + 1]) >= fy_thresh:
            return gp
    return Gap.empty()


def gap_dict(gp):
    """A Gap in GHz (JSON-ready)."""
    if gp is None:
        return None
    return dict(found=bool(gp.found), f_lower_GHz=float(gp.f_lower) * 1e-9,
                f_upper_GHz=float(gp.f_upper) * 1e-9,
                f_center_GHz=float(gp.f_center) * 1e-9,
                normalized_gap=float(gp.normalized_gap),
                lower_band=int(gp.lower_band))


def _close(a, b):
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(float(a), float(b), rel_tol=RTOL, abs_tol=1e-9)


def compare(recorded, recomputed):
    """recorded: dict(found?, normalized_gap, f_center [Hz], f_lower?, f_upper?)
    recomputed: Gap. Returns dict(recorded GHz, matches, differs)."""
    rec_found = recorded.get("found")
    if rec_found is None:
        rec_found = bool((recorded.get("normalized_gap") or 0.0) > 0.0)
    diffs = []
    if bool(rec_found) != bool(recomputed.found):
        diffs.append("found")
    if not _close(recorded.get("normalized_gap") or 0.0,
                  recomputed.normalized_gap if recomputed.found else 0.0):
        diffs.append("normalized_gap")
    if rec_found and recomputed.found:
        for key, val in (("f_center", recomputed.f_center),
                         ("f_lower", recomputed.f_lower),
                         ("f_upper", recomputed.f_upper)):
            if recorded.get(key) is not None and not _close(recorded[key], val):
                diffs.append(key)
    ghz = {k + "_GHz": (None if recorded.get(k) is None
                        else float(recorded[k]) * 1e-9)
           for k in ("f_lower", "f_upper", "f_center")}
    return dict(recorded=dict(found=bool(rec_found),
                              normalized_gap=recorded.get("normalized_gap"), **ghz),
                matches=not diffs, differs=diffs)


def rib_weights(rec, targets):
    """(weights, source) for objective_rib.mech_terms."""
    w = rec.get("score_weights")
    if isinstance(w, dict) and "mech_window_GHz" in w and "lambda_mech_window" in w:
        return w, "record score_weights"
    try:
        from objective_rib import rib_score_weights
        return rib_score_weights(targets), "configs/targets.yaml (current)"
    except Exception:
        return None, None


def analyze(entry, targets):
    rec = entry["rec"]
    cell = celltype_of(rec)
    f_t = target_hz(targets)
    res = dict(record_id=rec.get("id"), celltype=cell, status=rec.get("status"),
               source=entry["source"], ts=entry.get("ts"),
               notes=[n for n in [entry.get("note")] if n],
               mech_backend=rec.get("mech_backend"),
               optical_backend=rec.get("optical_backend"),
               params=rec.get("params") or {}, score=rec.get("score"),
               score_version=rec.get("score_version"),
               mech_gap_mode=rec.get("mech_gap_mode") if cell == "rib" else
               "breathing (sym study, fy >= %.1f)" % FY_THRESH,
               target_GHz=f_t * 1e-9, window_GHz=list(DEFAULT_WINDOW_GHZ),
               recorded=dict(mechanical_gap=rec.get("mechanical_gap"),
                             mechanical_center_GHz=(
                                 None if rec.get("mechanical_center_frequency") is None
                                 else rec["mechanical_center_frequency"] * 1e-9)),
               k=None, bands={}, fy=None, gaps={}, comparison={},
               scored=None, penalties=None, extra={})
    k, bands, fy, note = bands_of(rec)
    if note:
        res["notes"].append(note)
    res.update(k=k, bands=bands, fy=fy)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if cell == "rib":
            _analyze_rib(res, rec, bands, f_t, targets)
        else:
            _analyze_hole(res, rec, bands, fy, f_t)
    for s, b in bands.items():
        n_nan = int(np.isnan(b).sum())
        if n_nan:
            res["notes"].append("%s: %d NaN frequencies (mode not returned)"
                                % (s, n_nan))
    sc = res["scored"]
    if sc is not None and sc["found"]:
        lo, hi = res["window_GHz"]
        res["in_window"] = bool(lo <= sc["f_center_GHz"] <= hi)
    else:
        res["in_window"] = None
    # True / False / "legacy" (every difference explained by the legacy
    # breathing-gap definition) / None (nothing recorded to compare)
    cmp_ = list(res["comparison"].values())
    if not cmp_:
        res["matches"] = None
    elif all(c["matches"] for c in cmp_):
        res["matches"] = True
    elif all(c["matches"] or c.get("explained_by") for c in cmp_):
        res["matches"] = "legacy"
    else:
        res["matches"] = False
    return res


def _analyze_rib(res, rec, bands, f_t, targets):
    mode = rec.get("mech_gap_mode") or "complete"
    weights, wsrc = rib_weights(rec, targets)
    if weights is not None:
        res["window_GHz"] = [float(x) for x in weights["mech_window_GHz"]]
    if bands:
        gaps = {s: interval_gap([bands[s]], f_t, REL_TOL) for s in SECTORS
                if s in bands}
        if all(s in bands for s in SECTORS):
            gaps["complete"] = interval_gap([bands[s] for s in SECTORS], f_t,
                                            REL_TOL)
        res["gaps"] = {n: gap_dict(g) for n, g in gaps.items()}
        recorded = rec.get("mech_gaps") or {}
        for name, gp in gaps.items():
            if name in recorded:
                res["comparison"][name] = compare(recorded[name], gp)
        if mode in gaps:
            gp = gaps[mode]
            res["scored"] = dict(gap_dict(gp), name=mode)
            res["comparison"]["mechanical_gap"] = compare(
                dict(normalized_gap=rec.get("mechanical_gap"),
                     f_center=rec.get("mechanical_center_frequency")), gp)
        else:
            res["notes"].append("mech_gap_mode %r needs sectors the record "
                                "does not hold" % mode)
    else:
        rg = (rec.get("mech_gaps") or {}).get(mode)
        if rg:
            res["scored"] = dict(found=bool(rg.get("found")), name=mode,
                                 normalized_gap=rg.get("normalized_gap"),
                                 f_lower_GHz=rg.get("f_lower", 0) * 1e-9,
                                 f_upper_GHz=rg.get("f_upper", 0) * 1e-9,
                                 f_center_GHz=rg.get("f_center", 0) * 1e-9,
                                 lower_band=-1, recorded_only=True)
    tw = rec.get("mech_trusted_window_Hz")
    if tw:
        res["extra"]["trusted_window_GHz"] = [None if x is None else x * 1e-9
                                              for x in tw]
    sc = res["scored"]
    if weights is not None and sc is not None:
        from objective_rib import mech_terms
        G = sc["normalized_gap"] or 0.0
        fc = (sc["f_center_GHz"] or 0.0) * 1e9
        terms = mech_terms(G, fc, sc["found"], weights, f_t)
        ver = rec.get("score_version")
        res["penalties"] = dict(
            recomputed=terms, weights_source=wsrc,
            recorded={k: rec.get(k) for k in ("pen_size_m", "pen_window_m",
                                              "pen_center_m")})
        if rec.get("pen_window_m") is None:
            res["notes"].append(
                "pen_window_m not recorded (score version %s); shown as the "
                "weights of %s give it" % (ver, wsrc))
        elif not _close(rec["pen_window_m"], terms["pen_window_m"]):
            res["notes"].append("pen_window_m differs: recorded %.4g, "
                                "recomputed %.4g" % (rec["pen_window_m"],
                                                     terms["pen_window_m"]))


def _analyze_hole(res, rec, bands, fy, f_t):
    res["window_GHz"] = [MECH_F_MIN * 1e-9, MECH_F_MAX * 1e-9]
    if rec.get("fy_computed") is False:
        res["notes"].append("fy not computed in the run (fy = 1: every band "
                            "counts as breathing)")
    b = bands.get("sym")
    if b is None:
        if (rec.get("_freqs_m") is None and rec.get("_fy_arr") is None
                and rec.get("status") == "success"
                and str(rec.get("id", "")).startswith("cand_")):
            res["mech_gap_mode"] = "all modes (gap_near_frequency)"
            res["notes"].append(
                "hole records of src/objective.py (run_one.py / run_loop.py) "
                "store no band arrays; their gap is gap_near_frequency over "
                "all modes (centre within +/- 50 % of the target)")
        return
    if fy is None:
        fy = np.ones_like(b)
        res["notes"].append("no _fy_arr: fy = 1 assumed")
        res["fy"] = None
    window = (MECH_F_MIN * 0.5, MECH_F_MAX * 2.0)
    gb = breathing_gap(b, fy)
    gl = largest_gap(b, f_window=window)
    gi = interval_gap([b], f_t, REL_TOL)
    res["gaps"] = dict(breathing=gap_dict(gb), largest=gap_dict(gl),
                       interval=gap_dict(gi))
    res["scored"] = dict(gap_dict(gb), name="breathing")
    g_rec = rec.get("mechanical_gap", rec.get("G_m"))
    fc = rec.get("mechanical_center_frequency")
    if fc is None and rec.get("f_m_GHz") is not None:
        fc = rec["f_m_GHz"] * 1e9
    recorded = dict(normalized_gap=g_rec, f_center=fc or None)
    cmp_ = compare(recorded, gb)
    if not cmp_["matches"]:
        gold = legacy_breathing_gap(b, fy)
        if compare(recorded, gold)["matches"]:
            cmp_["explained_by"] = "legacy_breathing_gap"
            res["gaps"]["breathing_legacy"] = gap_dict(gold)
            res["notes"].append(
                "recorded G_m is the LEGACY breathing gap (both edges "
                "breathing, fy < 0.5 set to NaN point by point; "
                "run_opt_comsol.py before 2026-07-01 ~18:20), see "
                "gaps.breathing_legacy")
    res["comparison"]["mechanical_gap"] = cmp_
    if rec.get("G_m_largest") is not None:
        flc = rec.get("f_m_largest_GHz")
        res["comparison"]["G_m_largest"] = compare(
            dict(normalized_gap=rec["G_m_largest"],
                 f_center=None if not flc else flc * 1e9), gl)
    res["extra"]["fy_computed"] = rec.get("fy_computed")


# ---- outputs ----------------------------------------------------------------

def _f(x, fmt="%.2f", none="-"):
    if x is None:
        return none
    try:
        if not math.isfinite(float(x)):
            return none
    except (TypeError, ValueError):
        return none
    return fmt % x


def geometry_text(res):
    p = res["params"]
    nm = lambda key: _f(p.get(key) * 1e9 if p.get(key) is not None else None,  # noqa: E731
                        "%.0f")
    if res["celltype"] == "rib":
        keys = ("a", "w", "t", "s", "th", "r1", "r2")
        return "  ".join("%s=%s" % (k, nm(k)) for k in keys if k in p) + " nm"
    keys = ("a", "w", "t", "hx", "hy")
    return ("  ".join("%s=%s" % (k, nm(k)) for k in keys if k in p)
            + " nm (hx, hy semi-axes)")


def write_bands_csv(path, res):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sector", "k_a_over_pi", "band", "freq_GHz", "fy"])
        if res["k"] is None:
            return
        for s, b in res["bands"].items():
            fy = res["fy"] if s == "sym" else None
            for i, kk in enumerate(res["k"]):
                for j in range(b.shape[1]):
                    f = b[i, j]
                    w.writerow([s, "%.6g" % kk, j,
                                "" if not np.isfinite(f) else "%.9g" % (f * 1e-9),
                                "" if fy is None or not np.isfinite(fy[i, j])
                                else "%.6g" % fy[i, j]])


def summary_of(res):
    p = res["params"]
    return dict(
        record_id=res["record_id"], celltype=res["celltype"],
        status=res["status"], source=res["source"],
        saved=(datetime.datetime.fromtimestamp(res["ts"]).isoformat(
            timespec="seconds") if res["ts"] else None),
        mech_backend=res["mech_backend"], optical_backend=res["optical_backend"],
        params_nm={k: v * 1e9 for k, v in p.items() if isinstance(v, (int, float))},
        score=res["score"], score_version=res["score_version"],
        mech_gap_mode=res["mech_gap_mode"],
        target_GHz=res["target_GHz"], window_GHz=res["window_GHz"],
        has_bands=res["k"] is not None,
        n_k=None if res["k"] is None else int(res["k"].size),
        n_bands={s: int(b.shape[1]) for s, b in res["bands"].items()},
        recorded=res["recorded"], scored_gap=res["scored"],
        in_window=res["in_window"], gaps=res["gaps"],
        comparison=res["comparison"], matches_record=res["matches"],
        penalties=res["penalties"], extra=res["extra"], notes=res["notes"])


def plot(res, path, ymax=None):
    """Band diagram: GHz vs k a/pi, sectors by marker, gaps and window."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(10.0, 6.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    k, bands, fy = res["k"], res["bands"], res["fy"]
    lo_w, hi_w = res["window_GHz"]
    finite = [b[np.isfinite(b)] for b in bands.values()]
    f_top = max([float(x.max()) * 1e-9 for x in finite if x.size] + [hi_w * 1.2])
    top = ymax or f_top * 1.03
    handles = []

    ax.axhspan(lo_w, hi_w, color=CONE, lw=0, zorder=0)
    handles.append(Patch(color=CONE, label="centre window %g-%g GHz" % (lo_w, hi_w)))
    ax.axhline(res["target_GHz"], color=INK2, ls=(0, (4, 3)), lw=1.2, zorder=1)
    handles.append(Line2D([], [], color=INK2, ls=(0, (4, 3)), lw=1.2,
                          label="target %g GHz" % res["target_GHz"]))

    sc = res["scored"]
    if sc is not None and sc["found"]:
        ax.axhspan(sc["f_lower_GHz"], sc["f_upper_GHz"], color=GAP, alpha=0.28,
                   lw=0, zorder=1)
        name = {"complete": "complete gap", "sym": "y-even gap (scored)",
                "antisym": "y-odd gap (scored)",
                "breathing": "breathing gap (scored)"}.get(sc["name"], sc["name"])
        handles.append(Patch(color=GAP, alpha=0.4, label="%s %.1f %%"
                             % (name, 100 * sc["normalized_gap"])))

    # per-sector gaps: bars at the right edge, in the sector colour
    if res["celltype"] == "rib":
        for x, s in ((0.962, "sym"), (0.985, "antisym")):
            g = res["gaps"].get(s)
            if g and g["found"]:
                ax.plot([x, x], [g["f_lower_GHz"], g["f_upper_GHz"]],
                        color=SECTOR_COLOR[s], lw=5, solid_capstyle="butt",
                        zorder=6, transform=ax.get_yaxis_transform())
                handles.append(Line2D([], [], color=SECTOR_COLOR[s], lw=5,
                                      label="%s gap %.1f %% (bar)"
                                      % ({"sym": "y-even",
                                          "antisym": "y-odd"}[s],
                                         100 * g["normalized_gap"])))
    else:
        g = res["gaps"].get("largest")
        if g and g["found"]:
            ax.plot([0.985, 0.985], [g["f_lower_GHz"], g["f_upper_GHz"]],
                    color=INK2, lw=5, solid_capstyle="butt", zorder=6,
                    transform=ax.get_yaxis_transform())
            handles.append(Line2D([], [], color=INK2, lw=5,
                                  label="largest gap, any band %.1f %% (bar)"
                                  % (100 * g["normalized_gap"])))
    tw = res["extra"].get("trusted_window_GHz")
    if tw and tw[1] is not None and tw[1] < top:
        ax.axhline(tw[1], color=MUTED, ls=":", lw=1, zorder=1)
        handles.append(Line2D([], [], color=MUTED, ls=":", lw=1,
                              label="top of trusted window"))

    if k is not None:
        for s, b in bands.items():
            kk = np.repeat(k[:, None], b.shape[1], axis=1)
            ok = np.isfinite(b)
            col, mk = SECTOR_COLOR[s], SECTOR_MARKER[s]
            label = SECTOR_LABEL[s]
            if res["celltype"] == "hole":
                label = "study_mech bands"
            if res["celltype"] == "hole" and fy is not None:
                br = ok & (fy >= FY_THRESH)
                ax.scatter(kk[br], b[br] * 1e-9, s=26, marker=mk, color=col,
                           edgecolors=SURFACE, linewidths=0.6, zorder=4)
                ax.scatter(kk[ok & ~br], b[ok & ~br] * 1e-9, s=26, marker=mk,
                           facecolors="none", edgecolors=col, linewidths=1.0,
                           zorder=4)
                handles += [Line2D([], [], ls="", marker=mk, color=col, ms=6,
                                   label="mode with fy >= %.1f (breathing-like)"
                                   % FY_THRESH),
                            Line2D([], [], ls="", marker=mk, mfc="none",
                                   mec=col, ms=6, label="mode with fy < %.1f"
                                   % FY_THRESH)]
            elif s == "sym":
                ax.scatter(kk[ok], b[ok] * 1e-9, s=26, marker=mk, color=col,
                           edgecolors=SURFACE, linewidths=0.6, zorder=4)
                handles.append(Line2D([], [], ls="", marker=mk, color=col,
                                      ms=6, label=label + ", filled"))
            else:
                ax.scatter(kk[ok], b[ok] * 1e-9, s=30, marker=mk,
                           facecolors="none", edgecolors=col, linewidths=1.1,
                           zorder=5)
                handles.append(Line2D([], [], ls="", marker=mk, mfc="none",
                                      mec=col, ms=6, label=label + ", hollow"))
    else:
        ax.text(0.5, 0.5, "no band arrays in this record\n(summary only)",
                transform=ax.transAxes, ha="center", va="center", color=INK2)

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, top)
    ax.set_xlabel(r"wavevector $k_x a/\pi$   ($\Gamma$ at 0, X at 1)", color=INK2)
    ax.set_ylabel("frequency [GHz]", color=INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.grid(axis="y", color=CONE, lw=0.6, zorder=0)

    head = "%s cell  |  %s  |  %s" % (res["celltype"], res["record_id"],
                                      res["mech_gap_mode"])
    if res["score"] is not None:
        head += "  |  score %+.4f" % res["score"]
        if res["score_version"] is not None:
            head += " (v%s)" % res["score_version"]
    if sc is not None and sc["found"]:
        line = ("gap %.3f-%.3f GHz, centre %.3f GHz, %.2f %%, %s window"
                % (sc["f_lower_GHz"], sc["f_upper_GHz"], sc["f_center_GHz"],
                   100 * sc["normalized_gap"],
                   "inside" if res["in_window"] else "OUTSIDE"))
    else:
        line = "no scored gap"
    pen = (res["penalties"] or {}).get("recomputed")
    if pen:
        line += "; pen_window_m %.3g, pen_size_m %.3g" % (pen["pen_window_m"],
                                                         pen["pen_size_m"])
    if res["matches"] is False:
        line += "\nDIFFERS from the record: " + ", ".join(
            n for n, c in res["comparison"].items() if not c["matches"])
    elif res["matches"] == "legacy":
        g = res["gaps"].get("breathing_legacy") or {}
        line += ("\nrecorded G_m %.2f %% = legacy breathing-gap definition "
                 "(both edges breathing)" % (100 * (g.get("normalized_gap") or 0)))
    ax.set_title(head + "\n" + geometry_text(res) + "\n" + line, fontsize=9,
                 color=INK, loc="left")
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                    fontsize=7.5, frameon=False, borderaxespad=0.0)
    for text in leg.get_texts():
        text.set_color(INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def _gap_of(res):
    """(G, f_c GHz) to summarise: the recomputed scored gap, else recorded."""
    sc = res["scored"]
    if sc is not None and not sc.get("recorded_only"):
        return ((sc["normalized_gap"], sc["f_center_GHz"]) if sc["found"]
                else (0.0, None))
    r = res["recorded"]
    return r["mechanical_gap"], r["mechanical_center_GHz"] or None


def plot_summary(results, path):
    """G_m vs centre, and G_m vs each geometry variable, coloured by score."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    cells = sorted({r["celltype"] for r in results})
    keys = []
    if len(cells) == 1:
        keys = [k for k in ("a", "w", "t", "s", "th", "hx", "hy")
                if all(k in r["params"] for r in results)
                and np.ptp([r["params"][k] for r in results]) > 0]
    n = 1 + len(keys)
    ncol = min(n, 3)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.8 * nrow),
                             facecolor=SURFACE, squeeze=False)
    cmap = LinearSegmentedColormap.from_list("score", SEQ_RAMP)
    G = np.array([np.nan if _gap_of(r)[0] is None else 100 * _gap_of(r)[0]
                  for r in results], dtype=float)
    fc = np.array([np.nan if _gap_of(r)[1] is None else _gap_of(r)[1]
                   for r in results], dtype=float)
    S = np.array([np.nan if r["score"] is None else r["score"] for r in results],
                 dtype=float)
    bad = np.array([r["matches"] is False for r in results])
    fin = np.isfinite(S)
    # scores of old versions can be far below the rest: clip the colour
    # scale at the 20th percentile (lower scores saturate to the lightest)
    vmin, vmax = ((np.nanpercentile(S, 20), np.nanmax(S)) if fin.any()
                  else (0.0, 1.0))
    if not vmax > vmin:
        vmax = vmin + 1.0
    lo_w, hi_w = results[0]["window_GHz"]

    def draw(ax, x, xlabel):
        ax.set_facecolor(SURFACE)
        ok = np.isfinite(x) & np.isfinite(G)
        if (ok & ~fin).any():
            ax.scatter(x[ok & ~fin], G[ok & ~fin], s=40, facecolors="none",
                       edgecolors=MUTED, label="no score")
        sc = ax.scatter(x[ok & fin], G[ok & fin], s=44, c=S[ok & fin], cmap=cmap,
                        vmin=vmin, vmax=vmax, edgecolors=SURFACE, linewidths=0.8,
                        zorder=3)
        if (ok & bad).any():
            ax.scatter(x[ok & bad], G[ok & bad], s=90, marker="x", color=INK,
                       zorder=4, label="differs from record")
        ax.set_xlabel(xlabel, color=INK2, fontsize=9)
        ax.set_ylabel("mechanical gap [%]", color=INK2, fontsize=9)
        ax.tick_params(colors=INK2, labelsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
        ax.grid(color=CONE, lw=0.6, zorder=0)
        return sc

    ax0 = axes[0][0]
    ax0.axvspan(lo_w, hi_w, color=CONE, lw=0, zorder=0,
                label="window %g-%g GHz" % (lo_w, hi_w))
    sc = draw(ax0, fc, "gap centre [GHz]")
    ax0.legend(fontsize=7, loc="best", edgecolor=CONE)
    ax0.set_title("%d record(s), %s" % (len(results), "/".join(cells)),
                  fontsize=9, color=INK, loc="left")
    for i, key in enumerate(keys, start=1):
        ax = axes[i // ncol][i % ncol]
        x = np.array([r["params"][key] * 1e9 for r in results], dtype=float)
        draw(ax, x, "%s [nm]" % key)
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].set_visible(False)
    if fin.any():
        cb = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02,
                          extend="min")
        cb.set_label("score (each record's version; clipped at the 20th "
                     "percentile)", color=INK2, fontsize=8)
        cb.ax.tick_params(labelsize=7, colors=INK2)
        cb.outline.set_visible(False)
    fig.savefig(path, dpi=130, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


INDEX_COLS = ["record_id", "celltype", "status", "saved", "mech_backend",
              "score", "score_version", "mech_gap_mode", "has_bands", "n_k",
              "recorded_gap", "recorded_center_GHz", "gap", "f_lower_GHz",
              "f_upper_GHz", "f_center_GHz", "in_window", "pen_window_m",
              "sym_gap", "antisym_gap", "complete_gap", "largest_gap",
              "matches_record", "differs", "outdir"]


def index_row(s, outdir):
    g = s["scored_gap"] or {}
    gaps = s["gaps"] or {}
    pen = (s["penalties"] or {}).get("recomputed") or {}
    gv = lambda n: (gaps.get(n) or {}).get("normalized_gap")  # noqa: E731
    return dict(
        record_id=s["record_id"], celltype=s["celltype"], status=s["status"],
        saved=s["saved"], mech_backend=s["mech_backend"], score=s["score"],
        score_version=s["score_version"], mech_gap_mode=s["mech_gap_mode"],
        has_bands=int(s["has_bands"]), n_k=s["n_k"],
        recorded_gap=s["recorded"]["mechanical_gap"],
        recorded_center_GHz=s["recorded"]["mechanical_center_GHz"],
        gap=g.get("normalized_gap"), f_lower_GHz=g.get("f_lower_GHz"),
        f_upper_GHz=g.get("f_upper_GHz"), f_center_GHz=g.get("f_center_GHz"),
        in_window="" if s["in_window"] is None else int(s["in_window"]),
        pen_window_m=pen.get("pen_window_m"),
        sym_gap=gv("sym"), antisym_gap=gv("antisym"),
        complete_gap=gv("complete"), largest_gap=gv("largest"),
        matches_record=({None: "", True: 1, False: 0}.get(s["matches_record"],
                                                          s["matches_record"])),
        differs=";".join("%s:%s" % (n, "/".join(c["differs"]))
                         for n, c in s["comparison"].items() if not c["matches"]),
        outdir=outdir)


def write_index(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=INDEX_COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})


def print_table(summaries):
    head = ("%-18s %-4s %-8s %-9s %5s  %8s %8s  %-15s %6s %5s  %s"
            % ("id", "cell", "status", "mode", "bands", "G_m rec", "G_m now",
               "edges GHz", "f_c", "win", "record"))
    print(head)
    print("-" * len(head))
    for s in summaries:
        g = s["scored_gap"] or {}
        mode = s["mech_gap_mode"] or "-"
        mode = "breath" if mode.startswith("breathing") else mode
        now = g.get("normalized_gap") if not g.get("recorded_only") else None
        edges = ("%.2f-%.2f" % (g["f_lower_GHz"], g["f_upper_GHz"])
                 if g.get("found") else "-")
        match = {None: "n/a", True: "match", False: "DIFFERS",
                 "legacy": "legacy def"}[s["matches_record"]]
        print("%-18s %-4s %-8s %-9s %5s  %8s %8s  %-15s %6s %5s  %s" % (
            s["record_id"], s["celltype"], (s["status"] or "-")[:8], mode[:9],
            "yes" if s["has_bands"] else "no",
            _f(None if s["recorded"]["mechanical_gap"] is None
               else 100 * s["recorded"]["mechanical_gap"], "%.2f%%"),
            _f(None if now is None else 100 * now, "%.2f%%"), edges,
            _f(g.get("f_center_GHz") if g.get("found") else None, "%.2f"),
            {None: "-", True: "in", False: "OUT"}[s["in_window"]], match))


def print_notes(summaries, max_each=3):
    """Notes per record; a note shared by more than max_each records is
    printed once with its count (every note is in the summary.json files)."""
    by_text = {}
    for s in summaries:
        for n in s["notes"]:
            by_text.setdefault(n, []).append(s["record_id"])
    for text, ids in by_text.items():
        if len(ids) > max_each:
            print("  %d records (e.g. %s): %s" % (len(ids), ids[0], text))
        else:
            for cid in ids:
                print("  %s: %s" % (cid, text))


# ---- CLI --------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ids", nargs="*", help="record ids (or unique prefixes)")
    sel = ap.add_argument_group("selection")
    sel.add_argument("--latest", type=int, default=None, metavar="N",
                     help="the N most recently saved records that pass the filters")
    sel.add_argument("--all", action="store_true",
                     help="every record that passes the filters")
    sel.add_argument("--celltype", choices=("hole", "rib"), default=None,
                     help="read the database of this cell type [both]")
    sel.add_argument("--db", default=None,
                     help="SQLite file to read instead of the default one(s)")
    sel.add_argument("--json", default=None, metavar="PATH",
                     help="records from an out_json (no band arrays there: "
                          "they are looked up in SQLite by id)")
    flt = ap.add_argument_group("filters (for --latest / --all)")
    flt.add_argument("--status", default="success",
                     help="record status to keep, or 'any' [success]")
    flt.add_argument("--mech-backend", default="comsol",
                     help="mech_backend to keep, or 'any' [comsol; drops "
                          "surrogate_stub]")
    flt.add_argument("--since", default=None, metavar="DATE",
                     help="saved on/after this date (YYYY-MM-DD[ HH:MM], local)")
    flt.add_argument("--min-gap", type=float, default=None, metavar="G",
                     help="recorded mechanical_gap >= G (fraction, e.g. 0.1)")
    out = ap.add_argument_group("output")
    out.add_argument("--outdir", default=DEFAULT_OUTDIR,
                     help="output folder; one sub-folder per record "
                          "[python-scripts/results/mech_analysis]")
    out.add_argument("--ymax", type=float, default=None,
                     help="top of the frequency axis [GHz]")
    out.add_argument("--no-plot", action="store_true",
                     help="CSV / JSON / table only")
    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    if not (args.ids or args.all or args.latest is not None):
        ap.error("give record ids, --latest N or --all")
    if args.latest is not None and args.latest < 1:
        ap.error("--latest must be >= 1")
    try:
        with network_blocked():
            if args.db:
                db_paths = [args.db]
            elif args.celltype:
                db_paths = [default_db(args.celltype)]
            else:
                db_paths = [p for p in (default_db("rib"), default_db("hole"))
                            if os.path.isfile(p)]
            if args.json:
                entries = read_json(args.json, db_paths)
            else:
                if not db_paths:
                    raise AnalysisError("no database found")
                entries = []
                for p in db_paths:
                    entries.extend(read_db(p))
            jobs = select(entries, args)
            if not jobs:
                print("no record matches the selection")
                return 1
            if args.ids and (args.since or args.min_gap is not None):
                print("note: filters apply to --latest / --all, not to ids",
                      file=sys.stderr)
            targets = load_targets()
            os.makedirs(args.outdir, exist_ok=True)
            summaries, rows, results = [], [], []
            for entry in jobs:
                res = analyze(entry, targets)
                outdir = os.path.join(args.outdir, str(res["record_id"]))
                os.makedirs(outdir, exist_ok=True)
                write_bands_csv(os.path.join(outdir, "bands.csv"), res)
                summary = summary_of(res)
                with open(os.path.join(outdir, "summary.json"), "w") as fh:
                    json.dump(summary, fh, indent=2)
                if not args.no_plot:
                    plot(res, os.path.join(outdir, "band_diagram.png"),
                         ymax=args.ymax)
                summaries.append(summary)
                rows.append(index_row(summary, outdir))
                results.append(res)
            print_table(summaries)
            print_notes(summaries)
            if len(jobs) > 1:
                path = os.path.join(args.outdir, "mech_index.csv")
                write_index(path, rows)
                print("index: %s" % path)
                if not args.no_plot:
                    path = os.path.join(args.outdir, "mech_summary.png")
                    plot_summary(results, path)
                    print("summary figure: %s" % path)
            n_bands = sum(1 for s in summaries if s["has_bands"])
            n_diff = sum(1 for s in summaries if s["matches_record"] is False)
            n_leg = sum(1 for s in summaries if s["matches_record"] == "legacy")
            print("%d record(s): %d with band arrays, %d differ from the "
                  "record%s -> %s" % (
                      len(summaries), n_bands, n_diff,
                      " (+%d explained by the legacy breathing-gap definition)"
                      % n_leg if n_leg else "", args.outdir))
    except AnalysisError as err:
        print("plot_mech_bands: %s" % err, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
