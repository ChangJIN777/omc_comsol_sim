#!/usr/bin/env python3
"""Band diagrams of Tidy3D results, from the saved analysis (CSV) or the .hdf5.

What `optical: tidy3d` leaves in results/tidy3d/<celltype>/<record id>/:
  bands.csv        every resonance: k, frequency, Q (, amplitude, error)
  params.json      the BandParams of the run (geometry, index, window, filters)
  evaluation.json  the scored gap, the k points, the scoring settings
  plan.json        the simulation plan
  *.hdf5           the raw results (monitor time traces), one per k point.
                   With tidy3d.keep_hdf5: best | none they are deleted for
                   all cells but the best one, or for all.
This script plots the band structure of a record with the gap exactly as the
pipeline scores it.

The CSV is the normal input and needs NO tidy3d: any Python of the repo works.
The .hdf5 files are read, and ResonanceFinder is run again, only
  * with --reanalyze,
  * with a resonance-filter option (--minQ, --minamp, --maxerr, --freq-*-thz),
  * for .hdf5 files given directly, or
  * when a folder has no bands.csv;
that needs tidy3d (the `tidy` conda environment) and the files.

LOCAL FILES ONLY: nothing is uploaded or downloaded, the Tidy3D web API is
never called, and network connections are blocked while the script runs. The
result files and the databases are not modified (the sqlite record of a run
is read, read-only, for its optical_tol_thz); output goes to a new
`analysis/` folder (or --outdir).

Scoring is the repo's own (src/optical_tidy3d.score_resonances ->
comsol_optics.guided_interval_gaps / select_gap_by_value / optical_gap_tracked).
A rib record is scored as its run scored it: score version 4 (current: the
gap with the largest G_o - size - detuning - light-cone penalty, flat-top
detuning penalty of half-width optical_freq_tol_THz; weights of the run or
of configs/targets.yaml), version 3 (the same with a quadratic detuning
penalty) or version 2 (older records: Gaussian weight with optical_tol_thz);
--score-version overrides. With tidy3d, diamond_bandStructure_sim.py is reused for
the .hdf5 path (load_run, bloch_k, ResonanceFinder + filter_resonances,
BandParams) and for the exact unit-cell outline; without it the outline is
drawn from the saved dimensions with sharp corners.

Examples (from python-scripts/):
  python scripts/plot_tidy3d_bands.py cand_6551b0da7ff6          # from bands.csv
  python scripts/plot_tidy3d_bands.py rib/cand_6551b0da7ff6 --k-units ka/pi
  python scripts/plot_tidy3d_bands.py results/tidy3d/rib/cand_6551b0da7ff6/bands.csv
  python scripts/plot_tidy3d_bands.py some/bands.csv --a-um 0.57 --celltype rib
  python scripts/plot_tidy3d_bands.py --all                      # every record + index
  python scripts/plot_tidy3d_bands.py cand_3e6531b03ecc --reanalyze      # tidy env
  python scripts/plot_tidy3d_bands.py cand_3e6531b03ecc --minQ 1000      # tidy env
  python scripts/plot_tidy3d_bands.py path/a.hdf5 path/b.hdf5 --outdir out
  python scripts/plot_tidy3d_bands.py cand_6551b0da7ff6 --min-guided-q 1000

Inputs: a bands.csv, .hdf5 files, a run folder, a record id, or
<celltype>/<record id> under --results-dir. .hdf5 files are grouped by the
Bloch wavevector stored in them (the newest wins for a repeated k); file
names do not matter and non-result files such as batch.hdf5 are skipped.
--min-guided-q and --light-line-margin work on the CSV (they use Q and the
light line), the other filters do not.

Outputs per record: band_diagram.png, resonances.csv (k, f, Q, amplitude,
error, guided, in_window) and summary.json; with --all also an index CSV.
"""
import argparse
import contextlib
import csv
import dataclasses
import json
import os
import socket
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
DEFAULT_RESULTS = os.path.join(_ROOT, "results", "tidy3d")
C0 = 299_792_458.0
F_TARGET = C0 / 1550e-9                      # 193.4 THz, as targets.yaml
DEFAULT_TOL_THZ = 2.0                        # rib.optical_tol_thz
CURRENT_VERSION = 4                          # objective_rib.SCORE_VERSION

# colours: one sequential blue ramp for Q, the rest by role (shape and fill
# also differ, so nothing is told apart by colour alone)
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
SURFACE, CONE = "#fcfcfb", "#e9e8e4"
Q_RAMP = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
GAP, XGAP, TARGET = "#eda100", "#1baf7a", "#eb6834"


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


def load_scoring():
    """The repo's scoring (src/optical_tidy3d, numpy only: no tidy3d needed)."""
    sys.path.insert(0, os.path.join(_ROOT, "src"))
    import optical_tidy3d as ot
    return ot


class ModuleLoader:
    """diamond_bandStructure_sim (which imports tidy3d), loaded on first use.

    CSV / JSON inputs never need it. It is needed to read .hdf5 files
    (ResonanceFinder) and, if there, used for the exact unit-cell outline.
    """

    def __init__(self, ot):
        self.ot, self._bs = ot, None

    @property
    def available(self):
        return bool(self.ot.have_tidy3d())

    def get(self, why):
        if self._bs is None:
            if not self.available:
                raise AnalysisError(
                    f"tidy3d is needed to {why}. Run this script with the "
                    "Python of the `tidy` conda environment, e.g.\n"
                    "  conda run -n tidy python scripts/plot_tidy3d_bands.py ...\n"
                    "(records with a saved bands.csv are plotted without tidy3d: "
                    "drop --reanalyze and the resonance-filter options)")
            self._bs = self.ot.load_bs()
        return self._bs


# material names of diamond_bandStructure_sim.MATERIAL_INDEX
MATERIAL_INDEX = {"diamond_telecom": 2.386, "diamond_visible": 2.4064}
MODE_SYMMETRY = {"TE": (0, -1, 1), "TM": (0, 1, -1)}


class LiteParams:
    """The parameters the plot and the scoring need, from the saved JSON
    files of a record (no tidy3d). Lengths in um, as BandParams."""

    def __init__(self, d):
        nan = float("nan")
        self.geometry = d.get("geometry")
        self.mode = d.get("mode") or "TE"
        self.xsect = d.get("xsect") or "rect"
        idx = d.get("index")
        self.n = (float(idx) if idx is not None
                  else MATERIAL_INDEX.get(d.get("material", "diamond_telecom"), nan))
        self.a = float(d["a"])
        self.target_freq_thz = float(d.get("target_freq_thz") or F_TARGET * 1e-12)
        span = float(d.get("freq_span_thz") or 100.0)
        self.freq_range = ((self.target_freq_thz - span) * 1e12,
                           (self.target_freq_thz + span) * 1e12)
        self.minQ = float(d.get("minQ", 1.0))
        self.minamp = float(d.get("minamp", 1.0))
        self.maxerr = float(d.get("maxerr", 1000.0))
        for name in ("thickness", "t", "w", "spine_width", "beam_width",
                     "hole_x", "hole_y"):
            setattr(self, name, float(d[name]) if d.get(name) is not None else nan)
        self.fillet_r1 = float(d.get("fillet_r1") or 0.0)
        self.fillet_r2 = float(d.get("fillet_r2") or 0.0)

    @property
    def symmetry(self):
        sx, sy, sz = MODE_SYMMETRY.get(self.mode, (0, 0, 0))
        return (sx, sy, 0) if self.xsect == "wedge" else (sx, sy, sz)


# ---- inputs -----------------------------------------------------------------

def _has_hdf5(folder):
    return os.path.isdir(folder) and any(
        n.endswith(".hdf5") and n != "batch.hdf5" for n in os.listdir(folder))


def _has_cache(folder):
    """The saved analysis of a record: bands.csv (k, f, Q per resonance). It
    is what stays when the .hdf5 files are pruned (tidy3d.keep_hdf5)."""
    return os.path.isfile(os.path.join(folder, "bands.csv"))


def _has_results(folder):
    return os.path.isdir(folder) and (_has_hdf5(folder) or _has_cache(folder))


def rows_from_bands_csv(path):
    """Resonances saved by the pipeline, per k.

    bands.csv columns: k_2pi_over_a, freq_THz, Q, below_light_line and, in
    files written since 2026-11-24, amplitude, error. Older files have no
    amplitude / fit error: NaN. A file without k / frequency columns is refused.
    """
    def num(row, key):
        try:
            return float(row[key])
        except (KeyError, TypeError, ValueError):
            return float("nan")

    by_k = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"k_2pi_over_a", "freq_THz"} - set(reader.fieldnames or ())
        if missing:
            raise AnalysisError(f"{path} is not a bands.csv of the pipeline "
                                f"(no column {sorted(missing)})")
        for row in reader:
            by_k.setdefault(float(row["k_2pi_over_a"]), []).append(
                (float(row["freq_THz"]) * 1e12, num(row, "Q"),
                 num(row, "amplitude"), num(row, "error")))
    if not by_k:
        raise AnalysisError(f"{path} holds no resonances")
    rows = []
    for k in sorted(by_k):
        pairs = sorted(by_k[k])
        rows.append(dict(k=k, freq=np.array([x[0] for x in pairs]),
                         Q=np.array([x[1] for x in pairs]),
                         amplitude=np.array([x[2] for x in pairs]),
                         error=np.array([x[3] for x in pairs])))
    return rows


def record_folders(results_dir):
    """Every <results_dir>/<celltype>/<record id>/ that holds result files."""
    out = []
    if os.path.isdir(results_dir):
        for cell in sorted(os.listdir(results_dir)):
            sub = os.path.join(results_dir, cell)
            if os.path.isdir(sub):
                out += [os.path.join(sub, r) for r in sorted(os.listdir(sub))
                        if _has_results(os.path.join(sub, r))]
    return out


def resolve_inputs(inputs, results_dir, all_records=False):
    """-> list of jobs dict(folder=..., [files=[...]], [csv=path]).

    An input is a bands.csv file, .hdf5 file(s), a folder, a record id, or
    <celltype>/<record id> under results_dir. A bands.csv is a job of its own
    (its folder is read for the JSON files next to it); all .hdf5 files given
    directly form ONE job.
    """
    jobs, files = [], []
    if all_records:
        jobs += [dict(folder=f) for f in record_folders(results_dir)]
    for item in inputs:
        if os.path.isfile(item) and item.lower().endswith(".csv"):
            path = os.path.abspath(item)
            jobs.append(dict(folder=os.path.dirname(path), csv=path))
        elif os.path.isfile(item):
            files.append(os.path.abspath(item))
        elif os.path.isdir(item):
            jobs.append(dict(folder=os.path.abspath(item)))
        else:
            cands = [os.path.join(results_dir, item)]
            if os.path.isdir(results_dir):
                cands += [os.path.join(results_dir, c, item)
                          for c in sorted(os.listdir(results_dir))]
            hits = [c for c in cands if os.path.isdir(c)]
            if not hits:
                raise AnalysisError(
                    f"{item!r} is not a file, a folder, or a record under "
                    f"{results_dir} (records: "
                    f"{[os.path.basename(f) for f in record_folders(results_dir)]})")
            jobs.append(dict(folder=os.path.abspath(hits[0])))
    if files:
        jobs.append(dict(files=files, folder=os.path.dirname(files[0])))
    seen, unique = set(), []
    for job in jobs:
        key = (job["folder"], tuple(job.get("files", ())), job.get("csv"))
        if key not in seen:
            seen.add(key)
            unique.append(job)
    for job in unique:
        if ("files" not in job and "csv" not in job
                and not _has_results(job["folder"])):
            raise AnalysisError("no saved analysis (bands.csv) and no Tidy3D "
                                f"result files (*.hdf5) in {job['folder']}")
    if not unique:
        raise AnalysisError("nothing to analyze: give a bands.csv, .hdf5 files, a "
                            "run folder, a record id, or --all (see --help)")
    return unique


def load_files(bs, files):
    """As bs.load_run, for an explicit list of files: {"k=...": SimulationData}."""
    td = bs.td
    found = {}
    for path in files:
        try:
            sim_data = td.SimulationData.from_file(path)
        except Exception:
            print(f"  skipped (not a simulation result): {path}", file=sys.stderr)
            continue
        key = "k=%.6f" % bs.bloch_k(sim_data.simulation)
        mtime = os.path.getmtime(path)
        if key in found and found[key][0] > mtime:
            continue
        found[key] = (mtime, sim_data)
    if not found:
        raise AnalysisError("none of the files is a Tidy3D SimulationData")
    return {key: sd for key, (_, sd) in sorted(found.items())}


def _read_json(path):
    if os.path.exists(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            pass
    return None


def params_for(bs, folder, sim):
    """(BandParams, source). params.json of the run, else read from the
    simulation itself: period, slab thickness, index, mode, window."""
    p = bs._load_saved_params(folder) if folder else None
    if p is not None:
        return p, "params.json"
    a = float(sim.size[0])
    comp = sim.monitors[0].fields[0] if sim.monitors else "Ey"
    mode = {v[0]: k for k, v in bs.MODES.items()}.get(comp, "TE")
    eps = max([float(s.medium.permittivity) for s in sim.structures] + [1.0])
    src = sim.sources[0].source_time if sim.sources else None
    kw = dict(a=a, mode=mode, index=float(np.sqrt(eps)),
              xsect="rect" if sim.symmetry[2] != 0 else "wedge")
    if src is not None:
        kw.update(target_freq_thz=float(src.freq0) * 1e-12,
                  freq_span_thz=float(src.fwidth) / 0.6 * 1e-12)
    return bs.BandParams(**kw), "simulation (no params.json)"


def celltype_of(p, source, folder, evaluation):
    if evaluation and evaluation.get("celltype"):
        return evaluation["celltype"]
    if p.geometry and not source.startswith("simulation"):
        return "hole" if p.geometry == "holey" else "rib"
    parent = os.path.basename(os.path.dirname(folder or ""))
    return parent if parent in ("hole", "rib") else None


# ---- analysis ---------------------------------------------------------------

def resonances_of(bs, sim_datas, p):
    """Per k: dict(k, freq, Q, amplitude, error), sorted by frequency.

    The module's own steps: ResonanceFinder on the monitor time traces inside
    the frequency window, then filter_resonances (minQ, minamp, maxerr).
    """
    rows = []
    for sim_data in sim_datas.values():
        rd = bs.ResonanceFinder(freq_window=tuple(p.freq_range)).run(signals=sim_data.data)
        rd = bs.filter_resonances(rd, p.minQ, p.minamp, p.maxerr)
        cols = [np.atleast_1d(getattr(rd, name).to_numpy()).astype(float)
                for name in ("freq", "Q", "amplitude", "error")]
        order = np.argsort(cols[0])
        rows.append(dict(k=bs.bloch_k(sim_data.simulation),
                         freq=cols[0][order], Q=cols[1][order],
                         amplitude=cols[2][order], error=cols[3][order]))
    rows.sort(key=lambda d: d["k"])
    return rows


def _gap_dict(gp, penalty=None):
    out = dict(found=bool(gp.found), normalized_gap=float(gp.normalized_gap),
               f_lower_THz=float(gp.f_lower) * 1e-12,
               f_upper_THz=float(gp.f_upper) * 1e-12,
               f_center_THz=float(gp.f_center) * 1e-12)
    if penalty is not None:
        out.update(freq_penalty=float(penalty),
                   effective_gap=float(gp.normalized_gap * penalty))
    return out


def score(ot, celltype, rows, a_m, p, tol_thz, min_guided_q=None,
          light_line_margin=None, version=4, weights=None):
    """The pipeline's gap for these resonances (optical_tidy3d.score_resonances).

    rib, score version 3 or 4 (4 is current): among all guided gaps (sweep),
    or the gap at X (x_only), the one with the largest optical score
    contribution G_o - size - detuning - light-cone penalty
    (objective_rib.optical_terms with `weights`, which carry the form of the
    detuning penalty: flat-top for 4, quadratic for 3); the terms of the
    chosen gap are returned.
    rib, score version 2 (records of 2026-11-23/24): Gaussian-weighted choice,
    weight exp(-((f_c - f_t)/tol)^2). hole / unknown: the tracked
    fundamental gap with its +/- 30 % centre window. Two analysis knobs that
    are NOT in the pipeline drop suspect guided resonances first: min_guided_q
    (|Q| below it) and light_line_margin (within that fraction below the
    light line, where FDTD shows low-Q peaks that are not bands).
    """
    from comsol_optics import freq_penalty
    bands = []
    for d in rows:
        f = d["freq"]
        ll = d["k"] * C0 / a_m
        drop = np.zeros(f.size, dtype=bool)
        if min_guided_q:
            drop |= (f < ll) & (np.abs(d["Q"]) < min_guided_q)
        if light_line_margin:
            drop |= (f < ll) & (f > ll * (1.0 - light_line_margin))
        f = f[~drop]
        bands.append(dict(k=d["k"], freq=f))
    k_frac, freqs = ot.bands_to_arrays(bands)
    f_t = float(p.target_freq_thz) * 1e12
    k_mode = "x_only" if len(k_frac) == 1 else "sweep"
    out = {}
    v3 = celltype == "rib" and version >= 3
    tol_hz = tol_thz * 1e12 if (celltype == "rib" and not v3) else None
    gap_value = None
    if v3:
        from objective_rib import optical_terms

        def gap_value(gap):
            return optical_terms(gap.normalized_gap, gap.f_center, f_t, a_m,
                                 weights)["value"]
    try:
        gp, method = ot.score_resonances(
            celltype or "hole", k_frac, freqs, a_m, f_target=f_t, k_mode=k_mode,
            window_hz=tuple(p.freq_range), optical_tol_hz=tol_hz, out=out,
            gap_value=gap_value)
    except ValueError as err:
        return dict(found=False, not_scored=str(err), method=None, candidates=[])
    if v3:
        res = _gap_dict(gp)
        # a gap that was not found: centre taken as the target (rib_score)
        res.update(optical_terms(gp.normalized_gap, gp.f_center if gp.found else f_t,
                                 f_t, a_m, weights))
        res.update(method=method, k_mode=k_mode, candidates=out.get("candidates", []))
        return res
    pen = out.get("freq_penalty")
    if tol_hz is not None and pen is None:
        pen = freq_penalty(gp.f_center, f_t, tol_hz) if gp.found else 0.0
    res = _gap_dict(gp, pen)
    res.update(method=method, k_mode=k_mode, candidates=out.get("candidates", []))
    return res


def x_point_gap(rows, a_um):
    """The module's own figure (diamond_bandStructure_sim.band_edges): the two
    lowest resonances below the light line c/(2a) at X. None without X."""
    at_x = [d for d in rows if abs(d["k"] - 0.5) < 1e-9]
    if not at_x:
        return None
    f = np.sort(at_x[0]["freq"])
    f = f[f < C0 / (2 * a_um * 1e-6)] * 1e-12
    if f.size < 2:
        return dict(found=False)
    mid = 0.5 * (f[0] + f[1])
    return dict(found=True, f_lower_THz=float(f[0]), f_upper_THz=float(f[1]),
                f_center_THz=float(mid), normalized_gap=float((f[1] - f[0]) / mid))


def lite_params_for(folder, a_um=None):
    """(LiteParams, source) from the JSON files of a record, without tidy3d:
    params.json, else the band_params saved in evaluation.json or plan.json,
    else only the period given with --a-um."""
    d = _read_json(os.path.join(folder, "params.json"))
    if d and d.get("a") is not None:
        return LiteParams(d), "params.json"
    for name in ("evaluation.json", "plan.json"):
        j = _read_json(os.path.join(folder, name)) or {}
        detail = j.get("detail", j)
        bp = dict(detail.get("band_params") or {})
        if bp.get("a") is not None:
            win = detail.get("window_THz")
            if win and bp.get("target_freq_thz"):
                bp["freq_span_thz"] = float(win[1]) - float(bp["target_freq_thz"])
            return LiteParams(bp), f"{name} (no params.json)"
    if a_um:
        return LiteParams(dict(a=a_um)), "--a-um (no JSON file of the run)"
    raise AnalysisError(
        f"{folder} has no params.json, evaluation.json or plan.json to read the "
        "period from: give it with --a-um (and the cell type with --celltype)")


def recorded_scoring(folder, evaluation):
    """How the run that made this record scored it: dict(version, tol_thz,
    weights, source), each None when it was not recorded.

    evaluation.json (detail.score_version / optical_tol_thz), then the record
    of the same id in the sqlite file next to the results (read-only):
    score_version, optical_tol_thz, score_weights.
    """
    out = dict(version=None, tol_thz=None, weights=None, source=None)
    detail = (evaluation or {}).get("detail") or {}
    if detail.get("score_version") or detail.get("optical_tol_thz"):
        out.update(version=detail.get("score_version"),
                   tol_thz=detail.get("optical_tol_thz"), source="evaluation.json")
    record_id = os.path.basename(folder)
    results = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(folder))))
    for name in ("runs_rib.sqlite3", "runs.sqlite3"):
        db = os.path.join(results, name)
        if not os.path.isfile(db):
            continue
        try:
            import pathlib
            import sqlite3
            con = sqlite3.connect(pathlib.Path(db).resolve().as_uri() + "?mode=ro",
                                  uri=True)
            try:
                row = con.execute("SELECT record FROM candidates WHERE id = ?",
                                  (record_id,)).fetchone()
            finally:
                con.close()
            if not row:
                continue
            rec = json.loads(row[0])
            for key, field in (("version", "score_version"),
                               ("tol_thz", "optical_tol_thz"),
                               ("weights", "score_weights")):
                if out[key] is None and rec.get(field) is not None:
                    out[key] = rec[field]
                    out["source"] = out["source"] or name
            break
        except Exception:
            continue
    if out["tol_thz"] is not None:
        out["tol_thz"] = float(out["tol_thz"])
    return out


def recorded_tolerance(folder, evaluation):
    """(optical_tol_thz, where from) of the run, or (None, None)."""
    rec = recorded_scoring(folder, evaluation)
    return ((rec["tol_thz"], rec["source"]) if rec["tol_thz"] is not None
            else (None, None))


def rib_scoring_for(folder, evaluation, args):
    """(version, tol_thz, weights, note) to score a rib record with.

    Version: --score-version, else the one recorded with the run, else 2 if
    a Gaussian tolerance was recorded or given (--optical-tol-thz), else the
    current one (4). Version 2 needs the tolerance (recorded, given, or the
    old default 2 THz); versions 3 and 4 the weights (recorded with the run
    when it was scored with that version, else configs/targets.yaml as it is
    now, read for that version).
    """
    rec = recorded_scoring(folder, evaluation)
    version = args.score_version or rec["version"]
    if version is None:
        version = 2 if (args.optical_tol_thz or rec["tol_thz"]) else CURRENT_VERSION
        why = "not recorded with the run"
    else:
        why = "--score-version" if args.score_version else rec["source"]
    version = int(version)
    if version < 3:
        tol, src = args.optical_tol_thz, "--optical-tol-thz"
        if tol is None:
            tol, src = rec["tol_thz"], rec["source"]
        if tol is None:
            tol, src = DEFAULT_TOL_THZ, "default (not recorded with the run)"
        if not tol > 0:
            raise AnalysisError("--optical-tol-thz must be > 0")
        return 2, float(tol), None, "score version 2 (%s), optical_tol_thz = %g (%s)" % (
            why, tol, src)
    weights, src = rec["weights"], rec["source"]
    recorded = rec["version"] is not None and int(rec["version"]) == version
    if not (isinstance(weights, dict) and recorded):
        import yaml
        from objective_rib import rib_score_weights
        with open(os.path.join(_ROOT, "configs", "targets.yaml")) as fh:
            weights = rib_score_weights(yaml.safe_load(fh), version=version)
        src = "configs/targets.yaml"
    note = "score version %d (%s), weights from %s" % (version, why, src)
    if weights.get("optical_freq_tol_THz") is not None:
        note += "; flat-top detuning penalty, f_tol = %g THz, height %g" % (
            weights["optical_freq_tol_THz"], weights["lambda_optical_freq"])
    return version, None, weights, note


def compare_with_record(scored, evaluation):
    """Does the re-analysis agree with evaluation.json of the run?"""
    if not evaluation or "gap" not in evaluation:
        return None
    old = evaluation["gap"]
    same = (bool(old.get("found")) == bool(scored.get("found")))
    d_lo = d_up = None
    if same and old.get("found"):
        d_lo = scored["f_lower_THz"] - old["f_lower"] * 1e-12
        d_up = scored["f_upper_THz"] - old["f_upper"] * 1e-12
        same = abs(d_lo) < 0.05 and abs(d_up) < 0.05
    return dict(matches=bool(same), recorded=dict(
        found=bool(old.get("found")),
        normalized_gap=old.get("normalized_gap"),
        f_lower_THz=(old.get("f_lower") or 0.0) * 1e-12,
        f_upper_THz=(old.get("f_upper") or 0.0) * 1e-12),
        d_lower_THz=d_lo, d_upper_THz=d_up)


def analyze(loader, ot, job, args):
    """Everything about one record: resonances, gaps, comparison.

    The saved analysis (bands.csv) is used when it is there. The .hdf5 files
    are read, and ResonanceFinder is run again, only for explicit .hdf5
    inputs, with --reanalyze, with a resonance-filter / window option (they
    need the time traces), or when there is no bands.csv.
    """
    folder = job.get("folder")
    csv_path = job.get("csv") or os.path.join(folder, "bands.csv")
    has_csv = "files" not in job and os.path.isfile(csv_path)
    has_h5 = "files" in job or ("csv" not in job and _has_hdf5(folder))
    over = {k: getattr(args, k) for k in ("minQ", "minamp", "maxerr")
            if getattr(args, k) is not None}
    refilter = bool(over or args.freq_min_thz is not None
                    or args.freq_max_thz is not None)
    use_h5 = has_h5 and ("files" in job or args.reanalyze or refilter or not has_csv)
    if (args.reanalyze or refilter) and not has_h5:
        print(f"  note: {os.path.basename(folder)} has no .hdf5 files (pruned by "
              "tidy3d.keep_hdf5, or a bands.csv was given): --reanalyze / --minQ "
              "/ --minamp / --maxerr / --freq-*-thz are ignored, the resonances "
              "are those saved by the run", file=sys.stderr)
        over, refilter = {}, False
        args = argparse.Namespace(**dict(vars(args), freq_min_thz=None,
                                         freq_max_thz=None))
    evaluation = _read_json(os.path.join(folder, "evaluation.json")) if folder else None

    if use_h5:
        bs = loader.get("read .hdf5 files" if not (args.reanalyze or refilter)
                        or "files" in job else
                        "re-run ResonanceFinder (--reanalyze / filter options)")
        sim_datas = (load_files(bs, job["files"]) if "files" in job
                     else bs.load_run(folder))
        sim = next(iter(sim_datas.values())).simulation
        p, source = params_for(bs, folder, sim)
        a_um = float(sim.size[0])
        symmetry = [int(s) for s in sim.symmetry]
        if abs(a_um - p.a) > 1e-9:
            print(f"  note: the files have a = {a_um:.4f} um, params.json says "
                  f"{p.a:.4f}; using the files", file=sys.stderr)
            p = dataclasses.replace(p, a=a_um)
        p = dataclasses.replace(p, **over)
    else:
        p, source = lite_params_for(folder, args.a_um)
        a_um = float(p.a)
        symmetry = (((evaluation or {}).get("detail") or {}).get("symmetry")
                    or [int(s) for s in p.symmetry])
    if args.freq_min_thz is not None or args.freq_max_thz is not None:
        lo = args.freq_min_thz if args.freq_min_thz is not None else p.freq_range[0] * 1e-12
        hi = args.freq_max_thz if args.freq_max_thz is not None else p.freq_range[1] * 1e-12
        if not 0 < lo < hi:
            raise AnalysisError("need 0 < --freq-min-thz < --freq-max-thz")
        # the window of the analysis only; the scoring target stays the target
        window = (lo * 1e12, hi * 1e12)
    else:
        window = tuple(p.freq_range)
    p_win = _with_window(p, window)
    celltype = args.celltype or celltype_of(p, source, folder, evaluation)
    version, tol, weights, scoring_note = CURRENT_VERSION, None, None, None
    if celltype == "rib":
        version, tol, weights, scoring_note = rib_scoring_for(folder, evaluation, args)

    rows = (resonances_of(bs, sim_datas, p_win) if use_h5
            else rows_from_bands_csv(csv_path))
    a_m = a_um * 1e-6
    for d in rows:
        ll = d["k"] * C0 / a_m
        d["guided"] = d["freq"] < ll
        d["in_window"] = (d["freq"] >= window[0]) & (d["freq"] <= window[1])
    skw = dict(version=version, weights=weights)
    scored = score(ot, celltype, rows, a_m, p_win, tol, **skw)
    cleaned = bool(args.min_guided_q or args.light_line_margin)
    filtered = (score(ot, celltype, rows, a_m, p_win, tol, args.min_guided_q,
                      args.light_line_margin, **skw) if cleaned else None)
    defaults = not over and window == tuple(p.freq_range)
    record_id = os.path.basename(folder) if folder else "files"
    n_guided = int(sum(d["guided"].sum() for d in rows))
    low_q = int(sum((d["guided"] & (np.abs(d["Q"]) < 1000)).sum() for d in rows))
    near_ll = int(sum((d["guided"] & (d["freq"] > 0.97 * d["k"] * C0 / a_m)).sum()
                      for d in rows))
    if use_h5:
        data_source = "hdf5"
    else:
        has_amp = any(np.isfinite(d["amplitude"]).any() for d in rows)
        data_source = ("saved analysis (bands.csv%s%s)"
                       % ("" if has_amp else ", no amplitude / fit error",
                          "; .hdf5 present, use --reanalyze to read them"
                          if has_h5 else "; no .hdf5"))
    return dict(
        record_id=record_id, folder=folder, celltype=celltype, params=p,
        params_source=source, a_um=a_um, window_thz=[w * 1e-12 for w in window],
        rows=rows, scored=scored, scored_min_guided_q=filtered,
        x_point=x_point_gap(rows, a_um), tol_thz=tol, tol_source=scoring_note,
        score_version=(version if celltype == "rib" else None),
        score_weights=weights, scoring_note=scoring_note,
        symmetry=[int(s) for s in symmetry], data_source=data_source,
        from_hdf5=bool(use_h5),
        counts=dict(k_points=len(rows),
                    resonances=int(sum(len(d["freq"]) for d in rows)),
                    guided=n_guided, guided_low_q=low_q,
                    guided_near_light_line=near_ll,
                    negative_q=int(sum((d["Q"] < 0).sum() for d in rows)),
                    out_of_window=int(sum((~d["in_window"]).sum() for d in rows))),
        filters=dict(minQ=p.minQ, minamp=p.minamp, maxerr=p.maxerr,
                     min_guided_q=args.min_guided_q,
                     light_line_margin=args.light_line_margin),
        comparison=(compare_with_record(scored, evaluation) if defaults else None),
        evaluation=evaluation)


def simple_outline(p):
    """Top-view outline (um) of the unit cell from the saved dimensions, for
    the geometry panel when tidy3d is not there: sharp corners (the fillets
    are not drawn). None if the dimensions are not known."""
    a = p.a
    if p.geometry == "holey" and np.isfinite([p.beam_width, p.hole_x, p.hole_y]).all():
        bw = p.beam_width
        ang = np.linspace(0, 2 * np.pi, 96, endpoint=False)
        return [np.array([(-a / 2, -bw / 2), (a / 2, -bw / 2), (a / 2, bw / 2),
                          (-a / 2, bw / 2)]),
                np.stack([p.hole_x / 2 * np.cos(ang), p.hole_y / 2 * np.sin(ang)], axis=1)]
    if (p.geometry == "rib_block"
            and np.isfinite([p.t, p.w, p.spine_width]).all()):
        t, w, s = p.t, p.w, p.spine_width
        upper = [(-a / 2, s / 2), (-t / 2, s / 2), (-t / 2, w / 2), (t / 2, w / 2),
                 (t / 2, s / 2), (a / 2, s / 2)]
        return [np.array(upper + [(x, -y) for x, y in reversed(upper)])]
    return None


def _with_window(p, window_hz):
    """BandParams whose freq_range is window_hz (target unchanged if centred)."""
    lo, hi = window_hz
    if abs(lo - p.freq_range[0]) < 1 and abs(hi - p.freq_range[1]) < 1:
        return p

    class _Windowed:                     # only freq_range differs
        def __init__(self, base):
            self.__dict__["_base"] = base

        def __getattr__(self, name):
            return getattr(self._base, name)

        freq_range = (lo, hi)
    return _Windowed(p)


# ---- output -----------------------------------------------------------------

def geometry_text(res):
    p, cell = res["params"], res["celltype"]
    nm = lambda x: "%.0f" % (x * 1e3)                              # noqa: E731
    if (res["params_source"].startswith("simulation") or not p.geometry
            or not np.isfinite(getattr(p, "thickness", float("nan")))):
        return "a=%s nm (geometry not stored)" % nm(res["a_um"])
    if p.geometry == "holey":
        return ("a=%s  w=%s  t=%s  hx=%s  hy=%s nm (semi-axes)"
                % (nm(p.a), nm(p.beam_width), nm(p.thickness),
                   nm(p.hole_x / 2), nm(p.hole_y / 2)))
    return ("a=%s  w=%s  t=%s  s=%s  th=%s nm"
            % (nm(p.a), nm(p.w), nm(p.t), nm(p.spine_width), nm(p.thickness))
            + ("" if cell == "rib" else "  (%s)" % p.geometry))


def write_csv(path, res):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["k_2pi_over_a", "k_a_over_pi", "freq_THz", "Q", "amplitude",
                    "error", "guided", "in_window"])
        for d in res["rows"]:
            for i in range(len(d["freq"])):
                w.writerow([d["k"], 2 * d["k"], d["freq"][i] * 1e-12, d["Q"][i],
                            d["amplitude"][i], d["error"][i],
                            int(d["guided"][i]), int(d["in_window"][i])])


def summary_of(res):
    p = res["params"]
    return dict(
        record_id=res["record_id"], celltype=res["celltype"],
        xsect=getattr(p, "xsect", None), geometry=getattr(p, "geometry", None),
        params_source=res["params_source"], data_source=res["data_source"],
        a_um=res["a_um"],
        index=float(p.n), mode=p.mode, symmetry=res["symmetry"],
        target_THz=float(p.target_freq_thz), window_THz=res["window_thz"],
        light_line_at_X_THz=C0 / (2 * res["a_um"] * 1e-6) * 1e-12,
        score_version=res["score_version"], scoring=res["scoring_note"],
        score_weights=res["score_weights"],
        optical_tol_thz=res["tol_thz"],                  # score version 2 only
        optical_tol_source=(res["tol_source"] if res["tol_thz"] is not None else None),
        filters=res["filters"], counts=res["counts"],
        k_2pi_over_a=[float(d["k"]) for d in res["rows"]],
        scored_gap=res["scored"], scored_gap_min_guided_q=res["scored_min_guided_q"],
        x_point_gap=res["x_point"], comparison_with_evaluation_json=res["comparison"])


def plot(loader, res, path, k_units="2pi/a", show=False, geometry=True, ymax=None):
    """Band diagram (+ unit-cell outline) as a PNG."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    p, rows, scored = res["params"], res["rows"], res["scored"]
    a_m = res["a_um"] * 1e-6
    scale = 1.0 if k_units == "2pi/a" else 2.0           # k*a/pi = 2 * k
    k_x = 0.5 * scale
    outline, exact = None, False
    if geometry and res["folder"] and loader.available and os.path.isfile(
            os.path.join(res["folder"], "params.json")):
        try:                            # the exact outline, with fillets
            bs = loader.get("draw the unit cell")
            outline = bs.build_unit_cell(bs._load_saved_params(res["folder"])).outline
            exact = True
        except Exception:
            outline = None
    if (geometry and outline is None
            and not res["params_source"].startswith("simulation")):
        outline = simple_outline(p)
    fig, axes = plt.subplots(
        1, 2 if outline else 1, figsize=(10.5 if outline else 7.2, 6.0),
        gridspec_kw=dict(width_ratios=[3.2, 1]) if outline else None,
        facecolor=SURFACE)
    ax = axes[0] if outline else axes
    ax.set_facecolor(SURFACE)

    f_all = np.concatenate([d["freq"] for d in rows]) * 1e-12 if rows else np.array([])
    top = ymax or max(res["window_thz"][1],
                      float(f_all.max()) * 1.03 if f_all.size else 0.0)
    kk = np.linspace(0, 0.5, 100)
    light = kk * C0 / a_m * 1e-12
    ax.fill_between(kk * scale, np.minimum(light, top), top, color=CONE, lw=0,
                    label="light cone", zorder=0)
    ax.plot(kk * scale, light, color=MUTED, lw=1, zorder=1)

    # target and tolerance
    f_t = float(p.target_freq_thz)
    if res["tol_thz"] is not None:                       # score version 2
        ax.axhspan(f_t - res["tol_thz"], f_t + res["tol_thz"], color=TARGET,
                   alpha=0.16, lw=0, zorder=1,
                   label="target +/- %.3g THz (weight 1/e)" % res["tol_thz"])
    f_tol = (res.get("score_weights") or {}).get("optical_freq_tol_THz")
    if res["tol_thz"] is None and f_tol is not None:     # score version 4
        ax.axhspan(f_t - f_tol, f_t + f_tol, color=TARGET, alpha=0.16, lw=0,
                   zorder=1, label="target +/- %.3g THz (detuning penalty "
                                   "at half height)" % f_tol)
    ax.axhline(f_t, color=TARGET, ls=(0, (4, 3)), lw=1.4, zorder=2,
               label="target %.1f THz" % f_t)

    # gaps
    if scored.get("found"):
        ax.axhspan(scored["f_lower_THz"], scored["f_upper_THz"], color=GAP,
                   alpha=0.30, lw=0, zorder=1, label="scored gap %.1f %%"
                   % (100 * scored["normalized_gap"]))
    xp = res["x_point"]
    if xp and xp.get("found"):
        x0 = k_x * 0.955
        ax.plot([x0, x0], [xp["f_lower_THz"], xp["f_upper_THz"]], color=XGAP,
                lw=5, solid_capstyle="butt", zorder=3,
                label="gap at X (two lowest guided) %.1f %%"
                % (100 * xp["normalized_gap"]))

    # resonances: guided filled, coloured by log10|Q|; leaky hollow
    cmap = LinearSegmentedColormap.from_list("q", Q_RAMP)
    norm = Normalize(1.0, 6.0)
    gk, gf, gq, lk, lf = [], [], [], [], []
    for d in rows:
        for f, q, g in zip(d["freq"], d["Q"], d["guided"]):
            if g:
                gk.append(d["k"] * scale)
                gf.append(f * 1e-12)
                gq.append(np.log10(max(abs(q), 1.0)))
            else:
                lk.append(d["k"] * scale)
                lf.append(f * 1e-12)
    ax.scatter(lk, lf, s=34, facecolors="none", edgecolors=MUTED, linewidths=1.2,
               zorder=4, clip_on=False, label="resonance above the light line")
    sc = ax.scatter(gk, gf, s=46, c=gq, cmap=cmap, norm=norm, edgecolors=SURFACE,
                    linewidths=1.0, zorder=5, clip_on=False,
                    label="guided resonance (colour: Q)")
    if gk:
        cb = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.045)
        cb.set_label("log10 |Q| of guided resonances", color=INK2, fontsize=9)
        cb.ax.tick_params(labelsize=8, colors=INK2)
        cb.outline.set_visible(False)

    k_min = min([d["k"] for d in rows] + [0.5]) * scale
    left = 0.0 if len(rows) < 2 else max(0.0, k_min - 0.04 * scale)
    ax.set_xlim(left, k_x)
    ax.set_ylim(0, top)
    if k_units == "2pi/a":
        ax.set_xlabel(r"wavevector $k_x$  [$2\pi/a$]   (X at 0.5)", color=INK2)
    else:
        ax.set_xlabel(r"wavevector $k_x a/\pi$   (X at 1)", color=INK2)
    ax.set_ylabel("frequency [THz]", color=INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.grid(axis="y", color=CONE, lw=0.6, zorder=0)

    cell = res["celltype"] or "cell"
    head = ("%s cell, %s, %s, n = %.3f   |   %s"
            % (cell, getattr(p, "xsect", "?"), p.mode, p.n, res["record_id"]))
    if scored.get("found"):
        line = ("scored gap %.2f-%.2f THz, centre %.2f, %.1f %%"
                % (scored["f_lower_THz"], scored["f_upper_THz"],
                   scored["f_center_THz"], 100 * scored["normalized_gap"]))
        if "freq_penalty" in scored:
            line += ("; weight %.3g -> effective %.2f %%"
                     % (scored["freq_penalty"], 100 * scored["effective_gap"]))
        elif "pen_detune_o" in scored:
            line += ("\npenalties: size %.3g, detuning %.3g, light cone %.3g"
                     " -> optical term %.3f"
                     % (scored["pen_size_o"], scored["pen_detune_o"],
                        scored["pen_lightcone_o"], scored["value"]))
    else:
        line = "no scored gap" + (" (%s)" % scored["not_scored"]
                                  if scored.get("not_scored") else "")
    ax.set_title(head + "\n" + geometry_text(res) + "\n" + line, fontsize=9.5,
                 color=INK, loc="left")
    leg = ax.legend(loc="lower left" if len(rows) > 1 else "upper left",
                    fontsize=8, frameon=True, framealpha=0.92, edgecolor=CONE)
    for text in leg.get_texts():
        text.set_color(INK)

    if outline:
        ax2 = axes[1]
        ax2.set_facecolor(SURFACE)
        for i, poly in enumerate(outline):
            closed = np.vstack([poly, poly[:1]]) * 1e3
            ax2.fill(closed[:, 0], closed[:, 1],
                     color=Q_RAMP[1] if i == 0 else SURFACE,
                     alpha=0.55 if i == 0 else 1.0, lw=0)
            ax2.plot(closed[:, 0], closed[:, 1], color=Q_RAMP[3], lw=1)
        ax2.set_aspect("equal")
        ax2.set_title("unit cell (top view)" if exact else
                      "unit cell (top view,\nfillets not drawn)",
                      fontsize=9, color=INK, loc="left")
        ax2.set_xlabel("x [nm] (beam)", color=INK2, fontsize=8)
        ax2.set_ylabel("y [nm]", color=INK2, fontsize=8)
        ax2.tick_params(colors=INK2, labelsize=7)
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax2.spines[side].set_color(MUTED)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    if show:
        plt.show()
    plt.close(fig)


def report(res, outdir):
    c, s = res["counts"], res["scored"]
    print("%s  (%s, %s)" % (res["record_id"], res["celltype"] or "celltype unknown",
                            res["params_source"]))
    print("  data: %s" % res["data_source"])
    if res["scoring_note"]:
        print("  %s" % res["scoring_note"])
    print("  %d k point(s), %d resonances: %d guided (%d of them |Q| < 1000, "
          "%d within 3 %% of the light line), %d with Q < 0, %d outside the window"
          % (c["k_points"], c["resonances"], c["guided"], c["guided_low_q"],
             c["guided_near_light_line"], c["negative_q"], c["out_of_window"]))
    if s.get("found"):
        line = ("  scored gap (%s): %.2f-%.2f THz, centre %.2f THz, %.2f %%"
                % (s["method"], s["f_lower_THz"], s["f_upper_THz"],
                   s["f_center_THz"], 100 * s["normalized_gap"]))
        if "freq_penalty" in s:
            line += (", weight %.3g, effective %.3f %%"
                     % (s["freq_penalty"], 100 * s["effective_gap"]))
        print(line)
        if "pen_detune_o" in s:
            print("    optical score term %.4f = gap %.4f - size %.4g - detuning %.4g "
                  "- light cone %.4g   (%d candidate gap(s))"
                  % (s["value"], s["normalized_gap"], s["pen_size_o"],
                     s["pen_detune_o"], s["pen_lightcone_o"], len(s["candidates"])))
    else:
        print("  scored gap: none" + (" (%s)" % s["not_scored"]
                                      if s.get("not_scored") else ""))
    f = res["scored_min_guided_q"]
    if f is not None:
        knobs = ", ".join(x for x in (
            "guided |Q| >= %g" % res["filters"]["min_guided_q"]
            if res["filters"]["min_guided_q"] else "",
            "none within %g %% of the light line"
            % (100 * res["filters"]["light_line_margin"])
            if res["filters"]["light_line_margin"] else "") if x)
        print("  cleaned (%s): " % knobs
              + ("%.2f-%.2f THz, centre %.2f THz, %.2f %%, %s"
                 % (f["f_lower_THz"], f["f_upper_THz"], f["f_center_THz"],
                    100 * f["normalized_gap"],
                    "weight %.3g" % f["freq_penalty"] if "freq_penalty" in f
                    else "optical term %.4f" % f["value"] if "value" in f else "")
                 if f.get("found") else "none"))
    xp = res["x_point"]
    if xp is None:
        print("  gap at X: no file at the X point")
    elif xp.get("found"):
        print("  gap at X (two lowest guided): %.2f-%.2f THz, %.2f %%"
              % (xp["f_lower_THz"], xp["f_upper_THz"], 100 * xp["normalized_gap"]))
    else:
        print("  gap at X: fewer than two guided resonances")
    cmp_ = res["comparison"]
    if cmp_ is not None:
        r = cmp_["recorded"]
        print("  evaluation.json of the run: %s (recorded %.2f-%.2f THz, %s)"
              % ("MATCHES" if cmp_["matches"] else "DIFFERS",
                 r["f_lower_THz"], r["f_upper_THz"],
                 "%.2f %%" % (100 * r["normalized_gap"]) if r["found"] else "none"))
    print("  -> %s" % outdir)


def build_parser():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*",
                    help="a bands.csv, .hdf5 file(s), a run folder, a record id, "
                         "or <celltype>/<record id> under --results-dir")
    ap.add_argument("--reanalyze", action="store_true",
                    help="read the .hdf5 files and run ResonanceFinder again even "
                         "when a bands.csv is saved (needs tidy3d and the files)")
    ap.add_argument("--a-um", type=float, default=None,
                    help="period in um, only for a bands.csv without any JSON "
                         "file of the run next to it")
    ap.add_argument("--all", action="store_true",
                    help="every record folder under --results-dir, plus an index CSV")
    ap.add_argument("--results-dir", default=DEFAULT_RESULTS,
                    help="where the records are [python-scripts/results/tidy3d]")
    ap.add_argument("--outdir", default=None,
                    help="output folder [<record folder>/analysis]; with several "
                         "records each gets a sub-folder named after it")
    g = ap.add_argument_group("resonance filters (default: those of the run; "
                              "the first five need the .hdf5 files and tidy3d)")
    g.add_argument("--minQ", type=float, default=None, help="drop |Q| below this")
    g.add_argument("--minamp", type=float, default=None, help="drop amplitude below this")
    g.add_argument("--maxerr", type=float, default=None, help="drop fit error above this")
    g.add_argument("--freq-min-thz", type=float, default=None,
                   help="lower end of the harmonic-inversion window")
    g.add_argument("--freq-max-thz", type=float, default=None,
                   help="upper end of the harmonic-inversion window")
    g.add_argument("--min-guided-q", type=float, default=None,
                   help="ALSO score without guided resonances of |Q| below this "
                        "(a real guided mode has a very high Q; a low-Q peak "
                        "below the light line is usually spurious). Reported "
                        "next to the pipeline's gap, which is unchanged")
    g.add_argument("--light-line-margin", type=float, default=None, metavar="FRAC",
                   help="ALSO score without guided resonances within this "
                        "fraction below the light line (e.g. 0.03), where FDTD "
                        "shows peaks that are not bands. Combines with "
                        "--min-guided-q")
    s = ap.add_argument_group("scoring")
    s.add_argument("--score-version", type=int, choices=(2, 3, 4), default=None,
                   help="rib: score the gap as version 4 (gap minus size, "
                        "flat-top detuning and light-cone penalties; current), "
                        "version 3 (the same, quadratic detuning penalty) or "
                        "version 2 (Gaussian weight) [the version recorded "
                        "with the run, else 4]")
    s.add_argument("--optical-tol-thz", type=float, default=None,
                   help="rib, score version 2: Gaussian width of the gap weight "
                        "[the value recorded with the run, else %g]" % DEFAULT_TOL_THZ)
    s.add_argument("--celltype", choices=("hole", "rib"), default=None,
                   help="override the cell type read from the run")
    f = ap.add_argument_group("figure")
    f.add_argument("--k-units", choices=("2pi/a", "ka/pi"), default="2pi/a",
                   help="k axis: in 2*pi/a (X at 0.5, the module) or k*a/pi "
                        "(X at 1, the repo) [2pi/a]")
    f.add_argument("--ymax", type=float, default=None, help="top of the frequency axis [THz]")
    f.add_argument("--no-geometry", action="store_true", help="no unit-cell panel")
    f.add_argument("--show", action="store_true", help="also open the figure window")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        with network_blocked():
            jobs = resolve_inputs(args.inputs, args.results_dir, args.all)
            ot = load_scoring()
            loader = ModuleLoader(ot)
            index = []
            for job in jobs:
                res = analyze(loader, ot, job, args)
                if args.outdir:
                    outdir = (args.outdir if len(jobs) == 1
                              else os.path.join(args.outdir, res["record_id"]))
                else:
                    outdir = os.path.join(job["folder"], "analysis")
                os.makedirs(outdir, exist_ok=True)
                write_csv(os.path.join(outdir, "resonances.csv"), res)
                summary = summary_of(res)
                with open(os.path.join(outdir, "summary.json"), "w") as fh:
                    json.dump(summary, fh, indent=2)
                plot(loader, res, os.path.join(outdir, "band_diagram.png"),
                     k_units=args.k_units, show=args.show,
                     geometry=not args.no_geometry, ymax=args.ymax)
                report(res, outdir)
                index.append((summary, outdir))
            if args.all or len(jobs) > 1:
                path = os.path.join(args.outdir or args.results_dir,
                                    "analysis_index.csv")
                write_index(path, index)
                print("index: %s" % path)
    except AnalysisError as err:
        print("plot_tidy3d_bands: %s" % err, file=sys.stderr)
        return 2
    return 0


def write_index(path, index):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = ["record_id", "celltype", "xsect", "a_um", "k_points", "resonances",
            "guided", "guided_low_q", "gap_found", "f_lower_THz", "f_upper_THz",
            "f_center_THz", "normalized_gap", "freq_penalty", "effective_gap",
            "x_gap_lower_THz", "x_gap_upper_THz", "x_gap", "matches_evaluation",
            "outdir", "data_source", "optical_tol_thz", "score_version",
            "optical_term", "pen_size_o", "pen_detune_o", "pen_lightcone_o"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for s, outdir in index:
            g, x = s["scored_gap"], s["x_point_gap"] or {}
            cmp_ = s["comparison_with_evaluation_json"]
            w.writerow([s["record_id"], s["celltype"], s["xsect"], s["a_um"],
                        s["counts"]["k_points"], s["counts"]["resonances"],
                        s["counts"]["guided"], s["counts"]["guided_low_q"],
                        int(bool(g.get("found"))), g.get("f_lower_THz"),
                        g.get("f_upper_THz"), g.get("f_center_THz"),
                        g.get("normalized_gap"), g.get("freq_penalty"),
                        g.get("effective_gap"), x.get("f_lower_THz"),
                        x.get("f_upper_THz"), x.get("normalized_gap"),
                        "" if cmp_ is None else int(cmp_["matches"]), outdir,
                        s.get("data_source"), s.get("optical_tol_thz"),
                        s.get("score_version"), g.get("value"), g.get("pen_size_o"),
                        g.get("pen_detune_o"), g.get("pen_lightcone_o")])


if __name__ == "__main__":
    sys.exit(main())
