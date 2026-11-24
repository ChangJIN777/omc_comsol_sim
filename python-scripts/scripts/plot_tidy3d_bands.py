#!/usr/bin/env python3
"""Band diagrams from downloaded Tidy3D results (.hdf5), offline.

Re-analyzes the SimulationData files that `optical: tidy3d` leaves in
results/tidy3d/<celltype>/<record id>/ (one file per k point) and plots the
band structure with the gap exactly as the pipeline scores it.

LOCAL FILES ONLY: nothing is uploaded or downloaded, the Tidy3D web API is
never called, and network connections are blocked while the script runs. The
result files and the databases are not modified; output goes to a new
`analysis/` folder (or --outdir).

It reuses python-scripts/diamond_bandStructure_sim.py (load_run, bloch_k,
ResonanceFinder + filter_resonances, band_edges, BandParams with
params.json, build_unit_cell for the outline) and the repo's own scoring
(src/optical_tidy3d.score_resonances -> comsol_optics.guided_interval_gaps /
select_gap_by_penalized_size / optical_gap_tracked).

Needs tidy3d: run it with the Python of the `tidy` conda environment.

Examples (from python-scripts/):
  python scripts/plot_tidy3d_bands.py cand_f85eb9955630
  python scripts/plot_tidy3d_bands.py rib/cand_f85eb9955630 --k-units ka/pi
  python scripts/plot_tidy3d_bands.py results/tidy3d/hole/cand_808e4ba6dee3
  python scripts/plot_tidy3d_bands.py path/a.hdf5 path/b.hdf5 --outdir out
  python scripts/plot_tidy3d_bands.py --all
  python scripts/plot_tidy3d_bands.py cand_f85eb9955630 --min-guided-q 1000
  python scripts/plot_tidy3d_bands.py cand_f85eb9955630 --minQ 50 --optical-tol-thz 10

A record whose .hdf5 files were deleted (tidy3d.keep_hdf5: best | none) is
plotted from the resonances the run saved (bands.csv + params.json); the
resonance filters cannot be changed for it.

Inputs: .hdf5 files, a run folder, a record id, or <celltype>/<record id>
under --results-dir. Files are grouped by the Bloch wavevector stored in
them (the newest file wins for a repeated k); file names do not matter and
non-result files such as batch.hdf5 are skipped.

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


def load_module():
    """diamond_bandStructure_sim (imports tidy3d) and the repo scoring."""
    sys.path.insert(0, os.path.join(_ROOT, "src"))
    try:
        import optical_tidy3d as ot
        if not ot.have_tidy3d():
            raise ImportError("No module named 'tidy3d'")
        bs = ot.load_bs()
    except ImportError as err:
        raise AnalysisError(
            f"tidy3d is needed to read the result files ({err}). Run this "
            "script with the Python of the `tidy` conda environment, e.g.\n"
            "  conda run -n tidy python scripts/plot_tidy3d_bands.py ...")
    except Exception as err:
        if getattr(err, "tidy3d_stop", None) == "not_installed":
            raise AnalysisError(
                "tidy3d is needed to read the result files. Run this script "
                "with the Python of the `tidy` conda environment, e.g.\n"
                "  conda run -n tidy python scripts/plot_tidy3d_bands.py ...")
        raise
    return bs, ot


# ---- inputs -----------------------------------------------------------------

def _has_hdf5(folder):
    return os.path.isdir(folder) and any(
        n.endswith(".hdf5") and n != "batch.hdf5" for n in os.listdir(folder))


def _has_cache(folder):
    """A record whose .hdf5 files were pruned (tidy3d.keep_hdf5 best / none)
    but whose analysis is saved: bands.csv + params.json."""
    return all(os.path.isfile(os.path.join(folder, n))
               for n in ("bands.csv", "params.json"))


def _has_results(folder):
    return os.path.isdir(folder) and (_has_hdf5(folder) or _has_cache(folder))


def rows_from_bands_csv(path):
    """Resonances saved by the pipeline (bands.csv: k, freq, Q), per k.
    Amplitude and fit error are not in that file: NaN."""
    by_k = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            by_k.setdefault(float(row["k_2pi_over_a"]), []).append(
                (float(row["freq_THz"]) * 1e12, float(row["Q"])))
    rows = []
    for k in sorted(by_k):
        pairs = sorted(by_k[k])
        f = np.array([x[0] for x in pairs])
        rows.append(dict(k=k, freq=f, Q=np.array([x[1] for x in pairs]),
                         amplitude=np.full(f.size, np.nan),
                         error=np.full(f.size, np.nan)))
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
    """-> list of jobs dict(folder=... | files=[...], label=...).

    An input is a .hdf5 file, a folder, a record id, or <celltype>/<record id>
    under results_dir. All files given directly form ONE job.
    """
    jobs, files = [], []
    if all_records:
        jobs += [dict(folder=f) for f in record_folders(results_dir)]
    for item in inputs:
        if os.path.isfile(item):
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
        key = (job["folder"], tuple(job.get("files", ())))
        if key not in seen:
            seen.add(key)
            unique.append(job)
    for job in unique:
        if "files" not in job and not _has_results(job["folder"]):
            raise AnalysisError("no Tidy3D result files (*.hdf5) and no saved "
                                f"analysis (bands.csv + params.json) in {job['folder']}")
    if not unique:
        raise AnalysisError("nothing to analyze: give .hdf5 files, a run folder, "
                            "a record id, or --all (see --help)")
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
    if source == "params.json":
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
          light_line_margin=None):
    """The pipeline's gap for these resonances (optical_tidy3d.score_resonances).

    rib: Gaussian-weighted choice among all guided gaps (sweep) or the gap at
    X (x_only), weight exp(-((f_c - f_t)/tol)^2). hole / unknown: the tracked
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
    tol_hz = tol_thz * 1e12 if celltype == "rib" else None
    try:
        gp, method = ot.score_resonances(
            celltype or "hole", k_frac, freqs, a_m, f_target=f_t, k_mode=k_mode,
            window_hz=tuple(p.freq_range), optical_tol_hz=tol_hz, out=out)
    except ValueError as err:
        return dict(found=False, not_scored=str(err), method=None, candidates=[])
    pen = out.get("freq_penalty")
    if tol_hz is not None and pen is None:
        pen = freq_penalty(gp.f_center, f_t, tol_hz) if gp.found else 0.0
    res = _gap_dict(gp, pen)
    res.update(method=method, k_mode=k_mode, candidates=out.get("candidates", []))
    return res


def x_point_gap(bs, rows, p):
    """The module's own figure: two lowest guided resonances at X, or None."""
    at_x = [d for d in rows if abs(d["k"] - 0.5) < 1e-9]
    if not at_x:
        return None
    f_lo, f_up, gap, mid, fit = bs.band_edges(at_x[0]["freq"], p)
    if not np.isfinite(mid):
        return dict(found=False)
    return dict(found=True, f_lower_THz=float(f_lo), f_upper_THz=float(f_up),
                f_center_THz=float(mid), normalized_gap=float(gap / mid))


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


def analyze(bs, ot, job, args):
    """Everything about one record: resonances, gaps, comparison."""
    folder = job.get("folder")
    # A record whose .hdf5 files were pruned is plotted from its saved
    # resonances; ResonanceFinder cannot be re-run for it.
    from_cache = "files" not in job and not _has_hdf5(folder)
    if from_cache:
        p, source = bs._load_saved_params(folder), "params.json"
        sim_datas, a_um = None, float(p.a)
        symmetry = [int(s) for s in p.symmetry]
    else:
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
    over = {k: getattr(args, k) for k in ("minQ", "minamp", "maxerr")
            if getattr(args, k) is not None}
    if from_cache and (over or args.freq_min_thz is not None
                       or args.freq_max_thz is not None):
        print(f"  note: {os.path.basename(folder)} has no .hdf5 files (pruned by "
              "tidy3d.keep_hdf5): --minQ / --minamp / --maxerr / --freq-*-thz "
              "are ignored, the resonances are those saved by the run",
              file=sys.stderr)
        over = {}
        args = argparse.Namespace(**dict(vars(args), freq_min_thz=None,
                                         freq_max_thz=None))
    if args.freq_min_thz is not None or args.freq_max_thz is not None:
        lo = args.freq_min_thz if args.freq_min_thz is not None else p.freq_range[0] * 1e-12
        hi = args.freq_max_thz if args.freq_max_thz is not None else p.freq_range[1] * 1e-12
        if not 0 < lo < hi:
            raise AnalysisError("need 0 < --freq-min-thz < --freq-max-thz")
        # the window of the analysis only; the scoring target stays the target
        window = (lo * 1e12, hi * 1e12)
    else:
        window = tuple(p.freq_range)
    p = dataclasses.replace(p, **over)
    p_win = _with_window(p, window)
    evaluation = _read_json(os.path.join(folder, "evaluation.json")) if folder else None
    celltype = args.celltype or celltype_of(p, source, folder, evaluation)
    tol = args.optical_tol_thz
    if tol is None:
        tol = DEFAULT_TOL_THZ
    if not tol > 0:
        raise AnalysisError("--optical-tol-thz must be > 0")

    rows = (rows_from_bands_csv(os.path.join(folder, "bands.csv")) if from_cache
            else resonances_of(bs, sim_datas, p_win))
    a_m = a_um * 1e-6
    for d in rows:
        ll = d["k"] * C0 / a_m
        d["guided"] = d["freq"] < ll
        d["in_window"] = (d["freq"] >= window[0]) & (d["freq"] <= window[1])
    scored = score(ot, celltype, rows, a_m, p_win, tol)
    cleaned = bool(args.min_guided_q or args.light_line_margin)
    filtered = (score(ot, celltype, rows, a_m, p_win, tol, args.min_guided_q,
                      args.light_line_margin) if cleaned else None)
    defaults = not over and window == tuple(p.freq_range)
    record_id = os.path.basename(folder) if folder else "files"
    n_guided = int(sum(d["guided"].sum() for d in rows))
    low_q = int(sum((d["guided"] & (np.abs(d["Q"]) < 1000)).sum() for d in rows))
    near_ll = int(sum((d["guided"] & (d["freq"] > 0.97 * d["k"] * C0 / a_m)).sum()
                      for d in rows))
    return dict(
        record_id=record_id, folder=folder, celltype=celltype, params=p,
        params_source=source, a_um=a_um, window_thz=[w * 1e-12 for w in window],
        rows=rows, scored=scored, scored_min_guided_q=filtered,
        x_point=x_point_gap(bs, rows, p_win), tol_thz=tol,
        symmetry=symmetry,
        data_source=("saved analysis (bands.csv; .hdf5 pruned)" if from_cache
                     else "hdf5"),
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
    if res["params_source"] != "params.json":
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
        optical_tol_thz=(res["tol_thz"] if res["celltype"] == "rib" else None),
        filters=res["filters"], counts=res["counts"],
        k_2pi_over_a=[float(d["k"]) for d in res["rows"]],
        scored_gap=res["scored"], scored_gap_min_guided_q=res["scored_min_guided_q"],
        x_point_gap=res["x_point"], comparison_with_evaluation_json=res["comparison"])


def plot(bs, res, path, k_units="2pi/a", show=False, geometry=True, ymax=None):
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
    outline = None
    if geometry and res["params_source"] == "params.json":
        try:
            outline = bs.build_unit_cell(p).outline
        except Exception:
            outline = None
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
    if res["celltype"] == "rib":
        ax.axhspan(f_t - res["tol_thz"], f_t + res["tol_thz"], color=TARGET,
                   alpha=0.16, lw=0, zorder=1,
                   label="target +/- %.3g THz (weight 1/e)" % res["tol_thz"])
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
        ax2.set_title("unit cell (top view)", fontsize=9, color=INK, loc="left")
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
    if res["data_source"] != "hdf5":
        print("  from the %s: ResonanceFinder cannot be re-run, amplitude and "
              "fit error are not available" % res["data_source"])
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
              + ("%.2f-%.2f THz, centre %.2f THz, %.2f %%, weight %.3g"
                 % (f["f_lower_THz"], f["f_upper_THz"], f["f_center_THz"],
                    100 * f["normalized_gap"], f.get("freq_penalty", float("nan")))
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
                    help=".hdf5 file(s), a run folder, a record id, or "
                         "<celltype>/<record id> under --results-dir")
    ap.add_argument("--all", action="store_true",
                    help="every record folder under --results-dir, plus an index CSV")
    ap.add_argument("--results-dir", default=DEFAULT_RESULTS,
                    help="where the records are [python-scripts/results/tidy3d]")
    ap.add_argument("--outdir", default=None,
                    help="output folder [<record folder>/analysis]; with several "
                         "records each gets a sub-folder named after it")
    g = ap.add_argument_group("resonance filters (default: those of the run)")
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
    s.add_argument("--optical-tol-thz", type=float, default=None,
                   help="rib: Gaussian width of the gap weight [%g]" % DEFAULT_TOL_THZ)
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
            bs, ot = load_module()
            index = []
            for job in jobs:
                res = analyze(bs, ot, job, args)
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
                plot(bs, res, os.path.join(outdir, "band_diagram.png"),
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
            "outdir"]
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
                        "" if cmp_ is None else int(cmp_["matches"]), outdir])


if __name__ == "__main__":
    sys.exit(main())
