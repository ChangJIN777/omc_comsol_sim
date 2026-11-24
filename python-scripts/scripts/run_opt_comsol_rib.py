#!/usr/bin/env python3
"""Rib-cell branch of scripts/run_opt_comsol.py (celltype: rib).

Reached through the usual entry point:
    python scripts/run_opt_comsol.py --config configs/run_opt_comsol_rib.yaml
or directly, e.g. to redo the end-of-run characterization of the best cell:
    python scripts/run_opt_comsol_rib.py --config configs/run_opt_comsol_rib.yaml --characterize

Loop: ask -> evaluate (objective_rib.evaluate_rib_candidate: validity gate,
fabricability gate, optical 'opt TE' gap, complete mechanical gap) -> save to
the rib sqlite file and to out_json -> tell. Every candidate is saved,
including infeasible / unfabricable / failed ones.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src"))

from geometry import (load_bounds, n_dim, phys_to_u, rib_search_vars,  # noqa: E402
                      RibGeometry, resolve_rxtal)
from objective import check_xsect_backend                             # noqa: E402
from objective_rib import (evaluate_rib_candidate, score_version,     # noqa: E402
                           DEFAULT_OPTICAL_TOL_THZ, RIB_K_GRID_REL)
from optimizer import initialize_optimizer                            # noqa: E402
from database import save_result, db_path                             # noqa: E402

C0 = 299_792_458.0
_CFG_DIR = os.path.join(_HERE, "..", "configs")
_DEFAULT_CONFIG = os.path.join(_CFG_DIR, "run_opt_comsol_rib.yaml")
_RIB_OPTION_KEYS = ("mech_gap_mode", "n_k", "n_bands_mech", "optical_study",
                    "n_k_opt", "n_bands_opt", "optical_tol_thz", "rxtal_deg",
                    # deprecated, ignored with a FutureWarning by
                    # objective_rib.evaluate_rib_candidate:
                    "optical_rel_tol")


def _load_targets():
    with open(os.path.join(_CFG_DIR, "targets.yaml")) as fh:
        return yaml.safe_load(fh)


def _rib_options(cfg):
    rib = cfg.get("rib") or {}
    unknown = sorted(set(rib) - set(_RIB_OPTION_KEYS))
    if unknown:
        raise ValueError(f"unknown key(s) under `rib`: {unknown}; "
                         f"expected a subset of {_RIB_OPTION_KEYS}")
    return dict(rib)


def _f(x):
    return float("nan") if x is None else float(x)


def _fmt_time(seconds):
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    return f"{h}h{m:02d}m" if h else (f"{m}m{s:02d}s" if m else f"{s}s")


def _save_json(results, path):
    """All records, without the `_`-prefixed band arrays (those are in sqlite)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    slim = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    with open(path, "w") as fh:
        json.dump(slim, fh, indent=2)


# ---- drawing ----------------------------------------------------------------

def _rib_outline(g, n_periods):
    """Top-view (z, y) outline of n_periods rib cells, in nm (sharp corners)."""
    nm = 1e9
    a, w, t, s = g.a * nm, g.w * nm, g.t * nm, g.s * nm
    z0 = -n_periods * a / 2
    top_z, top_y = [z0], [s / 2]
    for i in range(n_periods):
        zc = z0 + (i + 0.5) * a
        top_z += [zc - t / 2, zc - t / 2, zc + t / 2, zc + t / 2]
        top_y += [s / 2, w / 2, w / 2, s / 2]
    top_z.append(-z0)
    top_y.append(s / 2)
    return top_z + top_z[::-1], top_y + [-y for y in top_y[::-1]]


def _draw_rib(ax, g, n_periods=3, fontsize=8):
    """Top view (z-y) of the rib cell with its dimensions."""
    nm = 1e9
    z, y = _rib_outline(g, n_periods)
    ax.fill(z, y, facecolor="#cfe8ff", edgecolor="#1f5fa6", lw=1.3)
    half_len, half_w = n_periods * g.a * nm / 2, g.w * nm / 2
    for zc in (-g.a * nm / 2, g.a * nm / 2):          # one period, dashed
        ax.plot([zc, zc], [-half_w, half_w], "k--", lw=0.6)
    ax.set_title(
        f"Rib cell, top view (z-y), fillets not drawn\n"
        f"a={g.a*nm:.0f}  w={g.w*nm:.0f}  t={g.t*nm:.0f}  s={g.s*nm:.0f}  "
        f"th={g.th*nm:.0f}  r1={g.r1*nm:.0f}  r2={g.r2*nm:.0f} nm\n"
        f"gap a-t = {(g.a - g.t)*nm:.0f} nm    arm (w-s)/2 = {(g.w - g.s)/2*nm:.0f} nm",
        fontsize=fontsize)
    ax.set_xlabel("z  [nm]  (beam / periodic axis)", fontsize=fontsize)
    ax.set_ylabel("y  [nm]", fontsize=fontsize)
    ax.set_aspect("equal")
    ax.set_xlim(-1.08 * half_len, 1.08 * half_len)
    ax.set_ylim(-1.15 * half_w, 1.15 * half_w)
    ax.tick_params(labelsize=fontsize - 1)


def _plot_progress(results, out_path):
    """Score history, gap history and the best cell so far."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    it = np.arange(len(results))
    score = np.array([_f(r.get("score")) for r in results])
    ok = np.array([r.get("status") == "success" for r in results], dtype=bool)

    ax = axes[0]
    if ok.any():
        ax.plot(it[ok], score[ok], "o", ms=4, color="#1f5fa6", label="success")
        best = np.maximum.accumulate(np.where(ok, score, -np.inf))
        ax.step(it, np.where(np.isfinite(best), best, np.nan), where="post",
                color="k", lw=1, label="best so far")
    if (~ok).any():
        ax.plot(it[~ok], score[~ok], "x", ms=4, color="0.6",
                label="infeasible / unfabricable / failed")
    ax.set_xlabel("iteration")
    ax.set_ylabel("score")
    ax.legend(fontsize=7)

    ax = axes[1]
    g_m = np.array([100 * _f(r.get("mechanical_gap")) for r in results])
    g_o = np.array([100 * _f(r.get("optical_gap")) for r in results])
    if ok.any():
        ax.plot(it[ok], g_m[ok], "s", ms=4, color="#c0392b", label="G_m [%]")
        ax.plot(it[ok], g_o[ok], "o", ms=4, color="#1f5fa6", label="G_o [%]")
    ax.axhline(20.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("iteration")
    ax.set_ylabel("normalized gap [%]")
    ax.legend(fontsize=7)

    ax = axes[2]
    succ = [r for r in results if r.get("status") == "success"]
    if succ:
        _draw_rib(ax, RibGeometry(**max(succ, key=lambda r: r["score"])["params"]))
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no successful evaluation yet", ha="center")

    n_status = {}
    for r in results:
        n_status[r.get("status")] = n_status.get(r.get("status"), 0) + 1
    fig.suptitle("Rib cell optimization  |  " +
                 "  ".join(f"{k}: {v}" for k, v in sorted(n_status.items())),
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---- optimization loop --------------------------------------------------------

def _same_scoring(rec, version, tol_thz):
    """Was this record scored the way the current run scores?

    Records without score_version are version 1 (the +/- optical_rel_tol
    window with quadratic penalties). For version 2 the Gaussian width must
    match too: another optical_tol_thz is another objective.
    """
    if int(rec.get("score_version") or 1) != int(version):
        return False
    if version == 1:
        return True
    old = rec.get("optical_tol_thz")
    return old is not None and abs(float(old) - float(tol_thz)) < 1e-12


def _warm_start(opt, prior, bounds, nd, rxtal=0.0, version=1,
                tol_thz=DEFAULT_OPTICAL_TOL_THZ):
    """Seed an Optuna study with prior successful records (u recomputed from
    the stored physical params with the CURRENT bounds). Only records of the
    same crystal rotation `rxtal` (deg; records without rxtal_deg count as 0)
    and of the same scoring (_same_scoring: score version and, for version 2,
    optical_tol_thz) are used, since scores of another definition would
    mislead the optimizer. Returns the count."""
    study = getattr(opt, "study", None)
    if study is None or not prior:
        return 0
    import optuna
    dists = {f"u{i}": optuna.distributions.FloatDistribution(0, 1)
             for i in range(nd)}
    n_added = 0
    for r in prior:
        if r.get("status") != "success" or not r.get("params"):
            continue
        if float(r.get("rxtal_deg") or 0.0) != rxtal:
            continue
        if not _same_scoring(r, version, tol_thz):
            continue
        u = phys_to_u(r["params"], bounds, celltype="rib")
        study.add_trial(optuna.trial.create_trial(
            params={f"u{i}": float(u[i]) for i in range(nd)},
            distributions=dists, value=float(r["score"])))
        n_added += 1
    return n_added


def main_rib(cfg, cfg_path):
    bounds = load_bounds()
    targets = _load_targets()
    names = rib_search_vars(bounds)
    nd = n_dim("rib", bounds)
    opts = _rib_options(cfg)
    optical, mech = cfg.get("optical", "comsol"), cfg.get("mech", "comsol")
    # optical: tidy3d -> the `tidy3d:` section (cloud opt-in; its max_credits
    # is the credit budget of this whole loop). See src/optical_tidy3d.py.
    tidy3d_options = cfg.get("tidy3d") if optical == "tidy3d" else None
    # xsect: rect (default) | wedge. Selects the template (optics + mechanics)
    # and the Tidy3D geometry; rect scores thickness-even optical modes only.
    xsect = check_xsect_backend(cfg.get("xsect"), optical)
    out_json, out_fig = cfg["out_json"], cfg["out_fig"]
    db = cfg.get("db_path") or db_path("rib")
    # crystal rotation about the slab normal (mechanics only): 0 or 45 deg
    rxtal = resolve_rxtal(opts.get("rxtal_deg"))

    if "comsol" in (optical, mech):
        from comsol_client import (require_rib_air_disk,
                                   require_rxtal_template, template_path)
        try:
            path = template_path("rib", xsect, must_exist=True)
            if optical == "comsol":            # the old air box: optics only
                require_rib_air_disk(path)
            if mech == "comsol":
                require_rxtal_template(path, rxtal)
        except (FileNotFoundError, RuntimeError) as err:   # missing or outdated
            print(f"[FAIL] {err}")
            sys.exit(1)

    prior = []
    if os.path.exists(out_json):
        if cfg.get("fresh_start", False):
            shutil.copy2(out_json, out_json + ".bak")
            print(f"Fresh start: backed up prior results to {out_json}.bak")
        else:
            with open(out_json) as fh:
                prior = json.load(fh)
            print(f"Auto-resumed: loaded {len(prior)} prior results from {out_json}")

    opt = initialize_optimizer(None, backend=cfg.get("backend", "auto"),
                               n_init=cfg["n_init"], n_dim=nd)
    # scoring of this run: version (2 as soon as optics is scored) and the
    # Gaussian width; prior records scored otherwise are not used as seeds
    version = score_version(optical)
    tol_thz = float(opts.get("optical_tol_thz", DEFAULT_OPTICAL_TOL_THZ))
    n_seed = _warm_start(opt, prior, bounds, nd, rxtal=rxtal, version=version,
                         tol_thz=tol_thz)
    n_other = sum(1 for r in prior if r.get("status") == "success"
                  and not _same_scoring(r, version, tol_thz))
    print(f"Rib cell | variables {names} | optimizer {type(opt).__name__} "
          f"(seeded with {n_seed}) | optical={optical} mech={mech} "
          f"xsect={xsect} {opts}")
    print(f"Scoring: version {version}"
          + (f", optical_tol_thz = {tol_thz:g}" if version != 1 else "")
          + (f" | {n_other} prior result(s) with another scoring are kept in "
             f"{out_json} but not used as seeds" if n_other else ""))
    print(f"Results: {out_json}  |  sqlite: {db}  |  figure: {out_fig}")

    # tidy3d.keep_hdf5: best -> only the best cell keeps its .hdf5 files.
    # Seeded with the best of earlier sessions (same scoring), which keeps its
    # files until a better cell with files appears. See optical_tidy3d.Hdf5Keeper.
    keeper = None
    if optical == "tidy3d":
        from optical_tidy3d import Hdf5Keeper
        keeper = Hdf5Keeper(tidy3d_options)
        keeper.seed([r for r in prior
                     if float(r.get("rxtal_deg") or 0.0) == rxtal
                     and _same_scoring(r, version, tol_thz)])
        if keeper.policy != "all":
            print(f"Tidy3D keep_hdf5: {keeper.policy} (params.json, plan.json, "
                  "bands.csv, evaluation.json are kept for every cell)")

    results = list(prior)
    best = None
    t_start = time.time()
    for it in range(cfg["n_iter"]):
        u = opt.ask()
        rec = evaluate_rib_candidate(u, bounds=bounds, targets=targets,
                                     optical_backend=optical,
                                     mech_backend=mech,
                                     tidy3d_options=tidy3d_options,
                                     xsect=xsect, **opts)
        if rec.get("optical_tidy3d_stop"):
            # Tidy3D did not run (cloud not allowed, dry run, credit budget
            # reached, ...) and would not for the next candidate either: save
            # everything and stop. Not a result, so the optimizer is not told.
            save_result(rec, db)
            results.append(rec)
            _save_json(results, out_json)
            _plot_progress(results, out_fig)
            print(f"[{it:3d}] STOP ({rec['optical_tidy3d_stop']}): {rec['error']}")
            break
        opt.tell(u, rec["score"])
        save_result(rec, db)
        results.append(rec)

        p = rec["params"]
        dims = "  ".join(f"{k}={p[k]*1e9:.0f}" for k in names)
        if rec["status"] == "success":
            if best is None or rec["score"] > best["score"]:
                best = rec
            print(f"[{it:3d}] {dims} nm | "
                  f"G_o={100*_f(rec['optical_gap']):.1f}%  "
                  f"G_m={100*_f(rec['mechanical_gap']):.1f}%  "
                  f"f_m={_f(rec['mechanical_center_frequency'])*1e-9:.2f} GHz  "
                  f"score={rec['score']:+.4f}"
                  + ("  *** BEST" if best is rec else "")
                  + f"  [{_fmt_time(time.time() - t_start)}]")
        else:
            print(f"[{it:3d}] {dims} nm | {rec['status'].upper()}  "
                  f"{rec.get('reasons') or rec.get('error', '')}")

        _save_json(results, out_json)
        _plot_progress(results, out_fig)
        if keeper is not None:
            # only now, with the record in sqlite and out_json: an interrupted
            # run never lacks files of a cell that was not yet superseded
            keeper.update(rec)

    print("\n" + "=" * 60 + "\nRIB OPTIMIZATION COMPLETE\n" + "=" * 60)
    print(f"Wall time: {_fmt_time(time.time() - t_start)}")
    if optical == "tidy3d":
        from optical_tidy3d import session_credits
        print(f"Tidy3D cost estimates of this run: {session_credits():.3f} FlexCredits")
        if keeper is not None and keeper.active:
            kept = keeper.best[1] if keeper.best else "none"
            print(f"Tidy3D keep_hdf5: best -> {keeper.freed_bytes / 1e6:.1f} MB "
                  f"deleted; .hdf5 kept for: {kept}")
    if best is not None:
        print(f"Best this run: score={best['score']:+.4f}  "
              f"G_o={100*_f(best['optical_gap']):.1f}%  "
              f"G_m={100*_f(best['mechanical_gap']):.1f}%")
        print("  params (nm):", {k: round(v * 1e9, 1)
                                 for k, v in best["params"].items()})
        print("  u:", [round(x, 4) for x in best["u"]])

    if cfg.get("skip_characterization", False) or "comsol" not in (optical, mech):
        print("(Characterization skipped.)")
        return
    # Fresh process -> clean COMSOL state, as for the hole cell.
    log = os.path.join(os.path.dirname(out_fig) or ".", "char_log_rib.txt")
    print(f"\n  [char] characterizing the best cell in a fresh process -> {log}")
    with open(log, "w") as fh:
        proc = subprocess.run([sys.executable, os.path.abspath(__file__),
                               "--characterize", "--config", cfg_path],
                              stdout=fh, stderr=subprocess.STDOUT)
    print("  [char] done" if proc.returncode == 0
          else f"  [char] FAILED (exit code {proc.returncode}); see {log}")


# ---- characterization of the best cell ---------------------------------------

def characterize_rib(cfg):
    """Full-zone bands of the best saved rib cell: opt TE + TM, mech sym +
    antisym. Saves best_rib_characterization.png and best_rib_char_data.npz
    next to out_fig. No mode profiles (see the proposal's TODO list)."""
    from acoustic_comsol import run_mechanical_comsol
    from bandgap import interval_gap
    from comsol_client import get_model, set_geometry, template_path
    from comsol_client import THICKNESS_PARAM, require_rib_air_disk
    from comsol_optics import (_full_opt_sweep, guided_interval_gaps,
                               opt_sweep_at, rib_k_grid, study_z_family,
                               select_gap_by_penalized_size)

    with open(cfg["out_json"]) as fh:
        results = json.load(fh)
    opts = _rib_options(cfg)
    # the crystal rotation of this config; records without rxtal_deg are 0
    rxtal = resolve_rxtal(opts.get("rxtal_deg"))
    # ... and its scoring: scores of another definition are not comparable
    version = score_version(cfg.get("optical", "comsol"))
    tol_thz = float(opts.get("optical_tol_thz", DEFAULT_OPTICAL_TOL_THZ))
    succ = [r for r in results if r.get("status") == "success"
            and float(r.get("rxtal_deg") or 0.0) == rxtal
            and _same_scoring(r, version, tol_thz)]
    if not succ:
        print("No successful rib result to characterize.")
        return None
    best = max(succ, key=lambda r: r["score"])
    g = RibGeometry(**best["params"])     # stored physical params, not u
    targets = _load_targets()
    f_o_t = C0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    f_m_t = targets["mechanical"]["target_frequency_GHz"] * 1e9
    tol_hz = float(opts.get("optical_tol_thz", DEFAULT_OPTICAL_TOL_THZ)) * 1e12
    study_o = opts.get("optical_study", "opt TE")
    out_dir = os.path.dirname(cfg["out_fig"]) or "."
    os.makedirs(out_dir, exist_ok=True)
    # Same cross-section as the optimization: its template, and for "rect"
    # the thickness-parity family of each optical study (None = all modes).
    xsect = check_xsect_backend(cfg.get("xsect"), "comsol")
    template = require_rib_air_disk(template_path("rib", xsect, must_exist=True))
    zkw = lambda study: dict(z_family=study_z_family(study, xsect),  # noqa: E731
                             thickness=THICKNESS_PARAM["rib"])
    n_k, n_bo = 15, 12
    n_bm = int(opts.get("n_bands_mech", 20))

    print("  [1/5] mech sym ...")
    sym = run_mechanical_comsol(g, n_k=n_k, n_bands=n_bm, template=template,
                                study="mech sym", rxtal_deg=rxtal)
    print("  [2/5] mech antisym ...")
    anti = run_mechanical_comsol(g, n_k=n_k, n_bands=n_bm, template=template,
                                 study="mech antisym", rxtal_deg=rxtal)
    k_m = np.asarray(sym["k_frac"], dtype=float)
    f_sym = np.asarray(sym["freqs_hz"], dtype=float)
    f_anti = np.asarray(anti["freqs_hz"], dtype=float)
    gap_m = interval_gap([f_sym, f_anti], f_m_t)

    model = get_model(template)
    set_geometry(model, g)
    print("  [3/5] opt TE (full zone) ...")
    k_te, f_te = _full_opt_sweep(model, g, "opt TE", n_k=n_k, n_bands=n_bo,
                                 **zkw("opt TE"))
    print("  [4/5] opt TM (full zone) ...")
    k_tm, f_tm = _full_opt_sweep(model, g, "opt TM", n_k=n_k, n_bands=n_bo,
                                 **zkw("opt TM"))
    print(f"  [5/5] {study_o} gap sweep (as in the optimization) ...")
    k_det = rib_k_grid(g.a, int(opts.get("n_k_opt", 13)), f_o_t, RIB_K_GRID_REL)
    k_det, f_det = opt_sweep_at(model, g, study_o, k_det,
                                int(opts.get("n_bands_opt", 10)), f_target=f_o_t,
                                **zkw(study_o))
    # as the optimization: every guided gap, the best gap * Gaussian weight
    gap_o, pen_o, _ = select_gap_by_penalized_size(
        guided_interval_gaps(k_det, f_det, g.a), f_o_t, tol_hz)

    np.savez(os.path.join(out_dir, "best_rib_char_data.npz"),
             params_names=np.array(sorted(best["params"])),
             params_values=np.array([best["params"][k] for k in sorted(best["params"])]),
             k_m=k_m, freqs_sym=f_sym, freqs_anti=f_anti,
             k_te=k_te, freqs_te=f_te, k_tm=k_tm, freqs_tm=f_tm,
             k_det=k_det, freqs_det=f_det)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    _draw_rib(axes[0, 0], g)

    ax = axes[0, 1]
    ax.axis("off")
    o_txt = (f"{100*gap_o.normalized_gap:.1f} %  ({gap_o.f_lower*1e-12:.1f}-"
             f"{gap_o.f_upper*1e-12:.1f} THz, {C0/gap_o.f_center*1e9:.0f} nm), "
             f"weight {pen_o:.3f}"
             if gap_o.found else "none")
    m_txt = (f"{100*gap_m.normalized_gap:.1f} %  ({gap_m.f_lower*1e-9:.2f}-"
             f"{gap_m.f_upper*1e-9:.2f} GHz)" if gap_m.found else "none")
    lines = ["Best rib cell", ""]
    lines += [f"  {k:3s} = {v*1e9:7.1f} nm" for k, v in best["params"].items()]
    lines += ["",
              f"  score (optimization) = {best['score']:+.4f}",
              f"  G_o (optimization)   = {100*_f(best.get('optical_gap')):.1f} %",
              f"  G_m (optimization)   = {100*_f(best.get('mechanical_gap')):.1f} %",
              "",
              f"  {study_o} guided gap : {o_txt}",
              f"  complete mech gap   : {m_txt}",
              f"  light line at X     : {C0/(2*g.a)*1e-12:.1f} THz"]
    ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes, va="top",
            family="monospace", fontsize=9)

    ax = axes[1, 0]
    ll = C0 * k_te / (2.0 * g.a) * 1e-12
    y_max = 1.6 * f_o_t * 1e-12
    ax.fill_between(k_te, ll, y_max, color="0.88", label="above light line")
    ax.plot(k_te, f_te * 1e-12, "o", ms=3, color="#1f5fa6")
    ax.plot(k_tm, f_tm * 1e-12, "o", ms=3, mfc="none", color="#c0392b")
    ax.plot([], [], "o", ms=3, color="#1f5fa6", label="opt TE (PEC at y=0)")
    ax.plot([], [], "o", ms=3, mfc="none", color="#c0392b", label="opt TM (PMC at y=0)")
    ax.axhline(f_o_t * 1e-12, color="k", ls=":", lw=0.8)
    if gap_o.found:
        ax.axhspan(gap_o.f_lower * 1e-12, gap_o.f_upper * 1e-12,
                   color="#f1c40f", alpha=0.35, label=f"{study_o} guided gap")
    ax.set_ylim(0, y_max)
    ax.set_xlabel(r"$k_z a / \pi$")
    ax.set_ylabel("Frequency [THz]")
    ax.set_title("Optical bands", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[1, 1]
    ax.plot(k_m, f_sym * 1e-9, "o", ms=3, color="#1f5fa6")
    ax.plot(k_m, f_anti * 1e-9, "o", ms=3, mfc="none", color="#c0392b")
    ax.plot([], [], "o", ms=3, color="#1f5fa6", label="y-even (mech sym)")
    ax.plot([], [], "o", ms=3, mfc="none", color="#c0392b", label="y-odd (mech antisym)")
    ax.axhline(f_m_t * 1e-9, color="k", ls=":", lw=0.8)
    if gap_m.found:
        ax.axhspan(gap_m.f_lower * 1e-9, gap_m.f_upper * 1e-9,
                   color="#f1c40f", alpha=0.35, label="complete gap")
    ax.set_xlabel(r"$k_z a / \pi$")
    ax.set_ylabel("Frequency [GHz]")
    ax.set_title("Mechanical bands", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")

    fig.tight_layout()
    out_path = os.path.join(out_dir, "best_rib_characterization.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Characterization -> {out_path}")
    return out_path


def _cli():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="YAML config (default: configs/run_opt_comsol_rib.yaml)")
    ap.add_argument("--characterize", action="store_true",
                    help="only characterize the best result in out_json")
    cli = ap.parse_args()
    with open(cli.config) as fh:
        cfg = yaml.safe_load(fh)
    if cfg.get("celltype") != "rib":
        ap.error("this script needs a config with `celltype: rib`")
    if cli.characterize:
        characterize_rib(cfg)
    else:
        main_rib(cfg, cli.config)


if __name__ == "__main__":
    _cli()
