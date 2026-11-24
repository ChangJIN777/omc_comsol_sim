"""Objective: geometry -> (optical gap, mechanical gap, scalar score).

Backend selection (config/CLI), feasibility gating, caching, and the penalized
score from configs/targets.yaml. Optical can be 'surrogate' (numpy, always
available), 'mpb', 'comsol', or 'tidy3d' (FDTD on the Tidy3D cloud, COSTS
FlexCredits; see optical_tidy3d.py). Mechanical is 'comsol' (or
'surrogate_stub' that returns NaN until the COMSOL template exists).

Score:
  S = min(G_o, G_m)
      - lo  * max(0, gmin_o - G_o)^2
      - lm  * max(0, gmin_m - G_m)^2
      - lof * ((f_o,c - f_o,target)/f_o,target)^2
      - lmf * ((f_m,c - f_m,target)/f_m,target)^2
Infeasible geometry -> large negative score, no solve.
"""
from __future__ import annotations

import os
import hashlib
import json
from typing import Optional

import yaml

from geometry import (u_to_geometry, check_feasibility, load_bounds,
                      resolve_xsect, diamond_index, z_mirror_symmetric)
from bandgap import gap_near_frequency

_CFG = os.path.join(os.path.dirname(__file__), "..", "configs")
C0 = 299_792_458.0


def _load(name):
    with open(os.path.join(_CFG, name)) as fh:
        return yaml.safe_load(fh)


def _hash_u(u, backends):
    key = json.dumps([round(float(x), 6) for x in u] + list(backends))
    return "cand_" + hashlib.sha1(key.encode()).hexdigest()[:12]


# Optical engines that model the plain rectangular beam only.
RECT_ONLY_OPTICS = ("surrogate", "mpb")


def check_xsect_backend(xsect, optical_backend):
    """Refuse, before any solve, an engine that cannot model the cross-section."""
    xsect = resolve_xsect(xsect)
    if xsect != "rect" and optical_backend in RECT_ONLY_OPTICS:
        raise ValueError(
            f"optical: {optical_backend} models the rectangular cross-section "
            f"only; xsect: {xsect} needs optical: comsol or tidy3d")
    return xsect


def physics_id_tag(optical_backend, xsect, z_family=None, tidy3d_tag=None):
    """Part of the record id that pins what was simulated, or None.

    The cross-section, the diamond index and (COMSOL optics) the scored
    thickness-parity family change the result for the same u, so they are in
    the id: records with another cross-section or index never replace each
    other, and neither do records made before these options existed.
    None (id unchanged) only for a rectangular cell without optics.
    """
    xsect = resolve_xsect(xsect)
    if optical_backend == "tidy3d":
        return tidy3d_tag               # holds engine, xsect, index, k sampling
    if optical_backend in (None, "none"):
        return None if xsect == "rect" else f"xsect={xsect}"
    tag = f"n={diamond_index():.4f}|xsect={xsect}"
    if optical_backend == "comsol":
        tag += f"|z={z_family or 'all'}"
    return tag


def evaluate_candidate(u, *, optical_backend="surrogate",
                       mech_backend="comsol", cache=None, bounds=None,
                       targets=None, on_stage=None, celltype="hole",
                       rib_options=None, tidy3d_options=None, xsect=None,
                       verbose=True):
    """Evaluate one normalized design vector u. Returns a result record dict.

    celltype="rib" evaluates the rib cell instead (objective_rib.py);
    rib_options is then a dict of its keyword options (mech_gap_mode, n_k,
    n_bands_mech, optical_study, n_k_opt, n_bands_opt, optical_rel_tol,
    rxtal_deg). For
    the hole cell, whose path below is unchanged, leave both at their defaults.

    tidy3d_options is the `tidy3d:` config section, used by
    optical_backend="tidy3d" only (either cell type). Without
    `allow_cloud: true` and a `max_credits` budget nothing is uploaded and the
    record comes back with status "optical_not_run" and the key
    `optical_tidy3d_stop`; loops should save it and stop.

    xsect: cross-section, "rect" (default) or "wedge" (geometry.XSECTS), for
    both cell types. It selects the COMSOL template (optics AND mechanics: one
    device) and the Tidy3D geometry; the surrogate and MPB are rect only and
    raise ValueError for "wedge" before any solve. For OPTICS, "rect" keeps
    only the modes that are even across the thickness (the rectangular beam
    has that mirror plane): Tidy3D by its symmetry plane, COMSOL by the parity
    of the solved field; "wedge" keeps all modes. Mechanics always uses all
    modes. Records carry `xsect`, `optical_index` and `optical_z_symmetry`,
    and xsect is part of the id.

    on_stage, if given, is called as on_stage(name, event, info=None) at each
    stage boundary (name in {"feasibility","optical","mechanical","score"},
    event in {"start","done","skip","fail"}) -- e.g. for CLI progress
    reporting. It is a no-op hook: evaluate_candidate never raises because
    of it and callers that omit it see no behavior change. The Tidy3D
    engine also sends its sub-steps as event "info" of stage "optical".

    verbose (optical_backend="tidy3d" only): False silences Tidy3D's web
    progress and the "info" sub-steps (`quiet: true` of the scripts). All
    of that output goes to stderr, never to stdout.
    """
    def _stage(name, event, info=None):
        if on_stage is not None:
            on_stage(name, event, info)

    bounds = bounds or load_bounds()
    targets = targets or _load("targets.yaml")
    xsect = check_xsect_backend(xsect, optical_backend)
    if celltype == "rib":
        from objective_rib import evaluate_rib_candidate
        return evaluate_rib_candidate(
            u, bounds=bounds, targets=targets,
            optical_backend=optical_backend, mech_backend=mech_backend,
            cache=cache, on_stage=on_stage, tidy3d_options=tidy3d_options,
            xsect=xsect, verbose=verbose, **(rib_options or {}))
    if celltype != "hole":
        raise ValueError(f"unknown celltype {celltype!r}")
    if rib_options:
        raise ValueError("rib_options given for celltype 'hole'")
    # Thickness-parity family scored by COMSOL optics ("even" = TE-like for
    # a rectangular beam, None = all modes for the wedge).
    z_family = None
    if optical_backend == "comsol":
        from comsol_optics import study_z_family
        z_family = study_z_family("opt TE", xsect)
    tidy3d_tag = None
    if optical_backend == "tidy3d":
        # engine, xsect, index, k sampling and numerics: a Tidy3D record never
        # replaces a COMSOL one.
        from optical_tidy3d import id_settings
        tidy3d_tag = id_settings(tidy3d_options, "hole", xsect=xsect)
    tag = physics_id_tag(optical_backend, xsect, z_family, tidy3d_tag)
    cid = _hash_u(u, (optical_backend, mech_backend) + ((tag,) if tag else ()))

    if cache is not None and cid in cache:
        return cache[cid]

    _stage("feasibility", "start")
    g = u_to_geometry(u, bounds)
    ok, reasons = check_feasibility(g, bounds, xsect)
    rec = dict(id=cid, u=[float(x) for x in u], params=g.as_dict(),
               optical_backend=optical_backend, mech_backend=mech_backend,
               optical_engine=optical_backend, xsect=xsect,
               optical_index=diamond_index(),
               optical_z_symmetry=z_mirror_symmetric(xsect),
               status="pending", reasons=reasons)

    if not ok:
        _stage("feasibility", "fail", reasons)
        rec.update(status="infeasible", optical_gap=0.0, mechanical_gap=0.0,
                   score=-10.0,
                   optical_center_frequency=0.0, mechanical_center_frequency=0.0)
        if cache is not None:
            cache[cid] = rec
        return rec
    _stage("feasibility", "done")

    f_o_target = C0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    f_m_target = targets["mechanical"]["target_frequency_GHz"] * 1e9

    # ---- optical ----
    _stage("optical", "start", optical_backend)
    try:
        if optical_backend == "surrogate":
            from optical_surrogate import optical_gap_surrogate
            o = optical_gap_surrogate(g)
            G_o, f_o_c = o["normalized_gap"], o["center_freq_Hz"]
        elif optical_backend == "mpb":
            from optical_mpb import run_optical_mpb
            d = run_optical_mpb(g)
            gp = gap_near_frequency(d["freqs_hz"], f_o_target)
            G_o, f_o_c = gp.normalized_gap, gp.f_center
        elif optical_backend == "comsol":
            from optical_comsol import run_optical_comsol
            d = run_optical_comsol(g, xsect=xsect)
            gp = gap_near_frequency(d["freqs_hz"], f_o_target)
            G_o, f_o_c = gp.normalized_gap, gp.f_center
            # per k: every returned mode and its thickness-parity fraction
            rec.update(optical_comsol_detail=dict(
                z_family=d.get("z_family"),
                freqs_all_hz=d.get("freqs_all_hz"),
                z_even_fraction=d.get("z_even_fraction")))
        elif optical_backend == "tidy3d":
            # Feasibility (min feature) was enforced above, before this solver.
            from optical_tidy3d import run_optical_tidy3d
            o = run_optical_tidy3d(g, celltype="hole", eval_id=cid,
                                   f_target=f_o_target, options=tidy3d_options,
                                   xsect=xsect, verbose=verbose,
                                   on_progress=lambda m: _stage(
                                       "optical", "info", m))
            G_o, f_o_c = o["gap"].normalized_gap, o["gap"].f_center
            rec.update(o["record"])
        else:
            raise ValueError(f"unknown optical_backend {optical_backend}")
    except Exception as e:  # solver failed -> record, don't crash loop
        _stage("optical", "fail", str(e))
        if getattr(e, "tidy3d_stop", None):
            # Nothing was run or charged (cloud not allowed, dry run, credit
            # budget, ...). Not cached: the next call may be allowed to run.
            from optical_tidy3d import not_run_record
            return not_run_record(rec, e)
        rec.update(status="optical_failed", error=str(e),
                   optical_gap=0.0, mechanical_gap=0.0, score=-5.0,
                   optical_center_frequency=0.0, mechanical_center_frequency=0.0)
        if cache is not None:
            cache[cid] = rec
        return rec
    _stage("optical", "done")

    # ---- mechanical ----
    _stage("mechanical", "start" if mech_backend != "surrogate_stub" else "skip",
           mech_backend)
    try:
        if mech_backend == "comsol":
            from acoustic_comsol import run_mechanical_comsol
            from comsol_client import template_path
            # same device as the optics: template by cross-section. Mechanics
            # is full thickness, all modes, whatever the cross-section.
            d = run_mechanical_comsol(
                g, template=template_path("hole", xsect, must_exist=True))
            gp = gap_near_frequency(d["freqs_hz"], f_m_target)
            G_m, f_m_c = gp.normalized_gap, gp.f_center
        elif mech_backend == "surrogate_stub":
            G_m, f_m_c = float("nan"), float("nan")
        else:
            raise ValueError(f"unknown mech_backend {mech_backend}")
    except Exception as e:
        _stage("mechanical", "fail", str(e))
        rec.update(status="mech_failed", error=str(e),
                   optical_gap=G_o, mechanical_gap=0.0, score=-2.0,
                   optical_center_frequency=f_o_c,
                   mechanical_center_frequency=0.0)
        if cache is not None:
            cache[cid] = rec
        return rec
    if mech_backend != "surrogate_stub":
        _stage("mechanical", "done")

    _stage("score", "start")
    score = _score(G_o, G_m, f_o_c, f_m_c, f_o_target, f_m_target, targets)
    rec.update(status="success", optical_gap=float(G_o),
               mechanical_gap=float(G_m) if G_m == G_m else None,
               optical_center_frequency=float(f_o_c),
               mechanical_center_frequency=(float(f_m_c) if f_m_c == f_m_c else None),
               score=float(score))
    _stage("score", "done")
    if cache is not None:
        cache[cid] = rec
    return rec


def _score(G_o, G_m, f_o_c, f_m_c, f_o_t, f_m_t, targets):
    s = targets["scoring"]
    gmin_o = targets["optical"]["min_gap"]
    gmin_m = targets["mechanical"]["min_gap"]
    # if mechanical not yet available (NaN), score on optical only (pre-screen)
    if G_m != G_m:
        base = G_o
        pen = s["lambda_optical"] * max(0.0, gmin_o - G_o) ** 2
        freqpen = s["lambda_optical_freq"] * ((f_o_c - f_o_t) / f_o_t) ** 2
        return base - pen - freqpen
    base = min(G_o, G_m)
    pen = (s["lambda_optical"] * max(0.0, gmin_o - G_o) ** 2
           + s["lambda_mech"] * max(0.0, gmin_m - G_m) ** 2)
    freqpen = (s["lambda_optical_freq"] * ((f_o_c - f_o_t) / f_o_t) ** 2
               + s["lambda_mech_freq"] * ((f_m_c - f_m_t) / f_m_t) ** 2)
    return base - pen - freqpen
