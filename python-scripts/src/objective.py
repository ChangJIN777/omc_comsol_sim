"""Objective: geometry -> (optical gap, mechanical gap, scalar score).

Backend selection (config/CLI), feasibility gating, caching, and the penalized
score from configs/targets.yaml. Optical can be 'surrogate' (numpy, always
available), 'mpb', or 'comsol'. Mechanical is 'comsol' (or 'surrogate_stub'
that returns NaN until the COMSOL template exists).

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

from geometry import u_to_geometry, check_feasibility, load_bounds
from bandgap import gap_near_frequency

_CFG = os.path.join(os.path.dirname(__file__), "..", "configs")
C0 = 299_792_458.0


def _load(name):
    with open(os.path.join(_CFG, name)) as fh:
        return yaml.safe_load(fh)


def _hash_u(u, backends):
    key = json.dumps([round(float(x), 6) for x in u] + list(backends))
    return "cand_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def evaluate_candidate(u, *, optical_backend="surrogate",
                       mech_backend="comsol", cache=None, bounds=None,
                       targets=None, on_stage=None):
    """Evaluate one normalized design vector u. Returns a result record dict.

    on_stage, if given, is called as on_stage(name, event, info=None) at each
    stage boundary (name in {"feasibility","optical","mechanical","score"},
    event in {"start","done","skip","fail"}) -- e.g. for CLI progress
    reporting. It is a no-op hook: evaluate_candidate never raises because
    of it and callers that omit it see no behavior change.
    """
    def _stage(name, event, info=None):
        if on_stage is not None:
            on_stage(name, event, info)

    bounds = bounds or load_bounds()
    targets = targets or _load("targets.yaml")
    cid = _hash_u(u, (optical_backend, mech_backend))

    if cache is not None and cid in cache:
        return cache[cid]

    _stage("feasibility", "start")
    g = u_to_geometry(u, bounds)
    ok, reasons = check_feasibility(g, bounds)
    rec = dict(id=cid, u=[float(x) for x in u], params=g.as_dict(),
               optical_backend=optical_backend, mech_backend=mech_backend,
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
            d = run_optical_comsol(g)
            gp = gap_near_frequency(d["freqs_hz"], f_o_target)
            G_o, f_o_c = gp.normalized_gap, gp.f_center
        else:
            raise ValueError(f"unknown optical_backend {optical_backend}")
    except Exception as e:  # solver failed -> record, don't crash loop
        _stage("optical", "fail", str(e))
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
            d = run_mechanical_comsol(g)
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
