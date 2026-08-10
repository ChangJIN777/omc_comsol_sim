"""Objective: 2D boomerang geometry -> (mechanical gap, [future] optical gap, score).

Adapted from python-scripts/src/objective.py (1D nanobeam project), with two
changes specific to this 2D phononic-crystal project:

1. Mechanical gap extraction is keyed on Z-PARITY (evenz/oddz), not a
   breathing-mode fy fraction -- see targets_2d.yaml's `gap_mode` and
   `_mechanical_gap` below. Both parities are always solved (the mechanical
   backend is cheap relative to optical) so "symmetry" vs "complete" gap_mode
   can be selected/re-scored after the fact without re-solving.
2. Optical is NOT YET IMPLEMENTED (src/optical_comsol_2d.py is a documented
   stub). `optical_backend="none"` (the default) simply excludes it from the
   score -- no solve attempted. Set `require_opt=True` and
   `optical_backend="comsol"` once that module is implemented; until then
   doing so will record a clean "optical_failed" status rather than crashing
   the loop (see evaluate_candidate's try/except below), so a run config that
   flips it on early won't silently corrupt results.

Score (mirrors python-scripts/src/objective.py's formula, extended with the
require_opt/require_mech toggle from python-scripts/scripts/run_opt_comsol.py
so a mechanical-only run -- the only backend implemented so far -- doesn't
need a fake optical term):
  S = (G_m if require_mech else 0)
      - lm  * max(0, gmin_m - G_m)^2                [if require_mech]
      - lmf * ((f_m,c - f_m,target)/f_m,target)^2    [if require_mech]
      + (G_o if require_opt else 0)
      - lo  * max(0, gmin_o - G_o)^2                 [if require_opt]
      - lof * ((f_o,c - f_o,target)/f_o,target)^2    [if require_opt]
Infeasible geometry -> large negative score, no solve.
"""
from __future__ import annotations

import os
import hashlib
import json

import numpy as np
import yaml

from geometry2d import u_to_geometry, check_feasibility, load_bounds
from bandgap import gap_near_frequency, largest_gap

_CFG = os.path.join(os.path.dirname(__file__), "..", "configs")
C0 = 299_792_458.0


def _load(name):
    with open(os.path.join(_CFG, name)) as fh:
        return yaml.safe_load(fh)


def _hash_u(u, backends):
    key = json.dumps([round(float(x), 6) for x in u] + list(backends))
    return "cand2d_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def _mechanical_gap(freqs_evenz, freqs_oddz, f_target, gap_mode):
    """Gap selection keyed on z-parity, per targets_2d.yaml's gap_mode.

    'symmetry': best single-family gap (evenz or oddz, whichever is larger
                near f_target) -- a quasi/symmetry-restricted gap.
    'complete': a gap with NO band of EITHER parity -- stack both families'
                bands together before searching.
    """
    if gap_mode == "complete":
        combined = np.hstack([freqs_evenz, freqs_oddz])
        return gap_near_frequency(combined, f_target)
    gp_evenz = gap_near_frequency(freqs_evenz, f_target)
    gp_oddz = gap_near_frequency(freqs_oddz, f_target)
    return gp_evenz if gp_evenz.normalized_gap >= gp_oddz.normalized_gap else gp_oddz


def evaluate_candidate(u, *, mech_backend="comsol", optical_backend="none",
                       require_mech=True, require_opt=False,
                       n_bands_mech=None, cache=None, bounds=None,
                       targets=None, on_stage=None):
    """Evaluate one normalized design vector u. Returns a result record dict.

    on_stage, if given, is called as on_stage(name, event, info=None) at each
    stage boundary (name in {"feasibility","mechanical","optical","score"},
    event in {"start","done","skip","fail"}) -- see scripts/run_one_2d.py for
    a stderr-reporting implementation. Never raises because of this hook.
    """
    def _stage(name, event, info=None):
        if on_stage is not None:
            on_stage(name, event, info)

    bounds = bounds or load_bounds()
    targets = targets or _load("targets_2d.yaml")
    cid = _hash_u(u, (mech_backend, optical_backend))

    if cache is not None and cid in cache:
        return cache[cid]

    _stage("feasibility", "start")
    g = u_to_geometry(u, bounds)
    ok, reasons = check_feasibility(g, bounds)
    rec = dict(id=cid, u=[float(x) for x in u], params=g.as_dict(),
               mech_backend=mech_backend, optical_backend=optical_backend,
               status="pending", reasons=reasons)

    if not ok:
        _stage("feasibility", "fail", reasons)
        rec.update(status="infeasible", mechanical_gap=0.0, optical_gap=0.0,
                   score=-10.0, mechanical_center_frequency=0.0,
                   optical_center_frequency=0.0)
        if cache is not None:
            cache[cid] = rec
        return rec
    _stage("feasibility", "done")

    f_m_target = targets["mechanical"]["target_frequency_GHz"] * 1e9
    f_o_target = C0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    gap_mode = targets["mechanical"].get("gap_mode", "symmetry")
    n_bands_mech = n_bands_mech or targets["mechanical"].get("n_bands", 10)

    # ---- mechanical (always solved when require_mech, both z-parities) ----
    _stage("mechanical", "start" if require_mech else "skip", mech_backend)
    G_m, f_m_c = 0.0, 0.0
    if require_mech:
        try:
            if mech_backend == "comsol":
                from acoustic_comsol_2d import run_mechanical_comsol_2d
                d = run_mechanical_comsol_2d(g, n_bands=n_bands_mech)
                gp = _mechanical_gap(d["freqs_evenz"], d["freqs_oddz"],
                                     f_m_target, gap_mode)
                G_m, f_m_c = gp.normalized_gap, gp.f_center
            else:
                raise ValueError(f"unknown mech_backend {mech_backend!r}")
        except Exception as e:  # solver failed -> record, don't crash loop
            _stage("mechanical", "fail", str(e))
            rec.update(status="mech_failed", error=str(e),
                       mechanical_gap=0.0, optical_gap=0.0, score=-5.0,
                       mechanical_center_frequency=0.0,
                       optical_center_frequency=0.0)
            if cache is not None:
                cache[cid] = rec
            return rec
        _stage("mechanical", "done")

    # ---- optical (not yet implemented; "none" = excluded from score) ------
    _stage("optical", "start" if require_opt else "skip", optical_backend)
    G_o, f_o_c = 0.0, f_o_target
    if require_opt:
        try:
            if optical_backend == "comsol":
                from optical_comsol_2d import run_optical_comsol_2d
                d = run_optical_comsol_2d(g)
                gp = gap_near_frequency(d["freqs_hz"], f_o_target)
                G_o, f_o_c = gp.normalized_gap, gp.f_center
            elif optical_backend == "none":
                raise ValueError(
                    "require_opt=True but optical_backend='none' -- set "
                    "optical_backend='comsol' (not yet implemented, see "
                    "src/optical_comsol_2d.py) or require_opt=False")
            else:
                raise ValueError(f"unknown optical_backend {optical_backend!r}")
        except Exception as e:
            _stage("optical", "fail", str(e))
            rec.update(status="optical_failed", error=str(e),
                       mechanical_gap=G_m, optical_gap=0.0,
                       score=-2.0 if not require_mech else None,
                       mechanical_center_frequency=f_m_c,
                       optical_center_frequency=0.0)
            if rec["score"] is None:
                # mechanical already succeeded -- keep it, just zero the score
                # contribution from optical rather than discarding the run.
                rec["score"] = _score(G_m, 0.0, f_m_c, f_o_target, f_m_target,
                                      f_o_target, targets,
                                      require_mech=require_mech,
                                      require_opt=False)
            if cache is not None:
                cache[cid] = rec
            return rec
        _stage("optical", "done")

    _stage("score", "start")
    score = _score(G_m, G_o, f_m_c, f_o_c, f_m_target, f_o_target, targets,
                   require_mech=require_mech, require_opt=require_opt)
    rec.update(status="success",
               mechanical_gap=float(G_m), optical_gap=float(G_o),
               mechanical_center_frequency=float(f_m_c),
               optical_center_frequency=float(f_o_c),
               require_mech=bool(require_mech), require_opt=bool(require_opt),
               score=float(score))
    _stage("score", "done")
    if cache is not None:
        cache[cid] = rec
    return rec


def _score(G_m, G_o, f_m_c, f_o_c, f_m_t, f_o_t, targets, *,
          require_mech, require_opt):
    s = targets["scoring"]
    base = 0.0
    pen = 0.0
    if require_mech:
        gmin_m = targets["mechanical"]["min_gap"]
        base += G_m
        pen += s["lambda_mech"] * max(0.0, gmin_m - G_m) ** 2
        pen += s["lambda_mech_freq"] * ((f_m_c - f_m_t) / f_m_t) ** 2
    if require_opt:
        gmin_o = targets["optical"]["min_gap"]
        base += G_o
        pen += s["lambda_optical"] * max(0.0, gmin_o - G_o) ** 2
        pen += s["lambda_optical_freq"] * ((f_o_c - f_o_t) / f_o_t) ** 2
    return base - pen


if __name__ == "__main__":
    # Smoke test with the surrogate-free path (feasibility only, no COMSOL):
    # exercises geometry + feasibility + score wiring without a solver.
    from geometry2d import u_to_geometry
    g = u_to_geometry([0.5, 0.5, 0.5, 0.3, 0.3, 0.5])
    ok, reasons = check_feasibility(g)
    print("geometry (nm):", {k: round(v * 1e9, 1) for k, v in g.as_dict().items()})
    print("feasible:", ok, reasons)
    print("(run acoustic_comsol_2d.py directly, with a COMSOL server up, "
          "for an actual band-structure solve)")
