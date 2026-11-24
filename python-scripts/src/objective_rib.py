"""Objective for the rib cell: geometry -> (optical gap, mechanical gap, score).

Rib counterpart of objective.evaluate_candidate, which dispatches here for
celltype="rib". Kept in its own module so the hole path is untouched.

Template: comsol/omc_unitcell_rib.mph (comsol/build_omc_unitcell_rib.m).
ANISOTROPIC cubic diamond (the hole template is isotropic), y >= 0 half model
with one study per y parity. The crystal is rotated by rxtal_deg about the
slab normal (template parameter rxt; mechanics only): 0 = beam along <100>
(default), 45 = beam along <110>. Only multiples of 45 deg keep the y mirror
of the half model (geometry.resolve_rxtal); 45 needs a template built with
rxt (comsol_client.require_rxtal_template, checked before any solve).

Gates before any solver, each saved as its own status:
  "infeasible"    the template cannot build the cell (geometry.check_feasibility)
  "unfabricable"  a fabrication limit fails (geometry.check_fabricability)

Mechanical gap, selected by mech_gap_mode:
  "complete"  no mode of either y parity (studies 'mech sym' + 'mech antisym',
              two solves per k)                             -- default
  "sym"       y-even modes only  (study 'mech sym')
  "antisym"   y-odd modes only   (study 'mech antisym')
The gap of every solved sector is recorded; mech_gap_mode picks the one that
is scored. Gaps are found by frequency interval (bandgap.interval_gap), since
each study returns the n_bands modes nearest its 8 GHz shift.

Optical: MPB and the TMM surrogate model the hole cell only and are refused.
  "none"    skip optics; the score uses the mechanical gap alone
  "comsol"  guided-mode gap of study optical_study (default 'opt TE': PEC on
            the y = 0 plane, i.e. Ey even in y). Every guided gap is a
            candidate (comsol_optics.guided_interval_gaps); the one with the
            largest gap * exp(-((f_c - f_target)/optical_tol)^2) is scored
            (score version 2; there is no window on the gap centre)
  "tidy3d"  the same gap choice on Tidy3D FDTD resonances (cloud, COSTS
            FlexCredits; needs the `tidy3d:` config section). The k sweep
            takes n_k_opt simulations per cell unless tidy3d.n_k says
            otherwise; n_bands_opt is not used. See optical_tidy3d.py for the
            geometry / polarization mapping and its limits.
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np

from geometry import (u_to_geometry, check_feasibility, check_fabricability,
                      fab_failure_text, resolve_xsect, diamond_index,
                      z_mirror_symmetric, wedge_depths, resolve_rxtal)
from bandgap import interval_gap, trusted_window

C0 = 299_792_458.0

MECH_STUDIES = {"sym": "mech sym", "antisym": "mech antisym"}
MECH_GAP_MODES = {"sym": ("sym",), "antisym": ("antisym",),
                  "complete": ("sym", "antisym")}
OPTICAL_BACKENDS = ("none", "comsol", "tidy3d")

# Scoring version of a rib record (key `score_version`).
#   1  optical gap centred within +/- optical_rel_tol of the target; score
#      min(G_o, G_m) minus quadratic size and centre-frequency penalties.
#      Also every record without optics (their score never changed).
#   2  no centre window; the optical gap enters as
#      G_o * exp(-((f_c - f_target)/optical_tol)^2)   (see _score_rib).
# Results of different versions (or tolerances) are not comparable: the
# version and the tolerance are in the id, and warm starts filter on them.
SCORE_VERSION = 2
DEFAULT_OPTICAL_TOL_THZ = 2.0
# Where the optical k sweep starts: at the light line of
# f_target*(1 - RIB_K_GRID_REL) (comsol_optics.rib_k_grid). This only sets the
# sampled range (and through it the trusted range of the gap search); it is
# not a rule on the gap centre. It was the default of the old optical_rel_tol
# and stays in the id tuple so that records without optics keep their ids.
RIB_K_GRID_REL = 0.30


def score_version(optical_backend) -> int:
    """Scoring version of a rib evaluation with this optical backend."""
    return 1 if optical_backend in (None, "none") else SCORE_VERSION
OPTICAL_STUDIES = ("opt TE", "opt TM")
MECH_BACKENDS = ("comsol", "surrogate_stub")


def _rib_id(u, settings):
    key = json.dumps([round(float(x), 6) for x in u] + ["rib"] + list(settings))
    return "cand_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def _num(x):
    """float, or None for NaN (so the record stays valid JSON)."""
    x = float(x)
    return x if x == x else None


def _bands_to_list(freqs):
    return [[_num(x) for x in row] for row in np.asarray(freqs, dtype=float)]


def _gap_to_dict(gp):
    return dict(found=bool(gp.found), normalized_gap=float(gp.normalized_gap),
                f_lower=float(gp.f_lower), f_upper=float(gp.f_upper),
                f_center=float(gp.f_center))


def _solve_optical_bands(g, study, k_frac, n_bands, f_target, xsect="rect",
                         info=None):
    """COMSOL sweep of one optical study of the rib template: [n_k, n_bands] Hz.

    The template follows xsect. For "rect" only the thickness-parity family of
    the study is returned (comsol_optics.study_z_family: even for 'opt TE');
    for "wedge" all modes. The diamond index of materials.yaml is set on the
    model by the sweep. `info` (a dict) receives, per k, every returned mode
    and its parity fraction.
    """
    from comsol_client import (get_model, set_geometry, template_path,
                               THICKNESS_PARAM)
    from comsol_optics import opt_sweep_at, study_z_family
    from comsol_client import require_rib_air_disk
    template = require_rib_air_disk(template_path("rib", xsect, must_exist=True))
    model = get_model(template)
    set_geometry(model, g)
    return opt_sweep_at(model, g, study, k_frac, n_bands, f_target=f_target,
                        z_family=study_z_family(study, xsect),
                        thickness=THICKNESS_PARAM["rib"], info=info)[1]


def _check_rxtal_template(xsect, rxtal):
    """Raise if the rib template of `xsect` cannot take the rotation `rxtal`."""
    from comsol_client import require_rxtal_template, template_path
    require_rxtal_template(template_path("rib", xsect, must_exist=True), rxtal)


def evaluate_rib_candidate(u, *, bounds, targets, optical_backend="none",
                           mech_backend="comsol", cache=None, on_stage=None,
                           mech_gap_mode="complete", n_k=15, n_bands_mech=20,
                           optical_study="opt TE", n_k_opt=13, n_bands_opt=10,
                           optical_tol_thz=DEFAULT_OPTICAL_TOL_THZ,
                           tidy3d_options=None, xsect=None, rxtal_deg=0.0,
                           verbose=True, optical_rel_tol=None):
    """Evaluate one normalized rib design vector u. Returns a result record.

    The record has every key of the hole record plus: celltype, the options
    above, fab_failures (status "unfabricable"), and after the solves
    mech_gaps / mech_trusted_window_Hz / _k_m / _freqs_m and
    optical_gap_detail / optical_min_resolvable_Hz / _k_o / _freqs_o.
    optical_gap is None when optics is skipped.

    optical_tol_thz: width [THz] of the Gaussian that weights the optical gap
    by its distance from the target frequency (1550 nm, targets.yaml):
        optical_freq_penalty = exp(-((f_c - f_target) / optical_tol)^2)
        optical_gap_effective = optical_gap * optical_freq_penalty
    There is no window on the gap centre any more: every guided gap in the
    trusted range that is wider than the resolvable step is a candidate
    (optical_gap_candidates), and the one with the largest effective gap is
    taken (comsol_optics.select_gap_by_penalized_size). The score uses the
    effective gap (see _score_rib). `optical_gap` stays the physical
    normalized gap of the chosen interval (= optical_gap_unpenalized).
    optical_rel_tol is DEPRECATED and ignored (a FutureWarning is issued), so
    older configs keep running.

    tidy3d_options: the `tidy3d:` config section (optical_backend="tidy3d").
    Such records add optical_engine and optical_tidy3d_detail; a solve that
    was not allowed to run has status "optical_not_run" and
    optical_tidy3d_stop (see optical_tidy3d.Tidy3DNotRun). verbose=False
    silences Tidy3D's web progress and the "info" sub-steps it sends to
    on_stage (all of it goes to stderr).

    xsect: "rect" (default) or "wedge", see objective.evaluate_candidate. It
    selects the template of optics and mechanics (omc_unitcell_rib.mph /
    omc_unitcell_rib_wedge.mph) and the Tidy3D geometry, sets the optical
    thickness-parity rule, and is part of the record and of the id.

    rxtal_deg: crystal rotation about the slab normal in degrees, 0 (beam
    along <100>, default) or 45 (<110>); any multiple of 45 is reduced mod 90
    (geometry.resolve_rxtal), anything else raises ValueError. Mechanics only
    (template parameter rxt, set by run_mechanical_comsol). Recorded as
    rxtal_deg, and part of the id only when nonzero, so the ids of unrotated
    cells are unchanged. A nonzero rotation with the COMSOL mechanics needs a
    template built with rxt; that is checked before any solve.
    """
    def _stage(name, event, info=None):
        if on_stage is not None:
            on_stage(name, event, info)

    if optical_backend not in OPTICAL_BACKENDS:
        raise ValueError(
            f"optical_backend {optical_backend!r} does not model the rib cell; "
            f"use one of {OPTICAL_BACKENDS}")
    if optical_study not in OPTICAL_STUDIES:
        raise ValueError(f"unknown optical_study {optical_study!r}; "
                         f"use one of {OPTICAL_STUDIES}")
    if mech_backend not in MECH_BACKENDS:
        raise ValueError(f"unknown mech_backend {mech_backend!r}; "
                         f"use one of {MECH_BACKENDS}")
    if mech_gap_mode not in MECH_GAP_MODES:
        raise ValueError(f"unknown mech_gap_mode {mech_gap_mode!r}; "
                         f"use one of {sorted(MECH_GAP_MODES)}")
    n_k, n_bands_mech = int(n_k), int(n_bands_mech)
    n_k_opt, n_bands_opt = int(n_k_opt), int(n_bands_opt)
    if optical_rel_tol is not None:
        import warnings
        warnings.warn(
            "rib option `optical_rel_tol` is deprecated and ignored: the gap "
            "centre is no longer restricted to a window. Use `optical_tol_thz` "
            "(Gaussian weight exp(-((f_c - f_target)/optical_tol)^2), default "
            f"{DEFAULT_OPTICAL_TOL_THZ} THz) and remove optical_rel_tol from "
            "the config.", FutureWarning, stacklevel=2)
    optical_tol_thz = float(optical_tol_thz)
    if not optical_tol_thz > 0:
        raise ValueError(f"optical_tol_thz must be > 0, got {optical_tol_thz}")
    optical_tol_hz = optical_tol_thz * 1e12
    version = score_version(optical_backend)

    xsect = resolve_xsect(xsect)
    rxtal = resolve_rxtal(rxtal_deg)
    if rxtal != 0.0 and mech_backend == "comsol":
        _check_rxtal_template(xsect, rxtal)
    z_family = None
    if optical_backend == "comsol":
        from comsol_optics import study_z_family
        z_family = study_z_family(optical_study, xsect)
    id_key = (optical_backend, mech_backend, mech_gap_mode, n_k,
              n_bands_mech, optical_study, n_k_opt, n_bands_opt,
              RIB_K_GRID_REL)
    tidy3d_tag = None
    if optical_backend == "tidy3d":
        # a Tidy3D record never replaces a COMSOL one, nor an x_only one a
        # k sweep, nor a wedge one a rectangular one.
        from optical_tidy3d import id_settings
        tidy3d_tag = id_settings(tidy3d_options, "rib", n_k_opt, xsect=xsect)
    # Cross-section, diamond index and parity family (objective.physics_id_tag).
    # A rectangular cell without optics keeps the id it always had.
    from objective import physics_id_tag
    tag = physics_id_tag(optical_backend, xsect, z_family, tidy3d_tag)
    if tag:
        id_key += (tag,)
    if rxtal != 0.0:     # unrotated cells keep the id they always had
        id_key += (f"rxtal={rxtal:g}",)
    if version != 1:     # optics scored: version and tolerance define the score
        id_key += (f"score=v{version}|otol={optical_tol_thz:g}THz",)
    cid = _rib_id(u, id_key)
    if cache is not None and cid in cache:
        return cache[cid]

    def _finish(rec):
        if cache is not None:
            cache[cid] = rec
        return rec

    _stage("feasibility", "start")
    g = u_to_geometry(u, bounds, celltype="rib")
    rec = dict(id=cid, u=[float(x) for x in u], params=g.as_dict(),
               celltype="rib",
               optical_backend=optical_backend, mech_backend=mech_backend,
               optical_engine=optical_backend, xsect=xsect,
               optical_index=(None if optical_backend == "none"
                              else diamond_index()),
               optical_z_symmetry=z_mirror_symmetric(xsect),
               status="pending", reasons=[],
               mech_gap_mode=mech_gap_mode, n_k=n_k, n_bands_mech=n_bands_mech,
               optical_study=optical_study, n_k_opt=n_k_opt,
               n_bands_opt=n_bands_opt,
               optical_rel_tol=None,    # deprecated key, kept in the schema
               optical_tol_thz=optical_tol_thz, score_version=version,
               rxtal_deg=rxtal)
    no_solve = dict(optical_gap=0.0, mechanical_gap=0.0, score=-10.0,
                    optical_center_frequency=0.0,
                    mechanical_center_frequency=0.0)

    ok, reasons = check_feasibility(g, bounds, xsect)
    if not ok:
        _stage("feasibility", "fail", reasons)
        rec.update(status="infeasible", reasons=reasons, **no_solve)
        return _finish(rec)
    if xsect == "wedge":
        rec.update(wedge_depths=wedge_depths(g))

    ok, failures = check_fabricability(g, bounds)
    if not ok:
        reasons = [fab_failure_text(f) for f in failures]
        _stage("feasibility", "fail", reasons)
        rec.update(status="unfabricable", reasons=reasons,
                   fab_failures=failures, **no_solve)
        return _finish(rec)
    _stage("feasibility", "done")

    f_o_target = C0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    f_m_target = targets["mechanical"]["target_frequency_GHz"] * 1e9
    nan = float("nan")

    # ---- optical ----
    G_o, f_o_c, pen_o = nan, nan, nan
    if optical_backend == "none":
        _stage("optical", "skip", optical_backend)
    elif optical_backend == "tidy3d":
        # Both gates above (buildable, fabricable) ran before this solver.
        _stage("optical", "start", f"{optical_backend}, {optical_study}")
        try:
            from optical_tidy3d import run_optical_tidy3d
            o = run_optical_tidy3d(
                g, celltype="rib", eval_id=cid, f_target=f_o_target,
                options=tidy3d_options, rel_tol=RIB_K_GRID_REL,
                optical_study=optical_study, n_k_default=n_k_opt,
                xsect=xsect, verbose=verbose,
                on_progress=lambda m: _stage("optical", "info", m),
                optical_tol_hz=optical_tol_hz)
            G_o, f_o_c = o["gap"].normalized_gap, o["gap"].f_center
            rec.update(o["record"])
            pen_o = o["record"].get("optical_freq_penalty")
            if pen_o is None:           # a runner without the gap choice
                from comsol_optics import freq_penalty
                pen_o = freq_penalty(f_o_c, f_o_target, optical_tol_hz)
                rec.update(optical_freq_penalty=pen_o)
            rec.update(
                optical_min_resolvable_Hz=o["record"]["optical_tidy3d_detail"]
                ["min_resolvable_Hz"],
                optical_light_line_edge_Hz=float(C0 / (2.0 * g.a)))
        except Exception as e:  # solver failed -> record, don't crash loop
            _stage("optical", "fail", str(e))
            if getattr(e, "tidy3d_stop", None):
                # Nothing was run or charged; not cached (see objective.py).
                from optical_tidy3d import not_run_record
                return not_run_record(rec, e)
            rec.update(status="optical_failed", error=str(e),
                       optical_gap=0.0, mechanical_gap=0.0, score=-5.0,
                       optical_center_frequency=0.0,
                       mechanical_center_frequency=0.0)
            return _finish(rec)
        _stage("optical", "done")
    else:
        _stage("optical", "start", f"{optical_backend}, {optical_study}")
        try:
            from comsol_optics import (guided_interval_gaps, resolvable_step_hz,
                                       rib_k_grid, select_gap_by_penalized_size)
            k_o = rib_k_grid(g.a, n_k_opt, f_o_target, RIB_K_GRID_REL)
            info = {}
            freqs_o = np.asarray(
                _solve_optical_bands(g, optical_study, k_o, n_bands_opt,
                                     f_o_target, xsect=xsect, info=info),
                dtype=float)
            # per k: every returned mode and its thickness-parity fraction
            rec.update(optical_comsol_detail=dict(
                z_family=z_family, freqs_all_hz=info.get("freqs_all_hz"),
                z_even_fraction=info.get("z_even_fraction")))
            # every guided gap, then the one with the largest gap * Gaussian
            gp, pen_o, candidates = select_gap_by_penalized_size(
                guided_interval_gaps(k_o, freqs_o, g.a), f_o_target,
                optical_tol_hz)
            G_o, f_o_c = gp.normalized_gap, gp.f_center
            rec.update(
                optical_gap_detail=_gap_to_dict(gp),
                optical_freq_penalty=float(pen_o),
                optical_gap_candidates=candidates,
                optical_min_resolvable_Hz=float(resolvable_step_hz(k_o, g.a)),
                optical_light_line_edge_Hz=float(C0 / (2.0 * g.a)),
                _k_o=[float(x) for x in k_o],
                _freqs_o=_bands_to_list(freqs_o))
        except Exception as e:  # solver failed -> record, don't crash loop
            _stage("optical", "fail", str(e))
            rec.update(status="optical_failed", error=str(e),
                       optical_gap=0.0, mechanical_gap=0.0, score=-5.0,
                       optical_center_frequency=0.0,
                       mechanical_center_frequency=0.0)
            return _finish(rec)
        _stage("optical", "done")

    # ---- mechanical ----
    G_m, f_m_c = nan, nan
    sectors = MECH_GAP_MODES[mech_gap_mode]
    if mech_backend == "surrogate_stub":
        _stage("mechanical", "skip", mech_backend)
    else:
        _stage("mechanical", "start", f"{mech_backend}, {mech_gap_mode}")
        try:
            from acoustic_comsol import run_mechanical_comsol
            from comsol_client import template_path
            bands, k_frac = {}, None
            for sector in sectors:
                d = run_mechanical_comsol(
                    g, n_k=n_k, n_bands=n_bands_mech,
                    template=template_path("rib", xsect, must_exist=True),
                    study=MECH_STUDIES[sector], rxtal_deg=rxtal)
                bands[sector] = np.asarray(d["freqs_hz"], dtype=float)
                k_frac = d["k_frac"]
            gaps = {sector: interval_gap([bands[sector]], f_m_target)
                    for sector in sectors}
            if mech_gap_mode == "complete":
                gaps["complete"] = interval_gap(
                    [bands[sector] for sector in sectors], f_m_target)
            window = trusted_window([bands[sector] for sector in sectors])
            gp = gaps[mech_gap_mode]
            G_m, f_m_c = gp.normalized_gap, gp.f_center
            rec.update(
                mech_gaps={name: _gap_to_dict(x) for name, x in gaps.items()},
                mech_trusted_window_Hz=[_num(window[0]), _num(window[1])],
                _k_m=[float(x) for x in k_frac],
                _freqs_m={sector: _bands_to_list(b) for sector, b in bands.items()})
        except Exception as e:
            _stage("mechanical", "fail", str(e))
            rec.update(status="mech_failed", error=str(e),
                       optical_gap=(_num(G_o) or 0.0), mechanical_gap=0.0,
                       score=-2.0,
                       optical_center_frequency=(_num(f_o_c) or 0.0),
                       mechanical_center_frequency=0.0)
            return _finish(rec)
        _stage("mechanical", "done")

    _stage("score", "start")
    solved = (G_o == G_o) or (G_m == G_m)
    # the optical gap enters the score weighted by its distance from the target
    G_o_eff = G_o * float(pen_o) if G_o == G_o else nan
    score = _score_rib(G_o_eff, G_m, f_m_c, f_m_target, targets)
    rec.update(status="success" if solved else "no_solver",
               optical_gap=_num(G_o), mechanical_gap=_num(G_m),
               optical_gap_unpenalized=_num(G_o),
               optical_gap_effective=_num(G_o_eff),
               optical_freq_penalty=(None if G_o != G_o else float(pen_o)),
               optical_center_frequency=_num(f_o_c),
               mechanical_center_frequency=_num(f_m_c),
               score=float(score))
    _stage("score", "done")
    return _finish(rec)


def _score_rib(G_o_eff, G_m, f_m_c, f_m_t, targets):
    """Rib score (score version 2), with either gap allowed to be absent (NaN).

    G_o_eff = G_o * exp(-((f_o,c - f_o,target) / optical_tol)^2) is the
    optical gap already weighted by its distance from the target. It replaces
    G_o everywhere, and the separate optical centre-frequency penalty
    (lambda_optical_freq) is gone:

      S = min(G_o_eff, G_m)
          - lo  * max(0, gmin_o - G_o_eff)^2
          - lm  * max(0, gmin_m - G_m)^2
          - lmf * ((f_m,c - f_m,target) / f_m,target)^2

    The mechanical terms are those of objective._score. One gap present: that
    gap minus its own penalties. None: 0.0. An optical gap that was searched
    for but not found has G_o_eff = 0 and takes the full size penalty; so
    does, nearly, a gap several tolerances away from the target. Since
    S <= min(G_o_eff, G_m), a mechanical gap cannot make up for an off-target
    optical one.
    """
    s = targets["scoring"]
    gaps, pen = [], 0.0
    if G_o_eff == G_o_eff:
        gaps.append(G_o_eff)
        pen += (s["lambda_optical"]
                * max(0.0, targets["optical"]["min_gap"] - G_o_eff) ** 2)
    if G_m == G_m:
        gaps.append(G_m)
        pen += (s["lambda_mech"]
                * max(0.0, targets["mechanical"]["min_gap"] - G_m) ** 2
                + s["lambda_mech_freq"] * ((f_m_c - f_m_t) / f_m_t) ** 2)
    if not gaps:
        return 0.0
    return min(gaps) - pen
