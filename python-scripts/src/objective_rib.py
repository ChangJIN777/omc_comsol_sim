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
            largest optical score contribution (optical_terms: gap minus its
            size, detuning and light-cone penalties) is scored. There is no
            window on the gap centre. Score: see rib_score (version 4).
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
#   2  no centre window; min(G_o * exp(-((f_c - f_t)/optical_tol)^2), G_m)
#      minus the size penalties and the mechanical centre penalty.
#   3  G_o + G_m minus size, detuning, light-cone, mechanical-window and
#      mechanical-centre penalties, the detuning penalty quadratic:
#      lof * ((f_o - f_t)/f_t)^2. Mechanics-only records are version 3 too:
#      the window penalty changed their score as well.
#   4  as 3 without the mechanical-centre (8 GHz) penalty, and with a
#      flat-top detuning penalty lof * [1 - 1/(1 + (|f_o - f_t|/f_tol)^8)]
#      (rib_score).
# Scores of different versions are not comparable: the version is in the id,
# and warm starts use records of the same version and weights only.
SCORE_VERSION = 4
# Exponent of the flat-top detuning penalty (recorded in score_weights as
# optical_freq_power, so another form is another objective).
OPTICAL_FREQ_POWER = 8.0
# Defaults of the rib-only weights (configs/targets.yaml, block `scoring_rib`).
# lambda_optical_freq and lambda_light_cone are LAMBDA_OPT_F and
# LAMBDA_LIGHTCONE of scripts/run_opt_comsol.py (the hole loop); the mechanical
# window is its MECH_F_MIN / MECH_F_MAX with a weight "like the light cone".
RIB_SCORING_DEFAULTS = dict(lambda_optical_freq=9.0, optical_freq_tol_THz=5.0,
                            lambda_light_cone=2000.0,
                            mech_window_GHz=(5.0, 10.0),
                            lambda_mech_window=2000.0)
# Where the optical k sweep starts: at the light line of
# f_target*(1 - RIB_K_GRID_REL) (comsol_optics.rib_k_grid). This only sets the
# sampled range (and through it the trusted range of the gap search); it is
# not a rule on the gap centre. It was the default of the old optical_rel_tol
# and stays in the id tuple.
RIB_K_GRID_REL = 0.30


def score_version(optical_backend=None) -> int:
    """Scoring version of a rib evaluation (the same for every backend)."""
    return SCORE_VERSION


def rib_score_weights(targets, version=None) -> dict:
    """Every number the rib score uses, from targets.yaml.

    Shared with the hole cell (unchanged keys): optical.min_gap,
    mechanical.min_gap and scoring.lambda_optical / lambda_mech. Rib only,
    from the optional block `scoring_rib` (RIB_SCORING_DEFAULTS when absent):
    lambda_optical_freq, optical_freq_tol_THz, lambda_light_cone,
    mech_window_GHz, lambda_mech_window. The hole cell's
    scoring.lambda_optical_freq and scoring.lambda_mech_freq are NOT used.

    version (default: the current one, 4) selects the set: version 3 has
    scoring.lambda_mech_freq (its 8 GHz centre term) and no
    optical_freq_tol_THz / optical_freq_power, which is how optical_terms and
    mech_terms tell the two forms apart (used to re-score older records).
    """
    version = SCORE_VERSION if version is None else int(version)
    s = targets["scoring"]
    rib = dict(RIB_SCORING_DEFAULTS, **(targets.get("scoring_rib") or {}))
    unknown = sorted(set(rib) - set(RIB_SCORING_DEFAULTS))
    if unknown:
        raise ValueError(f"unknown key(s) under `scoring_rib`: {unknown}; "
                         f"expected a subset of {sorted(RIB_SCORING_DEFAULTS)}")
    lo, hi = (float(x) for x in rib["mech_window_GHz"])
    if not 0 < lo < hi:
        raise ValueError(f"scoring_rib.mech_window_GHz must be 0 < low < high, "
                         f"got {rib['mech_window_GHz']}")
    w = dict(gmin_o=float(targets["optical"]["min_gap"]),
             gmin_m=float(targets["mechanical"]["min_gap"]),
             lambda_optical=float(s["lambda_optical"]),
             lambda_mech=float(s["lambda_mech"]),
             lambda_optical_freq=float(rib["lambda_optical_freq"]),
             lambda_light_cone=float(rib["lambda_light_cone"]),
             lambda_mech_window=float(rib["lambda_mech_window"]),
             mech_window_GHz=[lo, hi])
    if version >= 4:
        tol = float(rib["optical_freq_tol_THz"])
        if not tol > 0:
            raise ValueError(f"scoring_rib.optical_freq_tol_THz must be > 0, got {tol}")
        w.update(optical_freq_tol_THz=tol, optical_freq_power=OPTICAL_FREQ_POWER)
    else:
        w["lambda_mech_freq"] = float(s["lambda_mech_freq"])
    for key, val in w.items():
        if key.startswith("lambda") and not val >= 0:
            raise ValueError(f"{key} must be >= 0, got {val}")
    return w


def optical_terms(G_o, f_c, f_target, a, w) -> dict:
    """Optical part of the rib score for a gap G_o centred at f_c [Hz].

      pen_size_o      lo  * max(0, gmin_o - G_o)^2
      pen_detune_o    lof * [1 - 1/(1 + (|f_c - f_target| / f_tol)^8)]
      pen_lightcone_o llc * max(0, (f_c - c/2a) / (c/2a))^2    (a = period [m])
      value           G_o - the three penalties   (what the gap choice maximizes)
    The detuning penalty is flat-topped: ~0 within +/- f_tol
    (optical_freq_tol_THz) of the target, lof/2 at f_tol, lof far away
    (score version 4). Weights without optical_freq_tol_THz are those of a
    version-3 record, scored as it was: lof * ((f_c - f_target)/f_target)^2.
    The light-cone term is that of score_result in scripts/run_opt_comsol.py.
    No gap (G_o = 0): pass f_c = f_target, as that script does, so there is
    no detuning penalty, and a light-cone penalty only if the target itself
    lies above the zone-edge light line.
    """
    light_line = C0 / (2.0 * a)
    pen_size = w["lambda_optical"] * max(0.0, w["gmin_o"] - G_o) ** 2
    if w.get("optical_freq_tol_THz") is not None:
        r = (abs(f_c - f_target) / (w["optical_freq_tol_THz"] * 1e12)) ** w.get(
            "optical_freq_power", OPTICAL_FREQ_POWER)
        pen_detune = w["lambda_optical_freq"] * r / (1.0 + r)
    else:                                                  # score version 3
        pen_detune = w["lambda_optical_freq"] * ((f_c - f_target) / f_target) ** 2
    pen_lc = w["lambda_light_cone"] * max(0.0, (f_c - light_line) / light_line) ** 2
    return dict(pen_size_o=float(pen_size), pen_detune_o=float(pen_detune),
                pen_lightcone_o=float(pen_lc),
                value=float(G_o - pen_size - pen_detune - pen_lc))


def mech_terms(G_m, f_c, found, w, f_target=None) -> dict:
    """Mechanical part of the rib score for a gap G_m centred at f_c [Hz].

      pen_size_m    lm  * max(0, gmin_m - G_m)^2
      pen_window_m  lmw * max(0, (f_lo - f_c)/f_lo)^2   below the window
                    lmw * max(0, (f_c - f_hi)/f_hi)^2   above it
                    (0 inside [f_lo, f_hi] = mech_window_GHz, and 0 when no
                    gap was found: there is no centre to place)
      pen_center_m  0.0: score version 4 has no 8 GHz centre term. (Weights
                    of a version-3 record carry lambda_mech_freq, and then it
                    is lmf * ((f_c - f_target)/f_target)^2, as that run had.)
    """
    f_lo, f_hi = (x * 1e9 for x in w["mech_window_GHz"])
    pen_size = w["lambda_mech"] * max(0.0, w["gmin_m"] - G_m) ** 2
    pen_window = 0.0
    if found:
        pen_window = w["lambda_mech_window"] * (
            max(0.0, (f_lo - f_c) / f_lo) ** 2 + max(0.0, (f_c - f_hi) / f_hi) ** 2)
    pen_center = 0.0
    if w.get("lambda_mech_freq") is not None:              # score version 3
        pen_center = w["lambda_mech_freq"] * ((f_c - f_target) / f_target) ** 2
    return dict(pen_size_m=float(pen_size), pen_window_m=float(pen_window),
                pen_center_m=float(pen_center),
                value=float(G_m - pen_size - pen_window - pen_center))
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
                           tidy3d_options=None, xsect=None, rxtal_deg=0.0,
                           verbose=True, optical_rel_tol=None,
                           optical_tol_thz=None):
    """Evaluate one normalized rib design vector u. Returns a result record.

    The record has every key of the hole record plus: celltype, the options
    above, fab_failures (status "unfabricable"), and after the solves
    mech_gaps / mech_trusted_window_Hz / _k_m / _freqs_m and
    optical_gap_detail / optical_min_resolvable_Hz / _k_o / _freqs_o.
    optical_gap is None when optics is skipped.

    Optical gap and score (version 4, see rib_score): there is no window on
    the gap centre. Every guided gap in the trusted range that is wider than
    the resolvable step is a candidate (optical_gap_candidates); the one with
    the largest optical contribution to the score (optical_terms) is taken.
    The record carries the penalty terms pen_size_o, pen_detune_o,
    pen_lightcone_o, pen_size_m, pen_window_m, pen_center_m (always 0 from
    version 4 on, kept for the schema) and the weights used (score_weights). optical_rel_tol and optical_tol_thz are DEPRECATED
    and ignored (a FutureWarning is issued), so older configs keep running.

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
    for name, val, was in (("optical_rel_tol", optical_rel_tol,
                            "a +/- window on the gap centre (score version 1)"),
                           ("optical_tol_thz", optical_tol_thz,
                            "the width of a Gaussian weight (score version 2)")):
        if val is not None:
            import warnings
            warnings.warn(
                f"rib option `{name}` is deprecated and ignored: it was {was}. "
                "Score version 3 penalizes the detuning of the gap centre "
                "quadratically (targets.yaml: scoring_rib.lambda_optical_freq). "
                f"Remove {name} from the config.", FutureWarning, stacklevel=2)
    version = score_version(optical_backend)
    weights = rib_score_weights(targets)

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
    # the scoring version: a version-3 record never replaces a version-1 / 2
    # one of the same cell (those ids had no tag, or "score=v2|otol=...")
    id_key += (f"score=v{version}",)
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
               optical_rel_tol=None,    # deprecated keys, kept in the schema
               optical_tol_thz=None,
               score_version=version, score_weights=weights,
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
    G_o, f_o_c, opt_found = nan, nan, False

    def gap_value(gp):                  # what the optical gap choice maximizes
        return optical_terms(gp.normalized_gap, gp.f_center, f_o_target, g.a,
                             weights)["value"]
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
                gap_value=gap_value, score_version=version)
            G_o, f_o_c = o["gap"].normalized_gap, o["gap"].f_center
            opt_found = bool(o["gap"].found)
            rec.update(o["record"])
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
                                       rib_k_grid, select_gap_by_value)
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
            # every guided gap, then the one that contributes most to the score
            gp, candidates = select_gap_by_value(
                guided_interval_gaps(k_o, freqs_o, g.a), gap_value, f_o_target)
            G_o, f_o_c, opt_found = gp.normalized_gap, gp.f_center, bool(gp.found)
            rec.update(
                optical_gap_detail=_gap_to_dict(gp),
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
    G_m, f_m_c, mech_found = nan, nan, False
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
            G_m, f_m_c, mech_found = gp.normalized_gap, gp.f_center, bool(gp.found)
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
    score, terms = rib_score(G_o, f_o_c, opt_found, G_m, f_m_c, mech_found,
                             a=g.a, f_o_target=f_o_target, f_m_target=f_m_target,
                             weights=weights)
    rec.update(status="success" if solved else "no_solver",
               optical_gap=_num(G_o), mechanical_gap=_num(G_m),
               optical_gap_unpenalized=_num(G_o),
               # kept keys of version 2, in their version-3 meaning: there is
               # no multiplicative weight any more, so the "effective" gap is
               # the gap, and the frequency penalty is the detuning term that
               # is SUBTRACTED from the score (= pen_detune_o)
               optical_gap_effective=_num(G_o),
               optical_freq_penalty=terms["pen_detune_o"],
               optical_center_frequency=_num(f_o_c),
               mechanical_center_frequency=_num(f_m_c),
               score=float(score), **terms)
    _stage("score", "done")
    return _finish(rec)


def rib_score(G_o, f_o_c, opt_found, G_m, f_m_c, mech_found, *, a, f_o_target,
              f_m_target, weights):
    """Rib score, version 4. Returns (score, terms).

      S = G_o + G_m
          - lo  * max(0, gmin_o - G_o)^2                        pen_size_o
          - lof * [1 - 1/(1 + (|f_o,c - f_o,t| / f_tol)^8)]     pen_detune_o
          - llc * max(0, (f_o,c - c/2a) / (c/2a))^2             pen_lightcone_o
          - lm  * max(0, gmin_m - G_m)^2                        pen_size_m
          - lmw * (max(0, (f_lo - f_m,c)/f_lo)^2
                   + max(0, (f_m,c - f_hi)/f_hi)^2)             pen_window_m
    There is no penalty on the distance of the mechanical gap centre from
    8 GHz (pen_center_m is 0; f_m_target is not used), only the window.

    G_o, G_m are the normalized gaps (f_upper - f_lower)/f_center; NaN = that
    physics was not solved, and then neither the gap nor its penalties enter
    (optical: none, or mech: surrogate_stub). Both NaN: 0.0. The weights come
    from rib_score_weights. A gap that was searched for but not found
    (found = False, gap 0):
      optical     size penalty in full; the centre is taken as the target, so
                  no detuning penalty, and a light-cone penalty only if the
                  target lies above the zone-edge light line c/2a (as
                  score_result in scripts/run_opt_comsol.py);
      mechanical  size penalty in full, no window penalty.
    terms holds the six penalties (None for physics that was not solved).
    With the weights of a version-3 record (rib_score_weights(t, version=3))
    this returns the version-3 score.
    Unlike versions 1 and 2 the two gaps ADD: a wide mechanical gap raises
    the score of a cell with a poor optical gap, and the reverse.
    """
    terms = dict(pen_size_o=None, pen_detune_o=None, pen_lightcone_o=None,
                 pen_size_m=None, pen_window_m=None, pen_center_m=None)
    score, any_gap = 0.0, False
    if G_o == G_o:
        o = optical_terms(G_o, f_o_c if opt_found else f_o_target, f_o_target,
                          a, weights)
        score += o.pop("value")
        terms.update(o)
        any_gap = True
    if G_m == G_m:
        m = mech_terms(G_m, f_m_c, mech_found, weights, f_m_target)
        score += m.pop("value")
        terms.update(m)
        any_gap = True
    return (float(score) if any_gap else 0.0), terms
