#!/usr/bin/env python3
"""Co-optimize optical + mechanical bandgaps using COMSOL.

Optical target : TE bandgap >= 20%, center near 1550 nm
Mechanical target: breathing (fy > 0.5) gap >= 20%, center in 5-10 GHz

Runs all COMSOL solves through the already-open COMSOL model (cached by
comsol_client).  Uses Optuna (TPE) if installed, else Halton random search.

Usage:
    python scripts/run_opt_comsol.py --n-iter 20
    python scripts/run_opt_comsol.py --n-init 8 --n-iter 30 --n-k 9
    python scripts/run_opt_comsol.py --n-iter 5 --n-k 5   # quick test
"""
import argparse
import hashlib
import os
import re
import sys
import json
import subprocess
import time
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from geometry import u_to_geometry, check_feasibility, load_bounds, Geometry
from bandgap import largest_gap, all_gaps, Gap
from comsol_client import get_model
from database import save_result

C0 = 299_792_458.0
_TEMPLATE    = os.path.join(os.path.dirname(__file__), "..", "comsol",
                            "omc_unitcell_iso.mph")
_CHAR_SCRIPT = os.path.join(os.path.dirname(__file__), "characterize_best.py")

# ── Scoring constants ─────────────────────────────────────────────────────────
OPT_F_TARGET   = C0 / 1550e-9      # 193.4 THz
MECH_F_MIN     = 5e9               # lower mechanical center (informational only)
MECH_F_MAX     = 10e9              # upper mechanical center (informational only)
G_MIN_OPT      = 0.15               # minimum acceptable optical gap (relaxed from 0.20, 2026-07-01)
G_MIN_MECH     = 0.20               # minimum acceptable mechanical gap
LAMBDA_OPT     = 20.0              # penalty: optical gap < G_MIN_OPT
LAMBDA_MECH    = 5.0               # soft penalty: mech gap < G_MIN_MECH (reduced; G_m secondary)
LAMBDA_OPT_F   = 9.0               # penalty: optical center frequency (3x stronger, 2026-07-01)
LAMBDA_MECH_F  = 0.0               # mechanical frequency penalty DISABLED (f_m unconstrained)
# Thinness/narrowness are REWARD terms (2026-07-02), not threshold penalties:
# continuously additive across the full bounds range like G_o/G_m, so there's
# no value of t or w below which it "stops getting better" — smaller is
# always strictly preferred, all the way to the bounds minimum. Matches
# configs/bounds.yaml's t/w ranges; kept as constants here (not re-read from
# YAML per call) for speed, consistent with the other score constants.
# Must match configs/bounds.yaml's t/w min/max — updated 2026-07-02 when the
# floors were lowered (t: 220->150nm, w: 450->300nm) after the optical-only
# campaign kept pinning designs against the old floors.
T_MIN_BOUND, T_MAX_BOUND = 150e-9, 500e-9
W_MIN_BOUND, W_MAX_BOUND = 300e-9, 900e-9
LAMBDA_THICK   = 0.2                # weight on thinness reward (2026-07-02, user-chosen)
LAMBDA_WIDTH   = 0.2                # weight on narrowness reward (2026-07-02, user-chosen)
LAMBDA_LIGHTCONE = 2000.0          # penalty: gap center above the light line (2026-07-02).
                                    # Quadratic-in-fractional-excess is inherently weak near
                                    # the boundary (a 3.7% overshoot gives excess^2=0.00137),
                                    # so LAMBDA must be large for marginal violations (the
                                    # kind an optimizer pressing on the constraint actually
                                    # finds) to be genuinely severe. At LAMBDA=2000: 1% excess
                                    # -> pen=0.2, 3.7% -> pen=2.7, 5% -> pen=5.0 (score scale
                                    # here tops out around ~1.2).
FY_THRESH      = 0.5               # breathing mode classification threshold


# ── COMSOL helpers ────────────────────────────────────────────────────────────

def _study_dataset(model, study):
    matches = [d for d in model.datasets() if d.startswith(study + "//")]
    return matches[0] if matches else None


def _eval_list(model, dset, exprs, inner_id):
    result = model.evaluate(exprs, dataset=dset, inner=[inner_id])
    return np.column_stack([np.real(np.asarray(r)).ravel() for r in result])


def _fy_from_nodes(model, dset, inner_id):
    """Fractional u_y energy for one eigenmode (breathing ≈ 1)."""
    try:
        data = _eval_list(model, dset, ["real(u)", "real(v)", "real(w)"], inner_id)
        valid = ~np.isnan(data[:, 1])
        if valid.sum() == 0:
            return np.nan
        u2 = np.sum(data[valid, 0] ** 2)
        v2 = np.sum(data[valid, 1] ** 2)
        w2 = np.sum(data[valid, 2] ** 2)
        tot = u2 + v2 + w2
        return float(v2 / tot) if tot > 0 else np.nan
    except Exception:
        return np.nan


def _mech_sweep(model, g, study, n_k, n_bands, compute_fy=True):
    """Sweep Γ→X; return (k_frac, freqs_all, fy_arr) [n_k, n_bands].

    compute_fy=False: skip per-mode field evaluation (fy_arr set to 1.0),
    trusting that the study BCs (mech sym) already select breathing modes.
    This saves n_k × n_bands COMSOL field queries — ~30-50% of total time.
    """
    k_frac = np.linspace(1e-3, 1.0, n_k)
    freqs  = np.full((n_k, n_bands), np.nan)
    fy_arr = np.ones((n_k, n_bands))   # default: treat all as breathing
    for i, s in enumerate(k_frac):
        model.parameter("kF", f"{np.pi * s / g.a}[1/m]")
        model.solve(study)
        dset = _study_dataset(model, study)
        freqs_raw = np.real(np.asarray(model.evaluate("freq", dataset=dset)))
        sort_idx  = np.argsort(freqs_raw)
        n_av = min(n_bands, len(freqs_raw))
        freqs[i, :n_av] = freqs_raw[sort_idx[:n_av]]
        if compute_fy:
            inner_ids  = np.asarray(model.inner(dset)[0])
            sorted_ids = inner_ids[sort_idx[:n_av]]
            for j, iid in enumerate(sorted_ids):
                fy_arr[i, j] = _fy_from_nodes(model, dset, int(iid))
    return k_frac, freqs, fy_arr


def _study_tag(model, name):
    """Resolve a study's display name (e.g. 'opt TE') to its node tag ('std2')."""
    for tag in model.java.study().tags():
        if model.java.study(tag).label() == name:
            return tag
    raise ValueError(f"study {name!r} not found in model")


def _guided_shift_hz(kz, f_target=OPT_F_TARGET, margin=0.92):
    """Eigenvalue-search shift for one k-point: at/below the light line.

    The .mph template's "opt TE"/"opt TM" studies have a FIXED eigenfrequency
    search shift (c/1.55um, i.e. exactly the 1550 nm target) at every k. A
    shift-and-invert eigensolver returns the N eigenvalues closest to that
    shift — so whenever the fundamental (lowest-order) guided TE band's true
    frequency at a given k is far from the fixed target, it's simply never
    returned, and a higher-order band near the target gets scored instead
    (this is what produced the earlier "second-order gap" / "above light
    line" symptoms). Tracking the shift along the light line instead keeps
    the search anchored where genuinely guided modes actually live at each k,
    with a small margin below light line and never above the 1550 nm target.
    """
    f_lc = C0 * kz / (2.0 * np.pi)
    return min(f_target, f_lc * margin)


def _set_eig_search(model, study, f_shift_hz, neigs=None):
    """Point a study's Eigenfrequency feature at a new search shift (Hz),
    optionally widening how many eigenvalues it searches for. Best-effort:
    swallows exceptions so callers fall back to the model's baked-in
    shift/neigs if the Java property API is unavailable.

    Two COMSOL gotchas here (found by live introspection, not documented):
    - "neigs" is silently ignored unless "neigsactive" is explicitly "on"
      (it defaults to "off" in this template, so it was NEVER actually
      taking effect before this fix — always solving with the .mph's
      baked-in neigs=6, no matter what was passed here).
    - .set(name, <python int>) raises a Java overload-ambiguity TypeError
      (int vs bool overloads collide); must pass neigs as a string.
    """
    try:
        eig = model.java.study(_study_tag(model, study)).feature("eig")
        eig.set("shift", f"{f_shift_hz}[Hz]")
        if neigs is not None:
            eig.set("neigsactive", "on")
            eig.set("neigs", str(int(neigs)))
    except Exception:
        pass


def _opt_sweep(model, g, study, n_k, n_bands, min_THz=10.0):
    """Sweep zone-edge half of BZ for optical. Returns (k_frac, freqs) in Hz.

    Samples k ∈ [0.5, 1.0]×π/a only — the half nearest the zone edge.
    This is the range that matters for OMC Bragg reflection; modes near k=0
    are irrelevant.  With n_k points the same number of COMSOL solves is used
    but all in the relevant half — effectively giving 2× denser zone-edge
    sampling compared to a full-BZ sweep with the same n_k.

    At each k, the eigenvalue search shift tracks the light line (see
    `_guided_shift_hz`) rather than a fixed 1550 nm target, and `neigs` is
    widened so the fundamental guided band is included alongside whatever
    band actually forms the target-frequency gap.

    neigs floor is 6, not the wider margin used for plotting (`_full_opt_sweep`):
    during optimization we only need the lowest ~2 real bands near the gap
    plus room for the 1-3 near-zero scattering-boundary pseudo-modes (which
    get filtered out below by `min_THz`) — asking for more just costs solve
    time without changing the score.
    """
    k_frac = np.linspace(0.5, 1.0, n_k)   # zone-edge half only
    freqs  = np.full((n_k, n_bands), np.nan)
    neigs_req = max(n_bands, 6)
    for i, s in enumerate(k_frac):
        kz = np.pi * s / g.a
        model.parameter("kF", f"{kz}[1/m]")
        _set_eig_search(model, study, _guided_shift_hz(kz), neigs=neigs_req)
        model.solve(study)
        dset = _study_dataset(model, study)
        ev = (model.evaluate("freq", dataset=dset) if dset
              else model.evaluate("freq"))
        fs = np.sort(np.real(np.asarray(ev)))
        fs = fs[fs > min_THz * 1e12]
        n_av = min(n_bands, len(fs))
        freqs[i, :n_av] = fs[:n_av]
    return k_frac, freqs


# ── Gap detection for breathing family ───────────────────────────────────────

def breathing_gap(freqs_all, fy_arr, fy_thresh=FY_THRESH,
                  f_min=MECH_F_MIN, f_max=MECH_F_MAX):
    """Find the largest gap where AT LEAST ONE edge band is breathing-dominant.

    Only one side of the gap needs to be the mode family we actually care
    about (fy = mean breathing fraction over k >= fy_thresh) — the other
    edge can be any symmetry; we don't use it for anything besides marking
    where the gap closes. This intentionally does NOT require both edges to
    pass (an earlier version pointwise-NaN'd low-fy (k, band) entries and
    required both nanmax/nanmin to survive, which effectively demanded both
    edges be breathing AND could silently narrow a band's reported extent
    if its true frequency max fell at a low-fy k-point).

    Candidates are walked largest-gap-first (via `all_gaps`) so the first
    one clearing the fy check is the biggest breathing-adjacent gap, not
    necessarily the global largest gap.
    """
    freqs_all = np.asarray(freqs_all, dtype=float)
    fy_arr = np.asarray(fy_arr, dtype=float)
    n_bands = freqs_all.shape[1]
    fy_mean = np.nanmean(fy_arr, axis=0)   # one representative fy per band
    window = (f_min * 0.5, f_max * 2.0)
    for gp in all_gaps(freqs_all, f_window=window, nan_safe=True):
        i = gp.lower_band
        if not (0 <= i < n_bands - 1):
            continue
        if max(fy_mean[i], fy_mean[i + 1]) >= fy_thresh:
            return gp
    return Gap.empty()


# ── Optical gap detection: zone-edge gap (OMC-relevant metric) ───────────────

def optical_gap_tracked(freqs_o, f_target=OPT_F_TARGET, rel_tol=0.30):
    """Find the OMC-relevant optical bandgap from a zone-edge-half BZ sweep.

    Assumes freqs_o rows span k ∈ [0.5, 1.0]×π/a (from _opt_sweep).

    1. At zone edge (last row): evaluate SPECIFICALLY the fundamental
       (lowest-order) gap — between the two lowest surviving bands (the
       dielectric band and the first air band). This is a deliberate choice,
       not "whichever gap happens to be closest to f_target": a design only
       counts if its FUNDAMENTAL gap sits near 1550 nm. A higher-order gap
       that happens to land near the target is rejected (empty), even if
       nothing else qualifies — otherwise the optimizer can "win" by finding
       geometries whose 2nd/3rd-order gap coincidentally sits at 1550 nm
       while the true (lowest-order) guided gap is elsewhere, which isn't a
       real OMC mirror-cell design. Requires _opt_sweep/_full_opt_sweep to
       actually resolve the fundamental band (see `_guided_shift_hz`) —
       without that, this band may not even be present in freqs_o.
    2. Use gap center f_c as a tracer; track the gap width at every k-point
       provided.  All provided k-points are in the zone-edge half, so mode
       mixing near k=0 (which would falsely close the gap) is excluded.
    3. Return the minimum width across all tracked k-points.
    """
    freqs_o = np.asarray(freqs_o, dtype=float)
    n_k = freqs_o.shape[0]

    # ── Step 1: fundamental (lowest-order) gap at zone edge ────────────────────
    edge = freqs_o[-1, :]
    valid = np.sort(edge[np.isfinite(edge) & (edge > 10e12)])
    if len(valid) < 2:
        return Gap.empty()

    f_c    = 0.5 * (valid[0] + valid[1])
    g_edge = float(valid[1] - valid[0])
    if g_edge <= 0.0 or abs(f_c - f_target) >= f_target * rel_tol:
        return Gap.empty()
    f_c = float(f_c)

    # ── Step 2: track across all provided k-points ────────────────────────────
    gap_at_k = [g_edge]   # always include zone-edge value
    for i in range(n_k - 1):   # skip last row (already included above)
        row      = freqs_o[i, :]
        valid_row = np.sort(row[np.isfinite(row) & (row > 10e12)])
        if len(valid_row) < 2:
            continue
        below = valid_row[valid_row < f_c]
        above = valid_row[valid_row >= f_c]
        if len(below) > 0 and len(above) > 0:
            gap_at_k.append(above[0] - below[-1])

    min_gap = float(np.min(gap_at_k))
    if min_gap <= 0.0:
        return Gap.empty()

    ng = min_gap / f_c
    return Gap(f_c - min_gap / 2, f_c + min_gap / 2, f_c,
               min_gap, ng, -1, True)


# ── Score ─────────────────────────────────────────────────────────────────────

def score_result(G_o, G_m, f_o_c, f_m_c, t=None, w=None, a=None, *, require_opt=True,
                 require_mech=True, g_min_opt=G_MIN_OPT, g_min_mech=G_MIN_MECH):
    """Scalar score to maximise.  Higher is better.

    require_opt/require_mech: when False, that gap contributes NOTHING to
    the score — not as a reward (no +G term), not as a penalty (no gap-size
    penalty), and (for mechanical) not via the center-frequency penalty
    either. This is a hard "not in the cost function" switch, not just a
    relaxed threshold — use g_min_opt/g_min_mech for the latter.

    t/w: beam thickness/width [m]. Thinness/narrowness are added to `base`
    exactly like G_o/G_m — normalized linearly across the FULL bounds range
    (T_MIN_BOUND..T_MAX_BOUND / W_MIN_BOUND..W_MAX_BOUND) so smaller is
    ALWAYS strictly better, with no saturation point below which it stops
    improving. Only applied when require_mech=False (optical-only mode) —
    the co-optimized (require_mech=True) objective is intentionally left
    untouched, for backward compatibility with the seeded campaign's
    existing scores/behavior. Pass None to skip regardless of mode (e.g.
    when rescoring old records that didn't store params).

    a: lattice constant [m]. Used ONLY for the light-cone penalty below —
    unlike t/w, this applies whenever require_opt=True (both co-optimized
    AND optical-only runs), since a gap sitting above the light line isn't
    a real guided-mode gap regardless of which run produced it. Pass None
    to skip (e.g. old records without params).
    """
    if not (np.isfinite(G_o) and np.isfinite(G_m)):
        return -5.0
    base = 0.0
    pen_g = 0.0
    pen_of = 0.0
    pen_mf = 0.0
    pen_lc = 0.0
    if require_opt:
        base += G_o
        pen_g += LAMBDA_OPT * max(0.0, g_min_opt - G_o) ** 2
        pen_of = LAMBDA_OPT_F * ((f_o_c - OPT_F_TARGET) / OPT_F_TARGET) ** 2
        if a is not None:
            # Light line at the zone edge (k_z = pi/a), where the gap is
            # defined. Continuous: zero penalty at/below the light line,
            # growing quadratically with the fractional overshoot above it —
            # gives TPE a smooth signal back toward the guided region
            # instead of a step discontinuity.
            light_line = C0 / (2.0 * a)
            excess = max(0.0, (f_o_c - light_line) / light_line)
            pen_lc = LAMBDA_LIGHTCONE * excess ** 2
    if require_mech:
        base += G_m
        pen_g += LAMBDA_MECH * max(0.0, g_min_mech - G_m) ** 2
        if LAMBDA_MECH_F > 0:
            if f_m_c < MECH_F_MIN:
                pen_mf = LAMBDA_MECH_F * ((MECH_F_MIN - f_m_c) / MECH_F_MIN) ** 2
            elif f_m_c > MECH_F_MAX:
                pen_mf = LAMBDA_MECH_F * ((f_m_c - MECH_F_MAX) / MECH_F_MAX) ** 2
    else:
        if t is not None:
            thinness = (T_MAX_BOUND - t) / (T_MAX_BOUND - T_MIN_BOUND)
            base += LAMBDA_THICK * thinness
        if w is not None:
            narrowness = (W_MAX_BOUND - w) / (W_MAX_BOUND - W_MIN_BOUND)
            base += LAMBDA_WIDTH * narrowness
    return float(base - pen_g - pen_of - pen_mf - pen_lc)


# ── Full single-candidate evaluation ─────────────────────────────────────────

def _make_id(params: dict) -> str:
    """Hash physical geometry (nm-rounded) — bounds-independent unique ID."""
    key = json.dumps(
        {k: round(v * 1e9, 2) for k, v in sorted(params.items())},
        sort_keys=True,
    )
    return "opt_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def evaluate(model, u, g, *, n_k, n_bands_mech, n_bands_opt,
             study_mech, study_opt, compute_fy=True,
             require_opt=True, require_mech=True,
             g_min_opt=G_MIN_OPT, g_min_mech=G_MIN_MECH):
    """Set geometry, run both sweeps, return result dict.

    compute_fy=False: skip per-mode field queries in mechanical sweep
    (trusting mech sym BCs to have pre-selected breathing modes).

    require_opt/require_mech: whether each gap enters score_result() at
    all (see its docstring). The mechanical sweep is still ALWAYS run and
    G_m/f_m are still recorded even when require_mech=False — this only
    changes what the OPTIMIZER optimizes for; the data is kept so an
    optical-only run can still be inspected for "free" mechanical gaps
    later. require_opt/require_mech and the thresholds used are stored in
    the record itself so mixed-mode result files stay self-describing.
    """
    for name, val in g.as_dict().items():
        model.parameter(name, f"{val}[m]")

    # ── Mechanical ────────────────────────────────────────────────────────────
    k_m, freqs_m, fy_arr = _mech_sweep(model, g, study_mech, n_k, n_bands_mech,
                                        compute_fy=compute_fy)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        gp_m = breathing_gap(freqs_m, fy_arr)
        # Always also record the plain largest gap (no fy filtering) — cheap,
        # reuses freqs_m already computed above. When compute_fy=False this is
        # numerically identical to gp_m (fy_arr is all 1.0, no filtering
        # actually happens), but the field is kept distinct so runs made WITH
        # real fy have a genuine largest-vs-breathing comparison available.
        mech_window = (MECH_F_MIN * 0.5, MECH_F_MAX * 2.0)
        gp_m_largest = largest_gap(freqs_m, f_window=mech_window)
    G_m   = gp_m.normalized_gap if gp_m.found else 0.0
    f_m_c = gp_m.f_center        if gp_m.found else 0.0
    G_m_largest   = gp_m_largest.normalized_gap if gp_m_largest.found else 0.0
    f_m_largest_c = gp_m_largest.f_center        if gp_m_largest.found else 0.0

    # ── Optical ───────────────────────────────────────────────────────────────
    # Sweep zone-edge half of BZ (k ∈ [0.5,1.0]×π/a); detect gap nearest to
    # 1550 nm at zone edge then track across provided k-points.
    k_o, freqs_o = _opt_sweep(model, g, study_opt, n_k, n_bands_opt)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        gp_o = optical_gap_tracked(freqs_o)
    G_o   = gp_o.normalized_gap if gp_o.found else 0.0
    f_o_c = gp_o.f_center        if gp_o.found else OPT_F_TARGET

    sc = score_result(G_o, G_m, f_o_c, f_m_c, t=g.t, w=g.w, a=g.a, require_opt=require_opt,
                      require_mech=require_mech, g_min_opt=g_min_opt,
                      g_min_mech=g_min_mech)

    return dict(
        id=_make_id(g.as_dict()),
        u=[float(x) for x in u],
        params=g.as_dict(),
        optical_backend="comsol", mech_backend="comsol",
        G_o=G_o, G_m=G_m,
        f_o_THz=f_o_c * 1e-12,
        f_o_nm=C0 / f_o_c * 1e9 if f_o_c > 0 else float("nan"),
        f_m_GHz=f_m_c * 1e-9,
        G_m_largest=G_m_largest,
        f_m_largest_GHz=f_m_largest_c * 1e-9,
        fy_computed=bool(compute_fy),
        require_opt=bool(require_opt), require_mech=bool(require_mech),
        g_min_opt=float(g_min_opt), g_min_mech=float(g_min_mech),
        score=sc,
        status="success",
        # store band data for plotting
        _k_m=k_m.tolist(), _freqs_m=freqs_m.tolist(), _fy_arr=fy_arr.tolist(),
        _k_o=k_o.tolist(), _freqs_o=freqs_o.tolist(),
        optical_gap=G_o, mechanical_gap=G_m,
        optical_center_frequency=float(f_o_c),
        mechanical_center_frequency=float(f_m_c),
    )


# ── Optimizer (Optuna or Halton fallback) ─────────────────────────────────────

class HaltonSampler:
    def __init__(self, seed=0):
        self.rng = np.random.default_rng(seed)
        self.count = 0

    def _halton(self, i, base):
        f, r = 1.0, 0.0
        while i > 0:
            f /= base; r += f * (i % base); i //= base
        return r

    def ask(self):
        bases = [2, 3, 5, 7]
        u = [self._halton(self.count + 1, b) for b in bases]
        self.count += 1
        return u

    def tell(self, u, score):
        pass


def make_optimizer(n_init, seed=0, prior_results=None, bounds=None,
                   require_opt=True, require_mech=True,
                   g_min_opt=G_MIN_OPT, g_min_mech=G_MIN_MECH):
    try:
        import optuna
        from geometry import phys_to_u
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(n_startup_trials=n_init, seed=seed))
        # Warm-start: recompute u from PHYSICAL params using CURRENT bounds so
        # that Optuna seeds correctly regardless of any bounds changes.
        if prior_results:
            dists = {f"u{i}": optuna.distributions.FloatDistribution(0, 1)
                     for i in range(5)}
            n_added = 0
            for r in prior_results:
                params = r.get("params")
                if params is None or r.get("status") not in ("success", None):
                    continue
                # Derive u from physical params with current bounds — no migration needed
                u = phys_to_u(params, bounds)
                if len(u) != 5:
                    continue
                # Re-compute score with CURRENT constants AND this run's
                # require_opt/require_mech — not whatever mode each record was
                # originally scored under, so TPE's model is coherent with
                # what THIS run is actually optimizing for (files may mix
                # records from runs with different requirement settings).
                G_o  = float(r.get("G_o", 0) or 0)
                G_m  = float(r.get("G_m", 0) or 0)
                f_oc = float(r.get("optical_center_frequency", OPT_F_TARGET) or OPT_F_TARGET)
                f_mc = float(r.get("mechanical_center_frequency", 0) or 0)
                if f_oc <= 0:
                    f_oc = OPT_F_TARGET
                t_val = float(params.get("t")) if params.get("t") is not None else None
                w_val = float(params.get("w")) if params.get("w") is not None else None
                a_val = float(params.get("a")) if params.get("a") is not None else None
                sc = score_result(G_o, G_m, f_oc, f_mc, t=t_val, w=w_val, a=a_val,
                                  require_opt=require_opt,
                                  require_mech=require_mech, g_min_opt=g_min_opt,
                                  g_min_mech=g_min_mech)
                trial_params = {f"u{i}": float(u[i]) for i in range(5)}
                trial = optuna.trial.create_trial(params=trial_params,
                                                  distributions=dists,
                                                  value=float(sc))
                study.add_trial(trial)
                n_added += 1
            print(f"  [Optuna] seeded with {n_added} prior results (re-scored, u from physics)")
        return ("optuna", study)
    except ImportError:
        return ("halton", HaltonSampler(seed=seed))


def opt_ask(optimizer):
    kind, obj = optimizer
    if kind == "optuna":
        import optuna
        trial = obj.ask(
            {f"u{i}": optuna.distributions.FloatDistribution(0, 1)
             for i in range(5)})
        u = [trial.params[f"u{i}"] for i in range(5)]
        return u, trial
    else:
        u = obj.ask()           # Halton: 4-D; pad with 0.5 (t=mid-range)
        return u + [0.5], None


def opt_tell(optimizer, u, score, trial=None):
    kind, obj = optimizer
    if kind == "optuna":
        obj.tell(trial, score)
    else:
        obj.tell(u, score)


# ── End-of-run characterization ───────────────────────────────────────────────

# ── NPZ cache helpers (serialize/deserialize Gap + None arrays) ────────────────

def _gap_to_arr(gap):
    """Pack a Gap dataclass into a float64[7] array for NPZ storage."""
    return np.array([gap.f_lower, gap.f_upper, gap.f_center,
                     gap.gap_size, gap.normalized_gap,
                     float(gap.lower_band), float(gap.found)], dtype=np.float64)


def _arr_to_gap(arr):
    """Restore a Gap dataclass from a float64[7] array."""
    # Gap is already imported at module level via `from bandgap import Gap`
    return Gap(float(arr[0]), float(arr[1]), float(arr[2]),
               float(arr[3]), float(arr[4]),
               int(arr[5]), bool(arr[6]))


def _gaps_to_arr(gaps):
    """Pack a list[Gap] into a float64[n,7] array (empty n=0 if no gaps)."""
    if not gaps:
        return np.zeros((0, 7), dtype=np.float64)
    return np.stack([_gap_to_arr(gp) for gp in gaps])


def _arr_to_gaps(arr):
    """Restore a list[Gap] from a float64[n,7] array."""
    return [_arr_to_gap(row) for row in np.asarray(arr)]


def _n2e(arr):
    """None → empty 1-D array (for NPZ savez, which cannot store None)."""
    return np.array([], dtype=np.float64) if arr is None else np.asarray(arr)


def _e2n(arr):
    """Empty array → None (inverse of _n2e)."""
    return None if len(arr) == 0 else arr


def _full_opt_sweep(model, g, study, n_k=15, n_bands=6, min_THz=10.0):
    """Sweep the full Brillouin zone (k: 0→1×π/a) for optical bands.

    Distinct from _opt_sweep which covers only the zone-edge half; this version
    is used for visualization so the full dispersion is shown. Same
    light-line-tracking shift + widened neigs as `_opt_sweep` (see
    `_guided_shift_hz`) so the fundamental guided band isn't skipped here
    either — this is the sweep that feeds the characterization figure.
    """
    k_frac = np.linspace(1e-3, 1.0, n_k)
    freqs  = np.full((n_k, n_bands), np.nan)
    neigs_req = max(n_bands, 10)
    for i, s in enumerate(k_frac):
        kz = np.pi * s / g.a
        model.parameter("kF", f"{kz}[1/m]")
        _set_eig_search(model, study, _guided_shift_hz(kz), neigs=neigs_req)
        model.solve(study)
        dset = next((d for d in model.datasets() if d.startswith(study + '//')), None)
        ev = (model.evaluate("freq", dataset=dset) if dset
              else model.evaluate("freq"))
        fs = np.sort(np.real(np.asarray(ev)))
        fs = fs[fs > min_THz * 1e12]
        n_av = min(n_bands, len(fs))
        freqs[i, :n_av] = fs[:n_av]
    return k_frac, freqs


def _draw_geometry(ax, g, n_periods=3, fontsize=8):
    """Draw a top-view (z-y) unit cell sketch with annotated dimensions."""
    nm = 1e9
    a, w, hx, hy = g.a*nm, g.w*nm, g.hx*nm, g.hy*nm
    N = n_periods
    ax.add_patch(Rectangle((-N*a/2, -w/2), N*a, w,
                            facecolor="#cfe8ff", edgecolor="#1f5fa6", lw=1.5, zorder=1))
    for i in range(N):
        zc = -N*a/2 + (i + 0.5)*a
        ax.add_patch(Ellipse((zc, 0), 2*hx, 2*hy,
                              facecolor="white", edgecolor="#1f5fa6", lw=1.2, zorder=2))
    arrow = dict(arrowstyle="<->", color="k", lw=1)
    # period arrow
    ax.annotate("", xy=(-N*a/2,    -w/2-80), xytext=(-N*a/2+a, -w/2-80),
                arrowprops=arrow)
    ax.text(-N*a/2+a/2, -w/2-135, f"a = {a:.0f} nm", ha="center", fontsize=fontsize)
    # width arrow
    ax.annotate("", xy=(N*a/2+55, -w/2), xytext=(N*a/2+55, w/2), arrowprops=arrow)
    ax.text(N*a/2+90, 0, f"w = {w:.0f} nm", rotation=90, va="center", fontsize=fontsize)
    # hx arrow (one cell) — placed ABOVE the beam to avoid collision with 'a' arrow
    zc0 = -N*a/2 + 0.5*a
    ax.annotate("", xy=(zc0-hx, w/2+55), xytext=(zc0+hx, w/2+55),
                arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8))
    ax.text(zc0, w/2+100, f"2hx={2*hx:.0f}", ha="center", fontsize=fontsize-1,
            color="#555")
    # hy arrow
    ax.annotate("", xy=(-N*a/2-55, -hy), xytext=(-N*a/2-55, hy), arrowprops=arrow)
    ax.text(-N*a/2-100, 0, f"2hy\n{2*hy:.0f}", ha="center", va="center",
            fontsize=fontsize-1, color="#555")

    bridge   = a - 2*hx
    sidewall = w/2 - hy
    t_nm     = g.t * nm
    ax.set_title(f"Top view  (z-y)  |  t = {t_nm:.0f} nm\n"
                 f"bridge = {bridge:.0f} nm    sidewall = {sidewall:.0f} nm",
                 fontsize=fontsize)
    ax.set_xlabel("z  [nm]  (beam / periodic axis)", fontsize=fontsize)
    ax.set_ylabel("y  [nm]  (beam width)", fontsize=fontsize)
    ax.set_aspect("equal")
    ax.set_xlim(-N*a/2-160, N*a/2+230)
    ax.set_ylim(-w/2-230, w/2+160)
    ax.tick_params(labelsize=fontsize-1)


def _draw_cross_section(ax, g, fontsize=8):
    """Draw the beam cross-section (y-x plane, at z=0 through the hole center).

    Anisotropic (crystallographic) diamond etch leaves a triangular wedge
    below the rectangular slab: base flush with the slab's bottom face
    (x=-t/2, spanning y=-w/2..w/2), apex pointing in -x, base half-angles
    30 deg so the apex angle is 120 deg. Apex depth = (w/2)*tan(30deg)
    = w/(2*sqrt(3)) — matches comsol/omc_unitcell_iso.mph's Polygon 1
    vertices (-t/2,+-w/2), (-t/2-w/(2*sqrt(3)), 0). The elliptical hole is a
    constant-width (2*hy) prism cut along x, so it passes through both the
    rectangle and the triangle.
    """
    nm = 1e9
    w, t, hy = g.w*nm, g.t*nm, g.hy*nm
    tri_h = w / (2.0 * np.sqrt(3.0))          # triangle depth below the slab
    x_bot = -t/2 - tri_h                      # apex x (deepest point, -x direction)

    # Combined solid outline: rectangle top + triangle wedge below it.
    outline_y = [-w/2, w/2, w/2, 0.0, -w/2]
    outline_x = [t/2, t/2, -t/2, x_bot, -t/2]

    # Hole: constant-width (2*hy) prism cut through the FULL x-extent
    # (rectangle + triangle), matching the physical through-etch. Clip it to
    # the outline via a real polygon boolean difference (not just an
    # overlaid rectangle) so the drawn solid never appears to extend past
    # the true rect+triangle boundary — the outline determines the shape,
    # the hole only removes material where it actually overlaps solid.
    try:
        from shapely.geometry import Polygon
        outline_poly = Polygon(zip(outline_y, outline_x))
        hole_poly = Polygon([(-hy, x_bot), (hy, x_bot), (hy, t/2), (-hy, t/2)])
        solid = outline_poly.difference(hole_poly)
        pieces = list(solid.geoms) if solid.geom_type == "MultiPolygon" else [solid]
        for piece in pieces:
            if piece.is_empty:
                continue
            py, px = piece.exterior.xy
            ax.fill(py, px, facecolor="#cfe8ff", edgecolor="#1f5fa6", lw=1.5,
                    zorder=1)
        # Ghost contour of the ORIGINAL (pre-hole) rect+triangle outline, so
        # it's easy to see where material was removed vs. where the outline
        # was untouched by the hole.
        ax.plot(outline_y + [outline_y[0]], outline_x + [outline_x[0]],
                linestyle="--", color="#1f5fa6", lw=0.9, alpha=0.55, zorder=1.5)
    except Exception:
        # shapely unavailable: fall back to outline + overlaid hole rectangle
        # (may visually overshoot the taper — the boolean version above is
        # the correct rendering).
        ax.fill(outline_y, outline_x, facecolor="#cfe8ff", edgecolor="#1f5fa6",
                lw=1.5, zorder=1)
        ax.add_patch(Rectangle((-hy, x_bot), 2*hy, (t/2 - x_bot),
                                facecolor="white", edgecolor="#1f5fa6", lw=1.2,
                                zorder=2))

    arrow = dict(arrowstyle="<->", color="k", lw=1)
    # width (of the rectangular part)
    ax.annotate("", xy=(-w/2, t/2+30), xytext=(w/2, t/2+30), arrowprops=arrow)
    ax.text(0, t/2+80, f"w = {w:.0f} nm", ha="center", fontsize=fontsize)
    # thickness (rectangular slab only)
    ax.annotate("", xy=(w/2+30, -t/2), xytext=(w/2+30, t/2), arrowprops=arrow)
    ax.text(w/2+60, 0, f"t = {t:.0f} nm", rotation=90, va="center", fontsize=fontsize)
    # triangle depth
    ax.annotate("", xy=(w/2+30, -t/2), xytext=(w/2+30, x_bot), arrowprops=arrow)
    ax.text(w/2+60, (-t/2 + x_bot)/2, f"{tri_h:.0f}", rotation=90, va="center",
            fontsize=fontsize-1, color="#555")
    # base half-angle label
    ax.text(-w/2*0.55, -t/2 - tri_h*0.18, "30°", fontsize=fontsize-1, color="#555")
    # sidewall (placed above the width label to avoid overlap)
    sw = w/2 - hy
    ax.annotate("", xy=(hy, t/2+130), xytext=(w/2, t/2+130), arrowprops=arrow)
    ax.text((hy+w/2)/2, t/2+155, f"sw={sw:.0f}", ha="center", fontsize=fontsize-1,
            color="#555")

    ax.set_title(f"Cross-section  (y-x)  at hole centre\n"
                 f"t = {t:.0f} nm  |  triangle depth = {tri_h:.0f} nm",
                 fontsize=fontsize)
    ax.set_xlabel("y  [nm]  (width)", fontsize=fontsize)
    ax.set_ylabel("x  [nm]  (thickness, - = etch direction)", fontsize=fontsize)
    ax.set_aspect("equal")
    ax.set_xlim(-w/2-100, w/2+130)
    ax.set_ylim(x_bot - 100, t/2 + 210)
    ax.tick_params(labelsize=fontsize-1)


def _save_geometry_fig(g, out_path):
    """Standalone geometry figure: top view (3 periods) + cross-section."""
    fig, (ax_top, ax_xc) = plt.subplots(1, 2, figsize=(11, 4.5),
                                         gridspec_kw={"width_ratios": [3, 1]})
    _draw_geometry(ax_top, g, n_periods=3, fontsize=9)
    _draw_cross_section(ax_xc, g, fontsize=9)
    nm = 1e9
    fig.suptitle(
        f"Diamond OMC unit cell  |  "
        f"a={g.a*nm:.0f}  w={g.w*nm:.0f}  t={g.t*nm:.0f}  "
        f"hx={g.hx*nm:.0f}  hy={g.hy*nm:.0f} nm",
        fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Geometry  -> {out_path}")


def _save_optical_fig(k_TE, freqs_TE, k_TM, freqs_TM, gap_TE,
                      y_om, z_om, normE, wl_mode, g_nm, out_path):
    """Standalone optical figure: band structure + TE mode profile."""
    fig, (ax_bands, ax_mode) = plt.subplots(1, 2, figsize=(12, 5))
    _plot_opt_bands(ax_bands, k_TE, freqs_TE, k_TM, freqs_TM, gap_TE,
                    a_m=g_nm.a * 1e-9)
    ax_bands.tick_params(labelsize=9)
    ax_bands.set_xlabel(r"$k_z a / \pi$", fontsize=10)
    ax_bands.set_ylabel("Frequency [THz]", fontsize=10)
    ax_bands.set_title("Optical band structure  (TE + TM, below light cone = guided)",
                       fontsize=10)

    title = (f"TE dielectric band  |Ey|  (λ≈{wl_mode:.0f} nm)"
             if wl_mode > 0 else "TE dielectric band  |Ey|")
    _plot_mode_profile(ax_mode, y_om, z_om, normE, None, None, g_nm, title,
                       cmap="inferno")
    ax_mode.set_title(title, fontsize=10)

    gap_str = (f"G_o = {gap_TE.normalized_gap*100:.1f}%  "
               f"(center {C0/gap_TE.f_center*1e9:.0f} nm)"
               if gap_TE.found else "No TE gap found")
    fig.suptitle(f"Optical characterization  |  {gap_str}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Optical   -> {out_path}")


def _save_mech_fig(k_sym, freqs_sym, fy_sym, k_anti, freqs_anti, gap_mech,
                   y_mm, z_mm, mag_u, uy_arr, uz_arr, f_mmode, g_nm, out_path,
                   extra_gaps=None):
    """Standalone mechanical figure: band structure + breathing mode profile."""
    fig, (ax_bands, ax_mode) = plt.subplots(1, 2, figsize=(12, 5))
    _plot_mech_bands(ax_bands, k_sym, freqs_sym, fy_sym, k_anti, freqs_anti, gap_mech,
                     extra_gaps=extra_gaps)
    ax_bands.tick_params(labelsize=9)
    ax_bands.set_xlabel(r"$k_z a / \pi$", fontsize=10)
    ax_bands.set_ylabel("Frequency [GHz]", fontsize=10)
    ax_bands.set_title("Mechanical band structure  (sym + antisym)", fontsize=10)

    title = (f"Breathing mode  |uy|  ({f_mmode*1e-9:.2f} GHz)"
             if f_mmode > 0 else "Breathing mode  |uy|")
    _plot_mode_profile(ax_mode, y_mm, z_mm, mag_u, uy_arr, uz_arr, g_nm, title,
                       cmap="magma", mask_to_solid=True)
    ax_mode.set_title(title, fontsize=10)

    gap_str = (f"G_m = {gap_mech.normalized_gap*100:.1f}%  "
               f"(center {gap_mech.f_center*1e-9:.2f} GHz)"
               if gap_mech.found else "No breathing gap found")
    fig.suptitle(f"Mechanical characterization  |  {gap_str}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Mechanical -> {out_path}")


def _plot_opt_bands(ax, k_TE, freqs_TE, k_TM, freqs_TM, gap_TE, a_m=None):
    """Plot TE + TM optical bands on ax. Gap shaded. Light cone added if a_m given.

    Only shows the zone-boundary half of the BZ (k > 0.5) where the gap lives,
    and centers the y-axis on the gap so it is always visible.
    """
    C0_local = 299_792_458.0
    target_THz = C0_local / 1550e-9 * 1e-12

    # ── x-range: only show zone-boundary half (gap lives here) ──────────────
    k_lo = 0.5
    ax.set_xlim(k_lo, 1.0)

    # ── y-range: center on gap (or 1550 nm if no gap) ───────────────────────
    if gap_TE.found:
        gap_center = (gap_TE.f_lower + gap_TE.f_upper) / 2 * 1e-12   # THz
        gap_half   = (gap_TE.f_upper - gap_TE.f_lower) / 2 * 1e-12
        # Show 2× the half-gap on each side so adjacent bands are visible
        half_range = max(gap_half * 2.0, target_THz * 0.10)
    else:
        gap_center = target_THz
        half_range = target_THz * 0.20
    y_lo = gap_center - half_range
    y_hi = gap_center + half_range
    ax.set_ylim(y_lo, y_hi)

    # Light cone in zoomed region
    if a_m is not None and a_m > 0:
        k_lc = np.linspace(k_lo, 1.0, 200)
        f_lc = C0_local / (2.0 * a_m) * k_lc * 1e-12

        # If the light line sits entirely below the zoomed gap window, the
        # fill/line are invisible (off-screen) and the plot silently hides
        # the fact that these modes are ABOVE the light line (radiative, not
        # index-guided). Extend the window down (capped) so the boundary is
        # always visible, and flag it explicitly.
        if f_lc.max() < y_lo:
            y_lo = max(f_lc.max() * 0.85, y_hi - 4 * half_range)
            ax.set_ylim(y_lo, y_hi)
            ax.text(0.5, 0.95, "ABOVE LIGHT LINE — radiative, not index-guided",
                    transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                    fontweight="bold", color="#b00000",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8,
                              edgecolor="#b00000"))

        ax.fill_between(k_lc, f_lc, y_hi, color="gray", alpha=0.12, zorder=0,
                        label="Light cone")
        ax.plot(k_lc, f_lc, "-", color="gray", lw=0.9, alpha=0.5)

    # Filter to k > k_lo before plotting so lines don't clip awkwardly
    mask_TE = k_TE >= k_lo
    mask_TM = k_TM >= k_lo
    for b in range(freqs_TE.shape[1]):
        col = freqs_TE[:, b] * 1e-12
        valid = np.isfinite(col) & mask_TE
        if valid.any():
            ax.plot(k_TE[valid], col[valid], "-", color="#1f5fa6", lw=1.4,
                    label="opt TE study" if b == 0 else "")
    for b in range(freqs_TM.shape[1]):
        col = freqs_TM[:, b] * 1e-12
        valid = np.isfinite(col) & mask_TM
        if valid.any():
            ax.plot(k_TM[valid], col[valid], "--", color="#cc3333", lw=1.0,
                    alpha=0.7, label="opt TM study" if b == 0 else "")
    if gap_TE.found:
        lo, hi = gap_TE.f_lower*1e-12, gap_TE.f_upper*1e-12
        ax.axhspan(lo, hi, color="#ffd27f", alpha=0.55, zorder=0,
                   label=f"TE-study gap {gap_TE.normalized_gap*100:.1f}%")
    ax.axhline(target_THz, ls=":", color="#444", lw=1.0, label="1550 nm")
    ax.set_xlabel(r"$k_z a / \pi$", fontsize=8)
    ax.set_ylabel("Frequency [THz]", fontsize=8)
    ax.set_title("Optical band structure\n(guided region = below light cone)", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)


def _plot_mech_bands(ax, k_sym, freqs_sym, fy_sym, k_anti, freqs_anti, gap_mech,
                     extra_gaps=None):
    """Plot mech bands with continuous fy color coding via LineCollection.

    `gap_mech` (fy-verified breathing gap) is shaded green, as before.
    `extra_gaps` (optional list[Gap], unfiltered — may include non-breathing
    gaps) are outlined in gray so it's visible that OTHER candidate gaps
    exist in the raw "mech sym" spectrum besides the one chosen as "the"
    mechanical gap.
    """
    from matplotlib.collections import LineCollection

    cmap_fy = plt.cm.RdYlBu_r   # blue=low fy, red=high fy (breathing)
    norm_fy = plt.Normalize(0, 1)
    lc_last = None

    for b in range(freqs_sym.shape[1]):
        f_col = freqs_sym[:, b] * 1e-9
        fy_col = fy_sym[:, b]
        valid = np.isfinite(f_col) & (f_col > 0)
        if not valid.any():
            continue
        k_v = k_sym[valid]
        f_v = f_col[valid]
        fy_v = np.where(np.isfinite(fy_col[valid]), np.clip(fy_col[valid], 0, 1), 0.0)

        # Build segments so each segment can be colored by its midpoint fy
        pts  = np.array([k_v, f_v]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        fy_mid = (fy_v[:-1] + fy_v[1:]) / 2

        lc = LineCollection(segs, array=fy_mid, cmap=cmap_fy, norm=norm_fy,
                            linewidth=1.8, zorder=3)
        ax.add_collection(lc)
        lc_last = lc

    for b in range(freqs_anti.shape[1]):
        col = freqs_anti[:, b] * 1e-9
        valid = np.isfinite(col) & (col > 0)
        if valid.any():
            ax.plot(k_anti[valid], col[valid], "--", color="#999999", lw=0.9,
                    alpha=0.6, label="Antisym" if b == 0 else "")

    # Other candidate gaps (unfiltered by fy) — drawn as SHORT brackets at the
    # left margin only (not full-width axhspans), and capped to the largest
    # few, so they don't visually bury the actual dispersion curves.
    if extra_gaps:
        shown = [gp for gp in extra_gaps
                 if gp.normalized_gap >= 0.03
                 and not (gap_mech.found and abs(gp.f_lower*1e-9 - gap_mech.f_lower*1e-9) < 1e-6)]
        shown = shown[:5]
        for n, gp in enumerate(shown):
            lo, hi = gp.f_lower*1e-9, gp.f_upper*1e-9
            ax.axhspan(lo, hi, xmin=0.0, xmax=0.05, facecolor="none",
                       edgecolor="#888888", linestyle="--", linewidth=1.0,
                       alpha=0.8, zorder=1, label="Other sym gap" if n == 0 else "")
            ax.text(0.06, (lo + hi) / 2, f"{gp.normalized_gap*100:.0f}%",
                    fontsize=6, color="#666666", va="center", ha="left",
                    transform=ax.get_yaxis_transform())

    if gap_mech.found:
        lo, hi = gap_mech.f_lower*1e-9, gap_mech.f_upper*1e-9
        ax.axhspan(lo, hi, color="#c8f5c8", alpha=0.6, zorder=0,
                   label=f"Breathing gap {gap_mech.normalized_gap*100:.1f}%")

    # Colorbar for fy
    sm = plt.cm.ScalarMappable(cmap=cmap_fy, norm=norm_fy)
    sm.set_array([])
    try:
        cb = plt.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
        cb.set_label("fy  (breathing frac.)", fontsize=6)
        cb.ax.tick_params(labelsize=6)
    except Exception:
        pass

    ax.autoscale_view()
    ax.set_xlabel(r"$k_z a / \pi$", fontsize=8)
    ax.set_ylabel("Frequency [GHz]", fontsize=8)
    ax.set_title("Mechanical band structure  (sym, colored by fy)", fontsize=9)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)


def _last_outer_idx(model, dset):
    """Return the 1-based index of the last outer solution in the dataset.

    After a BZ sweep followed by a zone-edge re-solve, the dataset accumulates
    N_sweep + 1 outer solutions.  Pass this index to evaluate() so we select
    only the zone-edge re-solve data, not the entire sweep.
    """
    if not dset:
        return None
    try:
        n_all = len(np.asarray(model.evaluate("freq", dataset=dset)).ravel())
        n_one = len(np.asarray(model.evaluate("freq", dataset=dset, outer=1)).ravel())
        if n_one > 0:
            return max(1, n_all // n_one)
    except Exception:
        pass
    return None


def _mph_dset_to_java_tag(model, mph_path):
    """Convert an MPh dataset path ('Study//Solution 1') → COMSOL internal tag ('dset1').

    MPh uses double-slash paths; COMSOL labels use single slash.
    Searches model.java.result().dataset() for a matching label.
    """
    comsol_label = mph_path.replace("//", "/")
    rdsets = model.java.result().dataset()
    for tag in rdsets.tags():
        t = str(tag)
        try:
            if str(rdsets.get(t).label()) == comsol_label:
                return t
        except Exception:
            pass
    return None


_STABLE_CUTPLANE_CACHE = {}


def _get_stable_cutplane(model, dset):
    """Get (or create once) the x=0 (yz) CutPlane for this `dset`, cached
    and reused across calls instead of recreated fresh every time.

    Recreating a CutPlane fresh on every call was confirmed (by direct
    inspection) to give a DIFFERENT, non-deterministic point set each time
    — 305 vs 357 points for the identical setup — presumably from some
    internal auto-meshing/sampling choice that isn't seeded deterministically.
    That made comparisons ACROSS separate calls (e.g. "is mode 1's field
    similar to mode 2's?") meaningless, since the two calls weren't even
    sampling the same physical points. `dset` (e.g. "opt TE//Solution 1")
    is a SOLUTION dataset that COMSOL overwrites in place on every
    `model.solve()` — the CutPlane's own geometric definition (the x=0 cut)
    doesn't change across k-points/solves, so caching it here is safe and
    automatically reflects whatever the latest solve produced.
    """
    key = str(dset)
    tag = _STABLE_CUTPLANE_CACHE.get(key)
    if tag is not None:
        try:
            model.java.result().dataset(tag).label()  # verify still exists
            return tag
        except Exception:
            _STABLE_CUTPLANE_CACHE.pop(key, None)

    tag = f"omc_cp_{abs(hash(key)) % 1_000_000}"
    try:
        model.java.result().dataset().remove(tag)
    except Exception:
        pass
    internal_dset_tag = _mph_dset_to_java_tag(model, dset) if dset else None
    cp = model.java.result().dataset().create(tag, "CutPlane")
    if internal_dset_tag:
        cp.set("data", internal_dset_tag)
    cp.set("quickplane", "yz")   # plane normal to x; placed at x=0
    _STABLE_CUTPLANE_CACHE[key] = tag
    return tag


_STABLE_CUTPOINT_CACHE = {}


def _get_stable_cutpoint_grid(model, dset, g, ny=50, nz=50):
    """Get (or create once) an explicit regular CutPoint3D grid at x=0,
    spanning y in [0, w] and z in [-a/2, a/2], cached and reused across
    mode/solnum queries (same rationale as `_get_stable_cutplane`).

    Unlike CutPlane's auto-generated (sparse, non-deterministic point
    count/order) sampling, CutPoint3D takes EXPLICIT coordinates — verified
    by an exact point-by-point cross-check against COMSOL's own dense
    Multislice export (agreement to ~1%) to give the CORRECT field values.
    This is also the documented approach: COMSOL's mphinterp / Cut Point
    datasets are built for exactly this — interpolating a solution at
    arbitrary specified points.
    """
    key = (str(dset), ny, nz)
    tag = _STABLE_CUTPOINT_CACHE.get(key)
    if tag is not None:
        try:
            model.java.result().dataset(tag).label()
            return tag
        except Exception:
            _STABLE_CUTPOINT_CACHE.pop(key, None)

    tag = f"omc_cpt_{abs(hash(key)) % 1_000_000}"
    try:
        model.java.result().dataset().remove(tag)
    except Exception:
        pass
    internal_dset_tag = _mph_dset_to_java_tag(model, dset) if dset else None

    y = np.linspace(0.0, g.w, ny)
    z = np.linspace(-g.a / 2.0, g.a / 2.0, nz)
    ZZ, YY = np.meshgrid(z, y)   # shape [ny, nz]
    y_flat = YY.ravel()
    z_flat = ZZ.ravel()
    x_flat = np.zeros_like(y_flat)

    cpt = model.java.result().dataset().create(tag, "CutPoint3D")
    if internal_dset_tag:
        cpt.set("data", internal_dset_tag)
    cpt.set("pointx", [str(v) for v in x_flat])
    cpt.set("pointy", [str(v) for v in y_flat])
    cpt.set("pointz", [str(v) for v in z_flat])
    _STABLE_CUTPOINT_CACHE[key] = tag
    return tag


def _eval_mode_grid(model, dset, solnum, g, field_exprs, outer=None, abs_fields=True,
                    ny=50, nz=50):
    """Evaluate mode field on an EXPLICIT regular grid (default 50x50) at
    x=0, via CutPoint3D + Export/Data. Drop-in replacement for
    `_eval_mode_cutplane` with the same (y_nm, z_nm, fields) contract, but
    on a genuinely dense, deterministic regular grid instead of CutPlane's
    sparse (~300-900 point) auto-sampling.

    Writes the export to a permanent file under results/mode_exports/, same
    as `_eval_mode_cutplane`.
    """
    EXP_TAG = "omc_tmp_exp_grid"
    y_nm = z_nm = None
    fields = {k: None for k in field_exprs}
    out_path = None

    try:
        model.java.result().export().remove(EXP_TAG)
    except Exception:
        pass

    try:
        cpt_tag = _get_stable_cutpoint_grid(model, dset, g, ny=ny, nz=nz)

        labels = list(field_exprs.keys())
        exprs = list(field_exprs.values())

        export_dir = os.path.join(os.path.dirname(__file__), "..",
                                  "results", "mode_exports")
        os.makedirs(export_dir, exist_ok=True)
        study_slug = re.sub(r"[^A-Za-z0-9_-]+", "_", str(dset) or "study")
        kind = "abs" if abs_fields else "complex"
        fname = (f"{study_slug}_solnum{solnum}_outer{outer}_"
                 f"{'_'.join(labels)}_{kind}_grid{ny}x{nz}.txt")
        out_path = os.path.join(export_dir, fname)

        exp = model.java.result().export().create(EXP_TAG, "Data")
        exp.set("data", cpt_tag)
        exp.set("expr", exprs)
        exp.set("filename", out_path)
        if outer is not None:
            exp.set("outersolnum", str(int(outer)))
        exp.set("innerinput", "manual")
        exp.set("solnum", str(int(solnum)))
        exp.run()

        coords, vals = _parse_comsol_data_export(out_path, len(exprs))
        n = coords.shape[0]
        if n == 0:
            raise ValueError("CutPoint3D grid export returned 0 points")

        y_nm = coords[:, 1] * 1e9
        z_nm = coords[:, 2] * 1e9

        for label, v in zip(labels, vals):
            if len(v) != n:
                print(f"    [grid export size mismatch '{label}': {len(v)} vs {n}]")
                fields[label] = None
            else:
                fields[label] = np.abs(v) if abs_fields else v

    except Exception as e:
        print(f"    [grid eval failed: {e}]")
        y_nm = z_nm = None
    finally:
        try:
            model.java.result().export().remove(EXP_TAG)
        except Exception:
            pass

    return y_nm, z_nm, fields


def _java_eval_on_dataset(model, dataset_tag, expression, outer=None, inner=None):
    """Evaluate a COMSOL expression on any dataset using model.java directly.

    Bypasses MPh's model.evaluate() validation, which rejects derived datasets
    (CutPlane, CutLine, etc.) whose 'data' property is a dataset tag rather
    than a solution tag.  model.java.result().numerical() == model/'evaluations'
    in MPh nomenclature (confirmed in mph/node.py line 122).

    inner: list of 1-based eigenmode indices to select (e.g. [solnum]).
    outer: 1-based outer solution (k-point) index.

    Returns a 1D numpy array (complex or real).
    """
    TMP_TAG = "omc_byp_eval"
    try:
        try:
            model.java.result().numerical().remove(TMP_TAG)
        except Exception:
            pass

        ev = model.java.result().numerical().create(TMP_TAG, "Eval")
        ev.set("data", dataset_tag)
        ev.set("expr", expression)
        if outer is not None:
            ev.set("outersolnum", int(outer))   # JInt via JPype coercion

        # getData() triggers evaluation; isComplex() reflects last result.
        # MPh shape convention: [n_inner_solutions, n_pts] or [1, n_inner, n_pts].
        raw = np.array(ev.getData())
        is_complex = bool(ev.isComplex())
        if is_complex:
            raw = raw.astype(complex) + 1j * np.array(ev.getImagData())
        else:
            raw = raw.astype(float)

        # Flatten outer-expression dimension if present (mirrors MPh slicing at
        # model.py:599-606: results has shape [?, n_inner, n_pts]).
        if raw.ndim == 3:
            raw = raw[0]   # drop expression-batch dim → [n_inner, n_pts]

        # Select the requested inner eigenmode (1-based).
        if inner is not None and raw.ndim == 2:
            idx = max(0, int(np.array(inner, dtype=int).ravel()[0]) - 1)
            idx = min(idx, raw.shape[0] - 1)
            raw = raw[idx]   # → [n_pts]

        return raw.ravel()
    finally:
        try:
            model.java.result().numerical().remove(TMP_TAG)
        except Exception:
            pass


_COMPLEX_RE = re.compile(r'^([+-]?[\d.]+(?:[eE][+-]?\d+)?)([+-][\d.]+(?:[eE][+-]?\d+)?)i$')


def _parse_comsol_data_export(path, n_exprs):
    """Parse a COMSOL Results>Export>Data text file.

    Format: comment lines starting with '%', then data rows
    'X Y Z <expr1> <expr2> ...' where each value is either a plain real
    number or COMSOL's 'a+bi'/'a-bi' complex notation (no space, no '*').
    Returns (coords[n,3] float64, [values_expr1[n], values_expr2[n], ...]
    as complex128 arrays).
    """
    coords, vals = [], [[] for _ in range(n_exprs)]

    def _to_complex(s):
        m = _COMPLEX_RE.match(s)
        if m:
            return complex(float(m.group(1)), float(m.group(2)))
        return complex(float(s))

    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < 3 + n_exprs:
                continue
            coords.append((float(parts[0]), float(parts[1]), float(parts[2])))
            for i in range(n_exprs):
                vals[i].append(_to_complex(parts[3 + i]))
    coords = np.array(coords, dtype=np.float64) if coords else np.zeros((0, 3))
    return coords, [np.array(v, dtype=np.complex128) for v in vals]


def _eval_mode_cutplane(model, dset, solnum, g, field_exprs, outer=None, abs_fields=True):
    """Evaluate mode field on x=0 via a COMSOL CutPlane + Export/Data node.

    COMSOL's "Eval" numerical node cannot attach to a CutPlane dataset
    (throws "Invalid dataset type" — confirmed by direct introspection of
    the live model). Its Export/Data node CAN, and gives properly
    field-recovered values through COMSOL's own interpolation — unlike
    evaluating expressions at raw mesh nodes (`_eval_mode_on_midplane`),
    which is fine for nodal solid-mechanics DOFs (u,v,w) but NOT reliable
    for edge-element EM fields (ewfd.Ey etc. don't have simple per-node
    values the way nodal displacement does).

    IMPORTANT: the CutPlane dataset is created ONCE per `dset` and cached
    (`_get_stable_cutplane`), NOT recreated on every call. Recreating it
    fresh each time was confirmed to return a DIFFERENT, non-deterministic
    point set/count (305 vs 357 points for the identical setup) — which
    made unrelated eigenmodes' exported fields spuriously look ~95%
    "identical" to each other in one call ordering and completely
    uncorrelated in another, purely from mismatched point sampling, not
    any real difference in the underlying physics. Reusing the same
    CutPlane object gives an identical, stable point set across every
    mode/solnum query, making cross-mode comparisons actually valid.

    Writes the COMSOL export to a PERMANENT text file under
    results/mode_exports/ (not a deleted temp file) — same raw grid data
    COMSOL's own GUI export would give you, kept around so you can inspect
    or re-plot it yourself outside this pipeline.

    Returns (y_nm, z_nm, fields_dict) with fields as |value| (magnitude),
    matching this function's existing contract. All coords in nm.
    """
    EXP_TAG = "omc_tmp_exp"
    y_nm = z_nm = None
    fields = {k: None for k in field_exprs}
    out_path = None

    try:
        model.java.result().export().remove(EXP_TAG)
    except Exception:
        pass

    try:
        cp_tag = _get_stable_cutplane(model, dset)

        labels = list(field_exprs.keys())
        exprs  = list(field_exprs.values())

        export_dir = os.path.join(os.path.dirname(__file__), "..",
                                  "results", "mode_exports")
        os.makedirs(export_dir, exist_ok=True)
        study_slug = re.sub(r"[^A-Za-z0-9_-]+", "_", str(dset) or "study")
        kind = "abs" if abs_fields else "complex"
        fname = f"{study_slug}_solnum{solnum}_outer{outer}_{'_'.join(labels)}_{kind}.txt"
        out_path = os.path.join(export_dir, fname)

        exp = model.java.result().export().create(EXP_TAG, "Data")
        exp.set("data", cp_tag)
        exp.set("expr", exprs)
        exp.set("filename", out_path)
        if outer is not None:
            exp.set("outersolnum", str(int(outer)))   # must be str (int overload is ambiguous)
        exp.set("innerinput", "manual")
        exp.set("solnum", str(int(solnum)))
        exp.run()

        coords, vals = _parse_comsol_data_export(out_path, len(exprs))
        n = coords.shape[0]
        if n == 0:
            raise ValueError("CutPlane export returned 0 points")

        y_nm = coords[:, 1] * 1e9
        z_nm = coords[:, 2] * 1e9

        for label, v in zip(labels, vals):
            if len(v) != n:
                print(f"    [cutplane export size mismatch '{label}': {len(v)} vs {n}]")
                fields[label] = None
            else:
                fields[label] = np.abs(v) if abs_fields else v

    except Exception as e:
        print(f"    [cutplane eval failed: {e}]")
        y_nm = z_nm = None
    finally:
        try:
            model.java.result().export().remove(EXP_TAG)
        except Exception:
            pass

    return y_nm, z_nm, fields


def _eval_mode_on_midplane(model, dset, solnum, g, field_exprs, outer=None,
                            abs_fields=True):
    """Sample field expressions on mesh nodes near x=0, return (y_nm, z_nm, fields_dict).

    Evaluates x, y, z AND all field expressions in a SINGLE batch model.evaluate()
    call so COMSOL uses the same DOF grid for everything — no size mismatch.

    abs_fields: if True (default), apply np.abs() to each field array (good for
                amplitude heatmaps).  If False, return raw complex arrays so the
                caller can apply a phase correction before plotting arrows.
    """
    kw = {"inner": [solnum]}
    if dset:
        kw["dataset"] = dset
    if outer is not None and dset:
        kw["outer"] = outer

    labels  = list(field_exprs.keys())
    exprs   = list(field_exprs.values())

    all_exprs = ["x", "y", "z"] + exprs
    try:
        all_results = model.evaluate(all_exprs, **kw)
    except Exception as e:
        print(f"    [batch eval failed: {e}]")
        return None, None, {k: None for k in field_exprs}

    if not isinstance(all_results, (list, tuple)) or len(all_results) < 3:
        print("    [batch eval: unexpected result shape]")
        return None, None, {k: None for k in field_exprs}

    x_raw = np.real(np.asarray(all_results[0])).ravel()
    y_raw = np.real(np.asarray(all_results[1])).ravel()
    z_raw = np.real(np.asarray(all_results[2])).ravel()
    n     = len(x_raw)

    raw_fields = {}
    for i, label in enumerate(labels):
        try:
            arr = np.asarray(all_results[3 + i]).ravel()
            raw_fields[label] = np.abs(arr) if abs_fields else arr
        except Exception:
            raw_fields[label] = None

    for label, arr in raw_fields.items():
        if arr is not None and len(arr) != n:
            print(f"    [midplane: size mismatch '{label}': {len(arr)} vs {n}]")
            raw_fields[label] = None

    tol  = max(g.t * 0.10, 5e-9)
    mask = np.abs(x_raw) < tol
    if mask.sum() < 20:
        mask = np.abs(x_raw) < g.t * 0.30
    if mask.sum() < 10:
        print(f"    [midplane: only {mask.sum()} nodes near x=0]")
        return None, None, {k: None for k in field_exprs}

    y_nm = y_raw[mask] * 1e9
    z_nm = z_raw[mask] * 1e9

    masked_fields = {}
    for label, arr in raw_fields.items():
        if arr is None:
            masked_fields[label] = None
        else:
            sliced = arr[mask]
            if abs_fields:
                masked_fields[label] = np.where(np.isfinite(sliced), sliced, 0.0)
            else:
                masked_fields[label] = sliced  # keep complex; NaN left for caller

    return _merge_coincident_nodes(y_nm, z_nm, masked_fields)


def _merge_coincident_nodes(y_nm, z_nm, fields, bin_nm=3.0):
    """Average field values at coincident/near-coincident (y,z) mesh points.

    At every material interface (e.g. the elliptical hole boundary), COMSOL
    stores TWO mesh nodes at (numerically) the same location — one per
    adjoining domain — and the field's normal-ish component is genuinely
    discontinuous there (confirmed empirically: interface-node pairs <2nm
    apart differ by ~median 3.9e6 out of a ~3.5e7 total dynamic range for
    Ey, i.e. essentially double-valued). Feeding both conflicting values
    into `griddata` as separate points at ~the same spot corrupts the local
    interpolation, which is what produced the speckled/noisy TE mode plots
    (this does NOT affect the solid-mechanics case, which has no internal
    material interface — hence why that mode profile looked clean).
    Binning to a coarse (bin_nm) grid and averaging within each bin is a
    deliberately crude fix: it blurs the interface value rather than
    correctly modeling the jump, which is fine for a qualitative "does this
    look like a real mode" plot but would NOT be adequate for quantitative
    near-field analysis.
    """
    if y_nm is None or len(y_nm) == 0:
        return y_nm, z_nm, fields
    y_bin = np.round(y_nm / bin_nm).astype(np.int64)
    z_bin = np.round(z_nm / bin_nm).astype(np.int64)
    keys = y_bin.astype(np.int64) * 2_000_003 + z_bin  # unique-ish combined key
    _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    n_bins = len(counts)

    def _bin_average(vec):
        sums = np.zeros(n_bins, dtype=vec.dtype)
        np.add.at(sums, inverse, vec)
        return sums / counts

    y_avg = _bin_average(y_nm.astype(np.float64))
    z_avg = _bin_average(z_nm.astype(np.float64))
    fields_avg = {}
    for label, arr in fields.items():
        if arr is None:
            fields_avg[label] = None
        elif np.iscomplexobj(arr):
            fields_avg[label] = _bin_average(arr.real) + 1j * _bin_average(arr.imag)
        else:
            fields_avg[label] = _bin_average(arr.astype(np.float64))
    return y_avg, z_avg, fields_avg


def _plot_mode_profile(ax, y_nm, z_nm, mag, uy, uz, g_nm, title, cmap="magma",
                       uy_odd_y=False, mask_to_solid=False):
    """Render a 2D midplane mode profile as a smooth heatmap + displacement arrows.

    Scattered mesh-node data is interpolated onto a regular grid.  The half-model
    (y ≥ 0, PMC symmetry) is mirrored to show the full beam cross-section.
    The y-axis is extended by one beam half-width on each side to show the air clad.

    uy_odd_y: set True for breathing modes where uy is antisymmetric about y=0
              (top surface moves up, bottom moves down).  |mag| is always even.
    """
    from scipy.interpolate import griddata
    from matplotlib.patches import Rectangle, Ellipse

    if y_nm is None or mag is None:
        ax.text(0.5, 0.5, "Field not available\n(check COMSOL expressions)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="gray")
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        return

    # ── Mirror about y=0 ─────────────────────────────────────────────────────
    # |mag| is always even in y (same amplitude on both halves).
    # uy: odd for breathing (flipped sign), even for optical (no arrows anyway).
    # uz: even in y.
    y_half   = np.concatenate([y_nm,  -y_nm])
    z_half   = np.concatenate([z_nm,   z_nm])
    mag_half = np.concatenate([mag,    mag])
    uy_sign  = -1.0 if uy_odd_y else 1.0
    uy_half  = np.concatenate([uy, uy_sign * uy]) if uy is not None else None
    uz_half  = np.concatenate([uz, uz])            if uz is not None else None

    # ── Tile 3 unit cells along z (-a, 0, +a) ──────────────────────────────
    # At zone edge k=π/a:
    #   |field(z+a)| = |field(z)|   → amplitude tiles with the same sign
    #   u(z+a) = u(z)·exp(iπ) = -u(z) → displacement arrows flip sign each cell
    a = g_nm.a
    y_full   = np.concatenate([y_half,   y_half,   y_half])
    z_full   = np.concatenate([z_half - a, z_half, z_half + a])
    mag_full = np.concatenate([mag_half,  mag_half,  mag_half])
    # Adjacent cells carry opposite displacement sign (Bloch phase e^{iπ} = −1)
    uy_full  = np.concatenate([-uy_half, uy_half, -uy_half]) if uy_half is not None else None
    uz_full  = np.concatenate([-uz_half, uz_half, -uz_half]) if uz_half is not None else None

    # ── Grid bounds: 3 unit cells in z, full beam + one half-width air in y ──
    w = g_nm.w
    air_pad = w / 2
    y_min = -w / 2 - air_pad
    y_max =  w / 2 + air_pad
    z_min = -1.5 * a
    z_max =  1.5 * a

    nz, ny = 400, 300
    z_grid = np.linspace(z_min, z_max, nz)
    y_grid = np.linspace(y_min, y_max, ny)
    ZZ, YY = np.meshgrid(z_grid, y_grid)

    try:
        # method="linear" is justified now: the source data comes from an
        # explicit, dense 50x50 regular grid (via CutPoint3D — verified
        # exact to ~1% against COMSOL's own native evaluation), not
        # CutPlane's old sparse/scattered auto-sampling (~300-900 points,
        # non-deterministic point count between calls). With genuine dense
        # regular coverage, linear interpolation is filling in real
        # resolution, not fabricating it.
        mag_grid = griddata((z_full, y_full), mag_full.astype(float),
                            (ZZ, YY), method="linear", fill_value=0.0)

        if mask_to_solid:
            # "nearest" has no notion of "outside the sampled domain" the
            # way "linear"+fill_value=0 did — it always extrapolates to the
            # closest real sample, which incorrectly paints mechanical
            # displacement into the air/hole regions where no solid (and
            # hence no displacement field) exists. Zero those out using the
            # actual known geometry instead of relying on interpolation
            # side effects.
            hx, hy = g_nm.hx, g_nm.hy
            in_beam = np.abs(YY) <= w / 2
            in_hole = np.zeros_like(ZZ, dtype=bool)
            for i in (-1, 0, 1):
                in_hole |= (((ZZ - i * a) / hx) ** 2 + (YY / hy) ** 2) < 1.0
            mag_grid = np.where(in_beam & ~in_hole, mag_grid, 0.0)

        vmax = np.nanpercentile(mag_full, 98)
        pcm = ax.imshow(mag_grid, origin="lower", aspect="auto", cmap=cmap,
                        extent=[z_min, z_max, y_min, y_max],
                        vmin=0, vmax=vmax)
        plt.colorbar(pcm, ax=ax, label="|field| (a.u.)", pad=0.02)
    except Exception as e:
        print(f"    [griddata/imshow failed: {e}; falling back to scatter]")
        ax.scatter(z_full, y_full, c=mag_full, cmap=cmap, s=2)

    # In-plane displacement arrows: (uy, uz) = (v, w) projected onto yz cross-section.
    # uy_odd_y means breathing mode: uy points outward from y=0 on both sides.
    # uz is the z-component of displacement (along beam axis); uz is even in y.
    # Show arrows even if only uy is available (uz defaults to zero).
    if uy_full is not None:
        uy_r = np.real(uy_full)
        uz_r = np.real(uz_full) if uz_full is not None else np.zeros_like(uy_r)
        amp  = np.sqrt(uy_r**2 + uz_r**2)
        if amp.max() > 0:
            # Subsample to ~80 arrows across the full tiled domain
            threshold = 0.02 * amp.max()
            sig = amp > threshold
            if sig.sum() > 5:
                s = max(1, sig.sum() // 80)
                idx = np.where(sig)[0][::s]
                scale = amp[idx].max()
                ax.quiver(z_full[idx], y_full[idx],
                          uz_r[idx] / scale, uy_r[idx] / scale,
                          color="cyan", scale=10, width=0.004, alpha=0.85)

    # Beam + hole outlines for 3 unit cells
    for i in range(-1, 2):
        ax.add_patch(Rectangle((i * a - a/2, -w/2), a, w,
                                fill=False, edgecolor="white", lw=0.8, ls="--",
                                zorder=5))
        ax.add_patch(Ellipse((i * a, 0), 2 * g_nm.hx, 2 * g_nm.hy,
                              fill=False, edgecolor="white", lw=0.8, zorder=5))

    ax.set_xlabel("z [nm]  (beam axis)", fontsize=8)
    ax.set_ylabel("y [nm]  (width)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlim(z_min, z_max)
    ax.set_ylim(y_min, y_max)
    ax.tick_params(labelsize=7)


def _dump_all_modes(model, dset, ev, ev_solnum, g, study, outer, field_exprs, tag):
    """Export EVERY solved eigenmode's field profile plus a manifest of
    (solnum, frequency, exported file), for manual inspection.

    `ev_solnum[i]` is the RAW (unfiltered) 1-based solnum for `ev[i]` — do
    NOT assume "position i in ev" == "raw solnum i+1". After filtering out
    near-zero scattering-boundary artifacts (real modes, just not the ones
    we want — confirmed by cross-checking against COMSOL directly), the
    remaining frequencies' array positions no longer line up with COMSOL's
    own raw eigenmode numbering whenever any earlier raw entries got
    dropped. Passing the filtered position straight to COMSOL as "solnum"
    silently pulls the WRONG eigenmode's field — this was the actual
    root cause of the "speckled"/garbage mode profiles, not the extraction
    method. Always go through `ev_solnum` to get the correct raw index.

    NOTE: a near/far field-concentration filter was tried here to
    auto-detect "spurious" low-frequency modes, but a direct cross-check
    against dense (12473-node) COMSOL Multislice exports (done by hand)
    showed NO such sharp confined/unconfined split — the filter's earlier
    ~1.2 vs ~9.0 contrast was an artifact of undersampling (~305 points),
    not a real effect. Removed. This function is now purely a manual-
    inspection aid: dump every mode + frequency so a human can identify
    which one is physically real by comparing against COMSOL directly,
    with no automatic filtering.

    Returns the manifest path. Cheap: reuses the same per-mode CutPlane
    export `_eval_mode_cutplane` already does, just looped over every solnum.
    """
    try:
        export_dir = os.path.join(os.path.dirname(__file__), "..",
                                  "results", "mode_exports")
        os.makedirs(export_dir, exist_ok=True)
        study_slug = re.sub(r"[^A-Za-z0-9_-]+", "_", str(dset) or study)
        labels = list(field_exprs.keys())
        manifest_path = os.path.join(export_dir, f"{study_slug}_all_modes_manifest.txt")

        lines = ["# solnum  freq_Hz  freq_THz  exported_file"]
        for f_hz, solnum in zip(ev, ev_solnum):
            solnum = int(solnum)
            _eval_mode_grid(model, dset, solnum, g, field_exprs, outer=outer, ny=50, nz=50)
            fname = (f"{study_slug}_solnum{solnum}_outer{outer}_"
                     f"{'_'.join(labels)}_abs_grid50x50.txt")
            lines.append(f"{solnum}  {f_hz:.6e}  {f_hz*1e-12:.4f}  {fname}")
        with open(manifest_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"    [modes] dumped {len(ev)} modes -> {manifest_path}")
        return manifest_path
    except Exception as e:
        print(f"    [dump all modes failed: {e}]")
        return None


def _extract_opt_mode(model, g, study, gap_TE, solnum_override=None, dump_all=True):
    """Re-solve at zone edge and extract the TE dielectric band mode field.

    solnum_override: if given (1-based), force this eigenmode instead of the
    "highest mode below the gap" heuristic — use this once you've inspected
    the manifest from `dump_all` and identified which solnum is the real
    physical dielectric-band mode (spurious low-frequency scattering-BC
    modes can otherwise get picked by the heuristic).
    """
    model.parameter("kF", f"{np.pi / g.a}[1/m]")
    model.solve(study)
    dset = next((d for d in model.datasets() if d.startswith(study + '//')), None)

    # Determine which outer solution is the zone-edge re-solve (the last one).
    # The dataset accumulates one outer step per prior BZ sweep k-point plus
    # this new solve.  Using the wrong outer would mix all k-point field shapes.
    last_outer = _last_outer_idx(model, dset)

    # Pick solnum: highest mode below the gap lower edge (dielectric band).
    # Evaluate freq at the last outer step only to get the correct eigenvalues.
    #
    # IMPORTANT: COMSOL's "solnum" (the index the Export node needs) is the
    # RAW, unfiltered eigenmode position — but after a long sweep history,
    # some low RAW-solnum eigenvalues can be near-zero scattering-boundary
    # artifacts (confirmed real, not a solver bug — just modes we don't
    # want). Filtering those out with `ev > 10e12` shifts the array
    # positions, so "position i in the filtered array" no longer equals
    # "raw solnum i" whenever any earlier raw entries got dropped. Passing
    # that filtered position straight to `_eval_mode_grid` as if it were
    # the raw solnum silently extracts the WRONG (garbage) eigenmode's
    # field — this was the actual cause of the "speckled" mode profile,
    # not the extraction method itself. Track the raw index (`ev_solnum`)
    # alongside every filtered frequency so we always look up the correct
    # raw solnum before calling into COMSOL.
    ev = np.array([])
    ev_solnum = np.array([], dtype=int)
    try:
        freq_kw = ({"dataset": dset, "outer": last_outer} if (dset and last_outer)
                   else ({"dataset": dset} if dset else {}))
        raw_freq = np.real(np.asarray(model.evaluate("freq", **freq_kw)))
        order = np.argsort(raw_freq)             # ascending sort
        sorted_freq = raw_freq[order]
        sorted_solnum = order + 1                # 1-based RAW solnum per sorted entry
        mask = sorted_freq > 10e12                # optical range only
        ev = sorted_freq[mask]
        ev_solnum = sorted_solnum[mask]
        if dump_all and len(ev) > 0:
            _dump_all_modes(model, dset, ev, ev_solnum, g, study, last_outer,
                            {"Ey": "ewfd.Ey", "normE": "ewfd.normE"}, "opt")
        if solnum_override is not None:
            solnum = int(solnum_override)
        elif gap_TE.found:
            below = np.where(ev < gap_TE.f_lower)[0]
            solnum = int(ev_solnum[below[-1]]) if below.size > 0 else int(ev_solnum[0])
        else:
            solnum = int(ev_solnum[-1]) if len(ev_solnum) > 0 else 1   # highest available mode
    except Exception:
        solnum = 1

    # Plot the Ey component (the dominant field for this TE-like mode).
    # The earlier "speckled" appearance was NOT caused by plotting Ey vs
    # normE — it was the raw-solnum-vs-filtered-index bug now fixed above
    # (see project memory omc-gap-physics). With the correct mode selected,
    # Ey is smooth, matching COMSOL's own Ey export.
    # Look up f_hz via ev_solnum (raw index), NOT `ev[solnum-1]` — solnum is
    # now the RAW solnum, not a position in the filtered `ev` array.
    _pos = np.where(ev_solnum == solnum)[0]
    f_hz = float(ev[_pos[0]]) if len(_pos) > 0 else 0.0
    wl_nm = C0 / f_hz * 1e9 if f_hz > 0 else 0.0
    for Ey_expr, fallback_expr in (("ewfd.Ey", "ewfd.normE"),
                                    ("emw.Ey",  "emw.normE"),
                                    ("Ey",       "normE")):
        y_nm, z_nm, fields = _eval_mode_grid(
            model, dset, solnum, g,
            {"Ey": Ey_expr, "normE": fallback_expr},
            outer=last_outer, ny=50, nz=50)
        if y_nm is None:
            # 2nd attempt: mesh-node evaluation (with DOF ordering heuristic)
            y_nm, z_nm, fields = _eval_mode_on_midplane(
                model, dset, solnum, g,
                {"Ey": Ey_expr, "normE": fallback_expr},
                outer=last_outer)
        field_data = fields.get("Ey") if fields.get("Ey") is not None else fields.get("normE")
        if field_data is not None and y_nm is not None:
            # Keep ALL nodes including air — the evanescent optical field extends
            # beyond the beam walls, so we must not clip to beam width here.
            return y_nm, z_nm, field_data, wl_nm
    return None, None, None, wl_nm


def _extract_mech_mode(model, g, study, fy_sym, gap_mech=None, freqs_sym=None,
                        k_sym=None):
    """Re-solve at the gap's own k-point and extract the breathing mode field.

    If `gap_mech` is given (the fy-verified gap actually reported/shaded in
    the band plot), the displayed mode is picked from THAT gap's own edge
    bands — so the mode-profile figure always matches the gap it claims to
    illustrate, instead of independently picking whatever band has the
    globally highest fy (which can belong to a different, unrelated gap).

    `gap_mech.lower_band` is a direct column index into `fy_sym`/`freqs_sym`
    (band i; band i+1 is the other edge) — safe because the current
    `breathing_gap` sorts the already-ascending raw `freqs_sym` rows (a
    no-op) rather than scattering NaN into individual (k, band) points
    before sorting.

    IMPORTANT: the gap's f_lower/f_upper are extrema over the WHOLE swept
    k-range (`largest_gap`/`all_gaps` require max_k(band i) < min_k(band
    i+1) across ALL k, not specifically at the zone edge) — for many real
    geometries the gap closes somewhere in the BZ *interior*, not at k=pi/a.
    Always re-solving at the zone edge (as an earlier version did) can pick
    up a mode tens of GHz away from the actual gap if neither edge band is
    anywhere near it there. So: pick the edge band (i or i+1) with the
    higher MEAN fy over k, then re-solve at whichever already-swept k-point
    puts that band's frequency closest to the gap center — not necessarily
    the zone edge.
    """
    fy_mean = np.nanmean(fy_sym, axis=0)

    if (gap_mech is not None and gap_mech.found and freqs_sym is not None
            and k_sym is not None
            and 0 <= gap_mech.lower_band < fy_sym.shape[1] - 1):
        i = gap_mech.lower_band
        j = i if fy_mean[i] >= fy_mean[i + 1] else i + 1
        k_idx = int(np.nanargmin(np.abs(freqs_sym[:, j] - gap_mech.f_center)))
        k_frac = float(k_sym[k_idx])
        # Mechanical eigensolves have no scattering-boundary pseudo-modes
        # mixed into the raw solve order (unlike optical — see project
        # memory omc-gap-physics), so column position j == raw solnum here.
        solnum = j + 1
    else:
        k_frac = 1.0   # fallback: zone edge
        fy_edge = fy_sym[-1, :]
        fy_valid = np.where(np.isfinite(fy_edge), fy_edge, -1.0)
        solnum = int(np.argmax(fy_valid)) + 1 if fy_valid.max() > 0 else 1

    model.parameter("kF", f"{np.pi * k_frac / g.a}[1/m]")
    model.solve(study)
    dset = next((d for d in model.datasets() if d.startswith(study + '//')), None)

    # Determine last outer step (this re-solve) so we don't mix k-point data.
    last_outer = _last_outer_idx(model, dset)

    # Get frequency of this mode from the last outer step only.
    f_hz = 0.0
    try:
        freq_kw = ({"dataset": dset, "outer": last_outer} if (dset and last_outer)
                   else ({"dataset": dset} if dset else {}))
        ev = np.sort(np.real(np.asarray(model.evaluate("freq", **freq_kw))))
        f_hz = float(ev[solnum - 1]) if solnum <= len(ev) else 0.0
    except Exception:
        pass

    # Extract uy (heatmap = |v|) and complex displacements for phase-corrected arrows.
    exprs_to_try = [
        {"uy": "v", "uz": "w"},
        {"uy": "solid.uy", "uz": "solid.uz"},
    ]
    for exprs in exprs_to_try:
        # Amplitude heatmap (abs). Whichever method succeeds here (grid
        # vs. raw-mesh-node fallback) MUST also be used for the complex pass
        # below — mixing the two gives different point sets/order and the
        # arrow overlay silently drops (observed: cutplane 886 pts vs
        # midplane 1638 pts -> size-mismatch -> arrows None).
        y_nm, z_nm, fields = _eval_mode_grid(
            model, dset, solnum, g, exprs, outer=last_outer, ny=50, nz=50)
        used_grid = y_nm is not None
        if y_nm is None:
            y_nm, z_nm, fields = _eval_mode_on_midplane(
                model, dset, solnum, g, exprs, outer=last_outer)
        uy_data = fields.get("uy")   # |uy| for heatmap
        if uy_data is None or y_nm is None:
            continue

        # Restrict to beam nodes for heatmap (mech field is zero outside solid)
        w_nm = g.w * 1e9
        beam = np.abs(y_nm) <= w_nm / 2 * 1.1

        # ── Phase-corrected complex displacements for arrow directions ────────
        # COMSOL eigenfrequency modes have arbitrary complex phase.  We find the
        # global phase φ that maximises Re(Σ|uy_i| · uy_i), making the dominant
        # y-displacement real and positive — i.e. nodes near y=w/2 move outward.
        if used_grid:
            _, _, cfields = _eval_mode_grid(
                model, dset, solnum, g, exprs, outer=last_outer,
                abs_fields=False, ny=50, nz=50)
        else:
            _, _, cfields = _eval_mode_on_midplane(
                model, dset, solnum, g, exprs, outer=last_outer, abs_fields=False)
        uy_c = cfields.get("uy")
        uz_c = cfields.get("uz")

        uy_arr = uz_arr = None
        dbg_uy = "None" if uy_c is None else f"len={len(uy_c)} vs y_nm len={len(y_nm)}"
        if uy_c is not None and len(uy_c) == len(y_nm):
            uy_cb = uy_c[beam]
            # Use only finite (non-NaN) nodes for the phase reference; NaN appears
            # on nodes outside the solution domain (air) for solid mechanics fields.
            valid = np.isfinite(uy_cb)
            if valid.sum() > 5:
                uy_cb_v = uy_cb[valid]
                phi = np.angle(np.sum(uy_cb_v * np.abs(uy_cb_v)))  # amplitude-weighted
            else:
                phi = 0.0
            # Replace NaN with 0+0j so Re(…) stays finite everywhere
            def _clean(arr):
                return np.where(np.isfinite(arr), arr, 0.0 + 0.0j)
            uy_arr = np.real(_clean(uy_c[beam]) * np.exp(-1j * phi))
            dbg_uy = f"OK phi={np.degrees(phi):.1f}° valid={valid.sum()}/{len(uy_cb)} max={np.max(np.abs(uy_arr)):.3e}"
            if uz_c is not None and len(uz_c) == len(y_nm):
                uz_arr = np.real(_clean(uz_c[beam]) * np.exp(-1j * phi))
        print(f"    [mech arrows] uy_c: {dbg_uy}")
        print(f"    [mech arrows] uy_arr={'None' if uy_arr is None else f'max={np.max(np.abs(uy_arr)):.3e}'}  "
              f"uz_arr={'None' if uz_arr is None else f'max={np.max(np.abs(uz_arr)):.3e}'}")

        return (y_nm[beam], z_nm[beam], uy_data[beam], uy_arr, uz_arr, f_hz)

    return None, None, None, None, None, f_hz


def _run_characterization(model, best_rec, args, out_dir=None, no_cache=False):
    """Run full 4-study characterization and save a 2×3-panel figure.

    Top row   : geometry sketch | optical bands (TE+TM) | mechanical bands (sym+antisym)
    Bottom row: key numbers     | TE dielectric band mode |Ey| | breathing mode |u|

    All band sweeps use the full BZ (k: 0→1) with more k-points than the
    optimization pass for smooth dispersion curves.

    Band-structure and mode-profile data are cached to
    ``<out_dir>/best_char_data.npz``.  On a subsequent call with the same
    geometry (same ``u`` vector) the COMSOL sweeps are skipped entirely and
    the figure is regenerated from the cache.  Pass ``no_cache=True`` (or
    ``--no-cache`` in characterize_best.py) to force a fresh run.
    """
    if out_dir is None:
        out_dir = os.path.dirname(args.out_fig)
    out_path   = os.path.join(out_dir, "best_characterization.png")
    cache_path = os.path.join(out_dir, "best_char_data.npz")

    bounds = load_bounds()
    p = best_rec.get("params", {})
    # Use the stored PHYSICAL params directly, not u_to_geometry(best_rec["u"],
    # bounds) — u is normalized against whatever bounds.yaml said AT EVAL TIME,
    # so re-mapping it through the CURRENT bounds silently gives a different
    # geometry whenever bounds.yaml changes after the record was created
    # (e.g. widening t/w ranges). params are bounds-independent and exactly
    # what was actually simulated.
    g = Geometry(a=p["a"], w=p["w"], t=p.get("t", 0.22e-6), hx=p["hx"], hy=p["hy"])
    nm = 1e9

    print(f"\n{'='*60}")
    print("CHARACTERIZING BEST GEOMETRY (full 4-study sweep + mode profiles)")
    print(f"  a={p.get('a',0)*nm:.0f}  w={p.get('w',0)*nm:.0f}  "
          f"hx={p.get('hx',0)*nm:.0f}  hy={p.get('hy',0)*nm:.0f}  "
          f"t={p.get('t',0.22e-6)*nm:.0f} nm")
    print(f"  G_o={best_rec['G_o']*100:.1f}%  G_m={best_rec['G_m']*100:.1f}%  "
          f"score={best_rec['score']:+.4f}")
    print(f"{'='*60}")

    u_key = np.asarray(best_rec["u"], dtype=np.float64)

    # ── Try loading from NPZ cache ─────────────────────────────────────────────
    _use_cache = False
    if not no_cache and os.path.isfile(cache_path):
        try:
            c = np.load(cache_path, allow_pickle=False)
            if (c["u_key"].shape == u_key.shape and
                    np.allclose(c["u_key"], u_key, atol=1e-9)):
                print("  [cache] geometry matches — loading from NPZ, skipping COMSOL.")
                k_TE         = c["k_TE"]
                freqs_TE     = c["freqs_TE"]
                k_TM         = c["k_TM"]
                freqs_TM     = c["freqs_TM"]
                k_sym        = c["k_sym"]
                freqs_sym    = c["freqs_sym"]
                fy_sym       = _e2n(c["fy_sym"])
                k_anti       = c["k_anti"]
                freqs_anti   = c["freqs_anti"]
                k_TE_det     = c["k_TE_det"]
                freqs_TE_det = c["freqs_TE_det"]
                y_omode      = _e2n(c["y_omode"])
                z_omode      = _e2n(c["z_omode"])
                normE        = _e2n(c["normE"])
                wl_mode      = float(c["wl_mode"])
                y_mmode      = _e2n(c["y_mmode"])
                z_mmode      = _e2n(c["z_mmode"])
                mag_u        = _e2n(c["mag_u"])
                uy_arr       = _e2n(c["uy_arr"])
                uz_arr       = _e2n(c["uz_arr"])
                f_mmode      = float(c["f_mmode"])
                gap_TE       = _arr_to_gap(c["gap_TE_arr"])
                gap_mech     = _arr_to_gap(c["gap_mech_arr"])
                # Older caches (pre-fy-fix) won't have this key.
                gaps_mech_all = (_arr_to_gaps(c["gaps_mech_all_arr"])
                                 if "gaps_mech_all_arr" in c.files else [])
                # Sanity-check: cached G_o must be within 3 pp of the known value.
                # A large discrepancy means the cache was written after a corrupted
                # COMSOL run (e.g. wrong initial-guess modes after long optimization).
                G_o_cached = gap_TE.normalized_gap
                G_o_known  = best_rec.get("G_o", 0.0)
                if abs(G_o_cached - G_o_known) > 0.03:
                    print(f"  [cache] G_o mismatch: cached {G_o_cached*100:.1f}% "
                          f"vs record {G_o_known*100:.1f}% — invalidating cache.")
                else:
                    _use_cache = True
            else:
                print("  [cache] geometry changed — running fresh COMSOL sweeps.")
        except Exception as e:
            print(f"  [cache] load failed ({e}) — running fresh COMSOL sweeps.")

    # ── COMSOL sweeps (skipped when cache is valid) ────────────────────────────
    if not _use_cache:
        for name, val in g.as_dict().items():
            model.parameter(name, f"{val}[m]")

        N_K  = 15
        N_BO = 12   # bands for visualization (more than optimization for smooth curves)
        N_BM = 10

        print("  [1/4] opt TE  (full BZ, visualization)...")
        k_TE, freqs_TE = _full_opt_sweep(model, g, args.study_opt, n_k=N_K, n_bands=N_BO)

        print("  [2/4] opt TM  (full BZ, visualization)...")
        k_TM, freqs_TM = _full_opt_sweep(model, g, "opt TM", n_k=N_K, n_bands=N_BO)

        print("  [3/4] mech sym  (full BZ)...")
        k_sym, freqs_sym, fy_sym = _mech_sweep(model, g, args.study_mech,
                                                N_K, N_BM, compute_fy=True)

        print("  [4/4] mech antisym  (full BZ)...")
        k_anti, freqs_anti, _ = _mech_sweep(model, g, "mech antisym",
                                             N_K, N_BM, compute_fy=False)

        # Gap detection: re-run zone-edge-only sweep with EXACTLY the same
        # parameters as the optimization pass (same n_k, same n_bands_opt).
        print("  [gap] zone-edge sweep for gap detection (matching optimization)...")
        k_TE_det, freqs_TE_det = _opt_sweep(model, g, args.study_opt,
                                             args.n_k, args.n_bands_opt)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            gap_TE = optical_gap_tracked(freqs_TE_det)
            # fy_sym here is ALWAYS real (compute_fy=True above, regardless of
            # --no-fy at optimization time), so this is fy-VERIFIED: only bands
            # with breathing fraction >= FY_THRESH can define the gap. This can
            # legitimately differ from the fast --no-fy optimization-time G_m,
            # which trusts the "mech sym" BC alone without checking fy.
            mech_window = (MECH_F_MIN * 0.5, MECH_F_MAX * 2.0)
            gap_mech = breathing_gap(freqs_sym, fy_sym, f_min=MECH_F_MIN,
                                     f_max=MECH_F_MAX)
            # All raw gaps (unfiltered by fy) in the same window, for context —
            # shows whether other, non-breathing gaps exist that could be
            # mistaken for "the" mechanical gap.
            gaps_mech_all = all_gaps(freqs_sym, f_window=mech_window, nan_safe=True)

        # ── Mode profile extraction ──────────────────────────────────────────
        print("  [mode] optical TE profile at zone edge...")
        y_omode, z_omode, normE, wl_mode = _extract_opt_mode(
            model, g, args.study_opt, gap_TE,
            solnum_override=getattr(args, "opt_mode_solnum", None))

        print("  [mode] mechanical breathing mode profile at the gap's own k-point...")
        y_mmode, z_mmode, mag_u, uy_arr, uz_arr, f_mmode = _extract_mech_mode(
            model, g, args.study_mech, fy_sym, gap_mech=gap_mech, freqs_sym=freqs_sym,
            k_sym=k_sym)

        # ── Save to NPZ cache ────────────────────────────────────────────────
        try:
            os.makedirs(out_dir, exist_ok=True)
            np.savez(cache_path,
                     u_key       = u_key,
                     k_TE        = k_TE,
                     freqs_TE    = freqs_TE,
                     k_TM        = k_TM,
                     freqs_TM    = freqs_TM,
                     k_sym       = k_sym,
                     freqs_sym   = freqs_sym,
                     fy_sym      = _n2e(fy_sym),
                     k_anti      = k_anti,
                     freqs_anti  = freqs_anti,
                     k_TE_det    = k_TE_det,
                     freqs_TE_det= freqs_TE_det,
                     y_omode     = _n2e(y_omode),
                     z_omode     = _n2e(z_omode),
                     normE       = _n2e(normE),
                     wl_mode     = np.float64(wl_mode if wl_mode else 0.0),
                     y_mmode     = _n2e(y_mmode),
                     z_mmode     = _n2e(z_mmode),
                     mag_u       = _n2e(mag_u),
                     uy_arr      = _n2e(uy_arr),
                     uz_arr      = _n2e(uz_arr),
                     f_mmode     = np.float64(f_mmode if f_mmode else 0.0),
                     gap_TE_arr  = _gap_to_arr(gap_TE),
                     gap_mech_arr= _gap_to_arr(gap_mech),
                     gaps_mech_all_arr = _gaps_to_arr(gaps_mech_all))
            print(f"  [cache] saved → {cache_path}")
        except Exception as e:
            print(f"  [cache] save failed: {e}")

    # ── Build 2×3 figure ──────────────────────────────────────────────────────
    # Row 0: geometry | optical bands | mechanical bands
    # Row 1: summary  | TE mode |E|   | breathing mode |u|
    fig = plt.figure(figsize=(16, 9))
    gs  = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.36)
    ax_geo   = fig.add_subplot(gs[0, 0])
    ax_opt   = fig.add_subplot(gs[0, 1])
    ax_mech  = fig.add_subplot(gs[0, 2])
    ax_info  = fig.add_subplot(gs[1, 0])
    ax_omode = fig.add_subplot(gs[1, 1])
    ax_mmode = fig.add_subplot(gs[1, 2])

    _draw_geometry(ax_geo, g)
    _plot_opt_bands(ax_opt, k_TE, freqs_TE, k_TM, freqs_TM, gap_TE, a_m=g.a)
    _plot_mech_bands(ax_mech, k_sym, freqs_sym, fy_sym, k_anti, freqs_anti, gap_mech,
                     extra_gaps=gaps_mech_all)

    # Summary text panel (bottom-left)
    ax_info.axis("off")
    gap_o_str  = (f"{gap_TE.normalized_gap*100:.1f}%  "
                  f"({C0/gap_TE.f_center*1e9:.0f} nm)" if gap_TE.found else "none")
    gap_m_str  = (f"{gap_mech.normalized_gap*100:.1f}%  "
                  f"({gap_mech.f_center*1e-9:.2f} GHz)" if gap_mech.found else "none")
    lines = [
        "Best geometry summary",
        "",
        f"  a   = {g.a*nm:.0f} nm",
        f"  w   = {g.w*nm:.0f} nm",
        f"  t   = {g.t*nm:.0f} nm",
        f"  hx  = {g.hx*nm:.0f} nm",
        f"  hy  = {g.hy*nm:.0f} nm",
        f"  bridge   = {(g.a - 2*g.hx)*nm:.0f} nm",
        f"  sidewall = {(g.w/2 - g.hy)*nm:.0f} nm",
        "",
        f"  G_o   = {best_rec['G_o']*100:.1f}%   gap: {gap_o_str}",
        f"  G_m   = {best_rec['G_m']*100:.1f}%   gap: {gap_m_str}",
        f"  f_o   = {best_rec.get('f_o_nm', 0):.0f} nm",
        f"  f_m   = {best_rec.get('f_m_GHz', 0):.2f} GHz",
        "",
        f"  score = {best_rec['score']:+.4f}",
        f"  (target: 0.40)",
    ]
    ax_info.text(0.05, 0.95, "\n".join(lines), transform=ax_info.transAxes,
                 fontsize=8.5, va="top", family="monospace",
                 bbox=dict(boxstyle="round", facecolor="#f0f4ff", alpha=0.8))

    # Optical mode profile
    g_nm_scale = type('G', (), {
        'a': g.a*nm, 'w': g.w*nm, 'hx': g.hx*nm, 'hy': g.hy*nm})()
    opt_title = (f"TE dielectric band  |Ey|  "
                 f"(λ≈{wl_mode:.0f} nm)" if wl_mode > 0 else "TE mode  |Ey|")
    _plot_mode_profile(ax_omode, y_omode, z_omode, normE,
                       None, None, g_nm_scale, opt_title, cmap="inferno")

    # Mechanical mode profile
    print(f"  [arrows] uy_arr={'None' if uy_arr is None else f'shape={np.shape(uy_arr)} max={np.max(np.abs(uy_arr)):.3e}'}")
    print(f"  [arrows] uz_arr={'None' if uz_arr is None else f'shape={np.shape(uz_arr)} max={np.max(np.abs(uz_arr)):.3e}'}")
    mech_title = (f"Breathing mode  |uy|  "
                  f"({f_mmode*1e-9:.2f} GHz)" if f_mmode > 0 else "Breathing mode  |uy|")
    _plot_mode_profile(ax_mmode, y_mmode, z_mmode, mag_u,
                       uy_arr, uz_arr, g_nm_scale, mech_title, cmap="magma",
                       uy_odd_y=True, mask_to_solid=True)

    fig.suptitle(
        f"Best unit cell  |  score={best_rec['score']:+.4f}  "
        f"G_o={best_rec['G_o']*100:.1f}%  G_m={best_rec['G_m']*100:.1f}%  "
        f"f_o={best_rec.get('f_o_nm', 0):.0f} nm  "
        f"f_m={best_rec.get('f_m_GHz', 0):.2f} GHz",
        fontsize=11, y=0.995)

    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Composite -> {out_path}")

    # ── Three standalone figures ───────────────────────────────────────────────
    _save_geometry_fig(g,
                       os.path.join(out_dir, "best_geometry.png"))
    _save_optical_fig(k_TE, freqs_TE, k_TM, freqs_TM, gap_TE,
                      y_omode, z_omode, normE, wl_mode, g_nm_scale,
                      os.path.join(out_dir, "best_optical.png"))
    _save_mech_fig(k_sym, freqs_sym, fy_sym, k_anti, freqs_anti, gap_mech,
                   y_mmode, z_mmode, mag_u, uy_arr, uz_arr, f_mmode, g_nm_scale,
                   os.path.join(out_dir, "best_mechanical.png"),
                   extra_gaps=gaps_mech_all)

    return out_path


# ── Plot helpers ──────────────────────────────────────────────────────────────

def _scatter_gap_tw(ax, t_nm, w_nm, G, fy_mask, label, cmap, target=0.20,
                    required=True, beyond_lightcone=None):
    """Scatter G (G_o or G_m) vs. (t, w), marginalized over a/hx/hy.

    Points are optimizer-sampled, not on a grid, so no interpolation —
    each marker is one real evaluation from THIS run. Triangle = real fy,
    circle = --no-fy fallback (mixed only possible if resuming across runs
    with different --no-fy settings).

    required=False: this gap is NOT in the score function for this run
    (e.g. --no-require-mech). The target threshold is shown only as
    reference/context, clearly labeled "(not in cost function)" — it is
    NOT an active optimization target, so it must not look like one.

    beyond_lightcone: optional bool array (True = this gap's center sits
    above the light line at k=pi/a, i.e. NOT a real guided-mode gap). Those
    points are drawn faded (alpha=0.3) instead of full-opacity, so a large
    G_o that's actually a light-cone artifact doesn't visually read as a
    genuine result.
    """
    sc = None
    vmax = max(20.0, float(np.nanmax(G)) * 100) if len(G) else 20.0
    if beyond_lightcone is None:
        beyond_lightcone = np.zeros(len(G), dtype=bool)

    def _scat(mask, marker, size, edge, lw, alpha_full, lbl):
        nonlocal sc
        guided = mask & ~beyond_lightcone
        faded = mask & beyond_lightcone
        s_ = None
        if guided.any():
            s_ = ax.scatter(t_nm[guided], w_nm[guided], c=G[guided] * 100,
                            cmap=cmap, vmin=0, vmax=vmax, s=size, marker=marker,
                            edgecolors=edge, linewidths=lw, alpha=alpha_full, label=lbl)
        if faded.any():
            s2 = ax.scatter(t_nm[faded], w_nm[faded], c=G[faded] * 100,
                            cmap=cmap, vmin=0, vmax=vmax, s=size, marker=marker,
                            edgecolors=edge, linewidths=lw, alpha=0.3,
                            label=None if guided.any() else lbl)
            s_ = s_ if s_ is not None else s2
        if s_ is not None:
            sc = s_

    _scat(~fy_mask, "o", 40, "none", 0, 0.85, "no-fy (mech sym trusted)")
    _scat(fy_mask, "^", 55, "k", 0.5, 0.95, "real fy (breathing-verified)")

    if beyond_lightcone.any():
        ax.scatter([], [], marker="o", facecolors="gray", alpha=0.3,
                  label="faded = gap center above light line (not guided)")

    hit = G >= target
    hit_guided = hit & ~beyond_lightcone
    hit_faded = hit & beyond_lightcone
    tag = f"≥{target*100:.0f}% target" if required else \
          f"≥{target*100:.0f}% (ref. only — not in cost function)"
    if hit_guided.any():
        ax.scatter(t_nm[hit_guided], w_nm[hit_guided], facecolors="none",
                   edgecolors="red" if required else "gray",
                   linewidths=1.3, s=110, alpha=1.0, label=tag)
    if hit_faded.any():
        ax.scatter(t_nm[hit_faded], w_nm[hit_faded], facecolors="none",
                   edgecolors="red" if required else "gray",
                   linewidths=1.3, s=110, alpha=0.3,
                   label=None if hit_guided.any() else tag)
    if sc is not None:
        plt.colorbar(sc, ax=ax, label=f"{label} [%]")
    ax.set_xlabel("t  (thickness) [nm]")
    ax.set_ylabel("w  (width) [nm]")
    title = f"{label} vs. (t, w)  (marginalized over a, hx, hy)"
    if not required:
        title += "\n(NOT part of this run's score)"
    ax.set_title(title)
    ax.legend(fontsize=7, loc="best")
    ax.grid(alpha=0.3)


def _plot_progress(results, out_path, require_opt=True, require_mech=True,
                   g_min_opt=G_MIN_OPT, g_min_mech=G_MIN_MECH):
    """Save a 4-panel figure: score history, Pareto scatter, and G_o/G_m
    vs. (t, w) scatters (marginalized over a/hx/hy) for live geometry insight.

    Score range is dominated by infeasible results (score=-10) so we clip axes
    to the 5th-95th percentile of finite scores to make the interesting region
    visible. The running-best line (non-decreasing) shows actual progress.

    require_opt/require_mech/g_min_opt/g_min_mech: THIS run's actual scoring
    config (not necessarily the global G_MIN_OPT/G_MIN_MECH defaults) — reference
    lines/labels must reflect what's actually being optimized for, or a
    --no-require-mech run would misleadingly show a live 20% mechanical target.
    """
    all_scores = [r["score"] for r in results]
    G_o = [r["G_o"] * 100 for r in results]
    G_m = [r["G_m"] * 100 for r in results]

    # Filter by STATUS, not a numeric score threshold — "infeasible" (-10.0)
    # and "failed" (-5.0, a COMSOL exception) are both sentinel non-evaluations,
    # not real scores. An earlier `s > -9` threshold caught -10.0 but let -5.0
    # slip through, which dragged the 5th percentile down to exactly -5.0
    # whenever a few failed evals were present, clipping the y-axis/colorbar
    # to a near-empty range and squishing all the real scores into a sliver.
    finite_scores = [r["score"] for r in results
                     if r.get("status") == "success" and np.isfinite(r["score"])]
    if finite_scores:
        p5   = max(np.percentile(finite_scores, 5), -3.0)  # clip bad tail
        ymax = max(max(finite_scores) * 1.1, 0.05)         # always show best
    else:
        p5, ymax = -2.0, 0.5

    running_best = []
    curr = -np.inf
    for s in all_scores:
        curr = max(curr, s)
        running_best.append(curr)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9.6))

    # ── Top-left: score history ───────────────────────────────────────────────
    ax = axes[0, 0]
    iters = np.arange(len(all_scores))
    clipped = np.clip(all_scores, p5 - 0.1, ymax + 0.1)
    norm = plt.Normalize(vmin=p5, vmax=ymax)
    ax.scatter(iters, clipped, c=clipped, cmap="RdYlGn", norm=norm,
               s=22, alpha=0.65, zorder=2, label="Individual eval")
    ax.plot(iters, np.clip(running_best, p5 - 0.1, ymax + 0.1),
            lw=2.2, color="#1a7a3c", zorder=3, label="Running best")
    # "Practical target" depends on which gaps are actually in the cost
    # function — a hardcoded 0.40 (from the original G_MIN_OPT=G_MIN_MECH=
    # 0.20 co-optimized formula) is meaningless once thresholds change or a
    # gap is excluded, and doubly so once the thickness/width reward (only
    # active when require_mech=False) shifts the achievable score range.
    if require_opt and require_mech:
        practical_target = g_min_opt + g_min_mech
        target_label = f"Target ({practical_target:.2f} = g_min_opt+g_min_mech)"
    elif require_opt and not require_mech:
        # No mech term; thinness/narrowness rewards are unbounded improvers
        # (no natural "done" point), so treat "halfway there" on each as the
        # practical goalpost rather than their unreachable max of 1.0 each.
        practical_target = g_min_opt + 0.5 * (LAMBDA_THICK + LAMBDA_WIDTH)
        target_label = f"Target ({practical_target:.2f} = g_min_opt + half geom. reward)"
    elif require_mech and not require_opt:
        practical_target = g_min_mech
        target_label = f"Target ({practical_target:.2f} = g_min_mech)"
    else:
        practical_target = None
    if practical_target is not None:
        ax.axhline(practical_target, ls="--", color="#aa0000", lw=1.0, label=target_label)
    ax.axhline(0.0,  ls=":",  color="gray",    lw=0.7)
    ax.set_ylim(p5, ymax)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Score")
    ax.set_title("Score history  (y-axis: low tail clipped, best always visible)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # ── Top-right: Pareto scatter (color = clipped score) ─────────────────────
    ax = axes[0, 1]
    sc_arr = np.array(all_scores)
    cb_lo = max(np.percentile(finite_scores, 10) if finite_scores else -2.0, -2.0)
    cb_hi = max(max(finite_scores) if finite_scores else 0.5,
                cb_lo + 0.01)   # clip bad tail only; best results always visible
    sc = ax.scatter(G_o, G_m,
                    c=np.clip(sc_arr, cb_lo, cb_hi),
                    cmap="RdYlGn", vmin=cb_lo, vmax=cb_hi,
                    s=50, alpha=0.8, zorder=3, edgecolors="none")
    plt.colorbar(sc, ax=ax,
                 label=f"Score  (lower tail clipped at {cb_lo:.2f})")
    opt_label = (f"G_o target ({g_min_opt*100:.0f}%)" if require_opt
                else f"G_o ref. ({g_min_opt*100:.0f}%, not scored)")
    mech_label = (f"G_m target ({g_min_mech*100:.0f}%)" if require_mech
                 else f"G_m ref. ({g_min_mech*100:.0f}%, not scored)")
    ax.axvline(g_min_opt * 100, ls="--" if require_opt else ":",
              color="#cc3333" if require_opt else "#999999", lw=1.0, label=opt_label)
    ax.axhline(g_min_mech * 100, ls="--" if require_mech else ":",
              color="#1f5fa6" if require_mech else "#999999", lw=1.0, label=mech_label)
    ax.set_xlabel("Optical gap  G_o [%]")
    ax.set_ylabel("Mechanical breathing gap  G_m [%]")
    ax.set_title("Optical vs mechanical gap  (color = score)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # ── Bottom row: G_o / G_m vs (t, w), marginalized over a/hx/hy ────────────
    succ = [r for r in results if r.get("status") == "success" and r.get("params")]
    if succ:
        t_nm = np.array([r["params"]["t"] * 1e9 for r in succ])
        w_nm = np.array([r["params"]["w"] * 1e9 for r in succ])
        a_m  = np.array([r["params"]["a"] for r in succ])
        f_o_c = np.array([r.get("optical_center_frequency", 0.0) or 0.0 for r in succ])
        G_o_s = np.array([r.get("G_o", 0.0) for r in succ])
        G_m_s = np.array([r.get("G_m", 0.0) for r in succ])
        fy_mask = np.array([bool(r.get("fy_computed")) for r in succ])
        # Light line at k=pi/a, where the gap is defined — same criterion as
        # score_result()'s light-cone penalty (see LAMBDA_LIGHTCONE).
        light_line = C0 / (2.0 * a_m)
        beyond_lc = f_o_c > light_line
        _scatter_gap_tw(axes[1, 0], t_nm, w_nm, G_o_s, fy_mask,
                        "Optical gap G_o", "viridis", target=g_min_opt,
                        required=require_opt, beyond_lightcone=beyond_lc)
        _scatter_gap_tw(axes[1, 1], t_nm, w_nm, G_m_s, fy_mask,
                        "Mechanical gap G_m", "magma", target=g_min_mech,
                        required=require_mech)
    else:
        for ax in (axes[1, 0], axes[1, 1]):
            ax.set_axis_off()

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# ── Progress display ──────────────────────────────────────────────────────────

def _fmt_time(seconds):
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    if h > 0:
        return f"{h}h{m:02d}m"
    if m > 0:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def _progress_bar(done, total, width=24):
    frac = done / max(total, 1)
    filled = int(frac * width)
    bar = "█" * filled + "░" * (width - filled)
    return f"[{bar}] {done}/{total} ({frac*100:.0f}%)"


def _progress_status(it, n_iter, t_start, n_eval, best_rec, recent):
    """Return a dict with all progress fields (used by both terminal and file)."""
    elapsed = time.time() - t_start
    rate_per_min = n_eval / elapsed * 60 if elapsed > 0 else 0   # evals/min
    sec_per_eval = elapsed / n_eval if n_eval > 0 else float("inf")
    remaining = (n_iter - it - 1) * sec_per_eval
    return dict(
        it=it, n_iter=n_iter, n_eval=n_eval,
        elapsed=elapsed, rate_per_min=rate_per_min,
        sec_per_eval=sec_per_eval, remaining=remaining,
        best_rec=best_rec, recent=recent,
    )


def _print_progress(st):
    """Print one-line progress to the terminal."""
    bar = _progress_bar(st["it"] + 1, st["n_iter"])
    eta = f"~{_fmt_time(st['remaining'])} left" if st["remaining"] < 1e8 else "estimating..."
    best_str = "–"
    if st["best_rec"]:
        r = st["best_rec"]
        p = r.get("params", {})
        best_str = (f"G_o={r['G_o']*100:.1f}% G_m={r['G_m']*100:.1f}% "
                    f"a={p.get('a',0)*1e9:.0f} t={p.get('t',0)*1e9:.0f}nm "
                    f"f_o={r.get('f_o_nm',0):.0f}nm")
    rate_str = f"{st['rate_per_min']:.2f} evals/min" if st["rate_per_min"] > 0 else "–"
    print(f"  ↳ {bar}  elapsed {_fmt_time(st['elapsed'])}  {eta}  {rate_str}  │  best: {best_str}")


def _save_progress(st, path, out_json, require_opt=True, require_mech=True,
                   g_min_opt=G_MIN_OPT, g_min_mech=G_MIN_MECH):
    """Write a human-readable progress snapshot to *path* (overwrite each time)."""
    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = []
    lines.append(f"=== OMC Optimization Progress ===")
    lines.append(f"Updated : {now}")
    lines.append(f"Results : {out_json}")
    lines.append("")

    bar = _progress_bar(st["it"] + 1, st["n_iter"], width=30)
    pct = (st["it"] + 1) / st["n_iter"] * 100
    eta = f"~{_fmt_time(st['remaining'])}" if st["remaining"] < 1e8 else "?"
    lines.append(f"Progress: {bar}")
    lines.append(f"Elapsed : {_fmt_time(st['elapsed'])}  "
                 f"Rate: {st['rate_per_min']:.2f} evals/min  ETA: {eta}")
    lines.append("")

    if st["best_rec"]:
        r = st["best_rec"]
        p = r.get("params", {})
        lines.append("Best this run:")
        lines.append(f"  G_o = {r['G_o']*100:.1f}%   G_m = {r['G_m']*100:.1f}%   "
                     f"score = {r.get('score',0):+.4f}")
        lines.append(f"  a={p.get('a',0)*1e9:.0f}nm  w={p.get('w',0)*1e9:.0f}nm  "
                     f"hx={p.get('hx',0)*1e9:.0f}nm  hy={p.get('hy',0)*1e9:.0f}nm  "
                     f"t={p.get('t',0.22e-6)*1e9:.0f}nm")
        lines.append(f"  f_o = {r.get('f_o_nm',0):.0f} nm   "
                     f"f_m = {r.get('f_m_GHz',0):.2f} GHz")
    else:
        lines.append("Best this run: (none yet)")

    lines.append("")
    lines.append("Recent evaluations:")
    hdr = (f"  {'iter':>4}  {'a':>5} {'w':>5} {'t':>5} {'hx':>5} {'hy':>5}  "
          f"{'G_o':>6} {'G_m':>6}  {'f_o':>7} {'f_m':>8}  score")
    lines.append(hdr)
    lines.append("  " + "-" * (len(hdr) - 2))
    for rec in reversed(st["recent"]):
        p = rec.get("params", {})
        idx = rec.get("_it", "?")
        tag = "  *** BEST" if rec is st["best_rec"] else ""
        lines.append(
            f"  [{idx:>3}]  {p.get('a',0)*1e9:5.0f} {p.get('w',0)*1e9:5.0f} "
            f"{p.get('t',0.22e-6)*1e9:5.0f} "
            f"{p.get('hx',0)*1e9:5.0f} {p.get('hy',0)*1e9:5.0f}  "
            f"{rec.get('G_o',0)*100:5.1f}% {rec.get('G_m',0)*100:5.1f}%  "
            f"{rec.get('f_o_nm',0):7.0f} {rec.get('f_m_GHz',0):7.2f}GHz  "
            f"{rec.get('score',0):+.4f}{tag}"
        )

    lines.append("")
    if require_opt and require_mech:
        lines.append(f"Targets: G_o ≥ {g_min_opt*100:.0f}%   G_m ≥ {g_min_mech*100:.0f}%   "
                     "f_o ≈ 1550 nm   f_m ∈ 5-10 GHz")
    elif require_opt and not require_mech:
        lines.append(f"Targets: G_o ≥ {g_min_opt*100:.0f}%   f_o ≈ 1550 nm   "
                     f"(G_m not in cost function — thinness/narrowness reward "
                     f"active instead, t≤{T_MIN_BOUND*1e9:.0f}-{T_MAX_BOUND*1e9:.0f}nm "
                     f"w≤{W_MIN_BOUND*1e9:.0f}-{W_MAX_BOUND*1e9:.0f}nm range)")
    elif require_mech and not require_opt:
        lines.append(f"Targets: G_m ≥ {g_min_mech*100:.0f}%   f_m ∈ 5-10 GHz   "
                     "(G_o not in cost function)")
    else:
        lines.append("Targets: neither gap is in the cost function")
    lines.append("=" * 36)

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Co-optimize optical + mechanical OMC bandgaps via COMSOL.")
    ap.add_argument("--n-init",  type=int, default=6,
                    help="Space-filling / Halton trials before TPE kicks in")
    ap.add_argument("--n-iter",  type=int, default=14,
                    help="Total iterations (including n-init)")
    ap.add_argument("--n-k",     type=int, default=7,
                    help="k-points for optical zone-edge sweep (k∈[0.5,1.0]×π/a) "
                         "and full mechanical sweep; 7 is a good default")
    ap.add_argument("--n-bands-mech", type=int, default=8)
    ap.add_argument("--n-bands-opt",  type=int, default=6,
                    help="Optical bands kept per k-point during the OPTIMIZATION "
                         "loop (_opt_sweep; neigs floor is also 6 there). Only "
                         "the fundamental TE gap is scored, plus margin for the "
                         "1-3 near-zero scattering-boundary pseudo-modes that "
                         "get filtered out — 6 is enough; 10 was excessive. "
                         "The end-of-run/standalone characterization sweep "
                         "(_full_opt_sweep) keeps its own wider neigs floor of "
                         "10 for accurate plotting, independent of this flag.")
    ap.add_argument("--study-mech",   default="mech sym")
    ap.add_argument("--study-opt",    default="opt TE")
    ap.add_argument("--no-fy", action="store_true",
                    help="Skip per-mode fy field evaluation in mechanical sweep. "
                         "Use when mech sym BCs already select breathing modes — "
                         "saves ~30-50%% of total COMSOL time per evaluation.")
    ap.add_argument("--out-fig",      default="results/figures/opt_progress_iso.png")
    ap.add_argument("--out-json",     default="results/opt_results_t5d_iso.json")
    ap.add_argument("--progress-file", default="results/progress_iso.txt",
                    help="Human-readable status file updated after each iteration. "
                         "Monitor live with:  tail -f results/progress_iso.txt")
    ap.add_argument("--explore-every", type=int, default=0,
                    help="Force a random exploration sample every N iterations "
                         "(0=off). Prevents TPE from getting stuck in a local "
                         "optimum. Example: --explore-every 8 gives ~12%% random "
                         "samples in a 60-iteration run.")
    ap.add_argument("--resume", action="store_true",
                    help="(Legacy flag, now a no-op: prior results are always loaded "
                         "automatically if --out-json exists.)")
    ap.add_argument("--fresh-start", action="store_true",
                    help="Discard all prior results and start from scratch. "
                         "The existing --out-json is backed up to <file>.bak first. "
                         "Without this flag, prior results are always loaded automatically.")
    ap.add_argument("--skip-characterization", action="store_true",
                    help="Skip the end-of-run full characterization sweep "
                         "(opt TE+TM, mech sym+antisym over the full BZ). "
                         "Use for quick test runs where you don't need the figure.")
    ap.add_argument("--no-require-mech", action="store_true",
                    help="Exclude mechanical gap from the score entirely: no +G_m "
                         "reward, no gap-size penalty, no frequency penalty. The "
                         "mechanical sweep still runs and G_m/f_m are still "
                         "recorded (so you can inspect 'free' mechanical gaps "
                         "later) — this only changes what's optimized for.")
    ap.add_argument("--no-require-opt", action="store_true",
                    help="Symmetric to --no-require-mech, for optical gap. "
                         "Rarely used but supported for consistency.")
    ap.add_argument("--g-min-opt", type=float, default=G_MIN_OPT,
                    help=f"Minimum acceptable optical gap for scoring "
                         f"(default {G_MIN_OPT}). Ignored if --no-require-opt.")
    ap.add_argument("--g-min-mech", type=float, default=G_MIN_MECH,
                    help=f"Minimum acceptable mechanical gap for scoring "
                         f"(default {G_MIN_MECH}). Ignored if --no-require-mech.")
    args = ap.parse_args()
    require_opt  = not args.no_require_opt
    require_mech = not args.no_require_mech

    if not os.path.exists(_TEMPLATE):
        print(f"[FAIL] COMSOL template not found: {_TEMPLATE}")
        sys.exit(1)

    bounds = load_bounds()
    model  = get_model(_TEMPLATE)
    avail  = set(model.studies())
    print(f"Available studies: {sorted(avail)}")
    for s in [args.study_mech, args.study_opt]:
        if s not in avail:
            print(f"[FAIL] study '{s}' not in template — check --study-mech / --study-opt")
            sys.exit(1)

    prior = []
    if os.path.exists(args.out_json):
        if args.fresh_start:
            import shutil
            bak = args.out_json + ".bak"
            shutil.copy2(args.out_json, bak)
            print(f"Fresh start: backed up {len(json.load(open(bak)))} prior results "
                  f"to {bak}")
        else:
            with open(args.out_json) as fh:
                prior = json.load(fh)
            print(f"Auto-resumed: loaded {len(prior)} prior results from {args.out_json}")

    compute_fy = not args.no_fy
    optimizer = make_optimizer(args.n_init, prior_results=prior, bounds=bounds,
                               require_opt=require_opt, require_mech=require_mech,
                               g_min_opt=args.g_min_opt, g_min_mech=args.g_min_mech)
    print(f"Optimizer: {optimizer[0]}  |  n_init={args.n_init}  "
          f"n_iter={args.n_iter}  n_k={args.n_k}  compute_fy={compute_fy}")
    print(f"Scoring: require_opt={require_opt} (g_min={args.g_min_opt})  "
          f"require_mech={require_mech} (g_min={args.g_min_mech})")

    results = list(prior)   # carry over prior results for progress plot
    # Track best among NEW evaluations only (prior scores may be stale/fake)
    best_score = -np.inf
    best_rec   = None
    n_eval     = 0          # successful (non-infeasible) evaluations this run
    recent     = []         # rolling window of last 8 evaluated (not infeasible) records
    t_start    = time.time()
    print(f"Progress file: {args.progress_file}  (tail -f to monitor)")


    for it in range(args.n_iter):
        # ── Epsilon-greedy exploration ─────────────────────────────────────────
        # Every --explore-every iterations, enqueue a uniformly random point so
        # TPE is forced to evaluate a region far from its current optimum.
        # enqueue_trial() pushes the point into Optuna's queue; the very next
        # opt_ask() will return it, and opt_tell() will update the model normally.
        if (args.explore_every > 0 and it > 0 and it % args.explore_every == 0
                and optimizer[0] == "optuna"):
            explore_rng = np.random.default_rng(it * 997 + 13)
            random_params = {f"u{i}": float(explore_rng.random()) for i in range(5)}
            optimizer[1].enqueue_trial(random_params)
            print(f"[{it:3d}] EXPLORE: enqueued random point "
                  f"(every {args.explore_every} iters)")

        u, trial = opt_ask(optimizer)
        g = u_to_geometry(u, bounds)
        ok, reasons = check_feasibility(g, bounds)

        if not ok:
            rec = dict(id=_make_id(g.as_dict()), u=u, score=-10.0, G_o=0.0, G_m=0.0,
                       status="infeasible", reasons=reasons,
                       optical_backend="comsol", mech_backend="comsol",
                       optical_gap=0.0, mechanical_gap=0.0,
                       optical_center_frequency=0.0,
                       mechanical_center_frequency=0.0,
                       params=g.as_dict())
            opt_tell(optimizer, u, -10.0, trial)
            save_result(rec)
            results.append(rec)
            print(f"[{it:3d}] INFEASIBLE  {reasons}")
            continue

        print(f"\n[{it:3d}] u={[round(x,3) for x in u]}")
        print(f"       a={g.a*1e9:.0f}  w={g.w*1e9:.0f}  "
              f"hx={g.hx*1e9:.0f}  hy={g.hy*1e9:.0f}  t={g.t*1e9:.0f} nm")

        try:
            rec = evaluate(model, u, g,
                           n_k=args.n_k,
                           n_bands_mech=args.n_bands_mech,
                           n_bands_opt=args.n_bands_opt,
                           study_mech=args.study_mech,
                           study_opt=args.study_opt,
                           compute_fy=compute_fy,
                           require_opt=require_opt, require_mech=require_mech,
                           g_min_opt=args.g_min_opt, g_min_mech=args.g_min_mech)
        except Exception as e:
            rec = dict(id=_make_id(g.as_dict()), u=u, score=-5.0, G_o=0.0, G_m=0.0,
                       status="failed", error=str(e),
                       optical_backend="comsol", mech_backend="comsol",
                       optical_gap=0.0, mechanical_gap=0.0,
                       optical_center_frequency=0.0,
                       mechanical_center_frequency=0.0,
                       params=g.as_dict())
            print(f"       [FAIL] {e}")

        sc = rec["score"]
        opt_tell(optimizer, u, sc, trial)
        save_result(rec)
        results.append(rec)
        n_eval += 1

        is_best = sc > best_score
        if is_best:
            best_score = sc
            best_rec = rec

        rec["_it"] = it   # store iteration index for progress file display
        recent.append(rec)
        if len(recent) > 8:
            recent.pop(0)

        print(f"       G_o={rec['G_o']*100:.1f}%  "
              f"G_m={rec['G_m']*100:.1f}%  "
              f"f_o={rec.get('f_o_nm', float('nan')):.0f} nm  "
              f"f_m={rec.get('f_m_GHz', 0.0):.2f} GHz  "
              f"score={sc:+.4f}  "
              + ("*** BEST" if is_best else ""))

        st = _progress_status(it, args.n_iter, t_start, n_eval, best_rec, recent)
        _print_progress(st)
        _save_progress(st, args.progress_file, args.out_json, require_opt=require_opt,
                       require_mech=require_mech, g_min_opt=args.g_min_opt,
                       g_min_mech=args.g_min_mech)

        # Save JSON + progress plot after each iteration
        _save_json(results, args.out_json)
        _plot_progress(results, args.out_fig, require_opt=require_opt,
                       require_mech=require_mech, g_min_opt=args.g_min_opt,
                       g_min_mech=args.g_min_mech)

    # ── Final report ─────────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start
    print("\n" + "=" * 60)
    print("OPTIMIZATION COMPLETE")
    print("=" * 60)
    print(f"Total wall time : {_fmt_time(elapsed_total)}")
    print(f"Evaluations     : {n_eval} successful  ({args.n_iter} iterations total)")
    # Report best new result found this run
    if best_rec:
        best = best_rec
        print(f"\nBest NEW result (this run):")
        print(f"  Score      : {best['score']:+.4f}")
        print(f"  Optical gap: {best['G_o']*100:.1f}%  center={best.get('f_o_nm', float('nan')):.0f} nm")
        print(f"  Mech gap   : {best['G_m']*100:.1f}%  center={best.get('f_m_GHz', 0.0):.2f} GHz")
        print(f"  Params (nm): a={best['params']['a']*1e9:.0f}  w={best['params']['w']*1e9:.0f}  "
              f"hx={best['params']['hx']*1e9:.0f}  hy={best['params']['hy']*1e9:.0f}  "
              f"t={best['params'].get('t',0.22e-6)*1e9:.0f}")
        print(f"  u vector   : {[round(x,4) for x in best['u']]}")
    # Also show overall best across all runs
    all_succ = [r for r in results if r.get("G_o", 0) > 0 or r.get("G_m", 0) > 0]
    if all_succ:
        best_all = max(all_succ, key=lambda r: r.get("score", -99))
        if best_rec is None or best_all["id"] != best_rec["id"]:
            print(f"\nOverall best across all runs:")
            print(f"  Score={best_all['score']:+.4f}  "
                  f"G_o={best_all['G_o']*100:.1f}%  G_m={best_all['G_m']*100:.1f}%  "
                  f"a={best_all['params']['a']*1e9:.0f}  "
                  f"t={best_all['params'].get('t',0.22e-6)*1e9:.0f}nm")

    print(f"\nFigure -> {args.out_fig}")
    print(f"JSON   -> {args.out_json}")

    # ── End-of-run characterization (fresh subprocess → clean COMSOL state) ──
    # characterize_best.py is launched as a separate process so it gets its own
    # COMSOL instance, completely isolated from the optimizer's model state.
    out_dir  = os.path.dirname(args.out_fig)
    char_log = os.path.join(out_dir, "char_log.txt")
    if not args.skip_characterization:
        print(f"\n  [char] Launching characterize_best.py in a fresh process...")
        print(f"         Monitor: tail -f {char_log}")
        cmd = [sys.executable, _CHAR_SCRIPT,
               "--json",       args.out_json,
               "--out-dir",    out_dir,
               "--study-opt",  args.study_opt,
               "--study-mech", args.study_mech]
        if getattr(args, "no_fy", False):
            cmd.append("--no-fy")
        with open(char_log, "w") as char_log_fh:
            proc = subprocess.run(cmd, stdout=char_log_fh, stderr=subprocess.STDOUT)
        if proc.returncode == 0:
            print(f"  [char] Done → {os.path.join(out_dir, 'best_characterization.png')}")
        else:
            print(f"  [char] Characterization failed (exit code {proc.returncode})")
            print(f"         Check log: {char_log}")
            print(f"         Re-run manually: python scripts/characterize_best.py")
    else:
        print("(Characterization skipped — run  python scripts/characterize_best.py  when ready)")


def _save_json(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # Strip large band arrays from JSON to keep it readable
    slim = []
    for r in results:
        s = {k: v for k, v in r.items() if not k.startswith("_")}
        slim.append(s)
    with open(path, "w") as fh:
        json.dump(slim, fh, indent=2, default=lambda x: float(x) if hasattr(x, "__float__") else str(x))


if __name__ == "__main__":
    main()
