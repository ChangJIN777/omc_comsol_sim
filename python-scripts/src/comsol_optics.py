"""Optical eigen-search helpers and gap detection shared by both cell types.

The functions marked MOVED come unchanged from scripts/run_opt_comsol.py, which
re-imports them under the same names (hole behaviour identical). The functions
after them are new and used by the rib cell only.
"""
from __future__ import annotations

import numpy as np

from bandgap import Gap

C0 = 299_792_458.0
OPT_F_TARGET = C0 / 1550e-9      # 193.4 THz (same value as run_opt_comsol.py)


def _study_dataset(model, study):
    matches = [d for d in model.datasets() if d.startswith(study + "//")]
    return matches[0] if matches else None


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


# ---- rib cell (multimode): arbitrary k grid + guided-mode interval gap ------

def opt_sweep_at(model, g, study, k_frac, n_bands, min_THz=10.0,
                 f_target=OPT_F_TARGET):
    """_opt_sweep on a caller-supplied k grid. Returns (k_frac, freqs [Hz]).

    k_frac = k_z*a/pi. Same per-k search as _opt_sweep: shift on the light
    line (_guided_shift_hz), neigs = max(n_bands, 6), and frequencies below
    min_THz (scattering-boundary pseudo-modes) dropped.
    """
    k_frac = np.asarray(k_frac, dtype=float).ravel()
    freqs = np.full((k_frac.size, n_bands), np.nan)
    neigs_req = max(n_bands, 6)
    for i, s in enumerate(k_frac):
        kz = np.pi * s / g.a
        model.parameter("kF", f"{kz}[1/m]")
        _set_eig_search(model, study, _guided_shift_hz(kz, f_target=f_target),
                        neigs=neigs_req)
        model.solve(study)
        dset = _study_dataset(model, study)
        ev = (model.evaluate("freq", dataset=dset) if dset
              else model.evaluate("freq"))
        fs = np.sort(np.real(np.asarray(ev)).ravel())
        fs = fs[fs > min_THz * 1e12]
        n_av = min(n_bands, len(fs))
        freqs[i, :n_av] = fs[:n_av]
    return k_frac, freqs


def rib_k_grid(a, n_k, f_target=OPT_F_TARGET, rel_tol=0.30, safety=0.9):
    """k grid (k_frac = k_z*a/pi) for the rib optical sweep, ending at X.

    A guided mode of frequency f needs k_frac > 2*a*f/c (the light line), so
    modes at or above f_floor = f_target*(1 - rel_tol) cannot exist below
    k_start = 2*a*f_floor/c. The grid starts there (times `safety`), capped at
    0.5, the start of the hole cell's zone-edge sweep.
    """
    f_floor = f_target * (1.0 - rel_tol)
    k_start = safety * 2.0 * a * f_floor / C0
    k_start = min(max(k_start, 1e-3), 0.5)
    return np.linspace(k_start, 1.0, int(n_k))


def resolvable_step_hz(k_frac, a, vg_max_over_c=1.0):
    """Largest frequency change of one band between adjacent k samples.

    |df/dk_z| = v_g/(2*pi) <= vg_max*c/(2*pi) and dk_z = pi*dk_frac/a.
    """
    k_frac = np.asarray(k_frac, dtype=float).ravel()
    if k_frac.size < 2:
        return float("inf")
    return vg_max_over_c * C0 * float(np.max(np.abs(np.diff(k_frac)))) / (2.0 * a)


def guided_interval_gap(k_frac, freqs, a, f_target=OPT_F_TARGET, rel_tol=0.30,
                        min_THz=10.0, vg_max_over_c=1.0) -> Gap:
    """Guided-mode optical gap of one symmetry sector, by frequency interval.

    freqs: [n_k, n_bands] Hz from one optical study (NaN = not returned),
    rows matching k_frac. A mode is GUIDED at k if f < c*k_z/(2*pi). The gap
    is the widest interval that
      * is bounded below and above by guided modes,
      * contains no guided mode at any sampled k,
      * lies in the trusted range: above the highest "lowest returned mode"
        and the light line at the first k (smaller k is not sampled); below
        the light line at the last k and below the highest returned mode of
        any k where the solver stopped under the light line,
      * has its centre within +/- rel_tol of f_target,
      * is wider than resolvable_step_hz: a narrower interval could be a
        band crossing between two k samples. A returned gap can still be
        narrower than reported by up to that step.
    Band indices are not used (lower_band = -1).
    """
    k_frac = np.asarray(k_frac, dtype=float).ravel()
    freqs = np.asarray(freqs, dtype=float)
    if freqs.ndim != 2 or freqs.shape[0] != k_frac.size:
        raise ValueError("freqs must be [n_k, n_bands], one row per k_frac")
    if k_frac.size < 2:
        return Gap.empty()
    f_ll = C0 * k_frac / (2.0 * a)            # light line at each k
    f_lo, f_hi = float(np.min(f_ll)), float(np.max(f_ll))
    guided = []
    for row, ll in zip(freqs, f_ll):
        row = row[np.isfinite(row) & (row > min_THz * 1e12)]
        if row.size == 0:
            return Gap.empty()
        f_lo = max(f_lo, float(row.min()))
        if row.max() < ll:                    # solver stopped below the light line
            f_hi = min(f_hi, float(row.max()))
        guided.append(row[row < ll])
    if not f_hi > f_lo:
        return Gap.empty()
    f = np.concatenate(guided)
    f = np.unique(f[(f >= f_lo) & (f <= f_hi)])
    if f.size < 2:
        return Gap.empty()
    lo, hi = f[:-1], f[1:]
    width, center = hi - lo, 0.5 * (lo + hi)
    step = resolvable_step_hz(k_frac, a, vg_max_over_c)
    ok = (np.abs(center - f_target) < rel_tol * f_target) & (width > step)
    if not np.any(ok):
        return Gap.empty()
    ng = width / center
    i = int(np.argmax(np.where(ok, ng, -np.inf)))
    return Gap(float(lo[i]), float(hi[i]), float(center[i]),
               float(width[i]), float(ng[i]), -1, True)
