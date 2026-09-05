"""2D phononic band structure of the diamond "boomerang" unit cell via COMSOL,
driven from Python with `mph`.

Runs on YOUR Mac (COMSOL 6.2+ with a running COMSOL server), NOT in a sandbox.

This is the 2D-phononic analog of python-scripts/src/acoustic_comsol.py,
adapted from the MATLAB reference implementation in
omc-comsol-chang/buildBoomerangUnitCell.m (geometry) and
omc-comsol-chang/runBands_2D.m (physics, BCs, Brillouin-zone path). The two
things that are genuinely different from the 1D nanobeam project:

1. TRUE 2D Bloch periodicity. The unit cell is periodic in BOTH in-plane
   directions (x, y) via Floquet boundary conditions with wavevector
   (kx, ky), swept around the hexagonal Brillouin zone Gamma -> M -> K ->
   Gamma (see `_kxky_hexagonal` below, a direct port of the piecewise kx/ky
   expressions in runBands_2D.m). There is no separate in-plane symmetry
   reduction (unlike the 1D beam's y-even/y-odd split) -- the full unit cell
   is simulated in x and y.

2. Z-PARITY family split. The slab has finite thickness `th` and a z=0
   mirror plane. Mechanical modes split into Z-EVEN and Z-ODD families
   (symmetric / antisymmetric about the midplane) -- matching runBands_2D.m's
   `ZsymSel` + SymmetrySolid/Antisymmetry boundary-condition pair, driven by
   P.mbevenz there. A full/COMPLETE phononic gap must be free of bands from
   BOTH families; a family-restricted (symmetry) gap only needs one family
   clear -- see bandgap.py's module docstring and targets_2d.yaml's
   `gap_mode`.

Workflow (template-driven, same convention as the 1D project):
  1. Build ONE parametric model in the COMSOL GUI and save as
     comsol/trusty_boomerang.mph. See comsol/README_template_2d.md for the
     exact node-by-node recipe (geometry params a,w,r,r1,r2,th; Solid
     Mechanics; Floquet periodicity on BOTH x-face-pairs and y-face-pairs with
     wavevector (kx,ky); two Eigenfrequency studies "mech evenz" / "mech
     oddz" using the z=0 symmetry/antisymmetry boundary condition).
  2. This module loads that template, sets geometry + (kx,ky) at each
     Brillouin-zone point, solves each requested z-parity study, and returns
     eigenfrequencies so bandgap.py can find the gap(s).

Install once on the Mac:   pip install MPh        (talks to a COMSOL server)
Start server (separate terminal):  comsol mphserver   (or use mph.start())
"""
from __future__ import annotations

import os
import numpy as np

try:
    from comsol_client import have_mph as _have_mph_fn, get_model as _get_model
    _HAVE_MPH = _have_mph_fn()
except Exception:  # pragma: no cover
    _HAVE_MPH = False
    def _get_model(t):
        raise RuntimeError("comsol_client not available")

_TEMPLATE = os.environ.get("OMC2D_TEMPLATE") or os.path.join(
    os.path.dirname(__file__), "..", "comsol", "trusty_boomerang.mph")

# Study LABELS expected in the .mph template -- one Eigenfrequency study per
# z-parity family. MPh resolves study nodes by label, not by COMSOL tag.
#
# The current comsol/trusty_boomerang.mph does NOT use these labels: it has
# tag `study` labelled "Study_symmetric" (disables smech/asymBCs -> symBCs
# active -> z-even) and tag `std1` with no label set, so COMSOL reports the
# default "Study 1" (disables smech/symBCs -> asymBCs active -> z-odd). The
# parity mapping is right, only the names are wrong.
#
# Preferred fix is to relabel the two studies in the GUI to match the
# defaults below -- hard-coding "Study 1" is fragile because COMSOL re-derives
# default labels whenever another study is added. The env vars are an escape
# hatch for trying the existing template without a GUI edit:
#     OMC2D_STUDY_EVENZ="Study_symmetric" OMC2D_STUDY_ODDZ="Study 1" ...
STUDY_EVENZ = os.environ.get("OMC2D_STUDY_EVENZ", "mech evenz")
STUDY_ODDZ = os.environ.get("OMC2D_STUDY_ODDZ", "mech oddz")

# Geometry parameters the template must already define. Setting a parameter
# that does not exist SILENTLY CREATES it (MPh's model.parameter is a bare
# param().set with no existence check, and COMSOL's param.set creates on
# miss), so without _require_parameters below every candidate would solve the
# template's frozen geometry while appearing to succeed.
#
# Includes the derived helpers (hx, hy, selw, dsel) that the parameterized
# geometry expressions reference -- a template missing those would otherwise
# fail deep inside a COMSOL geometry rebuild instead of here with a clear
# message. Keep in sync with the model.param.set block in
# comsol/trusty_boomerang_script.m.
REQUIRED_PARAMS = ("a", "w", "r", "r1", "r2", "th",
                   "hx", "hy", "selw", "dsel", "k", "kx", "ky")


def _require_studies(model, studies, template):
    missing = set(studies) - set(model.studies())
    if missing:
        raise RuntimeError(
            f"template {template} has no study named {sorted(missing)}. "
            f"Available: {model.studies()}. Relabel the studies in the COMSOL "
            f"GUI, or set OMC2D_STUDY_EVENZ / OMC2D_STUDY_ODDZ.")


def _require_parameters(model, names, template):
    """Raise unless every name is ALREADY defined in the template.

    model.parameter(name) with no value raises if the parameter is undefined;
    that is the only cheap way to distinguish "set an existing parameter" from
    "silently create a new orphan one". See REQUIRED_PARAMS.
    """
    missing = []
    for name in names:
        try:
            model.parameter(name)
        except Exception:
            missing.append(name)
    if missing:
        raise RuntimeError(
            f"template {template} does not define parameter(s) {missing}. "
            f"Setting them would silently create unused parameters and every "
            f"candidate would solve the template's frozen geometry. The "
            f"template must be parameterized first -- see "
            f"comsol/README_template_2d.md.")


def _study_dataset(model, study):
    """Return the MPh dataset name for a study (e.g. 'mech evenz//Solution 1').

    MPh escapes '/' by doubling it, so a COMSOL label 'X//Solution 1' comes
    back from model.datasets() as 'X////Solution 1' -- hence the escaped
    prefix. Prefers the parametric dataset when a sweep is present. Raises
    rather than falling back to the default dataset: an unqualified
    model.evaluate("freq") resolves to whatever dataset is first in the model,
    which may belong to the OTHER parity -- a silent parity mix-up.
    """
    prefix = study.replace('/', '//') + '//'
    matches = [d for d in model.datasets() if d.startswith(prefix)]
    if not matches:
        raise RuntimeError(
            f"no dataset for study {study!r} (model solved?). "
            f"Available datasets: {model.datasets()}")
    parametric = [d for d in matches if 'Parametric' in d]
    return parametric[0] if parametric else matches[0]


def _eval_freq(model, study):
    return model.evaluate("freq", dataset=_study_dataset(model, study))


# ── Brillouin-zone path (ported from runBands_2D.m) ──────────────────────────
#
# k is a continuous parameter in [0, 3): segment [0,1) is Gamma->M (hexagonal)
# or Gamma->X (square), [1,2) is M->K / X->M, [2,3) is K->Gamma / M->Gamma.
# These are EXACT ports of the `model.param.set('kx', 'if(k<1, ...')`
# expressions in runBands_2D.m -- kept as plain Python instead of a COMSOL
# expression, since the Python driver solves one (kx,ky) point at a time and
# has no need for COMSOL's own parametric-sweep node (which the MATLAB
# pipeline uses because it drives the *native* GUI parametric batch solver).

def _kxky_hexagonal(k, a):
    if k < 1:
        kx = (-1.0 / np.sqrt(3.0)) * k * (np.pi / a)
        ky = (np.pi / a) * k
    elif k < 2:
        kx = (1.0 / np.sqrt(3.0)) * (np.pi / a) * (k - 2)
        ky = (k + 2) * np.pi / (3 * a)
    else:
        kx = 0.0
        ky = (3 - k) * 4 * np.pi / (3 * a)
    return kx, ky


def _kxky_square(k, a):
    if k < 1:
        kx = np.pi / a * k
        ky = 0.0
    elif k < 2:
        kx = np.pi / a
        ky = (k - 1) * np.pi / a
    else:
        kx = (3 - k) * np.pi / a
        ky = (3 - k) * np.pi / a
    return kx, ky


def bz_path(a, n_per_segment=15, unitcell="hexagonal"):
    """Return (k_norm, kx, ky) arrays sampling Gamma->M->K->Gamma (hexagonal)
    or Gamma->X->M->Gamma (square), matching runBands_2D.m's convention:
    3*n_per_segment points at k = i/n_per_segment for i=0..3*n_per_segment-1,
    PLUS one final point closing back at Gamma (k_norm=3) for plotting.
    """
    fn = _kxky_hexagonal if unitcell == "hexagonal" else _kxky_square
    n_tot = 3 * n_per_segment
    k_norm = np.array([i / n_per_segment for i in range(n_tot)] + [3.0])
    kx = np.zeros_like(k_norm)
    ky = np.zeros_like(k_norm)
    for i, k in enumerate(k_norm[:-1]):
        kx[i], ky[i] = fn(k, a)
    kx[-1], ky[-1] = fn(0.0, a)   # close the loop back at Gamma
    return k_norm, kx, ky


# ── Main driver ───────────────────────────────────────────────────────────────

def run_mechanical_comsol_2d(g, *, n_per_segment=15, n_bands=10,
                             parities=("evenz", "oddz"), unitcell="hexagonal",
                             template=_TEMPLATE):
    """Solve the 2D phononic band structure for one or both z-parities.

    Returns dict:
      k_norm, kx, ky : [n_k] Brillouin-zone path arrays (n_k = 3*n_per_segment+1)
      freqs_evenz, freqs_oddz : [n_k, n_bands] eigenfrequencies (Hz), or absent
        if that parity wasn't requested
      a : lattice constant [m] (for convenience downstream)

    KNOWN ISSUE (blocked on the template, see comsol/README_template_2d.md):
    this drives k-space by writing kx/ky as literals per point. In a template
    where kx/ky are COMSOL EXPRESSIONS of a swept scalar k (as
    comsol/trusty_boomerang.mph has them), that overwrite destroys the
    expressions and every point of the study's own k sweep solves the same
    wavevector. The right shape is to set only `k`, let COMSOL's Parametric
    step sweep it once per parity, and read back [n_k, n_bands] from the
    parametric dataset -- 54 eigensolves per candidate at kpts=9 instead of
    92. Do that as part of parameterizing the template; until then
    _require_parameters below raises before any of this can run.
    """
    if not _HAVE_MPH:
        raise RuntimeError("MPh not importable. On the Mac: pip install MPh, "
                           "and ensure a COMSOL server is reachable.")
    if not os.path.exists(template):
        raise FileNotFoundError(
            f"COMSOL template not found: {template}\n"
            "Build it once per comsol/README_template_2d.md.")

    study_by_parity = {"evenz": STUDY_EVENZ, "oddz": STUDY_ODDZ}
    unknown = set(parities) - set(study_by_parity)
    if unknown:
        raise ValueError(f"unknown parities {sorted(unknown)}, expected "
                         f"'evenz'/'oddz'")

    model = _get_model(template)
    _require_studies(model, [study_by_parity[p] for p in parities], template)
    _require_parameters(model, REQUIRED_PARAMS, template)

    for name, val in g.as_dict().items():
        model.parameter(name, f"{val}[m]")

    k_norm, kx_arr, ky_arr = bz_path(g.a, n_per_segment, unitcell)
    n_k = len(k_norm)

    out = dict(k_norm=k_norm, kx=kx_arr, ky=ky_arr, a=g.a)

    for parity in parities:
        study = study_by_parity[parity]
        freqs = np.full((n_k, n_bands), np.nan)
        for i in range(n_k):
            model.parameter("kx", f"{kx_arr[i]}[1/m]")
            model.parameter("ky", f"{ky_arr[i]}[1/m]")
            model.solve(study)
            ev = np.sort(np.real(np.atleast_1d(np.asarray(_eval_freq(model, study)))))
            if ev.size < n_bands:
                # Do NOT pad with NaN and carry on. gap_near_frequency runs
                # largest_gap with nan_safe=False, so an all-NaN column makes
                # np.max return NaN, the isfinite guard skips that band pair
                # AND every pair above it, and the candidate scores a clean
                # G_m = 0.0 with status="success" -- a good design silently
                # discarded. Fail loudly instead; objective2d records it as
                # mech_failed.
                raise RuntimeError(
                    f"{study}: eigensolver returned {ev.size} modes at k-point "
                    f"{i}/{n_k - 1} (k={k_norm[i]:.3f}), need {n_bands}. "
                    f"Raise neigs in the study, or lower n_bands.")
            freqs[i, :] = ev[:n_bands]
        out[f"freqs_{parity}"] = freqs

    return out


def save_bands(data, path):
    """Save a band-structure dict to .npz for later plotting. Returns `path`.

    Scalars are stored as 0-d arrays. The previous version filtered on
    `isinstance(v, np.ndarray)`, which silently dropped `a` -- the lattice
    constant is a plain Python float in the dict this module returns, and it is
    exactly what scripts/plot_bands_2d.py needs for the BZ annotation. A
    dropped key is invisible until something downstream KeyErrors on it, so
    coerce rather than filter.

    Anything that is neither array nor scalar is still skipped: np.savez would
    pickle it, and np.load would then need allow_pickle=True to read the file
    back, which turns a data file into a code-execution surface.
    """
    out = {}
    for key, val in data.items():
        if isinstance(val, np.ndarray):
            out[key] = val
        elif isinstance(val, (bool, int, float, np.bool_, np.integer, np.floating)):
            out[key] = np.asarray(val)
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)   # abspath so a bare filename works too
    np.savez(path, **out)
    return path


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry2d import u_to_geometry
    from bandgap import gap_near_frequency
    g = u_to_geometry([0.5, 0.5, 0.5, 0.3, 0.3, 0.5])
    data = run_mechanical_comsol_2d(g)
    gap_evenz = gap_near_frequency(data["freqs_evenz"], 8.0e9)
    gap_oddz = gap_near_frequency(data["freqs_oddz"], 8.0e9)
    print("evenz gap:", gap_evenz)
    print("oddz gap: ", gap_oddz)
