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
     comsol/omc2d_boomerang.mph. See comsol/README_template_2d.md for the
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

_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "omc2d_boomerang.mph")

# Study names expected in the .mph template -- one Eigenfrequency study per
# z-parity family. See comsol/README_template_2d.md section 4.
STUDY_EVENZ = "mech evenz"
STUDY_ODDZ = "mech oddz"


def _study_dataset(model, study):
    """Return the MPh dataset name for a study (e.g. 'mech evenz//Solution 1').

    MPh 1.3+ uses display names, not COMSOL node tags, for dataset lookup.
    Returns None if no matching dataset exists (model not yet solved).
    """
    matches = [d for d in model.datasets() if d.startswith(study + '//')]
    return matches[0] if matches else None


def _eval_freq(model, study):
    dset = _study_dataset(model, study)
    if dset is not None:
        return model.evaluate("freq", dataset=dset)
    return model.evaluate("freq")  # fallback: use default (last solved)


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
    """
    if not _HAVE_MPH:
        raise RuntimeError("MPh not importable. On the Mac: pip install MPh, "
                           "and ensure a COMSOL server is reachable.")
    if not os.path.exists(template):
        raise FileNotFoundError(
            f"COMSOL template not found: {template}\n"
            "Build it once per comsol/README_template_2d.md.")

    model = _get_model(template)

    # set geometry parameters (model parameters must be named a,w,r,r1,r2,th in m)
    for name, val in g.as_dict().items():
        model.parameter(name, f"{val}[m]")

    k_norm, kx_arr, ky_arr = bz_path(g.a, n_per_segment, unitcell)
    n_k = len(k_norm)

    study_by_parity = {"evenz": STUDY_EVENZ, "oddz": STUDY_ODDZ}
    out = dict(k_norm=k_norm, kx=kx_arr, ky=ky_arr, a=g.a)

    for parity in parities:
        if parity not in study_by_parity:
            raise ValueError(f"unknown parity {parity!r}, expected 'evenz'/'oddz'")
        study = study_by_parity[parity]
        freqs = np.full((n_k, n_bands), np.nan)
        for i in range(n_k):
            model.parameter("kx", f"{kx_arr[i]}[1/m]")
            model.parameter("ky", f"{ky_arr[i]}[1/m]")
            model.solve(study)
            ev = np.sort(np.real(np.asarray(_eval_freq(model, study))))
            n_av = min(n_bands, len(ev))
            freqs[i, :n_av] = ev[:n_av]
        out[f"freqs_{parity}"] = freqs

    return out


def save_bands(data, path):
    """Save a band-structure dict to .npz for later plotting."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **{k: v for k, v in data.items() if isinstance(v, np.ndarray)})


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
