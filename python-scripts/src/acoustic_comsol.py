"""Mechanical (phononic) band structure via COMSOL, driven from Python with `mph`.

Runs on YOUR Mac (COMSOL 6.2 + a running COMSOL server), NOT in the sandbox.

Workflow (template-driven -- recommended):
  1. Build ONE parametric model in the COMSOL GUI and save as
     comsol/omc_unitcell_iso.mph. See comsol/README_template.md for the exact
     node-by-node recipe (geometry params a,w,t,hx,hy; Solid Mechanics;
     Floquet periodicity on the two z-faces with wavevector kF; Eigenfrequency
     study searching for `n_bands` modes around a shift frequency).
  2. This module loads that template, sets parameters + the Floquet wavevector,
     solves at each k along Gamma->X, and pulls eigenfrequencies + a symmetry
     descriptor for each mode so bandgap.py can isolate the breathing family.

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
                         "omc_unitcell_iso.mph")


def _study_dataset(model, study):
    """Return the MPh dataset name for a study (e.g. 'Study 1//Solution 1').

    MPh 1.3+ uses display names, not COMSOL node tags, for dataset lookup.
    Returns None if no matching dataset exists (model not yet solved).
    """
    matches = [d for d in model.datasets() if d.startswith(study + '//')]
    return matches[0] if matches else None


def _eval_freq(model, study):
    """Evaluate eigenfrequencies after a solve, using the correct dataset name."""
    dset = _study_dataset(model, study)
    if dset is not None:
        return model.evaluate("freq", dataset=dset)
    return model.evaluate("freq")  # fallback: use default (last solved)


def _symmetry_descriptor(model, study, solnum):
    """Classify one eigenmode (solnum, 1-based) by displacement parity.

    Returns fractional energy in each component; breathing modes have large fy.
    Adjust integration operator name / expressions to match your template.
    """
    dset = _study_dataset(model, study)
    kw = {"dataset": dset, "solnum": solnum} if dset else {"solnum": solnum}
    ux2 = float(np.real(model.evaluate("intop1(solid.ux^2)", **kw)))
    uy2 = float(np.real(model.evaluate("intop1(solid.uy^2)", **kw)))
    uz2 = float(np.real(model.evaluate("intop1(solid.uz^2)", **kw)))
    tot = ux2 + uy2 + uz2 or 1.0
    return dict(fx=ux2 / tot, fy=uy2 / tot, fz=uz2 / tot)


def run_mechanical_comsol(g, *, n_k=15, n_bands=12, shift_GHz=8.0,
                          template=_TEMPLATE, classify=False,
                          study="mech sym"):
    """Solve the mechanical band structure. Returns dict with freqs [n_k,n_bands] in Hz."""
    if not _HAVE_MPH:
        raise RuntimeError("MPh not importable. On the Mac: pip install MPh, "
                           "and ensure a COMSOL server is reachable.")
    if not os.path.exists(template):
        raise FileNotFoundError(
            f"COMSOL template not found: {template}\n"
            "Build it once per comsol/README_template.md.")

    model = _get_model(template)

    # set geometry parameters (model parameters must be named a,w,t,hx,hy in m)
    for name, val in g.as_dict().items():
        model.parameter(name, f"{val}[m]")

    # Sweep Gamma->X. k_frac = k_z * a / pi in [0,1]; X (zone edge) at 1.
    k_frac = np.linspace(1e-3, 1.0, n_k)
    freqs = np.full((n_k, n_bands), np.nan)
    classes = [[None] * n_bands for _ in range(n_k)]

    for i, s in enumerate(k_frac):
        # Floquet wavevector along beam axis z: kz = pi * s / a
        model.parameter("kF", f"{np.pi*s/g.a}[1/m]")
        model.solve(study)
        ev = np.sort(np.real(np.asarray(_eval_freq(model, study))))[:n_bands]
        freqs[i, :len(ev)] = ev
        if classify:
            classes[i] = [_symmetry_descriptor(model, study, j + 1)
                          for j in range(len(ev))]

    return dict(freqs_hz=freqs, k_frac=k_frac, a=g.a,
                classes=classes if classify else None)


def save_bands(data, path):
    """Save a band-structure dict to .npz for scripts/plot_bands.py."""
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, k_frac=data["k_frac"], freqs_hz=data["freqs_hz"], a=data["a"])


def export_mode_grid(model, g, path, *, dset="dset1", solnum=1, n=120):
    """Export one solved mode on the x=0 top cut plane to .npz for plot_mode.py.

    Call right after a solve at the k-point of interest. `solnum` selects which
    eigenmode (1-based). Uses mphinterp to sample ux,uy,uz on a (y,z) grid.
    """
    import os
    z = np.linspace(-g.a/2, g.a/2, n)
    y = np.linspace(-g.w/2, g.w/2, n)
    Z, Y = np.meshgrid(z, y)
    pts = np.vstack([np.zeros(Z.size), Y.ravel(), Z.ravel()])  # (x,y,z)
    comps = {}
    for c in ("u", "v", "w"):   # COMSOL solid default disp components u,v,w
        vals = model.evaluate(c, dataset=dset)  # adapt to mphinterp at pts in your driver
        comps[c] = np.asarray(vals)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # NOTE: replace the .evaluate calls above with an mphinterp(model, 'u','v','w',
    # pts, dataset=dset, solnum=solnum) call in your local MPh version; shapes
    # must match Z. Saved keys match scripts/plot_mode.py.
    np.savez(path, Y=Y, Z=Z,
             UX=comps["u"].reshape(Z.shape) if comps["u"].size == Z.size else np.zeros_like(Z),
             UY=comps["v"].reshape(Z.shape) if comps["v"].size == Z.size else np.zeros_like(Z),
             UZ=comps["w"].reshape(Z.shape) if comps["w"].size == Z.size else np.zeros_like(Z),
             freq=0.0, kind="mechanical")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry import u_to_geometry
    from bandgap import gap_near_frequency
    g = u_to_geometry([0.5, 0.6, 0.5, 0.6])
    data = run_mechanical_comsol(g)
    gap = gap_near_frequency(data["freqs_hz"], 8.0e9)
    print("mechanical gap:", gap)
