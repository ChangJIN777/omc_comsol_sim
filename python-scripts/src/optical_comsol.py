"""Optical band structure via COMSOL Wave Optics, driven from Python with `mph`.

Use this (instead of MPB) when you need FABRICATION-AWARE optical bands on the
SAME geometry/mesh as the mechanical model: fillets, sidewall angle, etched
surface layers, etc. Runs on the Mac with COMSOL 6.2 + MPh.

Template recipe: comsol/README_template.md (Electromagnetic Waves, Frequency
Domain or Mode Analysis with Floquet periodicity along z). We search for guided
TE-like modes near the target frequency and sweep kz along Gamma->X.
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

C0 = 299_792_458.0
_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "omc_unitcell_iso.mph")


def _study_dataset(model, study):
    """Return the MPh dataset name for a study (e.g. 'Study 2//Solution 2').

    MPh 1.3+ uses display names, not COMSOL node tags, for dataset lookup.
    """
    matches = [d for d in model.datasets() if d.startswith(study + '//')]
    return matches[0] if matches else None


def run_optical_comsol(g, *, n_k=15, n_bands=6, template=_TEMPLATE,
                       study="opt TE"):
    """Return dict with optical band frequencies [n_k, n_bands] in Hz.

    Assumes the same template carries an EM physics + eigenfrequency/mode study
    named `study`. Set the diamond refractive index (2.40) in the template.
    """
    if not _HAVE_MPH:
        raise RuntimeError("MPh not importable. On the Mac: pip install MPh.")
    if not os.path.exists(template):
        raise FileNotFoundError(f"COMSOL template not found: {template}")

    model = _get_model(template)
    for name, val in g.as_dict().items():
        model.parameter(name, f"{val}[m]")

    # Sweep Gamma->X. k_frac = k_z * a / pi in [0,1]; X (zone edge) at 1.
    k_frac = np.linspace(1e-3, 1.0, n_k)
    freqs = np.full((n_k, n_bands), np.nan)
    for i, s in enumerate(k_frac):
        model.parameter("kF", f"{np.pi*s/g.a}[1/m]")
        model.solve(study)
        dset = _study_dataset(model, study)
        ev = (model.evaluate("freq", dataset=dset) if dset
              else model.evaluate("freq"))
        freqs[i, :len(ev)] = np.sort(np.real(np.asarray(ev)))[:n_bands]
    return dict(freqs_hz=freqs, k_frac=k_frac, a=g.a)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry import u_to_geometry
    from bandgap import gap_near_frequency
    g = u_to_geometry([0.5, 0.6, 0.5, 0.6])
    data = run_optical_comsol(g)
    print(gap_near_frequency(data["freqs_hz"], C0 / 1550e-9))
