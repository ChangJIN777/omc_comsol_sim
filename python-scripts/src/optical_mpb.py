"""Optical band structure of the 1D diamond nanobeam via MPB (MIT Photonic Bands).

This is the recommended FAST, accurate optical backend for *ideal* geometry
(sharp corners, vertical sidewalls). Run on a machine with pymeep installed:

    conda create -n omc -c conda-forge pymeep pymeep-extras
    conda activate omc
    python src/optical_mpb.py

For fabrication-aware optical bands (fillets, sidewall angle, surface oxide),
use the COMSOL backend (src/optical_comsol.py) instead -- MPB is the cheap,
high-throughput pre-screen; COMSOL is the fabrication-faithful check.

Model: a 3D supercell with the nanobeam along the z (periodic) axis, finite
extent in x (thickness) and y (width), surrounded by air + enough padding so the
guided modes decay. Diamond eps = n^2 ~ 5.76. We compute the lowest bands along
the irreducible BZ edge (Gamma -> X = pi/a) and report the TE-like (y-even)
photonic gap. Light-line filtering removes radiative (leaky) bands.
"""
from __future__ import annotations

import numpy as np

try:
    import meep as mp
    from meep import mpb
    _HAVE_MPB = True
except Exception:  # pragma: no cover - only importable where MPB is installed
    _HAVE_MPB = False


N_DIAMOND = 2.40
EPS_DIAMOND = N_DIAMOND ** 2


def run_optical_mpb(g, *, resolution_per_a=24, n_bands=8, n_k=21,
                    pad_y=1.5e-6, pad_x=1.2e-6, parity="EVEN_Y"):
    """Return dict with band data and the TE-like gap near the target.

    Parameters use SI meters from geometry.Geometry; MPB works in units of the
    lattice constant a, so everything is normalized by a internally.
    Frequencies are returned in Hz (f = c * f_mpb / a).
    """
    if not _HAVE_MPB:
        raise RuntimeError(
            "pymeep/MPB not importable. Install with:\n"
            "  conda create -n omc -c conda-forge pymeep pymeep-extras")

    a = g.a
    # normalized dimensions (units of a)
    w = g.w / a
    t = g.t / a
    rx = g.hx / a      # ellipse semi-axis along z (beam axis)
    ry = g.hy / a      # ellipse semi-axis along y (width)
    Ly = w + 2 * pad_y / a
    Lx = t + 2 * pad_x / a

    geometry_lattice = mp.Lattice(size=mp.Vector3(0, Ly, Lx))  # period along x-> use z? see note
    # NOTE: MPB convention here -- we make the *first* lattice vector the
    # periodic (beam) axis. We orient the beam along the x-axis of MPB and keep
    # y,z as the transverse supercell to match mpb's k along the first axis.
    geometry_lattice = mp.Lattice(size=mp.Vector3(1, Ly, Lx))

    beam = mp.Block(size=mp.Vector3(mp.inf, w, t),
                    material=mp.Medium(epsilon=EPS_DIAMOND))
    hole = mp.Ellipsoid(size=mp.Vector3(2 * rx, 2 * ry, mp.inf),
                        material=mp.air)
    geometry = [beam, hole]

    k_points = [mp.Vector3(kx, 0, 0) for kx in np.linspace(0, 0.5, n_k)]

    ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
                        geometry=geometry,
                        k_points=k_points,
                        resolution=resolution_per_a,
                        num_bands=n_bands)

    # TE-like (Ey-dominant) -> even about the y=0 mirror plane.
    if parity == "EVEN_Y":
        ms.run_yeven()
    elif parity == "ODD_Y":
        ms.run_yodd()
    else:
        ms.run()

    freqs_norm = np.array(ms.all_freqs)        # [n_k, n_bands], units c/a
    freqs_hz = freqs_norm * (299_792_458.0 / a)

    # crude light-line filter: guided modes have f < c*|k|/(2*pi*a) ... in
    # normalized units omega < k. Keep modes below the air light line.
    kx = np.linspace(0, 0.5, n_k)            # MPB units: fraction of 2*pi/a
    k_frac = 2 * kx                          # normalize so X (zone edge) = 1
    return dict(freqs_norm=freqs_norm, freqs_hz=freqs_hz,
                k_norm=kx, k_frac=k_frac, a=a, parity=parity)


def save_bands(data, path):
    """Save optical bands to .npz for scripts/plot_bands.py (--kind optical)."""
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, k_frac=data["k_frac"], freqs_hz=data["freqs_hz"], a=data["a"])


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry import u_to_geometry
    from bandgap import gap_near_frequency
    g = u_to_geometry([0.5, 0.6, 0.5, 0.6])
    data = run_optical_mpb(g)
    gap = gap_near_frequency(data["freqs_hz"], 299_792_458.0 / 1550e-9)
    print("TE-like gap:", gap)
