"""PLACEHOLDER -- 2D slab photonic band structure via COMSOL. NOT YET IMPLEMENTED.

See docs/optical_2d_plan.md for the full design writeup. This module exists
so src/objective2d.py already has a wired-up `optical_backend="comsol"`
branch to flip on once the physics/template below are built -- the plan is
kept here (next to the code it will become) rather than only in docs, so
whoever implements it has the contract already spelled out.

Why this is harder than the 1D project's optical backend
(python-scripts/src/optical_comsol.py):
  - The 1D nanobeam has a single guided-mode family swept along ONE Bloch
    direction (kz); "TE-like" there just means Ey-dominant, y-even.
  - Here the photonic crystal is 2D-periodic in-plane (kx,ky over the same
    hexagonal Brillouin zone as the mechanical problem -- see
    acoustic_comsol_2d.bz_path), and out-of-plane confinement in the finite
    slab means every band must ALSO be checked against the light line at
    each (kx,ky), not just filtered by a single symmetry family. Guided modes
    only exist below the light line; above it they're leaky/radiative
    (compare src/acoustic_comsol_2d.py's z-parity split, which has no
    analogous "light line" complication because mechanical waves in a free
    slab don't leak into vacuum the way photons do).

Planned interface (mirrors acoustic_comsol_2d.run_mechanical_comsol_2d):

    def run_optical_comsol_2d(g, *, n_per_segment=15, n_bands=8,
                              polarization="evenz", unitcell="hexagonal",
                              template=_TEMPLATE):
        '''Returns dict with k_norm, kx, ky, freqs_hz [n_k,n_bands], and
        (once implemented) a per-(k,band) boolean `below_light_line` mask so
        objective2d.py can discard leaky bands before calling bandgap.py.'''

Planned template additions (comsol/trusty_boomerang.mph), see
comsol/README_template_2d.md section 5 for the up-to-date recipe once
written:
  - Physics: Electromagnetic Waves, Frequency Domain (or Mode Analysis).
  - Floquet periodicity on BOTH x- and y-face-pairs, k=(kx,ky,0) -- same BZ
    path as the mechanical studies, so both physics can be swept together.
  - z-mirror symmetry BC (PEC/PMC on the z=0 plane) analogous to the
    mechanical evenz/oddz split, to classify "TE-like" (z-even, by the
    planned convention in targets_2d.yaml) vs "TM-like" modes.
  - Air cladding + scattering boundary condition or PML above/below the slab
    so leaky (above-light-line) modes are handled rather than spuriously
    confined by an artificial boundary.
  - Study named e.g. "opt evenz" / "opt oddz", searching `n_bands` modes near
    the target frequency (targets_2d.yaml: optical.target_wavelength_nm).

Until this is implemented, keep `require_opt: false` in run configs and
`optical_backend="none"` in objective2d.evaluate_candidate -- the mechanical
pipeline (src/acoustic_comsol_2d.py) does not depend on this module.
"""
from __future__ import annotations

import os

_TEMPLATE = os.environ.get("OMC2D_TEMPLATE") or os.path.join(
    os.path.dirname(__file__), "..", "comsol", "trusty_boomerang.mph")


def run_optical_comsol_2d(*args, **kwargs):
    raise NotImplementedError(
        "2D optical backend not yet implemented -- see the module docstring "
        "in src/optical_comsol_2d.py and docs/optical_2d_plan.md for the "
        "planned approach. Use optical_backend='none' (require_opt: false) "
        "in the meantime; the mechanical pipeline does not need this."
    )
