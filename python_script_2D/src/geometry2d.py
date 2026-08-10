"""Geometry mapping and feasibility checks for the 2D diamond "boomerang"
phononic-crystal unit cell (hexagonal lattice).

Adapted from python-scripts/src/geometry.py (1D nanobeam project) -- same
normalized-u <-> physical-meters pattern, same dependency-light design (stdlib
+ PyYAML only) so every backend can import it. The geometry itself mirrors
omc-comsol-chang/buildBoomerangUnitCell.m (MATLAB, celltype='boomerang',
unitcell='hexagonal'): a hexagonal-lattice unit cell of side `a` with one
boomerang-shaped air hole made of 3 rectangular "legs" (width w, length r) at
120 degrees, rounded with fillet radii r1 (tips) and r2 (junctions),
extruded through a slab of thickness `th`.

Design variables (normalized u in [0,1]) -> physical parameters (meters).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, asdict
from typing import Dict, List

import yaml

_CFG_DIR = os.path.join(os.path.dirname(__file__), "..", "configs")

# Order of the optimizer's normalized vector u.
VARS: List[str] = ["a", "w", "r", "r1", "r2", "th"]


def load_bounds(path: str | None = None) -> dict:
    path = path or os.path.join(_CFG_DIR, "bounds_2d.yaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


@dataclass
class Geometry2D:
    a: float    # hexagonal lattice constant [m]
    w: float    # boomerang leg width [m]
    r: float    # boomerang leg length [m]
    r1: float   # tip fillet radius [m]
    r2: float   # junction fillet radius [m]
    th: float   # slab thickness (out-of-plane, z) [m]

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


def phys_to_u(params: dict, bounds: dict | None = None) -> list:
    """Convert physical params dict (meters) -> normalized u in [0,1]^N.

    Uses the CURRENT bounds, so the result is always consistent with whatever
    bounds_2d.yaml says right now. Clips silently to [0,1] for out-of-range params.
    """
    bounds = bounds or load_bounds()
    v = bounds["variables"]
    u = []
    for name in VARS:
        lo, hi = v[name]["min"], v[name]["max"]
        val = float(params.get(name, lo))
        u.append(min(1.0, max(0.0, (val - lo) / (hi - lo))))
    return u


def u_to_geometry(u, bounds: dict | None = None) -> Geometry2D:
    """Map normalized u in [0,1]^6 -> Geometry2D (physical meters)."""
    bounds = bounds or load_bounds()
    v = bounds["variables"]

    u = list(u)
    if len(u) != len(VARS):
        raise ValueError(f"expected {len(VARS)} variables {VARS}, got {len(u)}")

    def lerp(name, uu):
        lo, hi = v[name]["min"], v[name]["max"]
        uu = min(1.0, max(0.0, float(uu)))
        return lo + uu * (hi - lo)

    phys = {name: lerp(name, uu) for name, uu in zip(VARS, u)}
    return Geometry2D(**phys)


def check_feasibility(g: Geometry2D, bounds: dict | None = None):
    """Return (ok: bool, reasons: list[str]). Cheap pre-simulation gate.

    CAVEAT (adaptation from the 1D project's exact sidewall/bridge formulas):
    the boomerang hole is a Boolean union of 3 rotated, filleted rectangles
    (see buildBoomerangUnitCell.m), so the exact minimum remaining solid width
    is geometry-dependent in a way that isn't a single closed-form expression
    worth re-deriving here. This function uses a CONSERVATIVE HEURISTIC:
      - leg width w and leg length r must each individually clear min_feature
        (a degenerate leg is not a fabricable feature).
      - fillet radii r1, r2 must not exceed half of whichever feature they
        round (a fillet larger than the feature it rounds is geometrically
        invalid / eats the whole feature).
      - the solid "web" between the hole and the hexagon cell boundary is
        approximated as (a*sqrt(3)/2 - r - w/2) -- the in-radius of the
        hexagonal cell minus the leg's radial extent -- and must clear
        min_feature.
    This has NOT been cross-validated against the true Boolean geometry via
    COMSOL's mphgeom at the bounds' extremes -- treat results near the edges
    of bounds_2d.yaml's ranges with caution, or tighten those bounds
    empirically once a few designs have been visually checked.
    """
    import math

    bounds = bounds or load_bounds()
    f = bounds["feasibility"]
    min_feature = f["min_feature"]
    reasons: List[str] = []

    if g.w < min_feature:
        reasons.append(f"leg width w {g.w*1e9:.0f} nm < min {min_feature*1e9:.0f} nm")
    if g.r < min_feature:
        reasons.append(f"leg length r {g.r*1e9:.0f} nm < min {min_feature*1e9:.0f} nm")
    if g.r1 > g.w / 2.0:
        reasons.append(f"tip fillet r1 {g.r1*1e9:.0f} nm exceeds w/2 {g.w/2*1e9:.0f} nm")
    if g.r2 > g.w / 2.0:
        reasons.append(f"junction fillet r2 {g.r2*1e9:.0f} nm exceeds w/2 {g.w/2*1e9:.0f} nm")

    web = (math.sqrt(3.0) / 2.0) * g.a - g.r - g.w / 2.0
    if web < min_feature:
        reasons.append(f"web (cell boundary to hole) {web*1e9:.0f} nm < min "
                       f"{min_feature*1e9:.0f} nm")

    return (len(reasons) == 0), reasons


if __name__ == "__main__":
    g = u_to_geometry([0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    ok, why = check_feasibility(g)
    print("geometry:", {k: f"{v*1e9:.1f} nm" for k, v in g.as_dict().items()})
    print("feasible:", ok, why)
