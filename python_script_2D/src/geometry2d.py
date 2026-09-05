"""Geometry mapping and feasibility checks for the 2D diamond "boomerang"
phononic-crystal unit cell (hexagonal lattice).

Adapted from python-scripts/src/geometry.py (1D nanobeam project) -- same
normalized-u <-> physical-meters pattern, same dependency-light design (stdlib
+ PyYAML only) so every backend can import it. The geometry itself mirrors
omc-comsol-chang/buildBoomerangUnitCell.m (MATLAB, celltype='boomerang',
unitcell='hexagonal'): a RHOMBIC primitive cell of the triangular lattice,
vertices (0,0), (a/2, a*sqrt(3)/2), (3a/2, a*sqrt(3)/2), (a,0) -- side `a`,
so `a` is the lattice constant -- with one boomerang-shaped air hole at the
cell centroid (3a/4, a*sqrt(3)/4), made of 3 rectangular "legs" (width w,
radial length r) at 120 degrees, extruded through a slab of thickness `th`.

FILLET CONVENTION (this was documented backwards before -- cross-check
against addFillet() at buildBoomerangUnitCell.m:119-149 before changing):
  r1 -> h_fil1, selected by h_disksel1, an annulus [w/(2*sqrt(2)), w] about
        the hole center. That catches the vertices at distance w/2 = the
        INNER (junction) corners where the three legs meet.
  r2 -> h_fil2, selected by h_disksel2, an annulus [r - sw, r + sw] about the
        hole center (sw = selection_width/2 = 25 nm in MATLAB). That catches
        the vertices at hypot(r, w/2) = the OUTER (tip) corners.
So r1 is the JUNCTION fillet and r2 is the TIP fillet.

Design variables (normalized u in [0,1]) -> physical parameters (meters).
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, asdict
from typing import Dict, List

import yaml

_CFG_DIR = os.path.join(os.path.dirname(__file__), "..", "configs")

# Order of the optimizer's normalized vector u. FREE variables only -- r1, r2
# and th are held fixed (bounds_2d.yaml:fixed) and are NOT part of u.
# Keep in sync with bounds_2d.yaml:variables and optimizer.N_DIM.
VARS: List[str] = ["a", "w", "r"]

# Geometry parameters that are held constant. Values live in
# bounds_2d.yaml:fixed so there is one source of truth.
FIXED_VARS: List[str] = ["r1", "r2", "th"]


def load_bounds(path: str | None = None) -> dict:
    path = path or os.path.join(_CFG_DIR, "bounds_2d.yaml")
    with open(path) as fh:
        b = yaml.safe_load(fh)
    missing_free = [n for n in VARS if n not in b.get("variables", {})]
    missing_fixed = [n for n in FIXED_VARS if n not in b.get("fixed", {})]
    if missing_free or missing_fixed:
        raise ValueError(
            f"{path}: `variables` is missing {missing_free} and/or `fixed` is "
            f"missing {missing_fixed}. VARS={VARS}, FIXED_VARS={FIXED_VARS}.")
    stray = set(b["variables"]) & set(b["fixed"])
    if stray:
        raise ValueError(f"{path}: {sorted(stray)} appear in BOTH `variables` "
                         f"and `fixed`; a parameter is one or the other.")
    return b


@dataclass
class Geometry2D:
    a: float    # lattice constant (rhombic primitive cell side) [m]
    w: float    # boomerang leg width (transverse) [m]
    r: float    # boomerang leg radial length from hole center [m]
    r1: float   # JUNCTION fillet radius (inner corners, at w/2) [m]
    r2: float   # TIP fillet radius (outer corners, at hypot(r, w/2)) [m]
    th: float   # slab thickness (out-of-plane, z; FULL thickness) [m]

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


def phys_to_u(params: dict, bounds: dict | None = None) -> list:
    """Convert physical params dict (meters) -> normalized u over the FREE vars.

    Only VARS are encoded; r1/r2/th are fixed and carry no u component. Uses
    the CURRENT bounds, so the result is always consistent with whatever
    bounds_2d.yaml says right now. Clips silently to [0,1] for out-of-range
    params. Raises if a FIXED param is supplied with a different value than the
    config -- silently dropping it would make the round-trip lossy.
    """
    bounds = bounds or load_bounds()
    v, f = bounds["variables"], bounds["fixed"]

    for name in FIXED_VARS:
        if name in params and not _close(float(params[name]), float(f[name])):
            raise ValueError(
                f"{name} is a FIXED parameter (bounds_2d.yaml:fixed = "
                f"{f[name]*1e9:.1f} nm) but was given {float(params[name])*1e9:.1f} "
                f"nm. It has no u component, so this value would be lost. Change "
                f"bounds_2d.yaml, or promote {name} to a free variable.")

    u = []
    for name in VARS:
        lo, hi = v[name]["min"], v[name]["max"]
        val = float(params.get(name, lo))
        u.append(min(1.0, max(0.0, (val - lo) / (hi - lo))))
    return u


def _close(x, y, rtol=1e-9):
    return abs(x - y) <= rtol * max(abs(x), abs(y), 1e-30)


def u_to_geometry(u, bounds: dict | None = None) -> Geometry2D:
    """Map normalized u over the FREE vars -> Geometry2D (physical meters).

    r1, r2 and th come from bounds_2d.yaml:fixed, not from u.
    """
    bounds = bounds or load_bounds()
    v, f = bounds["variables"], bounds["fixed"]

    u = list(u)
    if len(u) != len(VARS):
        raise ValueError(
            f"expected {len(VARS)} free variables {VARS}, got {len(u)}. "
            f"({FIXED_VARS} are fixed in bounds_2d.yaml and are not part of u.)")

    def lerp(name, uu):
        lo, hi = v[name]["min"], v[name]["max"]
        uu = min(1.0, max(0.0, float(uu)))
        return lo + uu * (hi - lo)

    phys = {name: lerp(name, uu) for name, uu in zip(VARS, u)}
    phys.update({name: float(f[name]) for name in FIXED_VARS})
    return Geometry2D(**phys)


def clearances(g: Geometry2D) -> Dict[str, float]:
    """Characteristic distances of the cell [m]. Pure geometry, no thresholds.

    Separated out from check_feasibility so tests can pin the FORMULAS against
    the MATLAB reference design independently of whatever thresholds
    bounds_2d.yaml happens to carry.

    in_radius : center-to-edge distance of the rhombic primitive cell. The
        rhombus has side `a` with 60/120 deg angles, so its height is
        a*sin(60) = a*sqrt(3)/2 and the center-to-edge distance is HALF that,
        a*sqrt(3)/4. (An earlier version used a*sqrt(3)/2 -- the in-radius of
        a hexagon of *side* a -- which is 2x too large. The hole sits at the
        centroid, so this is the distance a leg has to reach the wall.)
    d_tip_face : radial extent of a leg (the flat tip face).
    d_tip_corner : distance to the leg tip's sharp CORNER, which is the point
        of the hole actually closest to the cell wall. w is transverse to the
        leg (rect size [w r], base center, offset r/2 radially -- see
        buildBoomerangUnitCell.m:34-36), so the corner is at hypot(r, w/2),
        NOT at r + w/2.
    d_junction : distance to the inner corners where the three legs meet.
    web : minimum solid material between the hole and the cell wall.
    """
    in_radius = (math.sqrt(3.0) / 4.0) * g.a
    d_tip_corner = math.hypot(g.r, g.w / 2.0)
    return dict(
        in_radius=in_radius,
        d_tip_face=g.r,
        d_tip_corner=d_tip_corner,
        d_junction=g.w / 2.0,
        web=in_radius - d_tip_corner,
    )


def check_feasibility(g: Geometry2D, bounds: dict | None = None):
    """Return (ok: bool, reasons: list[str]). Cheap pre-simulation gate.

    Three groups of checks:

    1. Fabricability: leg width w and radial length r must each clear
       min_feature; fillet radii cannot exceed half the leg width.
    2. Web: solid material between the hole and the cell wall must clear
       min_web -- see clearances() for the geometry. web <= 0 means the leg
       punches through the wall, which is also the topology change that
       invalidates any hard-coded COMSOL boundary indices downstream.
    3. Fillet-selection validity: the MATLAB builder picks fillet vertices
       with two DISK ANNULI of fixed half-width (addFillet(),
       buildBoomerangUnitCell.m:119-149), and that scheme silently mis-selects
       over part of bounds_2d.yaml's box. Requiring the three conditions below
       keeps every candidate inside the region where the annuli resolve the
       junction and tip corners cleanly:
         d_tip_corner > w                  (else the junction annulus grabs
                                            the tips)
         d_tip_corner <= r + sw            (else the tip annulus misses the
                                            tips)
         w/2 + sqrt(3)*r1 < r - sw         (else the tip annulus grabs the
                                            junctions)
       The sqrt(3)*r1 term in the third check matters: h_fil1 runs BEFORE
       h_disksel2 is evaluated, and filleting a vertex whose edges meet at
       interior angle theta moves the tangent points outward by r1/tan(theta/2)
       along each edge -- ~sqrt(3)*r1 for the 60 deg junction. Comparing the
       pre-fillet distance w/2 would pass cases whose post-fillet vertices land
       inside the tip annulus. Verify the 60 deg assumption against mphgeom on
       the first rebuild.
       If you ever replace the annulus selection with a geometry-derived one,
       relax these together with it -- they encode the builder, not physics.

    STILL NOT cross-validated against COMSOL's actual Boolean geometry
    (mphgeom) at the bounds' extremes. min_web in particular is calibrated to
    the MATLAB reference design rather than derived -- see bounds_2d.yaml.
    """
    bounds = bounds or load_bounds()
    f = bounds["feasibility"]
    min_feature = f["min_feature"]
    min_web = f["min_web"]
    sw = f["fillet_select_halfwidth"]
    reasons: List[str] = []

    nm = lambda x: x * 1e9  # noqa: E731

    if g.w < min_feature:
        reasons.append(f"leg width w {nm(g.w):.0f} nm < min {nm(min_feature):.0f} nm")
    if g.r < min_feature:
        reasons.append(f"leg length r {nm(g.r):.0f} nm < min {nm(min_feature):.0f} nm")
    if g.r1 > g.w / 2.0:
        reasons.append(f"junction fillet r1 {nm(g.r1):.0f} nm exceeds w/2 "
                       f"{nm(g.w/2):.0f} nm")
    if g.r2 > g.w / 2.0:
        reasons.append(f"tip fillet r2 {nm(g.r2):.0f} nm exceeds w/2 "
                       f"{nm(g.w/2):.0f} nm")

    c = clearances(g)
    if c["web"] < min_web:
        reasons.append(f"web (cell wall to leg tip corner) {nm(c['web']):.1f} nm "
                       f"< min {nm(min_web):.1f} nm")

    if c["d_tip_corner"] <= g.w:
        reasons.append(f"fillet selection: tip corner {nm(c['d_tip_corner']):.0f} nm "
                       f"<= w {nm(g.w):.0f} nm -- junction annulus would also "
                       f"select the tips")
    if c["d_tip_corner"] > g.r + sw:
        reasons.append(f"fillet selection: tip corner {nm(c['d_tip_corner']):.0f} nm "
                       f"> r+sw {nm(g.r + sw):.0f} nm -- tip annulus would miss "
                       f"the tips")
    d_junction_filleted = g.w / 2.0 + math.sqrt(3.0) * g.r1
    if d_junction_filleted >= g.r - sw:
        reasons.append(f"fillet selection: post-fillet junction "
                       f"{nm(d_junction_filleted):.0f} nm >= r-sw "
                       f"{nm(g.r - sw):.0f} nm -- tip annulus would also select "
                       f"the junctions")

    return (len(reasons) == 0), reasons


# The MATLAB reference design (test_Boomerang.m:9-16). Kept here rather than in
# a test so both tests and ad-hoc debugging use the same anchor.
REFERENCE = dict(a=480e-9, w=140e-9, r=177e-9, r1=10e-9, r2=10e-9, th=220e-9)


if __name__ == "__main__":
    for label, params in (("u=0.5 midpoint", None), ("MATLAB reference", REFERENCE)):
        g = Geometry2D(**params) if params else u_to_geometry([0.5] * len(VARS))
        ok, why = check_feasibility(g)
        print(f"--- {label} ---")
        print("  geometry:", {k: f"{v*1e9:.1f} nm" for k, v in g.as_dict().items()})
        print("  clearances:", {k: f"{v*1e9:.1f} nm" for k, v in clearances(g).items()})
        print("  feasible:", ok, why)
