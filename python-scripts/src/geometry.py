"""Geometry mapping and feasibility checks for the 1D diamond OMC unit cell.

Design variables (normalized u in [0,1]) -> physical parameters (meters).
Rectangular cross-section nanobeam (width w, thickness t) with one elliptical
hole per period a (semi-axes hx along beam axis z, hy along width y).

Keep this module dependency-light (stdlib + PyYAML) so it can be imported by
every backend (MPB, COMSOL, surrogate) and by the optimizer.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, asdict
from typing import Dict, List

import yaml

_CFG_DIR = os.path.join(os.path.dirname(__file__), "..", "configs")

# Order of the optimizer's normalized vector u.
# u[4] = t (thickness) added in campaign 2; old 4-D results use u4=0 → t=t_min.
VARS: List[str] = ["a", "w", "hx", "hy", "t"]

# Cell types. "hole" is the original elliptical-hole cell; every function below
# defaults to it, so existing callers are unchanged. "rib" is the spine +
# centred-rib cell of comsol/build_omc_unitcell_rib.m.
CELLTYPES = ("hole", "rib")

# Rib parameters in canonical order. These ARE the global parameter names of
# comsol/omc_unitcell_rib.mph. NOTE the clash with the hole cell: for the rib,
# t is the RIB WIDTH along the beam and th is the slab thickness (hole cell: t
# is the thickness). The searched subset (bounds.yaml: rib.variables) keeps
# this order, which is the order of the optimizer's u for celltype "rib".
RIB_PARAMS: List[str] = ["a", "w", "t", "s", "th", "r1", "r2"]


def load_bounds(path: str | None = None) -> dict:
    path = path or os.path.join(_CFG_DIR, "bounds.yaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


@dataclass
class Geometry:
    a: float        # lattice constant [m]
    w: float        # beam width [m]
    t: float        # beam thickness [m]
    hx: float       # hole semi-axis along beam axis z [m]
    hy: float       # hole semi-axis along width y [m]

    celltype = "hole"   # class attribute (not a dataclass field, not in as_dict)

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


@dataclass
class RibGeometry:
    """Spine + centred-rib cell of comsol/build_omc_unitcell_rib.m (meters).

    Field names are the template's global parameter names, so as_dict() is sent
    to COMSOL as is. Unlike Geometry, t is the rib width and th the thickness.
    """
    a: float        # period along the beam axis z (cell width) [m]
    w: float        # rib tip-to-tip width along y, = 2L [m]
    t: float        # rib width along the beam axis z [m]  (NOT the thickness)
    s: float        # spine width along y, full width [m]
    th: float       # slab thickness along x [m]
    r1: float       # concave rib/spine (root) fillet radius, > 0 [m]
    r2: float       # convex rib-tip fillet radius, > 0 [m]

    celltype = "rib"    # class attribute (not a dataclass field, not in as_dict)

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


def _rib_cfg(bounds: dict) -> dict:
    if "rib" not in bounds:
        raise KeyError("bounds have no 'rib' section (see configs/bounds.yaml)")
    return bounds["rib"]


def rib_search_vars(bounds: dict | None = None) -> List[str]:
    """Names of the searched rib variables, in RIB_PARAMS order (= order of u).

    Every rib parameter must appear exactly once, either under rib.variables
    ({min, max}, searched) or under rib.fixed (a value in meters).
    """
    bounds = bounds or load_bounds()
    rib = _rib_cfg(bounds)
    variables = rib.get("variables") or {}
    fixed = rib.get("fixed") or {}
    unknown = [n for n in list(variables) + list(fixed) if n not in RIB_PARAMS]
    both = [n for n in RIB_PARAMS if n in variables and n in fixed]
    missing = [n for n in RIB_PARAMS if n not in variables and n not in fixed]
    if unknown or both or missing:
        raise ValueError(
            "bounds 'rib' section: each of "
            f"{RIB_PARAMS} must be in exactly one of variables/fixed "
            f"(unknown: {unknown}, in both: {both}, missing: {missing})")
    return [n for n in RIB_PARAMS if n in variables]


def n_dim(celltype: str = "hole", bounds: dict | None = None) -> int:
    """Length of the optimizer's u for a cell type."""
    if celltype == "hole":
        return len(VARS)
    if celltype == "rib":
        return len(rib_search_vars(bounds))
    raise ValueError(f"unknown celltype {celltype!r}; expected one of {CELLTYPES}")


def _rib_phys_to_u(params: dict, bounds: dict) -> list:
    v = _rib_cfg(bounds)["variables"]
    u = []
    for name in rib_search_vars(bounds):
        lo, hi = float(v[name]["min"]), float(v[name]["max"])
        val = float(params.get(name, lo))
        u.append(min(1.0, max(0.0, (val - lo) / (hi - lo))))
    return u


def _u_to_rib_geometry(u, bounds: dict) -> RibGeometry:
    rib = _rib_cfg(bounds)
    names = rib_search_vars(bounds)
    u = list(u)
    if len(u) != len(names):
        raise ValueError(f"expected {len(names)} rib variables {names}, got {len(u)}")
    # float(): PyYAML reads an exponent without a decimal point (20e-9) as str
    phys = {name: float(val) for name, val in (rib.get("fixed") or {}).items()}
    for name, uu in zip(names, u):
        lo = float(rib["variables"][name]["min"])
        hi = float(rib["variables"][name]["max"])
        uu = min(1.0, max(0.0, float(uu)))
        phys[name] = lo + uu * (hi - lo)
    return RibGeometry(**phys)


def _check_rib_feasibility(g: RibGeometry):
    """Can the template BUILD this cell? (ok, reasons), reasons 'GEOM_*: ...'.

    Geometric validity only: the rules in the header of
    comsol/build_omc_unitcell_rib.m (= checkRibGeometry in
    omc_mechanical_livelink_rib.m, validateRibParams.m for d = 0).
    Fabrication limits are separate: check_fabricability.
    """
    nm = 1e9
    reasons: List[str] = []
    for name, val in g.as_dict().items():
        if not val > 0:                       # also catches NaN
            reasons.append(f"GEOM_NONPOSITIVE: {name} = {val*nm:.1f} nm must be > 0")
    if reasons:
        return False, reasons
    if g.t >= g.a:
        reasons.append(f"GEOM_RIBS_MERGE: rib width t {g.t*nm:.0f} nm >= period a {g.a*nm:.0f} nm")
    if g.s >= g.w:
        reasons.append(f"GEOM_NO_RIB: spine s {g.s*nm:.0f} nm >= width w {g.w*nm:.0f} nm")
    if 2.0 * g.r2 >= g.t:
        reasons.append(f"GEOM_TIP_FILLET: 2*r2 {2*g.r2*nm:.0f} nm >= rib width t {g.t*nm:.0f} nm")
    if g.t + 2.0 * g.r1 >= g.a:
        reasons.append(f"GEOM_ROOT_FILLET: t + 2*r1 {(g.t + 2*g.r1)*nm:.0f} nm >= period a {g.a*nm:.0f} nm")
    if g.r1 + g.r2 >= (g.w - g.s) / 2.0:
        reasons.append(f"GEOM_FILLETS_MEET: r1 + r2 {(g.r1 + g.r2)*nm:.0f} nm >= rib arm "
                       f"(w - s)/2 {(g.w - g.s)/2*nm:.0f} nm")
    return (len(reasons) == 0), reasons


def check_fabricability(g, bounds: dict | None = None):
    """Can the cell be FABRICATED? Returns (ok, failures); no solver needed.

    Separate from check_feasibility (can the template build it). Thresholds
    come from bounds.yaml rib.fabrication; a missing or null threshold disables
    its rule. Each failure is a dict(code, rule, value, limit, kind), with
    kind "min" (value must be >= limit) or "max" (value must be <= limit);
    lengths in meters, ratios dimensionless. Hole cells have no rules here
    (their min_feature rule lives in check_feasibility) and always pass.
    """
    bounds = bounds or load_bounds()
    if getattr(g, "celltype", "hole") != "rib":
        return True, []
    fab = _rib_cfg(bounds).get("fabrication") or {}

    def ratio(num, den):
        return num / den if den > 0 else float("inf")

    arm = (g.w - g.s) / 2.0          # free-standing rib length beyond the spine
    gap = g.a - g.t                  # etched slot between neighbouring ribs
    rules = [
        # code, description, value, threshold key, kind
        ("FAB_MIN_RIB_WIDTH", "rib width t", g.t, "min_feature", "min"),
        ("FAB_MIN_SPINE_WIDTH", "spine width s", g.s, "min_feature", "min"),
        ("FAB_MIN_RIB_ARM", "rib arm (w - s)/2", arm, "min_feature", "min"),
        ("FAB_MIN_GAP", "air gap between ribs a - t", gap, "min_gap", "min"),
        ("FAB_MIN_THICKNESS", "slab thickness th", g.th, "min_thickness", "min"),
        ("FAB_ETCH_ASPECT", "etch aspect ratio th/(a - t)",
         ratio(g.th, gap), "max_etch_aspect", "max"),
        ("FAB_SOLID_ASPECT", "solid aspect ratio th/min(t, s)",
         ratio(g.th, min(g.t, g.s)), "max_solid_aspect", "max"),
        ("FAB_RIB_SLENDERNESS", "rib slenderness ((w - s)/2)/t",
         ratio(arm, g.t), "max_rib_slenderness", "max"),
        ("FAB_RIB_OVERHANG", "rib overhang ((w - s)/2)/th",
         ratio(arm, g.th), "max_rib_overhang", "max"),
        ("FAB_RELEASE_WIDTH", "widest strip to undercut max(t, s)",
         max(g.t, g.s), "max_release_width", "max"),
    ]
    failures = []
    for code, rule, value, key, kind in rules:
        limit = fab.get(key)
        if limit is None:
            continue
        limit = float(limit)
        bad = value < limit if kind == "min" else value > limit
        if bad or value != value:
            failures.append(dict(code=code, rule=rule, value=float(value),
                                 limit=limit, kind=kind))
    return (len(failures) == 0), failures


def fab_failure_text(failure: dict) -> str:
    """One-line text of a check_fabricability failure (for `reasons`)."""
    op = "<" if failure["kind"] == "min" else ">"
    return (f"{failure['code']}: {failure['rule']} = {failure['value']:.4g} "
            f"{op} limit {failure['limit']:.4g}")


def phys_to_u(params: dict, bounds: dict | None = None,
              celltype: str = "hole") -> list:
    """Convert physical params dict (meters) → normalized u in [0,1]^N.

    Uses the CURRENT bounds, so the result is always consistent with whatever
    bounds.yaml says right now — no migration needed when bounds change.
    Clips silently to [0,1] for params outside the current range.
    """
    bounds = bounds or load_bounds()
    if celltype == "rib":
        return _rib_phys_to_u(params, bounds)
    if celltype != "hole":
        raise ValueError(f"unknown celltype {celltype!r}; expected one of {CELLTYPES}")
    v = bounds["variables"]
    u = []
    for name in VARS:
        lo, hi = v[name]["min"], v[name]["max"]
        val = float(params.get(name, lo))
        u.append(min(1.0, max(0.0, (val - lo) / (hi - lo))))
    return u


def u_to_geometry(u, bounds: dict | None = None, celltype: str = "hole"):
    """Map normalized u in [0,1]^N -> Geometry (physical meters).

    Accepts 4-D (legacy, t fixed at t_min) or 5-D (t as free variable).

    celltype="rib" returns a RibGeometry instead; u then has one entry per
    searched rib variable (rib_search_vars), the rest come from rib.fixed.
    """
    bounds = bounds or load_bounds()
    if celltype == "rib":
        return _u_to_rib_geometry(u, bounds)
    if celltype != "hole":
        raise ValueError(f"unknown celltype {celltype!r}; expected one of {CELLTYPES}")
    v = bounds["variables"]

    u = list(u)
    if len(u) == 4 and "t" in v:
        u = u + [0.0]   # backward-compat: u4=0 → t = t_min (see bounds.yaml)
    if len(u) != len(VARS):
        raise ValueError(f"expected {len(VARS)} variables {VARS}, got {len(u)}")

    def lerp(name, uu):
        lo, hi = v[name]["min"], v[name]["max"]
        uu = min(1.0, max(0.0, float(uu)))
        return lo + uu * (hi - lo)

    phys = {name: lerp(name, uu) for name, uu in zip(VARS, u)}
    return Geometry(**phys)


def check_feasibility(g: Geometry, bounds: dict | None = None):
    """Return (ok: bool, reasons: list[str]). Cheap pre-simulation gate.

    Hole cell (Geometry): the min-feature rules below. Rib cell (RibGeometry):
    template validity only (_check_rib_feasibility); its fabrication limits
    are in check_fabricability.
    """
    if getattr(g, "celltype", "hole") == "rib":
        return _check_rib_feasibility(g)
    bounds = bounds or load_bounds()
    f = bounds["feasibility"]
    reasons: List[str] = []

    # Single absolute minimum-feature-size rule (2026-07-02): the remaining
    # solid material — sidewall (w/2 - hy) and bridge (a - 2*hx) — must each
    # be at least min_feature. Replaces the old mix of a generic f_min check,
    # separate min_sidewall/min_bridge values, and percentage-based
    # hy_frac_max/hx_frac_max caps, which could disagree with each other and
    # made it hard to reason about what the actual minimum feature was for a
    # given geometry (a % of a shrinking w is not a fixed nm value).
    sidewall = g.w / 2.0 - g.hy
    bridge = g.a - 2.0 * g.hx
    min_feature = f["min_feature"]
    if sidewall < min_feature:
        reasons.append(f"sidewall {sidewall*1e9:.0f} nm < min {min_feature*1e9:.0f} nm")
    if bridge < min_feature:
        reasons.append(f"bridge {bridge*1e9:.0f} nm < min {min_feature*1e9:.0f} nm")

    return (len(reasons) == 0), reasons


def filling_fraction(g: Geometry) -> float:
    """Air fraction of the hole within the cell footprint (rough, 2D in-plane)."""
    import math
    cell_area = g.a * g.w
    hole_area = math.pi * g.hx * g.hy
    return hole_area / cell_area


if __name__ == "__main__":
    g = u_to_geometry([0.5, 0.5, 0.5, 0.5])
    ok, why = check_feasibility(g)
    print("geometry:", {k: f"{v*1e9:.1f} nm" for k, v in g.as_dict().items()})
    print("feasible:", ok, why)
    print(f"air fill fraction: {filling_fraction(g):.3f}")
