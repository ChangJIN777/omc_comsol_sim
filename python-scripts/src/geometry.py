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

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


def phys_to_u(params: dict, bounds: dict | None = None) -> list:
    """Convert physical params dict (meters) → normalized u in [0,1]^N.

    Uses the CURRENT bounds, so the result is always consistent with whatever
    bounds.yaml says right now — no migration needed when bounds change.
    Clips silently to [0,1] for params outside the current range.
    """
    bounds = bounds or load_bounds()
    v = bounds["variables"]
    u = []
    for name in VARS:
        lo, hi = v[name]["min"], v[name]["max"]
        val = float(params.get(name, lo))
        u.append(min(1.0, max(0.0, (val - lo) / (hi - lo))))
    return u


def u_to_geometry(u, bounds: dict | None = None) -> Geometry:
    """Map normalized u in [0,1]^N -> Geometry (physical meters).

    Accepts 4-D (legacy, t fixed at t_min) or 5-D (t as free variable).
    """
    bounds = bounds or load_bounds()
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
    """Return (ok: bool, reasons: list[str]). Cheap pre-simulation gate."""
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
