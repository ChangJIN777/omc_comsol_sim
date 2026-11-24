"""Band structure of 1D diamond photonic-crystal unit cells with Tidy3D (FDTD, Bloch boundaries).

Self-contained: needs only tidy3d, numpy and matplotlib. It does not import anything else from this
repository, so the file can be copied to another machine and run there.

Method (same as the Tidy3D band-structure example and rib_v2_cavity.make_unitCell_sims in cavity_base.py):
one lattice period along x with a Bloch boundary, PML in y and z, a cloud of randomly placed / phased
point dipoles, point FieldTimeMonitors recorded after the pulse has decayed, and harmonic inversion
(ResonanceFinder) of the time traces. One simulation per k point, Gamma -> X (k = 0 ... 0.5 in units of
2*pi/a). The band edges are the two lowest guided resonances (below the light line) at X.

All lengths are in um, all frequencies given on the command line are in THz.

SAFE BY DEFAULT: without --run nothing is uploaded and no FlexCredits are spent.

Examples
--------
  python diamond_bandStructure_sim.py                        # build + check + geometry PNG, no cloud
  python diamond_bandStructure_sim.py --geometry holey       # another unit cell, still no cloud
  python diamond_bandStructure_sim.py --xsect wedge          # triangular underside (angled etch)
  python diamond_bandStructure_sim.py --estimate-cost        # uploads drafts, prints the cost, deletes them
  python diamond_bandStructure_sim.py --run --x-only         # one simulation (X point): gap only
  python diamond_bandStructure_sim.py --run                  # full Gamma -> X batch (Nk+1 simulations)
  python diamond_bandStructure_sim.py --load data/<run_name> # re-analyze downloaded results, no cloud
  python diamond_bandStructure_sim.py --sweep t 0.15 0.40 6            # sweep plan only, no cloud
  python diamond_bandStructure_sim.py --sweep t 0.15 0.40 6 --run --x-only --max-credits 2

Use as a module
---------------
  import diamond_bandStructure_sim as bs
  p = bs.BandParams(geometry="rib_block", t=0.3)
  sims = bs.build_simulations(p)        # {"sim_0": td.Simulation, ...}; nothing is run

Adding a geometry
-----------------
Write a function `my_cell(p: BandParams) -> UnitCell` and decorate it with @register_geometry("name").
"""

import argparse
import csv
import json
import os
import sys
from dataclasses import dataclass, field, fields, asdict, replace
from datetime import datetime
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
import tidy3d as td
import tidy3d.web as web  # explicit: `import tidy3d` alone does not load the web API (no td.web attribute)
from tidy3d.plugins.resonance import ResonanceFinder

C_0 = td.constants.C_0  # um/s
NM = 1e-3  # nm -> um

# refractive index of diamond (same values as `diamond_telecom` / `diamond` in cavity_base.py)
MATERIAL_INDEX = {
    "diamond_telecom": 2.386,
    "diamond_visible": 2.4064,
}

# mode -> (dipole / monitor field component, simulation symmetry)
# Tidy3D convention: -1 = odd = PEC plane, +1 = even = PMC plane, 0 = no symmetry.
#   TE: y = 0 PEC (Ey even in y), mid-thickness z = 0 PMC (Ex, Ey even in z: TE-like)
#   TM: y = 0 PMC (Ey odd in y),  mid-thickness z = 0 PEC (Ez even in z: TM-like)
MODES = {
    "TE": ("Ey", (0, -1, 1)),
    "TM": ("Ez", (0, 1, -1)),
}
XSECTS = ("rect", "wedge")


# --------------------------------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------------------------------
def _f(default, help, **kw):
    return field(default=default, metadata=dict(help=help, **kw))


@dataclass
class BandParams:
    """Every knob of the calculation. Each field is also a command line option (--<field name>)."""

    # what to simulate
    geometry: str = _f("rib_block", "unit cell type (see --list-geometries)")
    material: str = _f("diamond_telecom", "diamond_telecom (n = 2.386) or diamond_visible (n = 2.4064)",
                       choices=tuple(MATERIAL_INDEX))
    index: Optional[float] = _f(None, "refractive index; overrides --material when given", type=float)
    mode: str = _f("TE", "TE: Ey dipoles, symmetry (0,-1,1). TM: Ez dipoles, symmetry (0,1,-1)",
                   choices=tuple(MODES))
    # cross-section (rib_block and holey; rib_cos6 is rect only)
    xsect: str = _f("rect", "rect: slab only. wedge: slab + triangular underside of an angled etch; it has no "
                            "mid-thickness mirror, so the z symmetry is switched off (symmetry (0,sy,0))",
                    choices=XSECTS)
    wedge_angle: float = _f(30.0, "wedge: base angle (deg). Under a strip of width b the apex is b/2*tan(angle) "
                                  "below the slab (30 deg: b/(2*sqrt(3)), apex angle 120 deg)")
    wedge_layers: int = _f(12, "rib_block wedge: number of stacked layers that approximate the roof under the "
                               "filleted outline; the depth error is (deepest point)/(2*wedge_layers)")
    # shared dimensions (um)
    a: float = _f(0.55, "lattice constant (um)")
    thickness: float = _f(0.38, "slab thickness (um)")
    # rib_block / rib_cos6 (um)
    t: float = _f(0.275, "rib_block: block length along x (um), 0 <= t <= a")
    w: float = _f(0.9, "rib_block / rib_cos6: TOTAL tip-to-tip width in y, spine included (um)")
    spine_width: float = _f(0.2, "rib_block / rib_cos6: spine width (um)")
    fillet_r1: float = _f(0.0, "rib_block: concave fillet radius at the block/spine junction (um)")
    fillet_r2: float = _f(0.0, "rib_block: convex fillet radius at the block tip corners (um)")
    fillet_points: int = _f(8, "rib_block: number of vertices per fillet arc")
    # holey (um). hole_x / hole_y are FULL axes of the ellipse (hx / hy of holey_cavity in cavity_base.py)
    beam_width: float = _f(0.9, "holey: beam width in y (um)")
    hole_x: float = _f(0.275, "holey: full axis of the elliptical hole along x (um)")
    hole_y: float = _f(0.5, "holey: full axis of the elliptical hole along y (um)")
    hole_points: int = _f(128, "holey: number of vertices of the hole polygon")
    # frequency window and k sampling
    target_freq_thz: float = _f(193.4, "target frequency (THz); centre of the source and of the window")
    freq_span_thz: float = _f(100.0, "half-width of the frequency window (THz)")
    Nk: int = _f(8, "number of k intervals from Gamma to X (Nk+1 simulations)")
    # numerics
    min_steps_per_wvl: float = _f(20.0, "AutoGrid resolution")
    num_dipoles: int = _f(10, "number of random point dipoles")
    num_monitors: int = _f(5, "number of random point time monitors")
    dPML: float = _f(1.0, "distance from the structure to the PML, in units of the longest wavelength")
    runtime_fwidth: float = _f(100.0, "run time in units of 1/(source bandwidth)")
    t_start_fwidth: float = _f(2.0, "monitor start time in units of 1/(source bandwidth)")
    seed: int = _f(12345, "random seed of the dipole / monitor positions")
    # resonance filter and figure of merit
    minQ: float = _f(1.0, "drop resonances with |Q| below this")
    minamp: float = _f(1.0, "drop resonances with amplitude below this")
    maxerr: float = _f(1000.0, "drop resonances with fit error above this")
    freq_tolerance_thz: float = _f(3.0, "width of the Gaussian midgap penalty in the fitness (THz)")

    # derived quantities ---------------------------------------------------------------------------
    @property
    def n(self) -> float:
        return float(self.index) if self.index is not None else MATERIAL_INDEX[self.material]

    @property
    def medium(self) -> td.Medium:
        return td.Medium(permittivity=self.n**2)

    @property
    def polarization(self) -> str:
        return MODES[self.mode][0]

    @property
    def symmetry(self) -> Tuple[int, int, int]:
        sx, sy, sz = MODES[self.mode][1]
        # the wedge breaks the mid-thickness mirror: keep the y symmetry only (all thickness parities)
        return (sx, sy, 0) if self.xsect == "wedge" else (sx, sy, sz)

    @property
    def target_freq(self) -> float:  # Hz
        return _thz_to_hz(self.target_freq_thz)

    @property
    def freq_span(self) -> float:  # Hz
        return _thz_to_hz(self.freq_span_thz)

    @property
    def freq_range(self) -> Tuple[float, float]:  # Hz
        return (self.target_freq - self.freq_span, self.target_freq + self.freq_span)

    @property
    def f_light_X(self) -> float:  # Hz, light line at the X point
        return C_0 / (2 * self.a)


def _thz_to_hz(x: float) -> float:
    """exact decimal conversion: 193.4 -> 193.4e12 (193.4*1e12 can be off by one floating point ulp)"""
    return float(repr(float(x)) + "e12")


# --------------------------------------------------------------------------------------------------
# geometries
# --------------------------------------------------------------------------------------------------
@dataclass
class UnitCell:
    """What a geometry function returns."""

    structures: List[td.Structure]
    width_y: float  # extent of the structure in y (um); sets the simulation size
    src_lo: Sequence[float]  # (x, y, z) lower corner of the box where dipoles and monitors are placed
    src_hi: Sequence[float]  # upper corner; keep it inside the dielectric and in the y >= 0 half
    snap: List[tuple] = field(default_factory=list)  # AutoGrid snapping points (x, y, z); None = free
    fill_factor: float = float("nan")  # dielectric area / (a * width_y)
    outline: List[np.ndarray] = field(default_factory=list)  # polygons of one period, for plotting only
    depth_below: float = 0.0  # how far the structure reaches below the slab (wedge), um


GEOMETRIES: Dict[str, Callable[[BandParams], UnitCell]] = {}


def register_geometry(name: str):
    def deco(fn):
        GEOMETRIES[name] = fn
        return fn
    return deco


def _polygon_area(verts: np.ndarray) -> float:
    x, y = verts[:, 0], verts[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _slab(verts: np.ndarray, thickness: float, below: float = 0.0) -> td.PolySlab:
    """vertical prism of the slab; `below` extends it under the slab (a hole cut through a wedge)"""
    return td.PolySlab(vertices=verts, slab_bounds=(-thickness / 2 - below, thickness / 2), axis=2,
                       sidewall_angle=0, reference_plane="bottom")


def wedge_depth(width: float, angle_deg: float) -> float:
    """depth of the wedge apex below the slab for a strip of this width (um)"""
    return 0.5 * width * float(np.tan(np.radians(angle_deg)))


def _wedge_prism(axis: int, center: float, width: float, thickness: float, angle_deg: float,
                 bounds: Tuple[float, float]) -> td.PolySlab:
    """triangular prism under the slab: base = the slab's bottom face over `width`, apex wedge_depth below.

    axis 0: ridge along x (the beam), the triangle is drawn in (y, z).
    axis 1: ridge along y, the triangle is drawn in (x, z). PolySlab vertices are given in the plane normal
    to `axis`, remaining axes in ascending order. `center` is the ridge position along the first of them.
    """
    z0 = -thickness / 2
    verts = [(center - width / 2, z0), (center + width / 2, z0), (center, z0 - wedge_depth(width, angle_deg))]
    return td.PolySlab(vertices=verts, slab_bounds=bounds, axis=axis, sidewall_angle=0, reference_plane="bottom")


def rib_roof_distance(x, y, a, t, w, s, r1=0.0, r2=0.0):
    """in-plane distance d (um) from (x, y) to the filleted rib_block outline; > 0 inside. Vectorized.

    x is measured from the centre of a block (any period). The wedge of an angled etch is d*tan(wedge_angle)
    deep ("roof" of the outline): d = b/2 - |offset| under a straight strip of width b, a cone d = rho - r1
    around each concave root fillet (and over the whole block/spine junction, which is therefore the deepest
    place), a cone d = r2 - rho inside each convex tip fillet, and a hip plane d = w/2 - |y| at the tip.
    """
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    x = (x + a / 2) % a - a / 2  # nearest block
    u = t / 2 - np.abs(x)  # to the block flank line
    v = s / 2 - np.abs(y)  # to the spine edge line
    f = w / 2 - np.abs(y)  # to the tip line
    d = np.minimum(np.hypot(np.maximum(u + r1, 0), np.maximum(v + r1, 0)) - r1, f)
    corner = (u >= 0) & (u < r2) & (f < r2)  # inside the square that holds a tip fillet
    return np.where(corner, r2 - np.hypot(r2 - u, r2 - f), d)


def rib_roof_max_distance(t, w, s, r1=0.0) -> float:
    """largest d of rib_roof_distance: under the centre of the block/spine junction"""
    return float(min(w / 2, np.hypot(t / 2 + r1, s / 2 + r1) - r1))


def rib_roof_level_polygons(x0, a, t, w, s, r1, r2, delta, n_arc=8) -> List[np.ndarray]:
    """outline (N x 2 vertices, um) of the region rib_roof_distance >= delta of the cell centred at x = x0.

    This is the filleted outline offset inwards by delta: spine half-width s/2 - delta, block half-length
    t/2 - delta, root fillets of radius r1 + delta about the SAME centres, tip fillets of radius r2 - delta
    (sharp once delta >= r2), tip at w/2 - delta. Returns [] when nothing is left, one polygon that spans the
    period while the spine strip exists, otherwise one island around the block (the root arcs then meet in a
    cusp on the beam axis). Exact up to the n_arc chords per arc.
    """
    hs, ht, y_top = s / 2 - delta, t / 2 - delta, w / 2 - delta
    cx, cy, rho = t / 2 + r1, s / 2 + r1, r1 + delta  # root arc: centre (cx, cy) in the first quadrant
    if y_top <= 0:
        return []
    # first quadrant path, from the outer end (cell edge or beam axis) to the y axis; y never decreases
    th0 = -np.pi / 2 if hs > 0 else -np.pi + np.arcsin(min(cy / rho, 1.0))
    th1 = -np.pi if ht > 0 else -np.pi + np.arccos(min(cx / rho, 1.0))
    if th0 < th1 - 1e-12:  # the arc lies outside the quadrant: nothing is left
        return []
    path = [(a / 2, hs)] if hs > 0 else []
    ang = np.linspace(th0, th1, max(int(n_arc), 2) + 1) if rho > 0 else [th0]
    path += [(cx + rho * np.cos(q), cy + rho * np.sin(q)) for q in ang]
    if ht > 0:  # flank, tip fillet, tip
        rt = r2 - delta
        if rt > 0 and w / 2 - r2 >= cy:
            cq = np.linspace(0, np.pi / 2, max(int(n_arc), 2) + 1)
            path += [(t / 2 - r2 + rt * np.cos(q), w / 2 - r2 + rt * np.sin(q)) for q in cq]
        else:
            path.append((ht, max(y_top, cy)))
    pts = np.array(path, dtype=float)
    pts[:, 0] = np.clip(pts[:, 0], 0, None)
    pts[:, 1] = np.clip(pts[:, 1], 0 if hs <= 0 else hs, None)
    # cut at the tip plane y = y_top (short blocks: the tip can cut the junction arc), then close on x = 0
    if pts[:, 1].max() > y_top:
        k = int(np.argmax(pts[:, 1] > y_top))  # first point above the tip
        if k == 0:
            return []
        (xa, ya), (xb, yb) = pts[k - 1], pts[k]
        cut = (xa + (xb - xa) * (y_top - ya) / (yb - ya), y_top)
        pts = np.vstack([pts[:k], cut])
    if pts[-1, 0] > 1e-12:
        pts = np.vstack([pts, (0.0, pts[-1, 1] if ht <= 0 else y_top)])
    mirror_x = pts[::-1] * np.array([-1, 1])
    if hs > 0:  # spans the period: upper edge from -a/2 to a/2, then the lower edge back
        upper = np.vstack([pts * np.array([-1, 1]), pts[::-1]])
        poly = np.vstack([upper, upper[::-1] * np.array([1, -1])])
    else:  # island: counter-clockwise from the cusp on the +x axis
        poly = np.vstack([pts, mirror_x, pts * np.array([-1, -1]), pts[::-1] * np.array([1, -1])])
    keep = np.linalg.norm(poly - np.roll(poly, 1, axis=0), axis=1) > 1e-9  # drop duplicated vertices
    poly = poly[keep]
    if len(poly) < 3 or _polygon_area(poly) < 1e-12:
        return []
    poly[:, 0] += x0
    return [poly]


def _src_z(p: BandParams) -> float:
    """half-range in z of the dipole / monitor positions.

    rect: 0, everything on the mid-plane (which is a symmetry plane). wedge: there is no symmetry plane, and
    sources on one plane would excite modes with a node there only weakly, so they are spread over the slab.
    """
    return 0.3 * p.thickness if p.xsect == "wedge" else 0.0


def _check_common(p: BandParams):
    if p.a <= 0:
        raise ValueError("the lattice constant must be positive (a=%.4f)" % p.a)
    if p.thickness <= 0:
        raise ValueError("the slab thickness must be positive (thickness=%.4f)" % p.thickness)
    if p.xsect not in XSECTS:
        raise ValueError("unknown xsect '%s'; choose from %s" % (p.xsect, ", ".join(XSECTS)))
    if p.xsect == "wedge" and not 0 < p.wedge_angle < 90:
        raise ValueError("wedge_angle must be between 0 and 90 deg (wedge_angle=%.4f)" % p.wedge_angle)


def round_polygon_corners(vertices, radii, n_arc=8) -> np.ndarray:
    """replace every vertex i with radii[i] > 0 by a circular arc tangent to both adjacent edges"""
    verts = np.asarray(vertices, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n = len(verts)
    if not np.any(radii > 0):
        return verts
    d = np.zeros(n)  # tangent length of each fillet along its two edges
    geo = [None] * n
    for i in range(n):
        if radii[i] <= 0:
            continue
        P, A, B = verts[i], verts[i - 1], verts[(i + 1) % n]
        u1 = (A - P) / np.linalg.norm(A - P)
        u2 = (B - P) / np.linalg.norm(B - P)
        theta = np.arccos(np.clip(np.dot(u1, u2), -1, 1))  # angle between the two edges
        if theta < 1e-9 or np.pi - theta < 1e-9:  # spike or straight line: nothing to round
            continue
        d[i] = radii[i] / np.tan(theta / 2)
        bisector = (u1 + u2) / np.linalg.norm(u1 + u2)
        center = P + bisector * radii[i] / np.sin(theta / 2)
        geo[i] = (u1, u2, center, radii[i])
    for i in range(n):  # each edge must hold the tangent lengths of both of its end fillets
        j = (i + 1) % n
        L = np.linalg.norm(verts[j] - verts[i])
        if d[i] + d[j] > L + 1e-9:
            raise ValueError("fillet radii too large for the unit cell "
                             "(an edge of %.4f um would need %.4f um)" % (L, d[i] + d[j]))
    out = []
    for i in range(n):
        if geo[i] is None:
            out.append(verts[i])
            continue
        u1, u2, center, r = geo[i]
        T1 = verts[i] + u1 * d[i]
        T2 = verts[i] + u2 * d[i]
        ang1 = np.arctan2(T1[1] - center[1], T1[0] - center[0])
        ang2 = np.arctan2(T2[1] - center[1], T2[0] - center[0])
        dang = (ang2 - ang1 + np.pi) % (2 * np.pi) - np.pi  # the short way round
        angs = ang1 + dang * np.linspace(0, 1, max(int(n_arc), 2))
        out.extend(np.stack([center[0] + r * np.cos(angs), center[1] + r * np.sin(angs)], axis=1))
    out = np.array(out)
    keep = np.linalg.norm(out - np.roll(out, 1, axis=0), axis=1) > 1e-9  # drop duplicated vertices
    return out[keep]


def rib_block_vertices(x0, a, t, w, s, r1=0.0, r2=0.0, n_arc=8) -> np.ndarray:
    """vertices (N x 2, um) of one block-rib cell centred at x = x0.

    spine: width s across the full period a.  block: length t (x), TOTAL tip-to-tip width w (y).
    r1: concave fillet where the block meets the spine.  r2: convex fillet at the block tip corners.
    """
    tol = 1e-6  # um
    if a <= 0:
        raise ValueError("rib_block: the lattice constant must be positive (a=%.4f)" % a)
    if t < -tol or t > a + tol:
        raise ValueError("rib_block: need 0 <= t <= a (t=%.4f, a=%.4f)" % (t, a))
    if s <= 0:
        raise ValueError("rib_block: the spine width must be positive (spine_width=%.4f)" % s)
    if w < s - tol:
        raise ValueError("rib_block: w (total tip-to-tip width) must be >= spine_width "
                         "(w=%.4f, spine_width=%.4f)" % (w, s))
    if r1 < 0 or r2 < 0:
        raise ValueError("rib_block: fillet radii must be >= 0 (r1=%.4f, r2=%.4f)" % (r1, r2))
    t = min(max(t, 0), a)
    has_block = (t > tol) and (w > s + tol)
    if has_block and t >= a - tol:  # the block fills the period: plain strip of width w
        upper = [(-a / 2, w / 2), (a / 2, w / 2)]
        r_upper = [0, 0]
    elif has_block:
        upper = [(-a / 2, s / 2), (-t / 2, s / 2), (-t / 2, w / 2),
                 (t / 2, w / 2), (t / 2, s / 2), (a / 2, s / 2)]
        r_upper = [0, r1, r2, r2, r1, 0]
    else:  # no block: plain strip of width s
        upper = [(-a / 2, s / 2), (a / 2, s / 2)]
        r_upper = [0, 0]
    upper = np.array(upper, dtype=float)
    lower = upper[::-1] * np.array([1, -1])
    arr = np.concatenate([upper, lower])
    radii = np.concatenate([r_upper, r_upper[::-1]]).astype(float)
    try:
        arr = round_polygon_corners(arr, radii, n_arc)
    except ValueError as err:
        raise ValueError("rib_block: %s; need r1 + r2 <= (w - spine_width)/2 = %.4f, 2*r2 <= t = %.4f and "
                         "r1 <= (a - t)/2 = %.4f" % (err, (w - s) / 2, t, (a - t) / 2)) from None
    arr[:, 0] += x0
    return arr


@register_geometry("rib_block")
def rib_block_cell(p: BandParams) -> UnitCell:
    """spine + centred rectangular block with optional fillets (rib_v2_cavity in cavity_base.py)"""
    _check_common(p)
    a, t, w, s = p.a, p.t, p.w, p.spine_width
    # three cells, so no polygon edge lies on the Bloch boundary
    polys = [rib_block_vertices(x, a, t, w, s, p.fillet_r1, p.fillet_r2, p.fillet_points) for x in (-a, 0, a)]
    geometries = [_slab(v, p.thickness) for v in polys]
    has_block = t > 1e-6 and w > s + 1e-6
    depth = 0.0
    if p.xsect == "wedge" and not has_block:
        # plain strip: the roof is the triangular prism (exact)
        geometries.append(_wedge_prism(0, 0.0, s, p.thickness, p.wedge_angle, (-1.5 * a, 1.5 * a)))
        depth = wedge_depth(s, p.wedge_angle)
    elif p.xsect == "wedge":
        # Roof of the FILLETED outline: the wedge is rib_roof_distance*tan(wedge_angle) deep, so its edges
        # follow the root and tip fillets (hipped tip, cone under each fillet and under the junction).
        # Built as wedge_layers vertical PolySlabs, layer i being the outline offset inwards by the roof
        # distance of the MIDDLE of the layer (rib_roof_level_polygons, exact offsets, no sidewall_angle: the
        # offset outline changes topology when the spine or the block closes, which a slanted PolySlab cannot
        # do). Staircase error: at most half a layer, depth/(2*wedge_layers), in depth; none in plan view.
        if int(p.wedge_layers) < 1:
            raise ValueError("wedge_layers must be >= 1 (wedge_layers=%s)" % p.wedge_layers)
        n_layers, tan_a = int(p.wedge_layers), float(np.tan(np.radians(p.wedge_angle)))
        d_max = rib_roof_max_distance(t, w, s, p.fillet_r1)
        depth = d_max * tan_a
        z0 = -p.thickness / 2
        for i in range(n_layers):
            d_lo, d_hi = d_max * i / n_layers, d_max * (i + 1) / n_layers
            for x in (-a, 0, a):
                for v in rib_roof_level_polygons(x, a, t, w, s, p.fillet_r1, p.fillet_r2, (d_lo + d_hi) / 2,
                                                 p.fillet_points):
                    geometries.append(td.PolySlab(vertices=v, slab_bounds=(z0 - d_hi * tan_a, z0 - d_lo * tan_a),
                                                  axis=2, sidewall_angle=0, reference_plane="bottom"))
    structure = td.Structure(geometry=td.GeometryGroup(geometries=geometries), medium=p.medium, name="rib_v2")
    # dipoles / monitors inside the block, y > 0 half (the y < 0 half is the symmetry image)
    x_half = max(t, 0.2 * a) / 2
    y_top = w / 2 if has_block else s / 2  # no block: stay inside the spine
    lo, hi = [-x_half, 0.05 * y_top, -_src_z(p)], [x_half, 0.95 * y_top, _src_z(p)]
    # AutoGrid only snaps to bounding boxes, so snap grid lines to the block and spine faces. Fillets do
    # not move these faces; the arcs themselves are not refined (raise min_steps_per_wvl to resolve them).
    snap = [(sx * t / 2, sy * y, None) for sx in (-1, 1) for sy in (-1, 1) for y in (w / 2, s / 2)]
    return UnitCell([structure], w, lo, hi, snap, _polygon_area(polys[1]) / (a * w), [polys[1]], depth)


def rib_cos6_vertices(x0, a, hy, s, n_pts=501) -> np.ndarray:
    """vertices of one cos^6 rib cell (rib_cavity.gen_ribUnitCell in cavity_base.py); hy = bump amplitude"""
    x = np.linspace(-a / 2, a / 2, n_pts)
    upper = hy * np.cos(np.pi * x / a) ** 6 + s / 2
    lower = -hy * np.cos(np.pi * np.flip(x) / a) ** 6 - s / 2
    return np.array([np.concatenate([x, np.flip(x)]) + x0, np.concatenate([upper, lower])]).T


@register_geometry("rib_cos6")
def rib_cos6_cell(p: BandParams) -> UnitCell:
    """spine + cos^6 bumps; w is the total tip-to-tip width, so the bump amplitude is (w - spine_width)/2"""
    _check_common(p)
    if p.xsect != "rect":
        raise ValueError("rib_cos6: only xsect rect is available (xsect=%s)" % p.xsect)
    a, w, s = p.a, p.w, p.spine_width
    if s <= 0:
        raise ValueError("rib_cos6: the spine width must be positive (spine_width=%.4f)" % s)
    if w < s:
        raise ValueError("rib_cos6: w (total tip-to-tip width) must be >= spine_width "
                         "(w=%.4f, spine_width=%.4f)" % (w, s))
    hy = (w - s) / 2
    polys = [rib_cos6_vertices(x, a, hy, s) for x in (-a, 0, a)]
    structure = td.Structure(geometry=td.GeometryGroup(geometries=[_slab(v, p.thickness) for v in polys]),
                             medium=p.medium, name="rib_cos6")
    # at |x| = 0.1 a the bump is still 0.74 hy high, so this box is inside the dielectric
    y_top = s / 2 + 0.7 * hy
    lo, hi = [-0.1 * a, 0.05 * y_top, 0], [0.1 * a, 0.95 * y_top, 0]
    snap = [(0, sy * y, None) for sy in (-1, 1) for y in (w / 2, s / 2)]
    return UnitCell([structure], w, lo, hi, snap, _polygon_area(polys[1]) / (a * w), [polys[1]])


@register_geometry("holey")
def holey_cell(p: BandParams) -> UnitCell:
    """beam of width beam_width with one elliptical air hole per period (holey_cavity in cavity_base.py)"""
    _check_common(p)
    a, bw, hx, hy = p.a, p.beam_width, p.hole_x, p.hole_y
    if bw <= 0:
        raise ValueError("holey: the beam width must be positive (beam_width=%.4f)" % bw)
    if hx <= 0 or hy <= 0:
        raise ValueError("holey: the hole axes must be positive (hole_x=%.4f, hole_y=%.4f)" % (hx, hy))
    if hx >= a:
        raise ValueError("holey: need hole_x < a, otherwise neighbouring holes merge (hole_x=%.4f, a=%.4f)" % (hx, a))
    if hy >= bw:
        raise ValueError("holey: need hole_y < beam_width, otherwise the beam is cut "
                         "(hole_y=%.4f, beam_width=%.4f)" % (hy, bw))
    slab = td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, bw, p.thickness)), medium=p.medium, name="slab")
    structures = [slab]
    depth = 0.0
    if p.xsect == "wedge":
        # one prism along the beam under the slab, base = the beam width; the hole is cut through both
        depth = wedge_depth(bw, p.wedge_angle)
        structures.append(td.Structure(geometry=_wedge_prism(0, 0.0, bw, p.thickness, p.wedge_angle, (-1.5 * a, 1.5 * a)),
                                       medium=p.medium, name="wedge"))
    ang = np.linspace(0, 2 * np.pi, int(p.hole_points), endpoint=False)
    hole_verts = np.stack([hx / 2 * np.cos(ang), hy / 2 * np.sin(ang)], axis=1)
    # the hole comes last: it overrides the slab and the wedge
    structures.append(td.Structure(geometry=_slab(hole_verts, p.thickness, below=1.05 * depth), medium=td.Medium(), name="hole"))
    # whole period, y > 0 half. Some points fall inside the hole on purpose: the air band lives there.
    lo, hi = [-0.45 * a, 0.05 * bw / 2, -_src_z(p)], [0.45 * a, 0.95 * bw / 2, _src_z(p)]
    snap = [(sx * hx / 2, sy * y, None) for sx in (-1, 1) for sy in (-1, 1) for y in (bw / 2, hy / 2)]
    beam = np.array([(-a / 2, -bw / 2), (a / 2, -bw / 2), (a / 2, bw / 2), (-a / 2, bw / 2)])
    return UnitCell(structures, bw, lo, hi, snap, 1 - _polygon_area(hole_verts) / (a * bw), [beam, hole_verts], depth)


def build_unit_cell(p: BandParams) -> UnitCell:
    """validate the parameters and build the unit cell of p.geometry"""
    if p.geometry not in GEOMETRIES:
        raise ValueError("unknown geometry '%s'; choose from %s" % (p.geometry, ", ".join(GEOMETRIES)))
    if p.mode not in MODES:
        raise ValueError("unknown mode '%s'; choose from %s" % (p.mode, ", ".join(MODES)))
    if p.index is None and p.material not in MATERIAL_INDEX:
        raise ValueError("unknown material '%s'; choose from %s" % (p.material, ", ".join(MATERIAL_INDEX)))
    if p.n < 1:
        raise ValueError("the refractive index must be >= 1 (n=%.4f)" % p.n)
    cell = GEOMETRIES[p.geometry](p)
    if not all(l <= h for l, h in zip(cell.src_lo, cell.src_hi)):
        raise ValueError("%s: empty dipole / monitor region %s .. %s" % (p.geometry, cell.src_lo, cell.src_hi))
    return cell


# --------------------------------------------------------------------------------------------------
# simulations
# --------------------------------------------------------------------------------------------------
def build_simulations(p: BandParams, x_only: bool = False, cell: Optional[UnitCell] = None) -> Dict[str, td.Simulation]:
    """build (do NOT run) the Bloch simulations, Gamma -> X. Returns {"sim_i": td.Simulation}, k = i/(2 Nk).

    x_only=True returns only the X point ("sim_<Nk>"), which is enough to extract the band gap.
    """
    if int(p.Nk) < 1:
        raise ValueError("Nk must be >= 1 (Nk=%s)" % p.Nk)
    if p.num_dipoles < 1 or p.num_monitors < 1:
        raise ValueError("need at least one dipole and one monitor")
    Nk = int(p.Nk)
    cell = build_unit_cell(p) if cell is None else cell
    freq0 = p.target_freq
    freq_range = p.freq_range
    if freq_range[0] <= 0:
        raise ValueError("freq_span_thz must be smaller than target_freq_thz")
    lambda_max = C_0 / freq_range[0]  # um
    freqw = 2 * p.freq_span * 0.3  # Hz, width of the Gaussian pulse
    run_time = p.runtime_fwidth / freqw
    t_start = p.t_start_fwidth / freqw
    spacing = p.dPML * lambda_max
    # the domain stays centred on the slab; a wedge of depth d below it is kept `spacing` away from the PML
    sim_size = (p.a, cell.width_y + 2 * spacing, p.thickness + 2 * cell.depth_below + 2 * spacing)
    polarization = p.polarization

    rng = np.random.default_rng(p.seed)
    dipole_positions = rng.uniform(cell.src_lo, cell.src_hi, [p.num_dipoles, 3])
    dipole_phases = rng.uniform(0, 2 * np.pi, p.num_dipoles)
    dipoles = [td.PointDipole(source_time=td.GaussianPulse(freq0=freq0, fwidth=freqw, phase=dipole_phases[i]),
                              center=tuple(dipole_positions[i]), polarization=polarization, name="dipole_" + str(i))
               for i in range(p.num_dipoles)]
    monitor_positions = rng.uniform(cell.src_lo, cell.src_hi, [p.num_monitors, 3])
    monitors_time = [td.FieldTimeMonitor(fields=[polarization], center=tuple(monitor_positions[i]), size=(0, 0, 0),
                                         start=t_start, name="monitor_time_" + str(i))
                     for i in range(p.num_monitors)]
    grid_spec = td.GridSpec.auto(min_steps_per_wvl=p.min_steps_per_wvl, snapping_points=cell.snap)
    sims = {}
    for i in ([Nk] if x_only else range(Nk + 1)):
        sims[f"sim_{i}"] = td.Simulation(
            center=(0, 0, 0),
            size=sim_size,
            grid_spec=grid_spec,
            structures=cell.structures,
            sources=dipoles,
            monitors=monitors_time,
            run_time=run_time,
            shutoff=0,  # never stop early: the resonances must ring for the whole run time
            boundary_spec=td.BoundarySpec(x=td.Boundary.bloch((1 / 2) * i / Nk), y=td.Boundary.pml(), z=td.Boundary.pml()),
            normalize_index=None,
            symmetry=p.symmetry,
        )
    return sims


def bloch_k(sim: td.Simulation) -> float:
    """Bloch wavevector along x of a simulation, in units of 2*pi/a"""
    return float(getattr(sim.boundary_spec.x.plus, "bloch_vec", 0.0))


def grid_report(p: BandParams, cell: UnitCell, sims: Dict[str, td.Simulation]) -> str:
    """text summary of the geometry, the grid and the size of the job (no cloud access)"""
    sim = next(iter(sims.values()))
    b = sim.grid.boundaries
    xs, ys, zs = np.asarray(b.x), np.asarray(b.y), np.asarray(b.z)

    def steps(coords, half):  # grid steps (nm) inside |coordinate| <= half
        c = coords[np.abs(coords) <= half + 1e-9]
        d = np.diff(c) / NM
        return (d.min(), d.max()) if len(d) else (float("nan"), float("nan"))

    dx, dy, dz = steps(xs, p.a / 2), steps(ys, cell.width_y / 2), steps(zs, p.thickness / 2)
    worst = 0.0
    for pt in cell.snap:
        for coord, lines in zip(pt, (xs, ys, zs)):
            if coord is not None:
                worst = max(worst, float(np.min(np.abs(lines - coord))))
    lines = [
        "geometry        : %s (%s, n = %.4f, fill factor %.3f)" % (p.geometry, p.mode, p.n, cell.fill_factor),
        "cross-section   : %s" % ("rect" if p.xsect == "rect" else
                                  "wedge, base angle %.1f deg, %.1f nm deep" % (p.wedge_angle, cell.depth_below / NM)),
        "polarization    : %s dipoles, symmetry %s" % (p.polarization, str(p.symmetry)),
        "frequency window: %.1f .. %.1f THz (target %.1f THz); light line at X: %.1f THz"
        % (p.freq_range[0] * 1e-12, p.freq_range[1] * 1e-12, p.target_freq_thz, p.f_light_X * 1e-12),
        "domain (um)     : %.3f x %.3f x %.3f" % tuple(sim.size),
        "grid step (nm)  : dx %.1f..%.1f, dy %.1f..%.1f, dz %.1f..%.1f inside the structure"
        % (dx[0], dx[1], dy[0], dy[1], dz[0], dz[1]),
        "cells per period: %d along x" % int(np.sum(np.abs(xs) <= p.a / 2 + 1e-9) - 1),
        "snapping points : %d, worst distance to a grid line %.3g nm" % (len(cell.snap), worst / NM),
        "grid cells      : %.3g per simulation (full domain, before the symmetry reduction)" % sim.num_cells,
        "time steps      : %d (run time %.2f ps)" % (sim.num_time_steps, sim.run_time * 1e12),
        "simulations     : %d (k = %s)" % (len(sims), ", ".join("%.4g" % bloch_k(s) for s in sims.values())),
    ]
    if p.f_light_X >= p.freq_range[1]:
        lines.append("note            : the whole frequency window is below the light line at X")
    if p.f_light_X <= p.freq_range[0]:
        lines.append("WARNING         : the light line at X is below the frequency window: no guided band can be found")
    return "\n".join(lines)


def plot_geometry(p: BandParams, cell: UnitCell, sim: td.Simulation, path: str):
    """save a PNG with the unit cell outline, the z = 0 grid and the x = 0 cross-section"""
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), tight_layout=True)
    for i, poly in enumerate(cell.outline):
        closed = np.vstack([poly, poly[:1]]) / NM
        ax[0].fill(closed[:, 0], closed[:, 1], color="#4c78a8" if i == 0 else "white", alpha=0.6 if i == 0 else 1.0)
        ax[0].plot(closed[:, 0], closed[:, 1], color="#1f3b57", lw=1)
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x (nm)")
    ax[0].set_ylabel("y (nm)")
    ax[0].set_title("%s unit cell (fill factor %.3f)" % (p.geometry, cell.fill_factor))
    sim.plot(z=0.0, ax=ax[1])
    sim.plot_grid(z=0.0, ax=ax[1])
    ax[1].set_xlim(-p.a / 2, p.a / 2)
    ax[1].set_ylim(-0.75 * cell.width_y, 0.75 * cell.width_y)
    ax[1].set_title("z = 0: grid, dipoles, monitors")
    sim.plot(x=0.0, ax=ax[2])
    ax[2].set_title("x = 0 cross-section")
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# cloud (the only functions that talk to the Tidy3D server)
# --------------------------------------------------------------------------------------------------
class Budget:
    """running total of the estimated FlexCredits of everything started by this process"""

    def __init__(self, max_credits: Optional[float] = None):
        self.max_credits = max_credits
        self.estimated = 0.0
        self.last_estimate = 0.0  # estimate of the most recent batch, started or not


def run_cloud(sims: Dict[str, td.Simulation], run_dir: str, label: str, folder_name: str = "default",
              budget: Optional[Budget] = None, estimate_only: bool = False, verbose: bool = True,
              on_stage: Optional[Callable[[str], None]] = None):
    """upload the simulations as one Batch, print the maximum cost, then run and download to run_dir.

    Returns {"sim_i": td.SimulationData}, or None when estimate_only (the drafts are deleted again).
    COSTS FLEXCREDITS unless estimate_only. The estimate needs the upload, so it needs an API key.

    verbose: passed to web.Batch (Tidy3D's own upload / status / download progress) and switches the
    cost line printed here. on_stage: optional callable(message), told each step (uploading, estimated
    cost, running, loading); when given, the cost line goes to it instead of print. The defaults
    (verbose=True, on_stage=None) are the behaviour of the command line tool.
    """
    def say(message):
        if on_stage is not None:
            on_stage(message)

    budget = Budget() if budget is None else budget
    tasks = {"%s_%s" % (label, name): sim for name, sim in sims.items()}
    n = len(tasks)
    batch = web.Batch(simulations=tasks, folder_name=folder_name, verbose=bool(verbose))
    say("uploading %d %s to the Tidy3D cloud (folder %s)" % (n, "draft(s)" if estimate_only else "simulation(s)",
                                                           folder_name))
    batch.upload()
    cost = float(batch.estimate_cost(verbose=False))
    budget.last_estimate = cost
    if on_stage is not None:
        say("estimated maximum cost of %d simulation(s): %.3f FlexCredits" % (n, cost))
    elif verbose:
        print("estimated maximum cost of %d simulation(s): %.3f FlexCredits" % (n, cost))
    over = budget.max_credits is not None and budget.estimated + cost > budget.max_credits
    if estimate_only or over:
        try:
            batch.delete()
            say("drafts deleted, nothing was started")
        except Exception as err:  # the drafts cost nothing; just tell the user they are still there
            print("could not delete the uploaded drafts (%s); remove them in the web interface" % err)
        if over and not estimate_only:
            raise RuntimeError("stopping before start: %.3f (already started) + %.3f FlexCredits exceeds "
                               "--max-credits %.3f" % (budget.estimated, cost, budget.max_credits))
        return None
    budget.estimated += cost
    os.makedirs(run_dir, exist_ok=True)
    # one call starts, monitors and downloads every task; Tidy3D reports the status of each task itself
    say("running %d task(s) on the Tidy3D cloud, then downloading to %s" % (n, run_dir))
    batch_data = batch.run(path_dir=run_dir)
    out = {}
    for i, name in enumerate(sims):
        say("loading task %d/%d (%s_%s)" % (i + 1, n, label, name))
        out[name] = batch_data["%s_%s" % (label, name)]
    return out


def load_run(run_dir: str) -> Dict[str, td.SimulationData]:
    """load every SimulationData .hdf5 in run_dir (no cloud). Keys are "k=<bloch vector>".

    The k point is read from the Bloch boundary stored in each file, so the file names do not matter.
    If two files have the same k, the most recently modified one is used.
    """
    if not os.path.isdir(run_dir):
        raise FileNotFoundError("no such directory: %s" % run_dir)
    found = {}
    for name in sorted(os.listdir(run_dir)):
        path = os.path.join(run_dir, name)
        if not name.endswith(".hdf5"):
            continue
        try:
            sim_data = td.SimulationData.from_file(path)
        except Exception:
            continue  # not a simulation result (e.g. the batch.hdf5 bookkeeping file)
        key = "k=%.6f" % bloch_k(sim_data.simulation)
        mtime = os.path.getmtime(path)
        if key in found:
            print("warning: more than one result for %s in %s; using the newest" % (key, run_dir))
            if found[key][0] > mtime:
                continue
        found[key] = (mtime, sim_data)
    if not found:
        raise FileNotFoundError("no simulation results (*.hdf5) found in %s" % run_dir)
    return {key: sd for key, (_, sd) in sorted(found.items())}


# --------------------------------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------------------------------
@dataclass
class BandResult:
    # per k point: {"k": float, "freq": Hz array (sorted), "Q": array, "amplitude": array, "error": array}
    bands: List[dict]
    f_light_X: float  # Hz
    has_X: bool  # False if no simulation at k = 0.5 was available
    f_lower: float = float("nan")  # THz, dielectric band edge at X
    f_upper: float = float("nan")  # THz, air band edge at X
    bandGap: float = 0.0  # THz
    midGap: float = float("nan")  # THz
    fitness: float = 0.0  # -(gap/midgap) * exp(-((midgap - target)/tolerance)^2); more negative is better


def filter_resonances(resonance_data, minQ, minamp, maxerr):
    """drop unphysical / weak / poorly-fitted resonances from a ResonanceFinder dataset"""
    resonance_data = resonance_data.where(resonance_data.freq > 0, drop=True)
    resonance_data = resonance_data.where(abs(resonance_data.Q) > minQ, drop=True)
    resonance_data = resonance_data.where(resonance_data.amplitude > minamp, drop=True)
    resonance_data = resonance_data.where(resonance_data.error < maxerr, drop=True)
    return resonance_data


def find_resonances_full(signals, p: BandParams) -> Dict[str, np.ndarray]:
    """{"freq" (Hz, ascending), "Q", "amplitude", "error"} of the filtered resonances of FieldTimeData"""
    rd = ResonanceFinder(freq_window=tuple(p.freq_range)).run(signals=signals)
    rd = filter_resonances(rd, p.minQ, p.minamp, p.maxerr)
    cols = {name: np.atleast_1d(getattr(rd, name).to_numpy()).astype(float)
            for name in ("freq", "Q", "amplitude", "error")}
    order = np.argsort(cols["freq"])
    return {name: values[order] for name, values in cols.items()}


def find_resonances(signals, p: BandParams) -> Tuple[np.ndarray, np.ndarray]:
    """(freq in Hz sorted ascending, Q) of the filtered resonances of a tuple of FieldTimeData"""
    res = find_resonances_full(signals, p)
    return res["freq"], res["Q"]


def band_edges(freqs_X, p: BandParams) -> Tuple[float, float, float, float, float]:
    """(f_lower, f_upper, bandGap, midGap, fitness) in THz from the resonances at X (Hz).

    The band edges are the two lowest resonances below the light line c/(2a). Returns
    (nan, nan, 0, nan, 0) when fewer than two guided resonances exist.
    """
    fX = np.sort(np.asarray(freqs_X, dtype=float))
    fX = fX[fX < p.f_light_X]
    if len(fX) < 2:
        return float("nan"), float("nan"), 0.0, float("nan"), 0.0
    f_lower, f_upper = fX[0] * 1e-12, fX[1] * 1e-12
    bandGap = f_upper - f_lower
    midGap = (f_upper + f_lower) / 2
    fitness = -(bandGap / midGap) * np.exp(-(abs(midGap - p.target_freq_thz) / p.freq_tolerance_thz) ** 2)
    return f_lower, f_upper, bandGap, midGap, float(fitness)


def analyze_bands(sim_datas: Dict[str, td.SimulationData], p: BandParams) -> BandResult:
    """resonances of every k point plus the band gap at X. sim_datas: any dict of SimulationData."""
    bands = []
    for sim_data in sim_datas.values():
        res = find_resonances_full(sim_data.data, p)
        bands.append({"k": bloch_k(sim_data.simulation), "freq": res["freq"], "Q": res["Q"],
                      "amplitude": res["amplitude"], "error": res["error"]})
    bands.sort(key=lambda d: d["k"])
    result = BandResult(bands=bands, f_light_X=p.f_light_X, has_X=False)
    at_X = [d for d in bands if abs(d["k"] - 0.5) < 1e-9]
    if at_X:
        result.has_X = True
        result.f_lower, result.f_upper, result.bandGap, result.midGap, result.fitness = band_edges(at_X[0]["freq"], p)
    return result


def summarize(result: BandResult) -> str:
    if not result.has_X:
        return "no simulation at the X point (k = 0.5): band gap not extracted"
    if not np.isfinite(result.midGap):
        return "fewer than two guided resonances at X: no band gap extracted (bandGap = 0, midGap = nan, fitness = 0)"
    return ("band edges at X: %.2f / %.2f THz -> gap %.2f THz, midgap %.2f THz, fitness %.4f"
            % (result.f_lower, result.f_upper, result.bandGap, result.midGap, result.fitness))


def plot_bands(result: BandResult, p: BandParams, path: str, title: Optional[str] = None):
    """band diagram (frequency in THz against k in units of 2*pi/a) saved as a PNG"""
    fig, ax = plt.subplots(figsize=(5.5, 5), tight_layout=True)
    f_top = p.freq_range[1] * 1e-12
    k = np.linspace(0, 0.5, 100)
    light = k * C_0 / p.a * 1e-12
    ax.fill_between(k, light, max(f_top, light[-1]), color="0.85", lw=0, label="light cone")
    ax.plot(k, light, color="0.5", lw=1)
    guided_k, guided_f, leaky_k, leaky_f = [], [], [], []
    for d in result.bands:
        for f in d["freq"]:
            below = f < d["k"] * C_0 / p.a
            (guided_k if below else leaky_k).append(d["k"])
            (guided_f if below else leaky_f).append(f * 1e-12)
    ax.scatter(leaky_k, leaky_f, s=22, clip_on=False, zorder=3, facecolors="none", edgecolors="#4c78a8", label="resonance above the light line")
    ax.scatter(guided_k, guided_f, s=22, clip_on=False, zorder=3, color="#1f4e79", label="guided resonance")
    if np.isfinite(result.midGap):
        ax.axhspan(result.f_lower, result.f_upper, color="#f2c14e", alpha=0.35, lw=0,
                   label="gap at X: %.1f THz" % result.bandGap)
    ax.axhline(p.target_freq_thz, color="#c0392b", ls=":", lw=1.5, label="target %.1f THz" % p.target_freq_thz)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, f_top)
    ax.set_xticks([0, 0.5], [r"$\Gamma$", "X"])
    ax.set_xlabel(r"Wavevector $k_x$")
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(title or "Band diagram (%s, %s, a = %.0f nm)" % (p.geometry, p.mode, p.a / NM))
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# output files
# --------------------------------------------------------------------------------------------------
CSV_HEADER = ["timestamp", "geometry", "mode", "n", "a", "t", "w", "spine_width", "thickness", "fillet_r1",
              "fillet_r2", "beam_width", "hole_x", "hole_y", "filling_ratio", "f_lower_THz", "f_upper_THz",
              "bandGap_THz", "midGap_THz", "fitness", "target_freq_THz", "Nk", "min_steps_per_wvl", "run_dir"]


def save_bandgap_row(csv_path: str, p: BandParams, fill_factor: float, result: BandResult, run_dir: str = "") -> str:
    """append one row to csv_path (lengths in um); the header is written when the file is new.

    If an existing file has a different header, the row goes to a new time-stamped file instead, so a
    finished (paid) simulation never loses its result. Returns the path actually written.
    """
    row = [datetime.now().isoformat(timespec="seconds"), p.geometry, p.mode, p.n, p.a, p.t, p.w, p.spine_width,
           p.thickness, p.fillet_r1, p.fillet_r2, p.beam_width, p.hole_x, p.hole_y, fill_factor, result.f_lower,
           result.f_upper, result.bandGap, result.midGap, result.fitness, p.target_freq_thz, p.Nk,
           p.min_steps_per_wvl, run_dir]
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    new = not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0
    if not new:
        with open(csv_path, newline="") as f:
            if next(csv.reader(f), None) != CSV_HEADER:
                root, ext = os.path.splitext(csv_path)
                csv_path = "%s_%s%s" % (root, datetime.now().strftime("%Y%m%d_%H%M%S"), ext)
                print("warning: existing CSV has a different header; writing to %s" % csv_path)
                new = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        if new:
            writer.writerow(CSV_HEADER)
        writer.writerow(row)
    return csv_path


def save_bands_csv(path: str, result: BandResult, p: BandParams):
    """every resonance of every k point: k, frequency, Q, whether it is below the light line, and the
    amplitude and fit error of the harmonic inversion (the last two columns were added later: files without
    them, and bands without these arrays, are still valid)"""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["k_2pi_over_a", "freq_THz", "Q", "below_light_line", "amplitude", "error"])
        for d in result.bands:
            n = len(d["freq"])
            amp = d.get("amplitude", [""] * n)
            err = d.get("error", [""] * n)
            for i, (freq, q) in enumerate(zip(d["freq"], d["Q"])):
                writer.writerow([d["k"], freq * 1e-12, q, int(freq < d["k"] * C_0 / p.a), amp[i], err[i]])


# --------------------------------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------------------------------
def _field_type(f):
    return f.metadata.get("type") or type(f.default)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Band structure of a 1D diamond photonic-crystal unit cell with Tidy3D. "
                    "Lengths in um, frequencies in THz. Nothing is submitted without --run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples: see the docstring at the top of this file.")
    group = parser.add_argument_group("parameters (defaults in brackets)")
    for f in fields(BandParams):
        kw = dict(type=_field_type(f), default=argparse.SUPPRESS, help="%s [%s]" % (f.metadata["help"], f.default))
        if f.name == "geometry":
            kw["choices"] = tuple(GEOMETRIES)
        elif "choices" in f.metadata:
            kw["choices"] = f.metadata["choices"]
        names = ["--" + f.name]
        if "_" in f.name:
            names.append("--" + f.name.replace("_", "-"))
        group.add_argument(*names, dest=f.name, **kw)
    act = parser.add_argument_group("actions")
    act.add_argument("--run", action="store_true", help="submit to the Tidy3D cloud (COSTS FlexCredits)")
    act.add_argument("--x-only", action="store_true", help="only the X point (1 simulation): gives the band gap")
    act.add_argument("--estimate-cost", action="store_true",
                     help="upload drafts, print the maximum cost, delete the drafts; nothing is started")
    act.add_argument("--max-credits", type=float, default=None,
                     help="with --run: refuse to start once the summed cost estimates would exceed this")
    act.add_argument("--load", metavar="DIR", default=None, help="re-analyze downloaded results in DIR, no cloud")
    act.add_argument("--sweep", nargs=4, metavar=("PARAM", "START", "STOP", "NUM"), default=None,
                     help="sweep one numeric parameter (linspace); one CSV row per point. Plan only without --run")
    act.add_argument("--list-geometries", action="store_true", help="list the unit cell types and exit")
    out = parser.add_argument_group("output")
    out.add_argument("--data-dir", default="data", help="where results go [data]")
    out.add_argument("--run-name", default=None, help="sub-directory / task label [<geometry>_<mode>_<time>]")
    out.add_argument("--csv", default=None, help="band gap CSV [<data-dir>/<geometry>_bandGap.csv]")
    out.add_argument("--folder-name", default="default", help="Tidy3D web folder of the tasks [default]")
    out.add_argument("--no-plots", action="store_true", help="do not write PNG files")
    return parser


def params_from_args(args, base: Optional[BandParams] = None) -> BandParams:
    """BandParams from the options that were actually given on the command line, on top of `base`"""
    given = {f.name: getattr(args, f.name) for f in fields(BandParams) if hasattr(args, f.name)}
    return replace(base or BandParams(), **given)


def _sweep_points(args) -> Tuple[Optional[str], List[Optional[float]]]:
    if args.sweep is None:
        return None, [None]
    name, start, stop, num = args.sweep
    name = name.replace("-", "_")
    types = {f.name: _field_type(f) for f in fields(BandParams)}
    if types.get(name) is not float:
        raise ValueError("--sweep: '%s' is not a float parameter; choose from %s"
                         % (name, ", ".join(k for k, v in types.items() if v is float)))
    if int(num) < 1:
        raise ValueError("--sweep: NUM must be >= 1")
    return name, [float(v) for v in np.linspace(float(start), float(stop), int(num))]


def _load_saved_params(run_dir: str) -> Optional[BandParams]:
    path = os.path.join(run_dir, "params.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        saved = json.load(f)
    known = {f.name for f in fields(BandParams)}
    return BandParams(**{k: v for k, v in saved.items() if k in known})


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.list_geometries:
        for name, fn in GEOMETRIES.items():
            print("%-10s %s" % (name, (fn.__doc__ or "").strip()))
        return 0
    if args.run and (args.load or args.estimate_cost):
        raise SystemExit("--run cannot be combined with --load or --estimate-cost")
    if args.load and args.estimate_cost:
        raise SystemExit("--load cannot be combined with --estimate-cost")
    plt.switch_backend("Agg")  # files only; set here so that importing this module changes nothing
    try:
        td.config.logging.level = "ERROR"
    except AttributeError:  # tidy3d < 2.10
        td.config.logging_level = "ERROR"

    sweep_name, sweep_values = _sweep_points(args)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = params_from_args(args)
    run_name = args.run_name or "%s_%s_%s" % (base.geometry, base.mode, stamp)
    root_dir = args.load if args.load else os.path.join(args.data_dir, run_name)
    cloud = args.run or args.estimate_cost

    # 1) build everything first: a bad parameter must fail before any credit is spent
    plan = []
    for value in sweep_values:
        label = run_name if sweep_name is None else "%s_%s_%.6g" % (run_name, sweep_name, value)
        run_dir = root_dir if sweep_name is None else os.path.join(root_dir, "%s_%.6g" % (sweep_name, value))
        saved = _load_saved_params(run_dir) if args.load else None
        p = params_from_args(args, saved)
        if sweep_name is not None and saved is None:
            p = replace(p, **{sweep_name: value})
        try:
            cell = build_unit_cell(p)
            sims = build_simulations(p, x_only=args.x_only, cell=cell)
        except Exception as err:  # pydantic validation errors included
            where = "" if sweep_name is None else " at %s = %.6g" % (sweep_name, value)
            raise SystemExit("invalid parameters%s: %s" % (where, err))
        plan.append((label, run_dir, p, cell, sims))
    n_sims = sum(len(sims) for *_, sims in plan)

    # 2) report
    label, run_dir, p, cell, sims = plan[0]
    if not args.load:
        print(grid_report(p, cell, sims))
    if sweep_name is not None:
        print("sweep           : %s = %s (%d points)" % (sweep_name, ", ".join("%.6g" % v for v in sweep_values), len(plan)))
    if not args.load:
        print("TOTAL           : %d simulation(s) in %d batch(es)" % (n_sims, len(plan)))

    # 3) dry run: geometry figure, then stop
    if not cloud and not args.load:
        if not args.no_plots:
            preview_dir = os.path.join(args.data_dir, "preview")
            os.makedirs(preview_dir, exist_ok=True)
            png = os.path.join(preview_dir, "%s_%s_geometry.png" % (p.geometry, p.mode))
            plot_geometry(p, cell, next(iter(sims.values())), png)
            print("geometry figure : %s" % png)
        print("DRY RUN: nothing was uploaded or run. Add --run to submit (%d simulation(s)), "
              "--estimate-cost for the price, --x-only for the X point alone." % n_sims)
        return 0

    # 4) run / load, analyze, save
    budget = Budget(args.max_credits)
    csv_path = args.csv or os.path.join(args.data_dir, "%s_bandGap.csv" % base.geometry)
    if args.run:
        print("SUBMITTING %d simulation(s) to the Tidy3D cloud (FlexCredits will be charged)" % n_sims)
    for label, run_dir, p, cell, sims in plan:
        if sweep_name is not None:
            print("--- %s = %.6g ---" % (sweep_name, getattr(p, sweep_name)))
        if args.load:
            try:
                sim_datas = load_run(run_dir)
            except FileNotFoundError as err:
                if sweep_name is None:
                    raise SystemExit(str(err))
                print("skipped: %s" % err)
                continue
            a_data = float(next(iter(sim_datas.values())).simulation.size[0])
            if abs(a_data - p.a) > 1e-9:
                print("warning: the loaded simulations have a = %.4f um, not %.4f um; using the loaded value "
                      "(pass the original parameters, or keep params.json next to the data)" % (a_data, p.a))
                p = replace(p, a=a_data)
                try:
                    cell = build_unit_cell(p)
                except ValueError:
                    cell = replace(cell, fill_factor=float("nan"))
        else:
            sim_datas = run_cloud(sims, run_dir, label, args.folder_name, budget, estimate_only=args.estimate_cost)
            if sim_datas is None:
                continue
            with open(os.path.join(run_dir, "params.json"), "w") as f:
                json.dump(asdict(p), f, indent=2)
        result = analyze_bands(sim_datas, p)
        print(summarize(result))
        save_bands_csv(os.path.join(run_dir, "bands.csv"), result, p)
        if not args.no_plots:
            png = os.path.join(run_dir, "band_diagram.png")
            plot_bands(result, p, png)
            print("band diagram    : %s" % png)
        if result.has_X:
            print("band gap row    : %s" % save_bandgap_row(csv_path, p, cell.fill_factor, result, run_dir))
    if args.estimate_cost:
        print("ESTIMATE ONLY: nothing was started.")
    elif args.run:
        print("sum of the cost estimates of the started batches: %.3f FlexCredits" % budget.estimated)
    return 0


if __name__ == "__main__":
    sys.exit(main())
