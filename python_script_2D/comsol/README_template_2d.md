# Building the 2D boomerang unit-cell template (`trusty_boomerang.mph`)

> **STATUS: the template on disk does not satisfy this recipe yet.**
> `comsol/trusty_boomerang.mph` (exported as `trusty_boomerang_script.m`) is a
> *solved snapshot of one fixed design*, not a parametric template:
> - only `k`, `a`, `kx`, `ky` are real parameters; every geometry dimension is
>   a hard-coded literal, so `w,r,r1,r2,th` cannot be driven at all (and `a`
>   feeds only `kx`/`ky`, not the geometry);
> - the studies are labelled `Study_symmetric` / `Study 1`, not
>   `mech evenz` / `mech oddz`;
> - the periodic and parity BCs use absolute face indices (`[1 22]`, `[2 9]`,
>   `[3]`) rather than geometric selections, so they will attach to the wrong
>   faces as soon as the geometry changes.
>
> `src/acoustic_comsol_2d.py` now refuses to run against it (see
> `_require_parameters` / `_require_studies`) instead of silently solving the
> frozen cell for every candidate. See the to-do list for the plan.

Build this ONCE in the COMSOL 6.2+ GUI on your Mac, then
`src/acoustic_comsol_2d.py` just sets parameters and solves. The geometry and
physics are adapted from the MATLAB reference implementation:
`omc-comsol-chang/buildBoomerangUnitCell.m` (geometry) and
`omc-comsol-chang/runBands_2D.m` (boundary conditions, Brillouin-zone
convention) -- both driven by `omc-comsol-chang/test_Boomerang.m`.

**This is a genuinely different physical setup from `python-scripts/comsol/
README_template.md`** (the 1D nanobeam project): periodicity is TRUE 2D
(Bloch wavevector in BOTH in-plane directions), and there's no in-plane
mirror-symmetry reduction (no y-even/y-odd split) -- the only symmetry family
split is the slab's Z-PARITY (evenz/oddz about its midplane). Read this
whole file before building; don't reuse the 1D template's node names as-is.

## 1. Global parameters (Model > Parameters)

| name | expr      | meaning |
|------|-----------|---------|
| a    | 480[nm]   | lattice constant (rhombic primitive cell side) |
| w    | 140[nm]   | width of each boomerang leg rectangle (transverse) |
| r    | 177[nm]   | radial length of each leg, from the hole center |
| r1   | 10[nm]    | fillet radius at the leg **JUNCTIONS** (inner corners) |
| r2   | 10[nm]    | fillet radius at the leg **TIPS** (outer corners) |
| th   | 220[nm]   | slab thickness (out-of-plane, z), FULL thickness |
| k    | 0         | Brillouin-zone path parameter in [0,3), swept by the study |
| kx   | *expr*    | `if(k<1,(-1/sqrt(3))*k*(pi/a),if(k<2,(1/sqrt(3))*pi/a*(k-2),0))` |
| ky   | *expr*    | `if(k<1,(pi/a)*k,if(k<2,(k+2)*pi/(3*a),(3-k)*4*pi/(3*a)))` |

Values above match `test_Boomerang.m`'s defaults, for a sanity-check
cross-comparison against the MATLAB pipeline once both are built.

**All six geometry parameters must exist in the template, but only `a`, `w`
and `r` are swept.** `r1`, `r2` and `th` are fixed on the Python side
(`configs/bounds_2d.yaml:fixed`) and are not part of the optimizer's `u`
vector; the driver still writes all six so the template and Python can never
disagree about what was solved. To optimize one of them, move its entry into
`bounds_2d.yaml:variables`, add it to `geometry2d.VARS`, and bump
`optimizer.N_DIM` — stored `u` vectors are not comparable across that change.

**`r1` and `r2` were documented backwards in earlier revisions of this file.**
The authority is `addFillet()` at `buildBoomerangUnitCell.m:119-149`:
`h_disksel1` is the annulus `[w/(2*sqrt(2)), w]` about the hole center, which
catches the vertices at `w/2` -- the **junction** corners -- and `h_fil1`
applies `r1`. `h_disksel2` is `[r-25nm, r+25nm]`, catching the vertices at
`hypot(r, w/2)` -- the **tip** corners -- and `h_fil2` applies `r2`. The
mistake is invisible while `r1 == r2` and silently swaps two design variables
the moment the optimizer separates them.

**Keep the `k` parameter and the Parametric Sweep** (this file previously said
the opposite). `kx`/`ky` must stay COMSOL *expressions* of `k`, exactly as in
`runBands_2D.m:63-64`, for two reasons: it makes the sweep a single solve per
parity rather than one solve per k-point (54 eigensolves per candidate at
`kpts=9` vs 92 for a Python-driven point-by-point loop), and it keeps the
Python and MATLAB pipelines byte-comparable. `src/acoustic_comsol_2d.py`
still contains the older point-by-point loop that *overwrites* `kx`/`ky` with
literals -- that must be replaced with "set `k`, solve once, read the
parametric dataset" as part of this work. Writing literals into `kx`/`ky`
destroys the expressions for the rest of the process (`comsol_client.py`
caches the model), and every point of the study's own sweep then solves the
same wavevector.

## 2. Geometry (one hexagonal unit cell)

Mirrors `buildBoomerangUnitCell.m`'s Work Plane construction:

1. **Base plane** (Work Plane, z = -th/2): a Polygon tracing the rhombic
   hexagonal-lattice cell with vertices `(0,0)`, `(a/2, a*sqrt(3)/2)`,
   `(3a/2, a*sqrt(3)/2)`, `(a,0)`. (This rhombus, not a regular hexagon, is
   the primitive cell of the triangular/hexagonal lattice used here.)
2. **Boomerang hole** -- 3 rectangles, each `w` x `r`, positioned with their
   long axis pointing radially from a common center point at
   `(a*(1/2+1/4), a*sqrt(3)/4)`, rotated 0 deg / 120 deg / 240 deg from each
   other (`rec_1`/`rec_2`/`rec_3` in `buildBoomerangUnitCell.m`). Boolean
   **Difference** them from the base polygon.
3. **Fillets**: round the leg **junctions** with radius `r1` and the leg
   **tips** with radius `r2` (see `addFillet()` in `buildBoomerangUnitCell.m`
   for the exact selection logic -- two disk annuli around the hole center).
   The MATLAB annulus half-width is a fixed 25 nm, which mis-selects over part
   of `configs/bounds_2d.yaml`'s box; `check_feasibility()` now rejects that
   region, but a geometry-derived half-width (e.g.
   `min(r - w/2, hypot(r,w/2) - r)/2`) would be strictly better.
4. **Extrude** the 2D profile through thickness `th`, then subtract the lower
   half so only `z` in `[0, th/2]` is meshed, with the parity BC on the `z=0`
   face. (This is what the exported model does: work plane at `-th/2`, extrude
   `th`, then `symZComp = ext1 - symZPlaneExt`. Both work-plane `z` values and
   both extrude distances must become expressions in `th`.)

Build and visually check this against `mphgeom` output from the MATLAB
`buildBoomerangUnitCell.m` (run `test_Boomerang.m`'s commented-out debug
block) before proceeding -- the exact vertex/fillet-selection indices are
finicky enough that a visual cross-check is worth the five minutes.

## 3. Material

Diamond, anisotropic (cubic): `C11=1076[GPa]`, `C12=125[GPa]`, `C44=578[GPa]`,
density `3500[kg/m^3]`. The authority is
`omc-comsol-chang/LoadMaterialParams.m:17,24-26`, which the exported `.mph`
matches; `configs/materials.yaml` disagreed on both `C44` and density until
corrected, and is currently dead config that nothing in `src/` reads. Set the material
coordinate system's in-plane rotation to `rxtal = 45[deg]` about z (matches
`test_Boomerang.m`'s `P.anisoMat=1, P.rxtal=45` and
`omc-comsol-chang/RotateXtalTensor.m`). **Unlike the 1D project** (where only
the beam axis orientation matters), here BOTH in-plane lattice directions x
and y carry Bloch periodicity, so this single rotation angle is the only
orientation degree of freedom -- there is no separate "beam axis" choice.

## 4. Mechanical studies -- TWO, one per z-parity

Physics: **Solid Mechanics**, applied to the extruded domain.

- **Periodic Condition > Floquet periodicity**, TWO independent pairs (one
  per hexagonal lattice-translation direction), both using
  `kFloquet = (kx, ky, 0)`:
  - Pair 1: the two parallel side-wall faces corresponding to translation by
    lattice vector `(a, 0)` (matches `runBands_2D.m`'s `pbcX`, selections
    `P.xEnd1`/`P.xEnd2`).
  - Pair 2: the two parallel side-wall faces corresponding to translation by
    lattice vector `(a/2, a*sqrt(3)/2)` (matches `pbcY`,
    `P.yEnd1`/`P.yEnd2`).
  - **Caveat**: because the unit cell is a rhombus (not axis-aligned), these
    are NOT the coordinate planes `x=const`/`y=const` -- they're the
    rhombus's slanted side walls. `runBands_2D.m`'s in-code comments
    (`"periodic BCs for yz planes at x = +/- a/2"`) are stale copy-paste from
    the 1D nanobeam pipeline and describe the WRONG geometry for this cell;
    trust the actual `bndindex(...)` face selections in
    `buildBoomerangUnitCell.m`/`runBands_2D.m`, not those comments. In the
    GUI, identify the two pairs of parallel slanted faces visually (there are
    exactly 4 side walls beyond the top/bottom slab faces).
- **Face selections must be geometric, not index-based.** The exported model
  pins them to absolute indices (`pbcX [1 22]`, `pbcY [2 9]`, parity `[3]`),
  which COMSOL renumbers whenever the geometry changes -- and *every* design
  variable changes it (growing `r` past the wall splits a side face; changing
  `r1`/`r2` adds or removes fillet faces). The MATLAB pipeline instead
  resolves faces by coordinate and normal via `bndindex(...)`
  (`buildBoomerangUnitCell.m:106-112`); port that, or use named/box selections
  built inside the geometry sequence so COMSOL re-evaluates them per rebuild.
- **z=0 midplane selection**: a Box/plane selection picking out the z=0
  cross-section of the slab. **Do not copy `ZsymSel` from
  `buildBoomerangUnitCell.m:93-101`** -- it spans `x,y` in `[-a/2, +a/2]` with
  `condition='allvertices'`, but this cell spans `x` in `[0, 3a/2]` and `y` in
  `[0, a*sqrt(3)/2]`, so it selects nothing. It is dead code in MATLAB too
  (`runBands_2D.m` uses `P.zEnd` from `bndindex`, never `P.bndSel.Zsym`), which
  is why nobody noticed. The correct box is `x` in `[-d, 3a/2+d]`, `y` in
  `[-d, a*sqrt(3)/2+d]`, `z` in `[-d, +d]`.
  Apply, on that selection:
  - **Study "mech evenz"**: `Symmetry` boundary condition on that plane
    (z-even / symmetric modes).
  - **Study "mech oddz"**: `Antisymmetry` boundary condition on that same
    plane (z-odd / antisymmetric modes) -- a SEPARATE study, not a toggle on
    the same one, so `model.solve("mech evenz")` /
    `model.solve("mech oddz")` can each be called directly from Python
    without touching BC activation state between calls.
- Both studies: **Eigenfrequency**, searching for `n_bands` (e.g. 10)
  eigenfrequencies around a shift frequency (pick something near your
  expected mechanical target, e.g. 8[GHz] -- see `configs/targets_2d.yaml`).
- All OTHER outer faces (the slab's top/bottom z=+-th/2 faces): free
  (traction-free) by default -- do not add periodicity or symmetry there.

## 5. Optical study -- NOT YET BUILT

Skip this for now. See `docs/optical_2d_plan.md` and
`src/optical_comsol_2d.py`'s module docstring for the planned approach (a 2D
slab photonic band structure is a materially bigger lift than the 1D
project's optical study, mainly because of light-line filtering at every
in-plane `(kx,ky)` -- worth reading before starting).

## 6. Save as `comsol/trusty_boomerang.mph`

`src/acoustic_comsol_2d.py` requires parameters `a,w,r,r1,r2,th` (checked by
`_require_parameters`, plus `k,kx,ky`) and study **labels** `"mech evenz"` /
`"mech oddz"` (checked by `_require_studies`). MPh resolves studies by label,
not by COMSOL tag. Keep these names, or override at run time:

```sh
OMC2D_STUDY_EVENZ="Study_symmetric" OMC2D_STUDY_ODDZ="Study 1" \
OMC2D_TEMPLATE=/path/to/other.mph python scripts/run_one_2d.py
```

Relabelling in the GUI is safer than hard-coding `"Study 1"`, which COMSOL
re-derives whenever another study is added.

**Save a solution-free, mesh-free copy for Python.** The current file is
1.1 GB because it stores 27 sweep points x 10 modes x 2 studies of solution
vectors plus a fine mesh. In *File > Save As*, uncheck "Save solutions in
model file" and "Save mesh in model file"; that should land in the low
hundreds of KB, small enough to commit and fast enough to load per process.

Also consider reverting `geomRep('cadps')` (Parasolid) to the default COMSOL
kernel: Parasolid needs the CAD Import Module / Design Module / a LiveLink-CAD
product, and a headless `comsol mphserver` on a Structural-Mechanics-only
licence will fail at geometry rebuild. It doesn't fail today only because the
geometry is never rebuilt. `buildBoomerangUnitCell.m` never sets `geomRep`.

## Brillouin-zone sweep convention

The driver walks `k` continuously over `[0, 3)` in three straight segments of
the hexagonal BZ -- Gamma -> M -> K -> Gamma -- computing `(kx, ky)` from `k`
via the exact piecewise formulas in `runBands_2D.m` (ported to Python in
`src/acoustic_comsol_2d.py:_kxky_hexagonal`). `bz_path(a, n_per_segment)`
returns `3*n_per_segment` points plus one final point closing the loop back
at Gamma, matching `runBands_2D.m`'s `ds.k_norm` convention (`0` to `3`,
`M` at `k=1`, `K` at `k=2`, back to `Gamma` at `k=3`).
