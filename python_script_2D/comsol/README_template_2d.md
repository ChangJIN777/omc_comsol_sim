# Building the 2D boomerang unit-cell template (`omc2d_boomerang.mph`)

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
| a    | 480[nm]   | hexagonal lattice constant |
| w    | 140[nm]   | width of each boomerang leg rectangle |
| r    | 177[nm]   | length of each boomerang leg rectangle |
| r1   | 10[nm]    | fillet radius at leg tips |
| r2   | 10[nm]    | fillet radius at leg junctions (near hole center) |
| th   | 220[nm]   | slab thickness (out-of-plane, z) |
| kx   | 0[1/m]    | in-plane Floquet wavevector, x (swept by the Python driver) |
| ky   | 0[1/m]    | in-plane Floquet wavevector, y (swept by the Python driver) |

Values above match `test_Boomerang.m`'s defaults, for a sanity-check
cross-comparison against the MATLAB pipeline once both are built.

Note the departure from the MATLAB pipeline here: `runBands_2D.m` sweeps a
single scalar `k` and computes `kx`,`ky` from it via COMSOL expressions
(`if(k<1, ..., if(k<2, ..., ...))`) so that COMSOL's own GUI parametric-batch
solver can do the sweep natively. The Python/MPh driver instead computes
`kx`,`ky` directly in Python (see `src/acoustic_comsol_2d.py:bz_path`, an
exact port of those same piecewise expressions) and sets them per solve --
simpler here since MPh already drives point-by-point. **You do not need a
`k` parameter or COMSOL's Parametric Sweep/Batch node in this template.**

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
3. **Fillets**: round the leg tips with radius `r1`, and the junctions near
   the hole center with radius `r2` (see `addFillet()` in
   `buildBoomerangUnitCell.m` for the exact selection logic -- it uses a disk
   selection around the hole center to pick the right vertices).
4. **Extrude** the 2D profile through thickness `th` (centered at z=0, i.e.
   from `z=-th/2` to `z=+th/2`).

Build and visually check this against `mphgeom` output from the MATLAB
`buildBoomerangUnitCell.m` (run `test_Boomerang.m`'s commented-out debug
block) before proceeding -- the exact vertex/fillet-selection indices are
finicky enough that a visual cross-check is worth the five minutes.

## 3. Material

Diamond, anisotropic (cubic): `C11=1076[GPa]`, `C12=125[GPa]`, `C44=577[GPa]`,
density `3515[kg/m^3]` (see `configs/materials.yaml`). Set the material
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
- **z=0 midplane selection**: a Box/plane selection picking out the z=0
  cross-section of the slab (matches `ZsymSel` in `buildBoomerangUnitCell.m`).
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

## 6. Save as `comsol/omc2d_boomerang.mph`

`src/acoustic_comsol_2d.py` references parameter names `a,w,r,r1,r2,th,kx,ky`
and study names `"mech evenz"` / `"mech oddz"`. Keep these names, or update
`STUDY_EVENZ`/`STUDY_ODDZ` in that module to match whatever you name them.

## Brillouin-zone sweep convention

The driver walks `k` continuously over `[0, 3)` in three straight segments of
the hexagonal BZ -- Gamma -> M -> K -> Gamma -- computing `(kx, ky)` from `k`
via the exact piecewise formulas in `runBands_2D.m` (ported to Python in
`src/acoustic_comsol_2d.py:_kxky_hexagonal`). `bz_path(a, n_per_segment)`
returns `3*n_per_segment` points plus one final point closing the loop back
at Gamma, matching `runBands_2D.m`'s `ds.k_norm` convention (`0` to `3`,
`M` at `k=1`, `K` at `k=2`, back to `Gamma` at `k=3`).
