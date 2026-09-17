---
name: project-geom-runall-not-finalized
description: Latent repo-wide defect — builders ending in geom.runAll return empty P.xEnd*/P.yEnd*/P.zEnd/bndSel because runAll does not finalize the geometry
metadata:
  type: project
---

`geom.runAll` builds the sequence's geometry objects but does **not** finalize it
(it never runs the `fin` feature). Measured on COMSOL 6.3: immediately after
`runAll` the sequence reports `getNVertices = 0`, `getNBoundaries = 0`,
`getNDomains = 0`. So anything read at that point comes back empty —
`bndindex(...)` for all faces, cumulative selections via
`model.selection('geom1_..._bnd')`, and `P.bndSel.Zsym`. `geom.run` finalizes
and everything resolves.

**Why:** most builders in this repo end with `runAll` and then immediately do
their `bndindex` lookups (`buildCrossUnitCell.m`, `buildBoomerangUnitCell.m`,
`buildSnowflakeUnitCell.m`, `buildHoleStrip_3D.m`, ...). The bug stayed hidden
for years because essentially every `test_*.m` ships `P.fixed_bc = 0`, and that
is the only path that reads `P.xEnd*`/`P.yEnd*`; the *named* selections used for
the Floquet and symmetry BCs happen to work because `runBands.m` calls
`geom('geom1').run` itself before the physics. `buildHoleStrip_3D.m` has that
`run` call sitting commented out, which is someone hitting this before.
`buildCrossStrip.m` has been fixed (`run`, plus its `bndSel.Zsym` read moved
after it); the others have not.

**How to apply:** when writing or reviewing a builder, end with `geom.run`, not
`runAll`, if the builder returns entity indices or reads its own selections — and
read selections *after* that call, never right after a `runCurrent`. When a
`P.fixed_bc = 1` run or any `P.xEnd*`-consuming code fails with an empty index
list, suspect this before suspecting `bndindex`'s tolerance (its in-plane test is
relative, `max|d|*1e-10`, which is ~10^5 ulps at nm-scale coordinates — float
error is not the culprit). Warn the user if their run touches one of the
unfixed builders.

Related: [[feedback-preserve-legacy-geometry-behavior]],
[[feedback-scope-discipline]].
