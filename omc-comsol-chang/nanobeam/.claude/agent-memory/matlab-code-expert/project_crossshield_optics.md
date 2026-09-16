---
name: project-crossshield-optics
description: 2026-09-16 - cross-shield builder gained cavity-scoped optics for g_OM; the binding constraint is the shared mesh node and the mechanical DOF ceiling, not the geometry
metadata:
  type: project
---

On **2026-09-16** `BuildNanobeamCrossShieldFEM.m` was extended to build the
optical air cylinder so that optical and mechanical modes come out of ONE COMSOL
model and `CalcGOM` can run. `test_nanobeamCrossShieldFEM.m` now has
`P.solveOpt = 1`.

**Why:** `CalcGOM` requires `ds.ofem` and `ds.mfem` on the same model and Joins
their datasets, so g_OM cannot be obtained from a separate optical script. The
user was told this needed the geometry builder extended and chose that path over
the alternatives.

**How to apply:**

- The air cylinder is **scoped to the cavity**, not the structure:
  `P.airCylLen` defaults to `min(beamLenHalf, xS0 - P.lambda)`. A full-length
  cylinder genuinely does swallow the shield. x is the only separation between
  the cylinder and the shield - y and z overlap - so anything that shortens the
  beam or lengthens the clamp pad moves that clearance.
- The **real open question is not geometry, it is the mesh**. `SolveNanobeamFEM`
  makes one mesh node shared by both studies: the optical branch meshes globally
  at `lambda/5` with no shield refinement, then the mechanical branch drops that
  and re-meshes with `hmax = P.shieldHmax` on the shield. The two meshes differ
  by ~800x in shield element count. This is pre-existing behaviour that
  `test_nanobeamRectFEM.m` also exercises, but the divergence is far larger here.
  The clean fix is separate optical and mechanical mesh nodes, which needs a
  `SolveNanobeamFEM` edit - out of scope, flag it rather than doing it (see
  [[feedback-house-style]] rule 3).
- The mechanical DOF estimate at 6x5 cells and `shieldHmax = 33.3 nm` is ~4-6e6
  against `P.max_dof = 5e6`, and that was already true with optics off. Optics
  adds no mechanical DOF (the two physics are selected on different domains).
  `mAdjMesh` cannot rescue it because `applyMeshOverrides` re-applies the shield
  `hmax` on every coarsening pass.
- Still unvalidated by anyone: whether `P.freq = 7 GHz` is inside the cross
  lattice bandgap for `aShield = 1600 nm`. The shield cell dimensions in the test
  script are explicitly labelled placeholders pending an independent
  `buildCrossUnitCell` + `runBands` check.
