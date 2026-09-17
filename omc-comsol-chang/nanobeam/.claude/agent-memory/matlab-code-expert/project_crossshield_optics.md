---
name: project-crossshield-optics
description: 2026-09-16 - cross-shield builder now has two y build paths (half-y mirror, full-y flush-bottom one-sided shield) plus P.solveOpt optics; open risks are the shield mesh cost, the unverified bandgap, and stale in-file arithmetic
metadata:
  type: project
---

`BuildNanobeamCrossShieldFEM.m` carries TWO orthogonal gates, both added on
**2026-09-16**:

- `useOpt = isfield(P,'solveOpt') && P.solveOpt` - builds the cavity-scoped air
  cylinder so optical and mechanical modes come out of ONE model and `CalcGOM`
  can run. `test_nanobeamCrossShieldFEM.m` has a MASTER SWITCH block at the top
  and derives `P.calcG`, `P.plotOpt` and an optics tag in `datLoc` from it.
- `fullY = abs(P.meveny) < 1` - the asymmetric build. No y = 0 mirror plane, the
  full beam width is simulated, and the shield bottom is FLUSH with the beam
  bottom at `y = -P.w/2` with `P.nShieldY` cells running upward from there. The
  shield is one sided; `P.nShieldY` is the absolute row count, not a half-model
  count. `test_nanobeamCrossShieldFEM.m` now runs in this mode.

**Why:** `CalcGOM` requires `ds.ofem` and `ds.mfem` on the same model and Joins
their datasets, so g_OM cannot come from a separate optical script. The full-y
path is the user's explicit device change, not a modelling convenience.

**How to apply:**

- **`P.oeveny` and `P.meveny` must agree, and so must `P.oevenz`/`P.mevenz`.**
  There is one geometry and two physics. The builder hard-errors
  (`evenyMismatch`, `evenzMismatch`) inside the `useOpt` branch. In full-y the
  builder emits NO `beamYsym`/`PMLYsym`/`cylYsym` selections at all; that is
  safe because `SetupNanobeamFEM` reads them only from inside its
  `meveny/oeveny == +/-1` branches and `CalcGOM.m:92` is guarded the same way.
  No edit to those files was needed - confirm that again before assuming it.
- **Filename encodes the y convention, not just the count.**
  `CreateFileBase.m` (BOTH copies - root and `nanobeam/`, kept in sync per
  CLAUDE.md) writes `nsyH<n>_nsyP<2n>` for half-y and `nsyF<n>` for full-y.
  Without this a full-y and a half-y run at the same `P.nShieldY` collide on
  `P.fileBase` and `RunNanobeamFEM.m:64` hands back the wrong cached solve.
- **The air cylinder is scoped to the cavity**, not the structure:
  `P.airCylLen` defaults to `min(beamLenHalf, xS0 - P.lambda)`. x is the ONLY
  separation from the shield - y and z overlap in both builds - so anything that
  shortens the beam or lengthens the clamp pad moves that clearance. In full-y
  the Revolve sweeps 360 deg (`angle1 -180`, `angle2 +180`) instead of 180, so
  the air is a HALF cylinder (all y, z >= 0) and the optical mesh doubles.
- **The binding constraint is the shield mesh, and it is not close.** At the
  script's current cross cell (`a = 894`, `h = 863`, `w = 416 nm`) the default
  `shieldHmax = min((a-h)/3, th/3) = 10.3 nm`, giving ~5.2e6 shield tets and
  ~2.2e7 DOF. `test_nanobeamCrossShieldFEM.m` raised `P.max_dof` to `2.5e7` on
  2026-09-16 to accommodate that (other scripts still sit at `5e6`), so the
  budget is now marginal rather than 4x over - but it is untested at full size
  and a MUMPS factorisation there will exhaust a workstation. Shrink to
  `nShieldX = nShieldY = 2` and `solveMechPML = 0` for any mesh-only check.
  Switching to full-y adds under 1% (only the beam doubles;
  the shield and +x PML arm are the same `nShieldY*aShield = 4470 nm` tall
  either way). `mAdjMesh` cannot rescue it because `applyMeshOverrides`
  re-applies the shield `hmax` on every coarsening pass. Older comments quoting
  `33.3 nm` / `9.7e5 tets` were written for `a-h = 100 nm` and are wrong.
- **`SolveNanobeamFEM` makes ONE mesh node shared by both studies.** Optical
  meshes globally at `lambda/5` with no shield refinement, then mechanical drops
  that and re-meshes with the shield `hmax`. The clean fix is separate mesh
  nodes, which needs a `SolveNanobeamFEM` edit - out of scope, flag it (see
  [[feedback-house-style]] rule 3).
- **Still unvalidated by anyone:** whether `P.freq = 7 GHz` is inside the cross
  lattice bandgap. `v_shear/(2f) = 918 nm` vs `aShield = 894 nm` is encouraging
  but is not a bandgap calculation. Also `a - h = 31 nm` is below
  `P.shieldMinFeature = 50e-9`, so `isCrossFabricable` fails; widening the wall
  would fix both the fab warning and the mesh cost but changes the lattice, so
  it needs the independent `buildCrossUnitCell` + `runBands` check first.
- **Do not trust worked arithmetic in the test script's comments without
  re-deriving it.** Verified headless on R2025b at the 2026-09-16 parameters:
  `beamLenHalf = 11063.0 nm` (= `xpos(18) + a/2` = 10738.0 + 325.0),
  `airrad = 3500.0`, default `airCylLen = 9513.0`, hole 16 spans
  `[9266.5, 9609.5]`, `max(hy)/2 = 308.5` vs `w/2 = 400.0`. A figure of
  `beamLenHalf = 10999.4` was raised independently and could NOT be reproduced
  at these parameters - it corresponds to `maxdef ~ 0.193` (or `taperFunc`
  `'quadratic'` gives 10979.8), neither of which the script sets. Always rerun
  `LoadMaterialParams` + `CreateNanobeamGeom` under `matlab -batch` (stage S0b)
  before quoting a number.
- **The full-y clamp is worse and the audit says so.** 52.00% of the 800 nm beam
  end face lands on a 15.5 nm wall fin, vs 40.25% in half-y, because the flush
  bottom puts the beam axis near the middle of cell 1 rather than on a node
  centre. The fin band width is `wShield`, not a node-width difference, so the
  half-y closed form is the wrong shape and the full-y branch measures band
  overlap instead. `P.shieldPadLen ~ aShield/2 = 447 nm` is the documented fix
  and is deliberately left at 0 - that is the user's call, not the builder's.
