---
name: project-comsol-mesh-sequence
description: COMSOL LiveLink trap verified on 6.3 - adding any user feature to a mesh sequence flips it from physics-controlled to user-controlled, and it then meshes NOTHING unless an explicit meshing operation is present
metadata:
  type: project
---

A mesh sequence made by `model.mesh.create(tag, geom)` holds one global `size`
node and reports `isAutomatic = 1`: it is **physics-controlled** and COMSOL
generates the meshing operations itself at run time. Writing to `size`
(`set('hauto',N)`, `set('hmax',...)`) keeps it automatic. Creating **any other
feature - including a domain-scoped `Size` node -** flips it to user-controlled,
and a user-controlled sequence runs only what it holds. A `Size` node is not a
meshing operation, so the run yields zero elements everywhere and COMSOL emits
`No mesh on domains 1-N in the meshing sequence with tag mesh` for **every**
domain in the geometry.

Measured on COMSOL 6.3 via LiveLink, two-domain test block, `hauto = 5`:

| sequence | isAutomatic | elements |
|---|---|---|
| `size` only | 1 | 4877 |
| `size` + a top-level `Size` | 0 | **0** |
| + an explicit `FreeTet` | 0 | 3446, nested `hmax` honoured |

**Why it matters here:** this is what silently emptied the cross-shield mesh on
2026-09-16 (`P.PMLmesh`/`P.PMLmeshDiv`/`P.shieldHmax`/`P.beamHmax` each create a
`Size` node). The listed domain range being the *whole* geometry rather than a
subset is the tell that it is a sequence fault, not a selection fault.

**How to apply:**

- Never add a bare `Size` node to a sequence you have not also given a meshing
  operation. `SolveNanobeamFEM.m` now keeps the invariant with `ensureFreeTet` /
  `ensureMeshOperation`; nest sizing inside the operation the way
  `bands/optmechbands1D/RunOpticalBands.m:370-374` does. A top-level `Size`
  created *after* the operation would not apply to it anyway.
- A function that may add overrides must decide **all** its conditions before
  creating anything, and return untouched when none fire - otherwise it flips
  legacy callers off the physics-controlled mesh and changes their results (rule
  4 in [[feedback-house-style]]).
- Going user-controlled loses COMSOL's automatic **swept mesh in PML domains**,
  and on 2026-09-16 that was MEASURED and found not to matter for the
  cross-shield model. See [[project-pml-mesh-cost]] before spending any more
  effort there. `SolveNanobeamFEM.m` now has an opt-in swept path
  (`P.PMLmeshSwept`), default OFF, and the arithmetic for why it is off is in
  its `applyPMLSweep` doc block.
- `mesh.isAutomatic()` exists on 6.3 but prefer inspecting `mesh.feature.tags()`
  when you need version tolerance.
- Relevant to the "separate optical and mechanical mesh nodes" refactor flagged
  in [[project-crossshield-optics]]: any new mesh node inherits this trap.
