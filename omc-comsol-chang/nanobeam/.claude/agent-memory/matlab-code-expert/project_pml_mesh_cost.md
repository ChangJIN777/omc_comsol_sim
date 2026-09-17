---
name: project-pml-mesh-cost
description: Measured 2026-09-16 - the mechanical PML is under 1% of the cross-shield element count, a swept PML mesh cannot help, and "Out of memory on server" is the COMSOL server's 2 GB Java heap, not MUMPS
metadata:
  type: project
---

Two things that look like the same problem and are not. Both were investigated
on **2026-09-16** after `com.comsol.util.exceptions.UnexpectedServerException:
Out of memory on server. java.lang.OutOfMemoryError`.

## 1. The PML is NOT a mesh-cost problem, and a sweep cannot make it one

Recomputed headless at `test_nanobeamCrossShieldFEM.m`'s parameters
(`aShield 894`, `hShield 863`, `wShield 416`, `nShieldX 6`, `nShieldY 5`,
`th 250` with z symmetry so the slab is 125 nm, `PMLLen = PMLLenY = 2000`,
`PMLmeshDiv = 8`). All volumes are the z >= 0 half and are **independent of
beam length**, so they do not move with `P.a/hx/hy/maxdef`:

| | volume | hmax | elements |
|---|---|---|---|
| PML frame (x arm 1.117, y arm 1.341, corner 0.500) | 2.958 um^3 | 250 nm | 1.14e3 bulk |
| + graded shell down to the shield interface | | 10.33 -> 250 nm, 11 layers, 772 nm deep | ~4.9e4 |
| shield, SOLID only | 0.954 um^3 | 10.33 nm | 5.19e6 |

**The PML is 3.10x the shield BY VOLUME and 0.95% of it by element count.**
Its cost is not the bulk, it is the graded shell where it meets the shield.
A swept mesh cannot remove that shell, because the shell *is* the sweep source
face: `size_shield` is created after `size_pml` and so wins on the shared
`x = xS1` / `y = yS1` faces, which is correct - those faces carry the 31 nm
cross walls. Swept costs: 4 layers -> 2.96e4 prisms (1.68x better), 8 layers ->
5.90e4 (worse than the free tet, whose grading already coarsens outward).
Either way it moves the model by well under 1%.

So `P.PMLmeshSwept` exists in `SolveNanobeamFEM.m` and is **off by default on
purpose**. The levers that actually matter are `P.shieldHmax`, `P.nShieldX`,
`P.nShieldY` and widening `a - h` - see [[project-crossshield-optics]].

## 2. The exception is the Java heap, and it is 2 GB

`C:\Program Files\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver.ini`
ends with `-Xmx2g` / `-XX:MaxMetaspaceSize=1g`. That is the COMSOL **server**
heap, and `java.lang.OutOfMemoryError` can only come from there - MUMPS and
PARDISO allocate native memory and fail with their own "not enough memory for
LU factorization" text instead.

The allocation that overflows it is `mphxmeshinfo`, called at six places in
`SolveNanobeamFEM.m` purely to read `.ndofs`. It has **no option to skip the
bulk data** (checked the R2025b/6.3 `mli/mphxmeshinfo.m` header): it always
materialises `dofs.coords` (3 x ndofs doubles), `dofs.nodes/nameinds/
solvectorinds`, `nodes.coords` and `elements.tet.nodes/dofs`. At 5.2e6 tets /
2.2e7 DOF that is 1.5-2.5 GB on the Java side alone, before the solve.
`P.storeMPH = 1` with `P.mneigs = 30` is a second heap trigger - 30 complex
eigenvectors at 2.2e7 DOF is ~10 GB to serialise.

**How to apply:** when a run dies with `Out of memory on server`, ask first
*where in the log it died*. Before any solver output -> heap / `mphxmeshinfo`;
during `Eigenvalue solver` -> genuinely too large for MUMPS. Raising `-Xmx`
fixes only the first. Related: [[project-comsol-mesh-sequence]],
[[user-domain]].
