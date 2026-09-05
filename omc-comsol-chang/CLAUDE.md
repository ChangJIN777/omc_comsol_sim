# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running Simulations

There is no build system. Simulations are run by executing `test_*.m` scripts directly in MATLAB. COMSOL with LiveLink for MATLAB must be running before executing any script that calls `model.*` API methods.

**Band structure (root directory):**
```matlab
% In MATLAB, run a test script directly:
run('test_Boomerang.m')
```

**Nanobeam FEM (nanobeam/ directory):**
```matlab
run('nanobeam/test_nanobeamRectFEM.m')
```

To debug geometry interactively without a full solve, uncomment the block at the bottom of any test script:
```matlab
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.clear(); clear model
model = ModelUtil.create('model');
buildBoomerangUnitCellStrip(model,P);
mphlaunch(model);   % opens COMSOL GUI with the built geometry
```

Static analysis of a `.m` file (no COMSOL needed):
```matlab
checkcode('solveBands.m')
```

## Architecture

### Two Independent Pipelines

**Pipeline A — Photonic crystal band structure (root/)**
```
test_*.m  →  solveBands.m / solveOpticalBands.m
                ↑
            runBands.m calls buildXXXUnitCell.m first, then solveBands
```
Entry: a `test_*.m` sets `P` and calls `solveBands(P)` directly, or calls `runBands(P)` which builds the geometry internally. `solveBands` creates the COMSOL model, calls the appropriate `buildXXXUnitCell` function, applies Floquet–Bloch BCs, runs the eigenfrequency study, and calls `findGaps` on results.

**Pipeline B — Nanobeam FEM cavity (nanobeam/)**
```
test_*.m  →  RunNanobeamFEM(P, datLoc)
                ├── CreateNanobeamGeom*(P)      generates hole-dimension arrays
                ├── BuildNanobeamFEM(model,P)   COMSOL model + mesh
                ├── SetupNanobeamFEM(model,P)   study/solver configuration
                └── SolveNanobeamFEM(model,P)   solve + post-process
                        ├── CalcGOM             optomechanical coupling g_OM
                        └── CalcStrCplSiV       strain coupling to SiV
```
`RunNanobeamFEM` dispatches to snowflake or boomerang variants based on `P.celltype`. It skips the full solve if `.mat` and `.mph` output files already exist in `datLoc`.

### The `P` Struct

All geometry, material, physics, and I/O settings live in a single struct `P` passed through every function. `LoadMaterialParams(P)` must be called (or its logic present in the test script) before solving — it sets `P.E`, `P.rho`, `P.nu`, `P.D` (anisotropic stiffness in COMSOL Voigt ordering) from `P.beamMat`. The anisotropic stiffness tensor is only applied when `P.anisoMat = 1`; `P.rxtal` (degrees) then rotates it via `RotateXtalTensor`.

Symmetry is controlled by:
- `P.mbevenz` / `P.mbeveny` — mechanical even/odd symmetry (band structure pipeline)
- `P.mevenx` / `P.meveny` / `P.mevenz` — mechanical symmetry (nanobeam FEM, ±1 for even/odd, 0 for fixed BC)
- `P.oevenx` / `P.oeveny` / `P.oevenz` — optical symmetry (nanobeam FEM)

### Output File Naming

`CreateFileBase(P)` (called automatically if `P.fileBase` is unset) builds a filename string encoding all key geometry parameters, e.g. `boomerang_a_400nm_r_160nm_w_86nm_th_180nm_...`. Results are saved as `[datLoc, P.fileBase, '_mech.mat']` etc. All output goes under `datLoc`, which test scripts build as `[fullfile('.','test','<geometry>',currentDate),filesep]`.

**Never paste a path separator in by hand.** A literal `'.\test\...'` is only a path on Windows; on macOS and Linux the backslash is an ordinary filename character, so the script creates one file whose *name* contains backslashes rather than a directory tree. Git will happily track such a file, and it then makes the repository impossible to check out on Windows (`error: invalid path`). Always use `fullfile(...)` plus a trailing `filesep` where a trailing separator is needed.

### COMSOL Model Object

The `model` object (a Java object from the COMSOL LiveLink API) is created once with `ModelUtil.create('model')` and passed by reference through all Build/Setup/Solve functions. Geometry, physics, mesh, and study nodes are all added to this single object. `ModelUtil.clear()` must be called before creating a new model to avoid memory leaks.

## Key Conventions

- All physical dimensions in the `P` struct are in **SI units (meters)**. Display/filename strings convert to nm (`*1e9`).
- `P.meshSize` (band structure) and `P.mMesh` / `P.oMesh` (nanobeam) are integer quality levels 1–5 passed to COMSOL's mesh calibration; lower = coarser/faster.
- `P.max_dof` caps degrees of freedom. When `P.mAdjMesh = 1`, the solver automatically coarsens the mesh if DOF count exceeds this limit.
- Geometry builders always take `(model, P)` and return `[model, P]` — `P` may be augmented with computed geometry fields.
- `nanobeam/LoadMaterialParams.m` and `nanobeam/CreateFileBase.m` are duplicates of root-level files; keep them in sync if changing material constants or filename logic.
