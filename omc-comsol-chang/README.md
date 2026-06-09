# omc-comsol-chang

MATLAB/COMSOL simulation codebase for optomechanical crystal (OMC) cavities. Two parallel simulation pipelines: **photonic crystal band structure** (root directory) and **nanobeam FEM cavities** (`nanobeam/`).

---

## Repository Structure

```
omc-comsol-chang/
├── [geometry builders]       buildBoomerang*.m, buildHole*.m, buildSnowflake*.m, ...
├── [band structure solvers]  solveBands.m, solveOpticalBands.m, runBands*.m
├── [optimization scripts]    boomerang_optimize_*.m, cross_optimization.m, ...
├── [parameter sweeps]        sweep_*.m
├── [utilities]               LoadMaterialParams.m, CreateFileBase.m, findGaps.m, ...
├── [test scripts]            test_*.m
├── bandStruct_data/          Saved band structure data (.mat, .csv, .txt)
├── test/                     Output directory for test runs
└── nanobeam/
    ├── [FEM pipeline]        RunNanobeamFEM.m, CreateNanobeamGeom*.m, Build*.m, ...
    ├── [solvers]             SetupNanobeamFEM.m, SolveNanobeamFEM.m, ...
    ├── [coupling calc]       CalcGOM.m, CalcStrCplSiV.m
    ├── [optimization]        OptimizeNanobeamRect*.m, OptimizeNanobeamTri*.m, ...
    ├── [visualization]       PlotDefectCells.m, PlotEy.m, PlotDispStr.m, ...
    ├── [test scripts]        test_*.m
    └── test/                 Output directory for nanobeam test runs
```

---

## Pipeline A: Photonic Crystal Band Structure (Root)

### Simulation Flow

```
test_*.m  →  runBands.m  →  buildXXXUnitCell.m  →  solveBands.m  →  findGaps.m
                         →  solveOpticalBands.m  →  findGaps_optical.m
```

### Geometry Builders

All follow the signature `[model, P] = buildXXX(model, P)`.

| File | Description |
|------|-------------|
| `buildBoomerangUnitCell.m` | 3D boomerang-shaped hole unit cell (hexagonal lattice) |
| `buildBoomerangUnitCell_2D.m` | 2D version of the above |
| `buildLowerBoomerangUnitCell.m` | Boomerang variant with lowered symmetry (3D) |
| `buildLowerBoomerangUnitCell_2D.m` | 2D lowered boomerang |
| `buildBoomerangStrip_3D.m` | Boomerang unit cell embedded in a strip waveguide (3D) |
| `buildBoomerangUnitCellStrip.m` | Strip variant of boomerang unit cell |
| `buildBoomerangUnitCellStrip_v2.m` | Revised strip variant |
| `buildHoleUnitCell.m` | Cylindrical hole unit cell with optional sidewall angle |
| `buildHoleStrip_2D.m` | Hole array in strip waveguide (2D) |
| `buildHoleStrip_3D.m` | Hole array in strip waveguide (3D) |
| `buildHoleStrip_withWg_2D.m` | Hole strip with attached waveguide for coupling |
| `buildSnowflakeUnitCell.m` | 6-arm snowflake hole geometry |
| `buildSnowflakeStrip_2D.m` | Snowflake array in strip (2D) |
| `buildSnowflakeStrip_3D.m` | Snowflake array in strip (3D) |
| `buildRibUnitCell.m` | Rib waveguide unit cell |
| `buildRibUnitCell_LN.m` | Rib unit cell for lithium niobate |
| `buildCrossUnitCell.m` | Cross-shaped hole unit cell |

### Band Structure Solvers

| File | Description |
|------|-------------|
| `solveBands.m` | Core mechanical band structure solver (Floquet–Bloch BCs). Returns `ds` struct with frequencies, k-points, and gap info. |
| `solveBands_noSym.m` | Same as above but without symmetry exploiting (slower, used for validation) |
| `solveOpticalBands.m` | Photonic band structure solver using RF/EM physics |
| `runBands.m` | Wrapper: builds geometry then calls `solveBands` |
| `runBands_2D.m` | 2D version of `runBands` |
| `runBands_noSym.m` | Wrapper calling `solveBands_noSym` |

### Optimization & Sweep Scripts

| File | Description |
|------|-------------|
| `boomerang_optimize_sweep_diamond.m` | Multi-parameter optimization for boomerang geometry on diamond |
| `cross_optimization.m` | Cross geometry band gap optimization |
| `snowFlake_optimization.m` | Snowflake geometry band gap optimization |
| `rib_optimize_debug.m` | Debug/exploratory rib optimization |
| `rib_optimize_sweep.m` | Parameter sweep over rib geometry |
| `rib_optimize_sweep_diamond.m` | Rib sweep for diamond platform |
| `sweep_boomerang_code.m` | Sweep lattice constant, width, and hole dimensions (boomerang) |
| `sweep_boomerangLower_code.m` | Sweep for lowered-boomerang geometry |
| `sweep_boomerangLower_code_v2.m` | Revised sweep script |
| `sweep_boomerangLowerOptical_code.m` | Optical band structure sweep for lower boomerang |
| `sweep_Hole_optical_2d.m` | 2D optical sweep for hole geometry |
| `sweep_snowFlake_code.m` | Snowflake parameter sweep |
| `sweep_Snowflake_strip.m` | Snowflake strip sweep |

### Utility Functions

| File | Description |
|------|-------------|
| `LoadMaterialParams.m` | Sets elastic constants (E, ν, ρ, stiffness tensor C) for `diamond`, `silicon`, or `LN` onto `P` struct |
| `RotateXtalTensor.m` | Rotates a 6×6 elasticity tensor (Voigt notation) by a given crystal orientation matrix |
| `findGaps.m` | Finds mechanical band gaps from a dispersion relation |
| `findGaps_optical.m` | Finds photonic band gaps |
| `analyzeGapData.m` | Loads saved `.mat` band structure files and plots results |
| `bndindex.m` | Returns COMSOL boundary indices that lie on a given plane (used for symmetry BCs) |
| `edgeindex.m` | Returns COMSOL edge indices that lie on a given line |
| `CreateFileBase.m` | Generates a file basename string from the parameter struct `P` (for consistent naming of output files) |
| `geom.m` | Archived hardcoded COMSOL geometry — not part of active pipeline |

### Test Scripts (root)

All follow the pattern: set up `P` struct → call solver → optionally save.

| File | Parameters Varied |
|------|-------------------|
| `test_Boomerang.m` | Boomerang unit cell defaults |
| `test_BoomerangLower.m` | Lowered-boomerang variant |
| `test_BoomerangStrip.m` | Boomerang in strip waveguide |
| `test_Hole.m` | Cylindrical hole geometry |
| `test_Snowflake.m` | Snowflake geometry |
| `test_Rib.m` | Rib waveguide geometry |
| `test_Cross.m` | Cross-hole geometry |
| *(others)* | Optical/2D variants |

---

## Pipeline B: Nanobeam FEM Cavities (`nanobeam/`)

### Simulation Flow

```
test_nanobeam*.m
    └── RunNanobeamFEM.m
            ├── CreateNanobeamGeom.m   (cavity geometry)
            ├── BuildNanobeamFEM.m     (COMSOL model setup)
            ├── SetupNanobeamFEM.m     (study/solver config)
            └── SolveNanobeamFEM.m     (solve + post-process)
                    ├── CalcGOM.m          (optomechanical coupling g_om)
                    └── CalcStrCplSiV.m    (strain coupling to SiV centers)
```

### Geometry Creation

| File | Description |
|------|-------------|
| `CreateNanobeamGeom.m` | Main nanobeam cavity geometry: rectangular or triangular cross-section with defect taper |
| `CreateNanobeamGeom_asym.m` | Asymmetric cavity variant (for studying symmetry breaking) |
| `CreateNanobeamGeomSnowflake.m` | Nanobeam with snowflake mirror unit cells |
| `CreateNanobeamGeomBoomerang1D.m` | Nanobeam with 1D boomerang mirror lattice |
| `CreateNanobeamBlockTetGeom.m` | Block-with-tether geometry for FEM testing |

### COMSOL Model Building

| File | Description |
|------|-------------|
| `BuildNanobeamFEM.m` | Sets up COMSOL model with solid mechanics + RF physics, applies mesh settings, materials |
| `BuildNanobeamFEMGen.m` | Generalized geometry builder for nanobeam FEM. Handles multiple hole types (`hole`/ellipse, `tri`/triangular polygon, `anvil`/filleted triangle, `rib`/parabolic profile), rectangular and triangular beam cross-sections, optional air cylinder for optical simulations, optional mechanical PML on beam ends (for Q-factor studies), z-symmetry plane clipping, and both symmetric and asymmetric cavities. Populates `P.domSel` (domain selections) and `P.bndSel` (boundary selections) on return, used by downstream physics and solver setup. |
| `BuildNanobeamSnowflakeFEM.m` | Same for snowflake nanobeam |
| `BuildNanobeamBoomerang1DFEM.m` | Same for 1D boomerang nanobeam |
| `SetupNanobeamFEM.m` | Configures COMSOL study (eigenfrequency solvers for mechanical and optical) |
| `SetupNanobeamSnowflakeFEM.m` | Snowflake variant |

### Solvers

| File | Description |
|------|-------------|
| `SolveNanobeamFEM.m` | Runs eigenvalue solve, identifies mode indices, exports field data |
| `SolveNanobeamSnowflakeFEM.m` | Snowflake variant (near-identical to above) |
| `SolveNanobeamBoomerangFEM.m` | Boomerang variant (near-identical to above) |
| `RunNanobeamFEM.m` | Top-level orchestrator: calls Create → Build → Setup → Solve, handles saving |
| `RunNanobeamBands.m` | Band structure run for nanobeam mirror unit cell |

### Coupling Calculations

| File | Description |
|------|-------------|
| `CalcGOM.m` | Computes optomechanical coupling rate g_OM from overlap integral of optical and mechanical modes |
| `CalcStrCplSiV.m` | Computes strain coupling between mechanical mode and SiV color center spin; handles crystal orientation and Jahn–Teller distortion |

### Optimization Scripts

| File | Description |
|------|-------------|
| `OptimizeNanobeamRect.m` | Main optimization for rectangular cross-section nanobeam |
| `OptimizeNanobeamRect_20191016.m` | Historical snapshot (Oct 2019) |
| `OptimizeNanobeamRect_20200103.m` | Historical snapshot (Jan 2020) |
| `OptimizeNanobeamRect_20260222.m` | Historical snapshot (Feb 2026) |
| `OptimizeNanobeamTri.m` | Triangular cross-section optimization |
| `OptimizeNanobeamTri35_*.m` | Tri optimization at 35° sidewall angle, various dates |
| `OptimizeNanobeamTri45_*.m` | Tri optimization at 45° sidewall angle, various dates |
| `OptimizeNanobeamSnowflake.m` | Snowflake nanobeam optimization |
| `sweep_NanobeamFEM.m` | Parameter sweep over nanobeam geometry |
| `sweep_1DSnowflakeFEM.m` | Sweep for 1D snowflake nanobeam |

### Visualization

| File | Description |
|------|-------------|
| `PlotDefectCells.m` | Plots the defect taper cell geometry along the beam |
| `PlotEy.m` | Plots the Ey optical field component from FEM results |
| `PlotDispStr.m` | Plots mechanical displacement and strain distributions |
| `PlotStrCplSiV.m` | Plots strain coupling vs. SiV position along beam axis |
| `PlotStrCplSiVXY.m` | Strain coupling map in XY plane |
| `PlotStrCplSiVYZ.m` | Strain coupling map in YZ plane |
| `PlotStrCplSiV_v20200527.m` | Historical visualization snapshot |

### Utilities (nanobeam/)

| File | Description |
|------|-------------|
| `LoadMaterialParams.m` | Duplicate of root-level version (diamond/silicon/LN parameters) |
| `CreateFileBase.m` | Duplicate of root-level version (output filename generation) |

---

## Key Data Structures

### Parameter Struct `P`

Passed through the entire pipeline. Key fields:

| Field | Description |
|-------|-------------|
| `P.a` | Lattice constant (m) |
| `P.w` | Beam width (m) |
| `P.th` | Beam thickness (m) |
| `P.r`, `P.hx`, `P.hy` | Hole dimensions (m) |
| `P.celltype` | `'boomerang'`, `'hole'`, `'snowflake'`, `'rib'`, `'cross'` |
| `P.xsect` | Cross-section: `'rect'`, `'tri'`, `'isoFit'` |
| `P.beamMat` | Material: `'diamond'`, `'silicon'`, `'LN'` |
| `P.anisoMat` | `true` to use anisotropic stiffness tensor |
| `P.rxtal` | Crystal rotation matrix for anisotropic materials |
| `P.kpts` | Number of k-points along Γ–X |
| `P.nbands` | Number of eigenvalue bands to compute |
| `P.lambda` | Target optical wavelength (m) |
| `P.mbevenz`, `P.mbeveny` | Mechanical symmetry flags (even/odd) |
| `P.run_optical` | `true` to also solve optical band structure |
| `P.savedat` | `true` to save results to disk |
| `P.datLoc` | Output directory path |

### Results Struct `ds` / `bds`

Returned by solvers:

| Field | Description |
|-------|-------------|
| `ds.sym`, `ds.asym`, `ds.full` | Band structures by symmetry class |
| `ds.F` | Eigenfrequencies |
| `ds.k_norm` | Normalized k-points (0 to 0.5) |
| `ds.gapSize`, `ds.midGap` | Band gap size and center frequency |
| `ds.P` | Copy of input parameters |
| `ds.mfem`, `ds.ofem` | Mechanical / optical FEM result substructs (nanobeam) |
| `ds.cpl` | Coupling results (g_OM, strain coupling) |

---

## Potential Refactoring Opportunities

### 1. Collapse 20 Geometry Builders into a Dispatcher

Every `build*.m` file shares ~40–60% code: WorkPlane setup, symmetry plane creation, fillet application, and extrusion. The unique part is only the hole/polygon definition.

**Suggestion:** Create a single `buildUnitCell(model, P)` that dispatches on `P.celltype` and calls small private functions for the shape-specific part:

```matlab
function [model, P] = buildUnitCell(model, P)
    % shared setup ...
    switch P.celltype
        case 'boomerang', pts = boomerangPolygon(P);
        case 'hole',      pts = holePolygon(P);
        case 'snowflake', pts = snowflakePolygon(P);
    end
    % shared extrude, fillet, symmetry logic ...
end
```

This alone would eliminate ~1,000 lines of duplicated code.

### 2. Merge the Three Near-Identical Nanobeam Solvers

`SolveNanobeamFEM.m`, `SolveNanobeamSnowflakeFEM.m`, and `SolveNanobeamBoomerangFEM.m` are ~95% identical. The only differences are minor physics tag names.

**Suggestion:** Keep one `SolveNanobeamFEM.m` and pass the geometry variant as a field in `P` (e.g., `P.mirrorType`). Remove the two variant files.

### 3. Extract Shared Nanobeam Geometry Logic

`CreateNanobeamGeom.m`, `CreateNanobeamGeomSnowflake.m`, and `CreateNanobeamGeomBoomerang1D.m` share taper generation and parameter array setup. 

**Suggestion:** Extract a `CreateTaper(P)` helper and a `CreateHoleArray(P, taperParams)` helper. The three files then become thin wrappers.

### 4. Replace Replicated File I/O with a Logging Helper

Every optimization and sweep script repeats this exact pattern:
```matlab
datLoc = ['.\test\...', datestr(now,'mmddyyyy'), '\'];
if ~exist(datLoc, 'dir'); mkdir(datLoc); end
itr = fopen([datLoc, 'results.csv'], 'at+');
fprintf(itr, '%.4e\t%.4e\t...\r\n', P.a, P.th, ...);
fclose(itr);
```

**Suggestion:** Create `LogResult(datLoc, P, result)` that handles directory creation, file opening, formatting, and closing. Reduces ~10 lines per call site to 1.

### 5. Centralize Material Constants

`LoadMaterialParams.m` is duplicated verbatim in the root and `nanobeam/` directories. Material values are hardcoded as literals.

**Suggestion:** Keep one copy in root (or a `lib/` folder), have `nanobeam/LoadMaterialParams.m` call it, or better, use a `materials.json` or `materials.mat` config file so the constants can be updated in one place without touching code.

### 6. Add a Parameter Template System

All 30+ test files manually construct `P` from scratch with many repeated defaults. When a default changes, every test file must be updated.

**Suggestion:** Create `DefaultParams(celltype)` that returns a `P` struct with sensible defaults. Test files then only override what differs:
```matlab
P = DefaultParams('boomerang');
P.a = 450e-9;  % only the thing being tested
```

### 7. Delete Version-Dated Files — Use Git

Files like `OptimizeNanobeamRect_20191016.m`, `OptimizeNanobeamRect_20200103.m`, `PlotStrCplSiV_v20200527.m` are point-in-time snapshots that now live in the working tree. They create confusion about which version is current and inflate the codebase by ~15 files.

**Suggestion:** Verify these are committed in git history, then delete them. `git log -- OptimizeNanobeamRect.m` shows all historical states.

### 8. Name Magic Numbers as Constants

Scattered literals with no explanation:
- `50e-9` — fillet radius (appears in 8+ files)
- `10e9` — target mechanical frequency (GHz)
- `1640e-9` — target optical wavelength (nm)
- Mesh quality integers `1`–`5`

**Suggestion:** A single `constants.m` (or top-of-file `const` struct) documents intent and makes global changes easy.

### 9. Standardize Naming Conventions

Mixed conventions make the codebase harder to navigate:
- Struct fields: `P.beamMat` (camelCase) vs. `P.dat_loc` (snake_case) — pick one
- Result structs: `ds`, `bds`, `mfem`, `ofem` — no consistent prefix
- Function names: `buildBoomerangUnitCell` vs. `RunNanobeamFEM` vs. `sweep_boomerang_code` — three conventions

### 10. Break Up Long Functions

- `OptimizeNanobeamRect.m` (~700 lines): separate geometry setup, COMSOL run, result analysis, and file I/O into named subfunctions or separate files
- `CalcStrCplSiV.m` (~400 lines): the inner loop over SiV positions and orientations can be extracted
- `RunNanobeamFEM.m` (~330 lines): the try-catch branching for different geometry types could be its own dispatch function

---

## Dependencies

- **MATLAB** R2019b or later (uses `string`, `logical`, table operations)
- **COMSOL Multiphysics** with LiveLink for MATLAB — required for all `model.*` API calls
  - Solid Mechanics module (structural eigenfrequency studies)
  - RF / Wave Optics module (optical eigenfrequency studies)
- No external MATLAB toolboxes required beyond base MATLAB

---

## Quickstart

```matlab
% Band structure — boomerang geometry on diamond
P.xsect    = 'rect';
P.beamMat  = 'diamond';
P.celltype = 'boomerang';
P.a        = 400e-9;
P.w        = 86e-9;
P.th       = 200e-9;
P.savedat  = false;
bds = solveBands(P);

% Nanobeam FEM cavity
datLoc = '.\test\nanobeam\';
P = struct('xsect','rect','beamMat','diamond','a',400e-9,'w',500e-9,'th',200e-9);
ds = RunNanobeamFEM(P, datLoc);
```
