# `rib` parameter reference

This page is the reference for every `P` field the `'rib'` celltype uses. `'rib'` merges two earlier devices:

- the visible-rib design of `comsol_templates/visible_rib_omc.m`: spine + centred rib, fillets, cavity taper;
- the lithium-niobate rib cell: sidewall angle `P.d`, explicit beam width `P.w`.

It is used in both pipelines:

- **Unit cell** (band structure): `solveBands` / `solveBands_noSym` → `runBands` / `runBands_noSym` → `buildRibUnitCell_LN`. Optical: `solveOpticalBands` → `runOpticalBand_1D` / `_3D` → `buildRibUnitCell_LN`.
- **Cavity** (nanobeam FEM): `RunNanobeamFEM` → `createRibCavityGeom` → `buildRibCavity` → `SetupNanobeamFEM` / `SolveNanobeamFEM`.

All defaults live in one place, `ribDefaultParams.m`. Every builder, file-naming branch and solver branch for this celltype calls it, and it fills in **only the fields you have not set**. For where this fits in the repo, see the [Rib section of README.md](README.md#rib-celltype-rib).

> `buildRibUnitCell_LN.m` keeps its historical name so that every existing caller keeps working. It now contains the merged builder, and its default material is diamond. `buildRibUnitCell.m` is an unrelated device (`'2D_ribs'`, a slab with T-shaped voids).

---

## 1. Geometry overview

**Units.** Every length in `P` is in SI metres (e.g. `P.a = 300e-9`). The sidewall angle `P.d` is in **radians**. This page quotes nm and degrees for readability.

The beam runs along **x**. A continuous **spine** of width `s` carries **one rib per cell**. Each rib is centred on its cell, is `t` wide in x, and reaches `|y| = L` from the beam axis (tip to tip `w = 2L`). The slab is `th` thick in **z**, centred on z = 0. All widths are measured at the bottom face, z = −th/2.

```
 y
 ^            |<-- t -->|
 |   L  ---   +---------+        r2 : convex fillet on each rib-tip corner
 |            |         |
 |            |         |        r1 : concave fillet where rib meets spine
 |            |         |
 | s/2 -- +---+         +---+
 |        |      spine      |
 +--------+--------+--------+----> x      (cell centre at x_n; one cell shown)
 |        |<------ a ------>|
          (mirror image below y = 0)
```

There are two ways the cell gets built:

- **Vertical walls (`d = 0`, default).** The template's 2D feature chain with fillets (`addRibCell.m`), extruded by `th`. Each cell is first drawn off to the side at `y0 = 4L` and then moved into place, so each fillet's selection box catches only that cell's corners.
- **Sidewall angle (`d > 0`).** The former LN construction (`addRibPrisms.m`): trapezoid cross-sections for the spine, the rib flanks and the rib tip. Every wall leans in by `th·tan(d)` toward the top face. **In-plane fillets cannot be combined with this, so `r1 = r2 = 0` is required, and it is the default when `d ≠ 0`.** The slab is not mirror-symmetric in z, so z-symmetry is rejected.

In the unit cell the spine spans exactly `x ∈ [-a/2, a/2]`, so the Floquet faces at `x = ±a/2` cut only the spine. The old LN cell put half-ribs on `x = ±a/2`. That is the same lattice shifted by a/2, so the band structure is identical, but the faces differ.

---

## 2. Parameter tables

Column key:

- **Pipe**: U = unit cell, C = cavity, U+C = both.
- **I/O**: *in* = input with a default; *req* = required input with no default; *out* = computed and written back.
- **Src**: where the default comes from (see the footnotes below the tables).

### 2.1 Cell geometry

| Field | Meaning | Default | Pipe | I/O | Src |
|---|---|---|---|---|---|
| `celltype` | Selects this device in every dispatcher | `'rib'` | U+C | in | |
| `xsect` | Cross-section label; only `'rect'` is supported (use `d` for sloped walls). Not validated | `'rect'` | U+C | in | |
| `a` | Lattice constant. Unit cell: the cell width. Cavity: the **mirror** lattice constant | 300 nm | U+C | in | T1 |
| `t` | Rib width along **x**, bottom face | 100 nm | U+C | in | T2, D |
| `s` | Spine width along **y** (full width), bottom face | 100 nm | U+C | in | D, T3 |
| `d` | Sidewall angle in **radians**; 0 = vertical walls | 0 | U+C | in | LN |

### 2.2 Rib length (the rib's extent along y)

| Field | Meaning | Default | Pipe | I/O | Src |
|---|---|---|---|---|---|
| `L` | Explicit rib half-length. **Highest precedence** | unset | U+C | in / **out** | |
| `w` | Explicit tip-to-tip width, = 2L. Used if `L` is unset | unset | U+C | in / **out** | LN |
| `L_ref` | Rule: L at `a = a_ref` | 400 nm | U+C | in | T4 |
| `a_ref` | Rule: anchor point. **Not** the lattice constant | 240 nm | U+C | in | T4 |
| `dLda` | Rule: slope dL/da | −0.4998 | U+C | in | T5 |
| `ribDerived` | Records the `L`/`w` values a builder wrote back (§3.1) | — | U+C | **out** | |

### 2.3 Fillets

| Field | Meaning | Default | Pipe | I/O | Src |
|---|---|---|---|---|---|
| `r1` | Concave fillet where rib meets spine; 0 = none | 40 nm if `d = 0`, else 0 | U+C | in | T6, TODO |
| `r2` | Convex fillet at the rib tip; 0 = none | 40 nm if `d = 0`, else 0 | U+C | in | T6, TODO |

Both radii are global, applied to every cell.

### 2.4 Cavity taper (cavity only)

| Field | Meaning | Default | I/O |
|---|---|---|---|
| `ncellHalf` | Cells in the x ≥ 0 half; the centre cell counts once | 14 | in |
| `ribatctr` | 1 = rib on x = 0; 0 = spine gap on x = 0 | 1 | in |
| `ndef` | Number of defect (tapered) cells: cells k ≤ `ndef` | 5 | in |
| `maxdef` | Centre a = `a·(1 − maxdef)` | 0.2 (→ 240 nm) | in |
| `taperFunc` | `'quadratic'`, `'cubic'`, `'linear'` or `'none'` (as `CreateNanobeamGeom`) | `'quadratic'` | in |
| `taperConvention` | `'centerspacing'` or `'cellwidth'` (§3.2) | `'centerspacing'` | in |
| `ribGeomHalf`, `ribGeom`, `ribGeomCols` | Per-cell table `[x a t L s xL xR]` for the half / full beam, and its column indices | — | **out** |
| `a_cell`, `L_cell`, `sp_cell`, `w_hole` | Half-beam a_k, L_k, centre spacings, and `2·L_k` | — | **out** |
| `beamLenHalf`, `beamLen` | Beam length (half / full) | — | **out** |

Sources: the defaults reproduce the template's cell positions (T1); `taperConvention` was the user's choice.

### 2.5 Thickness, air and z-offset

| Field | Meaning | Default | Pipe | I/O | Src |
|---|---|---|---|---|---|
| `th` | Slab thickness | **110 nm** (nominal design 120 nm) | U+C | in | T7 |
| `airrad` | Cavity air-cylinder radius | 3.1 µm | C | in | T8 |
| `airZOffset` | Height of the cavity air-cylinder axis above the beam mid-plane | 60 nm | C | in | T7 |
| `airradUC` | Unit-cell air radius; used only with `run_optical = 1` | `5·a` | U | in | LN |
| `run_optical` | Unit cell: 1 = add the air half-cylinder (optical runners only) | unset (0) | U | in | LN |
| `lbd0` | Model parameter `lbd0`, read by the LN refractive index in `runOpticalBand_1D`/`_3D` | 1.55 µm | U | in | LN |

### 2.6 Material and crystal

| Field | Meaning | Default | Pipe | I/O |
|---|---|---|---|---|
| `beamMat` | `'diamond'`, `'LN'`, … (`LoadMaterialParams`) | `'diamond'` | U+C | in |
| `anisoMat` | 1 = anisotropic stiffness tensor `P.D` | 1 | U+C | in |
| `rxtal` | In-plane crystal rotation in degrees. The LN scripts set 45 | 0 | U+C | in |
| `nbeam` | Refractive index (cavity optics; set by `LoadMaterialParams` if absent) | 2.4028 | C | in |

### 2.7 Symmetry and boundary conditions

| Field | Meaning | Pipe | I/O |
|---|---|---|---|
| `mbeveny`, `mbevenz` | ≠ 0 builds only y ≥ 0 / z ≥ 0. Missing = 0. `solveBands` sets them per symmetry sector; `solveBands_noSym` does not use them | U | in |
| `TwoSymPlanes`, `zSymCondition`, `solveasym` | How `solveBands` steps through sectors. With `d ≠ 0` use `TwoSymPlanes = 0`, `zSymCondition = 0`, or `solveBands_noSym` | U | **req** (solveBands) |
| `mevenx`, `meveny`, `mevenz`, `oevenx`, `oeveny`, `oevenz` | Cavity symmetries. `oevenx` is overridden when `ribatctr = 1` (§3.3). With `d ≠ 0` z symmetry is rejected | C | **req** |
| `mEndBC` | Cavity mechanical end face: `'lowreflecting'` or `'fixed'` (§3.5) | C | in (`'lowreflecting'`) |
| `fixed_bc` | Must stay 0. There is no flat outer y face, so no `yEnd2` | U | **req** (runBands) |
| `solveMechPML`, `asymCav`, `asym` | Must be 0 (`asym` is read by `CreateFileBase`) | C | **req** |

### 2.8 Solver, mesh and I/O fields

These are the standard fields for each runner; the test scripts show working values.

- **Unit cell, mechanics:** `kpts`, `nbands`, `freq`, `meshSize`, `max_dof`, `holeatedge` (defaulted), `plotgeom`, `savedat`, `savebndplot`, `saveplots`, `saveMPH`, `bandStruct_2D`, `completeBandGaps`, `datLoc`, `fileBase`/`prefname`.
- **Unit cell, optics:** adds `optical_freq`, `bandStructureDim` (the LN scripts use 1 → `runOpticalBand_1D`), and **`unitcell`**. `unitcell` is read by the optical runners; the LN scripts set `'rectrangular'` (sic); anything except `'hexagonal'` gives the square-lattice path.
- **Cavity:** see `nanobeam/test_nanobeamRibFEM.m`: `solveMech`/`solveOpt`, `lambda`, `freq`, mesh fields, `calcG`/`calcS`, plot/store flags.

**Builder outputs** (all **out**):

| Field | Pipe | Content |
|---|---|---|
| `L`, `w` | U | Rib length actually used (+ `ribDerived`) |
| `w` | C | `2·max(L_k)` (+ `ribDerived.w` when computed) |
| `domSel.beam` | U+C | Beam domains |
| `domSel.air` | U | Air domains |
| `domSel.cyl` | C | Air-cylinder domains |
| `bndSel.Zsym` | U | z = 0 face(s) |
| `bndSel.airOuter` | U | Curved outer air surface |
| `zEnd2` | U | = `bndSel.airOuter`; the 3D optical runner's scattering boundary |
| `xEnd1`, `xEnd2`, `yEnd1`, `zEnd` | U | Boundary-index lookups |
| `bndSel.beam*`, `bndSel.cyl*`, `geomname`, `xc` | C | As `BuildNanobeamFEMGen` |

Named selections:

- Unit cell: `geom1_xboundaries_bnd`, `geom1_yboundaries_bnd`, `geom1_ZsymSel`, `geom1_beamDomSel`, `geom1_airDomSel`, `geom1_airOuter_bnd`.
- Cavity: `beam_beamSel`, which `CalcGOM` uses by name.

**Footnotes (sources).**
- **T** = `comsol_templates/visible_rib_omc.m`.
  - T1: final cell set's tag positions; mirror pitch 300 nm.
  - T2: `r1_0_2` size (line 23175).
  - T3: IR set half-spine height of 50 nm.
  - T4: centre cell, L = 400 nm at a = 240 nm.
  - T5: the IR generator's rule L = 600 − 0.4998·(a − 600) nm.
  - T6: IR set fillets.
  - T7: z-offset defect (§3.4).
  - T8: `cyl_cut`.
- **D** = diamond rib sweep (`rib_optimize_sweep_diamond.m`, `test/diamond_ribUnitCell*/`).
- **LN** = the former `buildRibUnitCell_LN` and its scripts.
- **TODO** = the visible set's value is not in the export.

---

## 3. Notes on the non-obvious parameters

### 3.1 Rib-length precedence and the L(a) rule
`ribLength.m` decides L, in this order:

1. **`P.L`** (explicit half-length) → L = `P.L` for every cell.
2. **`P.w`** (explicit tip-to-tip width, the LN scripts' field) → L = `P.w/2` for every cell.
3. Otherwise the **rule** `L(a) = L_ref + dLda·(a − a_ref)`.

Setting both `P.L` and `P.w` with `P.w ≠ 2·P.L` is an error.

- **Rule example:** a = 300 nm → L = 400 − 0.4998·60 = **370.012 nm** (740 nm tip to tip). `a_ref` is only the anchor point of the line. `dLda = 0` gives L = `L_ref` everywhere.
- **Write-back.** The builders write the length they used back to `P.L` and `P.w` (cavity: `P.w = 2·max(L_k)`). Fields they *computed* are recorded in `P.ribDerived`, and a value equal to its recorded derived value is not treated as explicit. Re-running a builder on its own output, or after changing `a`, therefore still follows the rule. Setting a new `P.w` makes it explicit again.
- **Cavity.** An explicit `P.L`/`P.w` gives every cell the same rib length (no taper in L).

### 3.2 `taperConvention`
Rib dimensions are the same in both conventions:

- cells k ≤ `ndef`: `a_k = a·(1 − maxdef·f((k−1)/ndef))`;
- cells beyond: `a_k = a`;
- `L_k = ribLength(P, a_k)`.

The conventions differ only in how cells are placed:

- **`'centerspacing'` (default).** Exactly `CreateNanobeamGeom.m`.
  - `x_{n+1} = x_n + a_{n+1}`; the first spacing out of the centre is `a·(1 − maxdef)`.
  - Spine segments run between midpoints of neighbouring centres.
  - Defaults give 0, 240, 483.75, 738.75, 1012.5, 1312.5, …, 3712.5 nm (7.725 µm beam).
- **`'cellwidth'` (the template's).**
  - `x_{k+1} = x_k + (a_k + a_{k+1})/2`, with the rib centred in a spine segment of width `a_k`.
  - Reproduces the template: 0, 241, 487, 743, 1013, 1302, …, 3702 nm (7.704 µm beam).
  - The default lands up to 10.5 nm further out than this.

### 3.3 `ribatctr` and optical parity
`SetupNanobeamFEM` uses `oevenx = P.oevenx^(1 − P.ribatctr)`. A rib on the centre puts dielectric at x = 0, the analogue of `holeatctr = 0`, so `oevenx` is forced to 1. With `ribatctr = 0`, `P.oevenx` is used as given.

### 3.4 Why `th` = 110 nm and `airZOffset` = 60 nm
In the template, the work plane stayed at z = −125 nm after the extrusion was changed to 120 nm, while its z-symmetry cut sits at z = −60 nm. It therefore solves an **effective 110 nm beam**, with the air-cylinder axis 60 nm above the symmetry plane. The defaults reproduce that model. **For the intended design use `th = 120e-9` and `airZOffset = 0`.** Both values are irrelevant to the LN scripts, which set their own `th` and do not use the cavity.

### 3.5 `mEndBC`
This is the cavity's mechanical end-face condition.

- `'lowreflecting'` (default for `'rib'`, as in the template) adds a LowReflectingBoundary. Expect complex eigenfrequencies.
- `'fixed'` clamps the face, the behaviour every other celltype keeps by default.

### 3.6 `r1` and `r2` naming, and fillets with a sidewall angle
In `P`, `r1` is the concave **root** fillet and `r2` the convex **tip** fillet. The template's feature tags run the other way (`fil2_*` = tip, `fil1_*` = root), and `addRibCell` keeps those tags.

Fillets are drawn in the x–y plane and extruded, which only works with vertical walls. With `d > 0` they are therefore not supported: they default to 0, and setting them raises `rib:filletWithSidewall`. The old LN builder had no active fillets either (its `fil1` was disabled and referenced an undefined parameter).

---

## 4. Validation rules

Most checks are in `validateRibParams.m`; the others are in the file named in the table. Define:

- a_min = `a` (unit cell) or `a·(1 − maxdef)` (cavity);
- L_min, L_max = the rib half-length at the two ends of [a_min, a] (constant if `L`/`w` is explicit).

| Condition | Result |
|---|---|
| `a`, `t`, `s`, `th`, `L_ref` or `a_ref` not a positive finite scalar; `dLda` not finite; `r1`/`r2` negative; cavity counts or ranges invalid | error, `MATLAB:validateRibParams:*` |
| `d` not in [0, π/2) (e.g. passed in degrees) | error `rib:badSidewallAngle` |
| `taperConvention` / `mEndBC` unknown (cavity) | error `rib:badTaperConvention` / `rib:badEndBC` (also `SetupNanobeamFEM:badEndBC`) |
| L ≤ 0; or both `L` and `w` given with `w ≠ 2L` (`ribLength.m`) | error `rib:nonPositiveL` / `rib:inconsistentLW` |
| `t` ≥ a_min | error `rib:ribTooWide` |
| L_min ≤ `s`/2 | error `rib:ribShorterThanSpine` |
| `d > 0` and (`r1` > 0 or `r2` > 0) | error `rib:filletWithSidewall` |
| `d > 0` and `t − 2·th·tan d ≤ 0`, or `s − 2·th·tan d ≤ 0` | error `rib:sidewallTooSteep` |
| `d > 0` with z symmetry (unit cell: `mbevenz ≠ 0`; cavity: the z cut would be on, in `buildRibCavity.m`) | error `rib:zSymWithSidewall` |
| `d = 0`: 2·`r2` ≥ `t` | error `rib:tipFilletTooLarge` |
| `d = 0`: `r1` ≥ (a_min − `t`)/2 | error `rib:rootFilletTooLarge` |
| `d = 0`: `r1` + `r2` ≥ L_min − `s`/2 | warning `rib:filletsMeet` |
| `d = 0`: 4·L_min + `s`/2 − 10 nm ≤ L_max (a cell's off-to-the-side drawing position overlaps placed cells) | warning `rib:stagingOverlap` |
| Air does not enclose the beam (cavity `airrad` − \|`airZOffset`\|; unit-cell `airradUC`) | error `rib:airTooSmall` |
| Cavity: spine segment ≤ `t` / unknown `taperFunc` (`createRibCavityGeom.m`) | error `rib:segmentTooNarrow` / `rib:badTaperFunc` |
| Cavity: missing fields, PML or asymmetric cavity requested (`buildRibCavity.m`) | error `rib:missingFields` / `rib:noPML` / `rib:noAsym` |
| Root helpers not on the path (nanobeam side) | error `rib:pathMissing` |
| `run_optical = 1` passed to the mechanics runners | error `runBands:ribOptical` / `runBands_noSym:ribOptical` |
| A selection or domain set resolves empty after a build | warning `rib:selectionEmpty` / `:emptyBndSel` / `:emptyDomSel` |

---

## 5. Output file names (`fileBase`)

| Solver | Name | Notes |
|---|---|---|
| `solveBands` (mechanics, symmetric sectors) | `rib_a_…nm_s_…nm_w_…nm_t_…nm_th_…nm_d_…deg_r1_…nm_r2_…nm_rxt_…_` | New. Previously `solveBands` had no rib branch at all |
| `solveBands_noSym` (mechanics, LN scripts) | `hole_unitCella_…nm_s_…nmw_…nm_t_…nm_th_…nm_d_…deg` + `_r1_…nm_r2_…nm` **only if a fillet is nonzero** | Legacy format kept, quirks included, so existing LN data under `test/` still matches |
| `solveOpticalBands` | `optical_rib_a_…nm_s_…nmw_…nm_t_…nm_th_…nm_d_…deg` + the same conditional fillet suffix | Legacy format kept |
| `CreateFileBase` (cavity, both copies) | `rib_a…nm_t…nm_s…nm_<L>_r1_…nm_r2_…nm_sw…deg_nc…_nd…_d…_<taper>_<cs\|cw>_ctr…_az…nm_<lr\|fx>`, then `_X…_Y…_Z…`. `<L>` is `Lr…nm_ar…nm_dL…` (rule) or `Lc…nm` (explicit). The shared tail prepends `th…nm_` (and `rxt…o_` when `rxtalInFilename = 1`) | |

- `w` is always the width actually built (explicit, or `2·L(a)`).
- `d` is printed in whole degrees; lengths in whole nm; `dLda` and `maxdef` to 4 decimals.
- **Cached-result collisions.** Runs that differ only in fields a name omits, or by less than the rounding, reuse each other's cached results. Examples: the legacy names omit `beamMat`, `rxtal` and the symmetry settings; the cavity name omits `airrad`, `lambda` and mesh settings. Use `P.prefname`, a fresh `datLoc`, or a hand-set `P.fileBase`.

---

## 6. Usage

```matlab
% Unit cell, mechanical bands, visible design (root folder): test_rib_unitCell.m
P = ribDefaultParams();
P.mbevenz = 1;  P.mbeveny = 1;  P.TwoSymPlanes = 1;  P.zSymCondition = 1;  P.solveasym = 1;
P.kpts = 10;  P.nbands = 20;  P.freq = 0;  P.meshSize = 4;  P.max_dof = 3e6;
P.fixed_bc = 0;  P.run_optical = 0;
P.savedat = 1;  P.savebndplot = 1;  P.saveplots = 0;  P.saveMPH = 0;
P.bandStruct_2D = 0;  P.completeBandGaps = 1;
P.datLoc = [fullfile('.','test','rib',datestr(now,'mmddyyyy')),filesep];
bds = solveBands(P);

% LN-style cell with sidewalls (as rib_optimize_sweep*.m / test_rib_unitCell_LN.m):
%   P.beamMat = 'LN'; P.w = 1747e-9; P.d = 15*pi/180; P.mbevenz = 0; ...
%   bds = solveBands_noSym(P);                       % mechanics
%   P.run_optical = 1; bds = solveOpticalBands(P);   % optics (bandStructureDim = 1)

% Cavity (nanobeam folder): nanobeam/test_nanobeamRibFEM.m, which also adds the root
% folder with addpath(<omc-comsol-chang>,'-end')
[ds, model] = RunNanobeamFEM(P, datLoc);
```

- Run scripts:
  - `test_rib_unitCell.m`: diamond visible design, L(a) rule, mechanical bands (`solveBands`, all four symmetry sectors), output in `test/ribUnitCell/`.
  - `test_rib_unitCell_explicitW.m`: explicit `P.w`, d = 0, 20 nm fillets, mechanical (`solveBands`, `zSymCondition = 0`) or optical bands, output in `test/ribUnitCell[_optical]/`.
  - `test_rib_unitCell_LN.m`: LiNbO3, d = 15°, explicit `P.w`, no fillets (the old LN shape), optical (`solveOpticalBands`) or mechanical (`solveBands_noSym`) bands, output in `test/LN_ribUnitCell[_optical]/`, which reuses the legacy LN results.
  - `nanobeam/test_nanobeamRibFEM.m`: the cavity.
  - `optimize_rib_unitCell.m`: geometry optimization of the unit cell for a target gap (§8).
- Tests: `runtests('RibUnitCellTest','Tag','PureMATLAB')` (no COMSOL) and `runtests('RibUnitCellTest','Tag','COMSOL')` (builds geometry; skipped without LiveLink).
- Both test scripts end in a commented-out geometry block that finishes with `mphlaunch(model)`.

---

## 7. Compatibility with the old LN rib scripts, and limitations

**Unchanged for the LN scripts:**
- the celltype `'rib'` and the entry point `buildRibUnitCell_LN`;
- `P.w`, and `P.d` in radians;
- the sidewall construction and the air half-cylinder (`5a`, y ≥ 0);
- the `solveBands_noSym` and `solveOpticalBands` file names;
- with `d ≠ 0` and no fillet fields, the geometry is the LN shape.

**Different:**
- The rib is centred in the cell (band-equivalent; the Floquet faces cut only the spine).
- The builder no longer sets model parameters a/s/w/t/th/d/airrad/k/kx (the runners own those) or creates result tables tbl1–4. Those tables used to collide with `runOpticalBand_3D`'s own.
- The optical runners now assign air and beam materials and the scattering boundary from the builder's named selections, instead of the hard-coded domains `[1]`, `[2..10]` and boundaries `[3 12]`.
- An LN call **with `d = 0` and no fillet fields** now gets 40 nm fillets (the visible default). Set `r1 = r2 = 0` for the legacy shape.

**Limitations / TODOs:**
- The visible set's fillet radii are unknown (TODO).
- None of this has been run in COMSOL; check new geometry with `mphlaunch`.
- `runOpticalBand_3D` still puts a Floquet periodic condition on the single y = 0 plane for `'rib'` (as before). A symmetry plane would be the right condition.
- Fillets and an explicit L are global (no per-cell override in the cavity).
- No mechanical PML or asymmetric cavity.

---

## 8. Optimization (`optimize_rib_unitCell.m`)

`optimize_rib_unitCell.m` is a function modelled on `optimize_hole_unitCell.m`. It runs a Nelder–Mead (`fminsearch`) search over selected rib dimensions to maximize a fitness computed from the `solveBands` band structure. Edit the **TUNABLE KNOBS** at the top, or pass overrides:

```matlab
[xBest, fitnessBest, dsBest, PBest, OPT] = optimize_rib_unitCell();                 % defaults
optimize_rib_unitCell(struct('vars',{{'t','s','L','a'}}, 'targetFreq', 10e9));      % OPT overrides
optimize_rib_unitCell(struct(), struct('rxtal',45));                                % P overrides
optimize_rib_unitCell(struct('dryRun',true));                                       % no COMSOL
```

Overrides are merged over the defaults. Nested structs (`OPT.init`, `OPT.lb`, `OPT.ub`, `OPT.fminOpts`) merge field by field.

### Design variables and bounds

`OPT.vars` is any subset of `{'t','s','L','a','d'}`; the default is **`{'t','s','L'}`** with `a` fixed. `L` is applied as an explicit `P.L`, which wins over the L(a) rule. Variables not listed keep their `P` value.

| Var | Init | Bounds | Source |
|---|---|---|---|
| `t` | 100 nm | [90, 200] nm | `ribDefaultParams`; `t0` in `rib_optimize_sweep_diamond.m`; logged runs in `test/diamond_ribUnitCell*/` (100–105 nm). Lower bound keeps 2·r2 < t (40 nm fillets); upper bound keeps r1 < (a − t)/2 at a = 300 nm |
| `s` | 100 nm | [50, 200] nm | as `t`. 50 nm is the fabrication floor used in `rib_optimize_sweep_diamond.m` |
| `L` | 370 nm | [250, 500] nm | L(a = 300 nm) from the rule; the logged runs used w = 2L = 600–840 nm |
| `a` | 300 nm | [240, 360] nm | `ribDefaultParams` (template mirror pitch); the logged runs used 250–262 nm |
| `d` | 0 | [0, 15°] (radians) | vertical walls; the LN scripts use 15°. Requires `r1 = r2 = 0` (passed as P overrides) |

- **Normalization.** The search runs in normalized coordinates v ∈ [1, 2], mapped linearly onto [lb, ub]. With the offset, `fminsearch`'s default 5 % initial simplex becomes a uniform 5–10 % step of every range. Proposals outside the box are clamped (as in the hole optimizer) and penalized by `OPT.boxPenaltySlope` per unit of normalized distance.
- **Optional constraint (off by default).** `OPT.fixFootprint = true` holds the footprint area A = a·s + t·(2L − s) at its initial value by solving for `L`, the analogue of the hole optimizer's fixed ellipse area. `L` must then not be in `OPT.vars`.

### Target, fitness and symmetry
- **`OPT.targetFreq` = 11.2 GHz:** the visible template's mechanical eigen-shift (`comsol_templates/visible_rib_omc.m:23488`, unit GHz at `:1262`).
- **`OPT.fitnessFcn`** takes a handle `@(ds, targetFreq)` or one of these names:
  - `'closestGap_sym'` (default): size of the symmetric gap (`ds.sym`) nearest the target;
  - `'closestGapRatio_sym'`: gap-to-midgap ratio of that gap. It is dimensionless, so lower `TolFun` to about 1e-3;
  - `'closestGap'`: complete gap (`ds.full`). This is complete only over the sectors actually solved, so set `P.TwoSymPlanes = 1` (d = 0 only) or `P.solveasym = 1`.
- **Symmetry defaults:** `TwoSymPlanes = 0`, `zSymCondition = 1`, `mbeveny = 1`, `solveasym = 0`. `solveBands` then solves one sector into `ds.sym`: y-even **and** z-even, i.e. the optomechanical breathing modes, on a quarter model. If `d` is optimized or nonzero, the script switches to `zSymCondition = 0`, `mbevenz = 0` (y-even, both z parities), because angled walls are not z-symmetric.

### Feasibility
- Before any COMSOL call, each candidate goes through `ribLength` and `validateRibParams(P, 'unitcell')`. The checks use `mbevenz = 1` when `zSymCondition = 1`, because that is the symmetry `solveBands` will impose.
- An invalid candidate is not solved. It scores `OPT.noGapPenalty·(1 + violation)` and is logged as infeasible with its error ID (`rib:ribTooWide`, `rib:rootFilletTooLarge`, …).
- The violation is a sum of relative constraint excesses: t/a, s/(2L), 2·r2/t, (t + 2·r1)/a, and for d > 0 the top-face widths. It grows monotonically outside the feasible region, so the simplex is pushed back rather than sitting on a flat plateau.
- A candidate with no usable gap scores `OPT.noGapPenalty`. If nothing feasible is found, `fitnessBest` is NaN (warning `optimize_rib_unitCell:noUsableDesign`) and the final re-run is skipped.

### Output
- `OPT.rootDir` defaults to `test/ribUnitCell_opt/opt_<mmddyyyy_HHMMSS>/`. The **time stamp** keeps a second run on the same day from reusing the first run's cached results; the hole optimizer has that problem.
- Each feasible evaluation solves in its own `eval_NNNN/` folder. The best design is re-run with plots in `BEST/`.
- **`opt_log.csv`** columns: `eval`, each variable (`t_nm`, …, `d_deg`), its normalized value (`v_t`, …), `feasible`, `errorID`, `boxDist`, `fitness`, `J`. Rows are appended per evaluation, so the log survives a crash.
- All paths use `fullfile`/`filesep`. On macOS/Linux, `solveBands` appends a literal `'\'` to `datLoc`, which creates an empty sub-folder named `\`. The script removes it after each solve. Its cache lookup never hits there, which is harmless because each evaluation uses a fresh folder.

### Dry run
`OPT.dryRun = true` replaces `solveBands` with an analytic surrogate: a smooth bump with its optimum at t = 130, s = 80, L = 400 nm and a 2 GHz gap. **This is not physics.** It exercises normalization, clamping, the feasibility penalty, the CSV log and the BEST path without COMSOL. `RibUnitCellTest` runs it under the `PureMATLAB` tag.

### Caveats
- **`P.rxtal` defaults to 0** (beam along ⟨100⟩, as in the template and the cavity) and is a `TODO(confirm)` knob. At 45° (⟨110⟩) the in-plane shear velocity drops about 9 % and the longitudinal velocity rises about 5 %, so the optimized gap moves by several percent. Use the same value as the cavity runs.
- Each evaluation is a full `solveBands` run (one symmetry sector at the defaults). Budget accordingly with `MaxFunEvals`.
