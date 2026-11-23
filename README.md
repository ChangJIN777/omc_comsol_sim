# omc_comsol_sim

FEM simulation code for **diamond optomechanical crystals (OMCs)** — photonic and
phononic band structures of periodic unit cells, and full nanobeam cavity
simulations with optomechanical / SiV strain coupling. Everything here targets
COMSOL Multiphysics, driven either from MATLAB (LiveLink) or from Python (MPh).
The repo holds **several parallel pipelines at very different maturity levels**:
one working MATLAB pipeline that is the physics ground truth, one mature Python
optimization loop for the 1D nanobeam cells (elliptical hole and rib), one in-progress Python port for the
2D "boomerang" cell, plus vendored upstream scripts and a large pile of
committed simulation output. Read the table below before running anything.

---

## Directory map

| Path | What it is | Language | Status |
|---|---|---|---|
| [`omc-comsol-chang/`](omc-comsol-chang/) | The working pipeline: photonic + phononic band structures (root) and nanobeam FEM cavities (`nanobeam/`), including the `'rib'` device (spine + centred rib; see [`rib_README.md`](omc-comsol-chang/rib_README.md)). Derived from `omc-comsol-master/`. **Physics ground truth.** | MATLAB + COMSOL LiveLink (194 `.m`) | **Active** |
| [`python_script_2D/`](python_script_2D/) | Python/MPh closed-loop optimizer for the 2D hexagonal-lattice "boomerang" phononic cell. | Python (12 `.py`) | **Active, not yet runnable end-to-end** |
| [`python-scripts/`](python-scripts/) | Python closed-loop optimizer for the 1D nanobeam unit cell (optical @1550 nm + mechanical @8 GHz): the elliptical-hole cell and, via `celltype: rib`, the rib cell. Has real results for the hole cell. | Python (35 `.py`) | **Active** |
| [`comsol_templates/`](comsol_templates/) | COMSOL GUI-session exports used as sources for builders: `visible_rib_omc.m` (origin of the visible-wavelength rib design) and `BoomerangWithAirDisk.m`. | MATLAB | Reference |
| [`omc-comsol-master/`](omc-comsol-master/) | Vendored copy of the Loncar-group shared `omc-comsol` scripts. Untouched since 2024-07-23. | MATLAB + a little Python/MPB (82 `.m`) | Reference / upstream |
| `bands/` | Earlier standalone band-structure code: 2D "Cross" unit cell, 1D opt+mech nanobeam. Carries its own duplicate copies of `bndindex.m`, `LoadMaterialParams.m`, `RotateXtalTensor.m`, `findGaps.m`. | MATLAB (25 `.m`) | Legacy — superseded by `omc-comsol-chang/` |
| `cavities/` | Empty. The only tracked file is a macOS `.DS_Store`. | — | Placeholder, delete |
| `test/` | Committed output of 2024-07 MATLAB runs (BlockTet, boomerang, cross, Holely, rib): 853 `.png`, 19 `.fig`, 23 `.mat`. 84 MB, 895 tracked files. | data | Output data — not code |
| `photonic_bandStructure_code.m` | Orphaned COMSOL 6.1 GUI export (Jul 2024) of a **2D photonic** band structure (`ElectromagneticWavesFrequencyDomain`, Floquet on `kx`/`ky`, rhombic triangular-lattice cell, `a=360 nm`). Nothing references it. | MATLAB | Reference — useful starting point for the unbuilt 2D optical model |

### How `omc-comsol-chang/` relates to `omc-comsol-master/`

`omc-comsol-master/` is the upstream: its README says "avoid editing the master
branch," and credits the original COMSOL v3.5 scripts to Michael Burek (from
Oskar Painter's group), the v5.x port to Cleaven Chia, and the MPB interface to
Graham Joe. `omc-comsol-chang/` was added a fortnight after it and is a
working fork — 50 shared `.m` filenames, with `bndindex.m`, `edgeindex.m` and
`RotateXtalTensor.m` still byte-identical, while `LoadMaterialParams.m` and
`findGaps.m` have diverged. **Edit `omc-comsol-chang/`, read
`omc-comsol-master/` when you need the original behaviour.**

---

## Which pipeline to use

### `omc-comsol-chang/` — MATLAB, the one that works

The reference implementation for geometry, boundary conditions, and material
physics. Every Python module here is a port of something in this directory; when
Python and MATLAB disagree, MATLAB is right until proven otherwise.

No build system. Run a front-panel `test_*.m` script directly in MATLAB, with
COMSOL LiveLink already connected:

```matlab
run('test_Boomerang.m')                      % 2D boomerang phononic band structure
run('nanobeam/test_nanobeamRectFEM.m')       % full nanobeam cavity FEM
```

Flow: `test_*.m` → `solveBands.m` → `runBands.m` (1D cells, including `'rib'`)
or `runBands_2D.m` (2D cells) → `build*UnitCell.m` → `findGaps.m`. See
[`omc-comsol-chang/README.md`](omc-comsol-chang/README.md) and
[`omc-comsol-chang/CLAUDE.md`](omc-comsol-chang/CLAUDE.md).

The **`'rib'` device** (spine + centred rib, optional sidewall angle `P.d`) runs
through both the band-structure and nanobeam-cavity pipelines. Its parameters,
defaults and validation rules are in
[`omc-comsol-chang/rib_README.md`](omc-comsol-chang/rib_README.md). Related tools:

- `optimize_rib_unitCell.m`: `fminsearch` optimizer for the rib unit cell, with
  a dry-run mode.
- `trackRibBands.m`: a standalone tool that follows the two gap-edge mechanical
  bands at X across an ordered list of designs by mode shape, bridging ambiguous
  matches with midpoint solves.
- `RibUnitCellTest.m` and `BandTrackingTest.m`: `matlab.unittest` classes. Run
  them with `runtests(...,'Tag','PureMATLAB')`, no COMSOL needed, or
  `'Tag','COMSOL'`.

### `python_script_2D/` — 2D boomerang optimizer, **not runnable end-to-end yet**

The 42 numpy-only tests pass (`conda run -n omc python3 tests/test_pipeline_2d.py`),
and geometry mapping, feasibility gating, band-gap extraction, scoring and a
Bayesian-optimization loop are all implemented and tested. The design space is 3 free variables `(a, w, r)`, with
`r1`, `r2` and `th` fixed in `configs/bounds_2d.yaml:fixed`.

`comsol/trusty_boomerang.mph` **has been rebuilt** from the parameterized export
script and is now parametric — it carries all 13 parameters
(`a,w,r,r1,r2,th,hx,hy,selw,dsel,k,kx,ky`) and the `leg1/leg2/leg3` geometry
tags, saved solution-free at 2.7 MB instead of the original 1.1 GB. Three things
still block a sweep:

- **Study labels don't match.** The template has `Study_symmetric` and an
  unlabelled study (COMSOL reports `Study 1`); `src/acoustic_comsol_2d.py`
  resolves studies by label and expects `mech evenz` / `mech oddz`, so
  `_require_studies` raises. Relabel in the GUI, or override for a smoke test:
  `OMC2D_STUDY_EVENZ="Study_symmetric" OMC2D_STUDY_ODDZ="Study 1"`.
- **The BCs are pinned to absolute COMSOL face indices** (`pbcX [1 22]`,
  `pbcY [2 9]`, parity `[3]`), valid only for the default geometry. COMSOL
  renumbers faces when topology changes, and the solve then succeeds against the
  wrong walls without complaining. Now that the geometry *does* rebuild, this is
  live rather than dormant — **the single highest risk in the project.** Fix by
  porting `bndindex.m`; the already-correct `ZsymSel` selection can take over
  for the parity plane immediately.
- **The k-sweep is still wrong.** The driver writes `kx`/`ky` as literals per
  point, destroying the COMSOL expressions in `k` that the study sweeps. Set `k`
  and read the parametric dataset instead.

[`python_script_2D/docs/TODO.md`](python_script_2D/docs/TODO.md) is the
authoritative status and roadmap — start there, not with the README.

### `python-scripts/` — 1D nanobeam optimizer, working

Closed-loop optimization of a 1D diamond OMC unit cell for a simultaneous optical
gap near 1550 nm and a mechanical gap near 8 GHz. Two cell types, selected with
`celltype:` in the run config (default `hole`):

- **`hole`** (rectangular nanobeam + elliptical hole): 5 design variables
  `(a, w, hx, hy, t)`, where `t` = slab **thickness**. It has three optical
  backends (numpy TMM surrogate, MPB, COMSOL) and a COMSOL mechanical backend.
  Results are committed under `results/`: SQLite run stores, Optuna/TPE result
  JSON, band and mode figures, and several `_archive_*` generations.
- **`rib`** (spine + centred rib, template `comsol/omc_unitcell_rib.mph`, built by
  `comsol/build_omc_unitcell_rib.m`):
  - **Search space:** `a, w, t, s` are searched; `th = 250 nm` and
    `r1 = r2 = 40 nm` are fixed.
  - **Naming:** here `t` = **rib width** and `th` = thickness, deliberately
    different from the hole cell.
  - **Rejected before any solve:** designs that fail the geometric-validity
    rules or the separate fabricability check (`FAB_*` reason codes, thresholds
    in `configs/bounds.yaml`). They are saved, not discarded.
  - **Mechanical score:** the complete gap by default, from the
    `mech sym` + `mech antisym` studies.
  - **Optics:** the `opt TE` study via `src/comsol_optics.py`. It is not yet
    validated on the rib template; run the optical check solve first.
  - **Where it runs and stores results:** `run_one.py` / `run_loop.py` and
    `run_opt_comsol.py` (which dispatches to `scripts/run_opt_comsol_rib.py`),
    with results in a separate `results/runs_rib.sqlite3`. No rib results are
    committed yet.

15 numpy-only tests (`python tests/test_pipeline.py`).

`python_script_2D/` reuses this project's architecture: `comsol_client.py` is
near-verbatim (5 lines differ) and `bandgap.py`, `database.py`, `optimizer.py`
are lightly adapted (mostly docstrings), because that logic is genuinely
geometry-agnostic. Everything physics-specific was rewritten. See
[`python-scripts/README.md`](python-scripts/README.md) and
[`python-scripts/CLAUDE.md`](python-scripts/CLAUDE.md).

---

## Getting started

### Prerequisites

| Need | For |
|---|---|
| COMSOL Multiphysics 6.2+ (6.3 for the 2D template) with **LiveLink for MATLAB** | everything |
| MATLAB R2016b+ | `omc-comsol-chang/`, `bands/`, `omc-comsol-master/` |
| Python 3.10+ with `numpy`, `pyyaml`, `optuna` (+ `pandas`, `matplotlib` for the 1D project) | both Python pipelines |
| `pip install MPh` + a reachable COMSOL server (`comsol mphserver`) | Python → COMSOL |
| conda env `omc` (`conda create -n omc -c conda-forge pymeep pymeep-extras`) | MPB optical backend; also the env this repo's Python is run in |

COMSOL preferences: set the MATLAB installation folder under *LiveLink
Connections*, and *Geometry → Geometry representation in new models* to **COMSOL
kernel** (see `omc-comsol-master/README.md`). Note that
`python_script_2D/comsol/trusty_boomerang_script.m` currently overrides this with
`geomRep('cadps')` (Parasolid), which needs the CAD Import or Design Module.

### Minimal commands

```bash
# MATLAB pipeline — start MATLAB from "COMSOL with MATLAB", then:
#   run('omc-comsol-chang/test_Boomerang.m')

# 1D nanobeam (Python) — surrogate backends, no solver needed
cd python-scripts && pip install -r requirements.txt
python tests/test_pipeline.py
python scripts/run_one.py --config configs/run_one.yaml

# 2D boomerang (Python) — tests only; the COMSOL path is blocked, see docs/TODO.md
cd python_script_2D && conda run -n omc python3 tests/test_pipeline_2d.py
```

The project convention for Python is `conda run -n omc python3 …`.

---

## Conventions

Stated once here because several of these have been documented inconsistently
elsewhere in the repo.

- **SI units everywhere.** MATLAB `P` structs and COMSOL parameters are in
  **meters**, not nm (`P.a = 480e-9`). The Python geometry dataclasses are also
  in meters; only display strings use nm.
- **Fillet convention: `r1` = JUNCTION fillet (inner corners, at `w/2` from the
  hole centre), `r2` = TIP fillet (outer corners, at `hypot(r, w/2)`).** The
  authority is `addFillet()` in `omc-comsol-chang/buildBoomerangUnitCell.m:119-149`.
  This was documented backwards in at least four files and is invisible while
  `r1 == r2`; it silently swaps two design variables as soon as an optimizer
  separates them. Check before you trust any comment naming these.
- **Diamond material constants** — the authority is
  `omc-comsol-chang/LoadMaterialParams.m`:
  `rho = 3500 kg/m^3`, cubic `C11 = 1076`, `C12 = 125`, `C44 = 578 GPa`
  (isotropic approximation `E = 1050 GPa`, `nu = 0.2`). In-plane crystal
  rotation is `rxtal = 45°` about z for the boomerang cell and **`rxtal = 0`
  (beam along ⟨100⟩) for every rib script and the rib cavity**. Any YAML or
  `.mph` carrying different numbers is wrong, not the MATLAB. Known exceptions
  and inconsistencies:
  - `nu = 0.2` does not match the cubic constants, which give ν⟨100⟩ ≈ 0.104.
    It only matters for `anisoMat = 0` runs.
  - The Python hole templates `python-scripts/comsol/omc_unitcell*.mph` are
    isotropic (`E = 1050 GPa`, `nu = 0.1`).
  - The rib `.mph` uses the anisotropic cubic tensor.
- **Brillouin-zone sweep parameterization**: a single scalar `k ∈ [0,3)` with
  `kx`/`ky` as COMSOL *expressions* of `k`. Segments are Γ→M (`k<1`), M→K
  (`1≤k<2`), K→Γ (`k≥2`) for the hexagonal cell. Never overwrite `kx`/`ky` with
  literals — that destroys the expressions and every sweep point then solves the
  same wavevector.
- **Unit cells are half-slabs.** Geometry is built for `z ∈ [0, th/2]` with a
  symmetry (z-even) or antisymmetry (z-odd) BC at `z = 0`; `th` always means the
  **full** slab thickness. Exceptions for the rib:
  - With a sidewall angle (`P.d ≠ 0`) the MATLAB rib cell is not mirror-symmetric
    in z, so it is built at full thickness and z-symmetry is rejected.
  - The Python/MPh templates (`omc_unitcell*.mph`, `omc_unitcell_rib.mph`) use a
    different frame: the beam runs along **z**, the thickness along **x**, and
    the model is the y ≥ 0 half with an air box.

---

## Repo hygiene

`.gitignore` covers `*.mph`, `*.png`, `*.fig`, `*.mat`, `*.mph.lock`,
`*.mph.recovery*`, `.DS_Store`.

**Onboarding trap: `.mph` files are not in version control.** COMSOL templates
must be rebuilt locally from the exported `.m` scripts — e.g.
`python_script_2D/comsol/trusty_boomerang_script.m` regenerates that project's
template via `model = trusty_boomerang_script(); mphsave(model, 'trusty_boomerang.mph')`,
and `python-scripts/comsol/build_omc_unitcell_rib.m` builds
`omc_unitcell_rib.mph` for the Python rib cell.
Treat the exported `.m` as the source of truth and re-export it whenever you
change a model in the GUI, or the change is lost to everyone else.

Known issues worth fixing (the ignore rules were added *after* these files were
committed, and `.gitignore` does not untrack anything):

- `test/` — 895 tracked files / 84 MB of simulation output. `.git` is 370 MB.
- `python-scripts/results/` — 115 tracked files including SQLite run stores.
- 51 tracked `__pycache__` files (mostly `python-scripts/src/__pycache__/`).
- `.claude/settings.local.json` is tracked, and so is now public on GitHub. It
  holds local permission rules with machine paths (no secrets). It should be
  untracked and ignored.
- 8 tracked `.mph` files under `bands/` and `omc-comsol-master/`, plus a tracked
  `bands/mechbands2D/Cross.zip`.
- 12 tracked `.DS_Store` files, including the only file in `cavities/`.
- A stale **1.1 GB** `python-scripts/comsol/trusty_boomerang.mph` — a leftover
  copy of the 2D template that does not belong in the 1D project at all. (The
  live one, `python_script_2D/comsol/trusty_boomerang.mph`, is now 2.7 MB after
  being re-saved solution-free.) Both are correctly git-ignored.
- `python-scripts/` has no `.gitignore` of its own, which is why its `results/`
  tree is tracked; `python_script_2D/.gitignore` already excludes `results/`.
- Four diverging duplicate copies of `bndindex.m`, `LoadMaterialParams.m`,
  `RotateXtalTensor.m` and `findGaps.m` across `bands/`,
  `omc-comsol-master/mechbands1D/`, `omc-comsol-master/mechbands2D/` and
  `omc-comsol-chang/`. `LoadMaterialParams.m` and `findGaps.m` have already
  diverged between master and chang — decide which copy is canonical before
  more drift accumulates.

---

## Further reading

- **2D pipeline status/roadmap (authoritative):** [`python_script_2D/docs/TODO.md`](python_script_2D/docs/TODO.md)
- **2D COMSOL template recipe:** [`python_script_2D/comsol/README_template_2d.md`](python_script_2D/comsol/README_template_2d.md)
- **MATLAB pipeline:** [`omc-comsol-chang/README.md`](omc-comsol-chang/README.md) (structure, builders, solver flow) · [`CLAUDE.md`](omc-comsol-chang/CLAUDE.md) (running & debugging)
- **Rib device (parameters, validation, file names, optimization):** [`omc-comsol-chang/rib_README.md`](omc-comsol-chang/rib_README.md)
- **Upstream scripts, COMSOL/MPB setup, provenance:** [`omc-comsol-master/README.md`](omc-comsol-master/README.md)
- **1D nanobeam project:** [`python-scripts/README.md`](python-scripts/README.md) · [`CLAUDE.md`](python-scripts/CLAUDE.md)
- **2D project overview and 1D↔2D differences:** [`python_script_2D/README.md`](python_script_2D/README.md) · [`CLAUDE.md`](python_script_2D/CLAUDE.md)
