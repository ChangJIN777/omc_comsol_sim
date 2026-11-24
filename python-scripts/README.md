# 1D Diamond OMC Unit-Cell Bandgap Optimization

Closed-loop optimization of a 1D diamond optomechanical-crystal **unit cell**
for a simultaneous optical bandgap near **1550 nm** and a mechanical bandgap
near **8 GHz**. Two cell types, selected with `celltype:` in the run config:

- `hole` (default): rectangular cross-section nanobeam + elliptical hole.
- `rib`: spine + centred rib (see "Rib cell" below).

Everything below the Rib section describes the hole cell unless stated otherwise.

See `docs/strategy.md` for the full strategy, field review, and physics notes.

## Quick start (numpy-only smoke test, runs anywhere)
```bash
pip install -r requirements.txt           # numpy/pyyaml/pandas/matplotlib (+optuna)
python tests/test_pipeline.py             # all 15 tests pass with numpy only (hole + rib)
python tests/test_optical_tidy3d.py       # Tidy3D engine, offline (tidy3d mocked)
python tests/test_xsect_index.py          # xsect switch, parity filter, n = 2.386
python scripts/run_one.py  --config configs/run_one.yaml   # surrogate optical + mech by default
python scripts/run_loop.py --config configs/run_loop.yaml  # surrogate pre-screen by default
```
The surrogate is a coarse numpy TMM model for loop debugging only -- not physically accurate.

Each script takes a single `--config PATH` flag (defaulting to its own YAML in
`configs/`), so `python scripts/run_one.py` runs with `configs/run_one.yaml`.
Edit the YAML to change backends, iteration counts, the `u` vector, etc.

The `u` key in `configs/run_one.yaml` takes 5 normalized values in `[0,1]` --
`ua uw uhx uhy ut` -- mapping to `(a, w, hx, hy, t)` via
`src/geometry.py:u_to_geometry` (physical ranges in `configs/bounds.yaml`). `t`
(beam thickness) is a free design variable optimized by `run_loop.py`
(`src/optimizer.py:N_DIM = 5`). A legacy 4-value `u` is still accepted by
`run_one.py` and pins `t` to its minimum.

`run_one.py` prints per-stage progress (feasibility / optical / mechanical / score,
with elapsed time) to **stderr** as it runs, so `stdout` stays pure JSON even when
redirected to a file. Set `quiet: true` in the config to suppress the progress lines.

## Visualization (runs anywhere)
```bash
python scripts/plot_geometry.py   # unit cell: top view + cross-section
python scripts/plot_bands.py      # band diagram + gap (kind: optical-surrogate by default)
python scripts/plot_mode.py       # mode-field rendering (demo: true by default)
```
Like `run_one.py`/`run_loop.py`, each script takes a single `--config PATH` flag
defaulting to its own YAML in `configs/`, so the commands above run with no
external solver needed. Edit the YAML to change settings that used to be CLI
flags: `u`, `periods` in `configs/plot_geometry.yaml`; `u`, `kind`,
`target_nm`/`target_GHz` in `configs/plot_bands.yaml`; `demo` in
`configs/plot_mode.yaml`; and `npz`/`out` in the relevant config.

For real band/mode figures, the MPB and COMSOL backends expose `save_bands(...)`
and `export_mode_grid(...)`; instead of a `--npz` flag, set
`npz: path/to/file.npz` in `configs/plot_bands.yaml` / `configs/plot_mode.yaml`
(and for `plot_bands.py`, set `kind: optical` or `kind: mechanical`).

## Real physics (on the Mac)
COMSOL runs on your Mac, not in any sandbox. **Setup + diagnostic:** see
`docs/comsol_setup.md`, then `python scripts/check_comsol.py`.

Optical (fast, ideal geometry) -- MPB via conda:
```bash
conda create -n omc -c conda-forge pymeep pymeep-extras && conda activate omc
# set optical: mpb, mech: surrogate_stub in configs/run_one.yaml, then:
python scripts/run_one.py --config configs/run_one.yaml
```
Mechanical + fabrication-aware optical -- COMSOL 6.2:
1. Build `comsol/omc_unitcell.mph` once (see `comsol/README_template.md`).
2. Python driver: `pip install MPh`, start a COMSOL server, then set
   `optical: comsol`, `mech: comsol` in `configs/run_loop.yaml` and run
   `python scripts/run_loop.py --config configs/run_loop.yaml`.
   Or MATLAB LiveLink: `comsol/omc_mechanical_livelink.m`.

## Rib cell (`celltype: rib`)
A spine + centred-rib cell, the same geometry as the MATLAB `'rib'` device
(`../omc-comsol-chang/rib_README.md`). It is set up with `celltype: rib` in the
run config, and the hole cell is unaffected.

**Template.** It is built once in MATLAB with LiveLink
(`.mph` files are not in git):
```matlab
cd python-scripts/comsol
build_omc_unitcell_rib();      % writes comsol/omc_unitcell_rib.mph
```
- **Frame:** beam along z, thickness along x, modelled as the y >= 0 half.
- **Air:** a half-disk (half-cylinder along the beam, y >= 0) centred on the
  mid-depth of the beam, of radius `air_r` = distance to the farthest beam
  corner + `air_gap`; the curved face carries the scattering condition.
  `air_gap` (default 1.55 um, one wavelength) is a global parameter of the
  model: change it in the GUI or with `-p air_gap=2.5[um]` in
  `scripts/run_comsol_batch.sh`. Templates built before this have an air box
  (no parameter `air_r`) whose faces sit inside the evanescent tail of the
  guided modes; the COMSOL optical paths refuse them (delete and rebuild).
  Mechanics does not use the air and is unaffected.
- **Material:** anisotropic cubic diamond (C11/C12/C44). `rib.rxtal_deg`
  rotates the crystal about the slab normal (`P.rxtal` of
  `../omc-comsol-chang`; template parameter `rxt`): `0` (default) puts the
  beam along <100>, `45` along <110>. Only multiples of 45 deg are allowed
  (others break the y mirror of the half model; 90 is the same as 0), and it
  acts on mechanics only. It is saved in every record as `rxtal_deg` and is
  part of the candidate id only when nonzero, so existing ids are unchanged.
  A template built before `rxt` is refused for 45 (delete and rebuild). The
  hole templates are isotropic.
- **Studies:** `mech sym`, `mech antisym`, `opt TE` and `opt TM`.

The MATLAB driver `comsol/omc_mechanical_livelink_rib.m` gives a reference band
structure. Its run of the template's default cell is saved in `comsol/rib_bands/`.

**Parameters.** The names are those of the template, and **`t` means something different here**:

| Name | Meaning | Searched range / fixed value |
|---|---|---|
| `a` | period along the beam | 300-1000 nm |
| `w` | rib tip-to-tip width (= 2L) | 50-3000 nm |
| `t` | **rib width** along the beam (NOT the thickness) | 50-1000 nm |
| `s` | spine width | 50-1000 nm |
| `th` | slab thickness | fixed 250 nm |
| `r1`, `r2` | root / tip fillet radii (must be > 0) | fixed 40 nm |

The search space is `u` in [0,1]^4 over `(a, w, t, s)`, set in the `rib:`
section of `configs/bounds.yaml`. Each parameter sits under `variables`
(searched) or `fixed`, so moving `th` to `variables` needs no code change.

**Before any solver,** a candidate must pass two checks. Failures are saved, not discarded:
1. **Geometric validity** (status `infeasible`, `GEOM_*` codes): `t < a`,
   `s < w`, `2*r2 < t`, `t + 2*r1 < a`, `r1 + r2 < (w - s)/2`.
2. **Fabricability** (status `unfabricable`, `FAB_*` codes in `fab_failures`):
   - minimum rib width, spine width, rib arm `(w - s)/2`, air gap `a - t` and
     thickness, with a 50 nm default;
   - etch and solid aspect ratios;
   - rib slenderness and overhang;
   - an optional maximum release width.

   The thresholds live in `bounds.yaml` (`rib.fabrication`). The ratio limits
   have no repo source yet; set them from the real process.

About 40 % of the box passes both checks with the defaults. For `a` > ~775 nm,
1550 nm lies above the zone-edge light line, so no guided optical gap is possible.

**Scoring.**
- **Mechanical:** `rib.mech_gap_mode` = `complete` (default; runs `mech sym`
  and `mech antisym` at every k), `sym` or `antisym`. Gaps are found by
  frequency interval (`bandgap.interval_gap`), not band index, because each
  study returns the modes nearest 8 GHz.
- **Optical:** `optical: comsol` solves `opt TE` on the rib template and keeps
  only guided modes, using the light-line tracking and pseudo-mode filter in
  `src/comsol_optics.py`. `optical: none` scores the mechanical gap only.
  There is no window on the gap centre (the old `rib.optical_rel_tol` is
  deprecated: accepted with a warning and ignored). Every guided gap that lies
  below the zone-edge light line, in the trusted range, and is wider than the
  resolvable step is a candidate (`optical_gap_candidates`), and the gap centre
  `f_c` enters through a Gaussian weight of width `rib.optical_tol_thz`
  (default 2 THz) about the target `f_t` = 193.4 THz (1550 nm):
  `optical_freq_penalty = exp(-((f_c - f_t)/optical_tol)^2)`,
  `optical_gap_effective = optical_gap * optical_freq_penalty`.
  The candidate with the largest effective gap is taken.
- **Score (version 2, with optics):**
  `S = min(G_o_eff, G_m) - lo*max(0, gmin_o - G_o_eff)^2
       - lm*max(0, gmin_m - G_m)^2 - lmf*((f_m - f_m,t)/f_m,t)^2`,
  with `G_o_eff = optical_gap_effective`. The weight is 1 on target, 1/e at
  2 THz and 0.02 at 4 THz from it, so an off-target optical gap scores like no
  gap. `optical_gap` in the record stays the physical gap. Records carry
  `score_version` (2 with optics, 1 without, where nothing changed) and
  `optical_tol_thz`; both are in the id, and `run_opt_comsol.py` seeds the
  optimizer only with prior results of the same version and tolerance.
  `optical: tidy3d` scores the same gap from Tidy3D FDTD (see "Optical engine:
  COMSOL vs Tidy3D"). `surrogate` and `mpb` are not available for the rib.
  With the default `xsect: rect` only thickness-even modes are scored.
- **Validate optics first.** Before any optimization with `optical: comsol`,
  solve `opt TE` and `opt TM` once at k = pi/a on the default cell and check
  that the modes near 193 THz are guided and confined to the rib:
  `scripts/run_comsol_batch.sh --study std2 --study std4` (see "Getting a
  COMSOL seat"). The run of 2026-11-23 is in `results/rib_optical_validation/`:
  both studies return modes below the light line (c/2a = 263 THz) on either
  side of 193 THz (TE 163.6, 185.2, 218.3 THz; TM 183.4, 206.3 THz) with 76-79 %
  of the electric energy in the beam, but with Q of only 200-340. A guided
  mode at X should not radiate, so that loss is probably the scattering
  boundary sitting inside the evanescent tail (air box 2.5*th in x); look at
  the fields in `solved_std2.mph` / `solved_std4.mph` before trusting the gap.
  That run used the OLD air box: repeat it after rebuilding the template with
  the air half-disk (the old file can be kept, renamed, for the comparison:
  `--template ... --allow-old-air`). TODO(confirm): Q of the guided modes with
  the half-disk, and that it grows with `air_gap` (e.g. `-p air_gap=2.5[um]`).

**Running.**
```bash
python scripts/run_one.py      --config configs/run_one_rib.yaml        # one cell (optical: none)
python scripts/run_loop.py     --config configs/run_loop_rib.yaml       # ask/tell loop
python scripts/run_opt_comsol.py --config configs/run_opt_comsol_rib.yaml
#   run_opt_comsol.py dispatches celltype: rib to scripts/run_opt_comsol_rib.py
#   (loop, JSON, progress figure with rib drawing, characterization)
```
For a dry run without COMSOL, set `mech: surrogate_stub` and `optical: none`.

**Results** go to a separate database, `results/runs_rib.sqlite3` (override with
`OMC_DB_PATH_RIB`), so the hole database `results/runs.sqlite3` is never touched.
Rib records carry `celltype: "rib"`, plus `mech_gaps`, `mech_gap_mode`,
`fab_failures` and the other rib fields. Rib and hole results are not directly
comparable, because one uses an anisotropic material and the other an isotropic one.

## Getting a COMSOL seat
The licence server (`27005@license.seas.harvard.edu`) has 3 `WAVEOPTICS` seats
and they are often all taken. Rule for every tool here: **a seat is never held
while idle.** Nothing is checked out while waiting, and a seat is given back as
soon as the work that needed it ends.

`scripts/check_comsol_license.sh` asks the server what is free without checking
anything out. `-p` prints `FEATURE TOTAL USED FREE STATUS` lines for scripts.

**1. Batch solve on the batch licence pool:** `scripts/run_comsol_batch.sh`
```bash
scripts/run_comsol_batch.sh                              # rib template, opt TE at X
scripts/run_comsol_batch.sh --study std2 --study std4    # opt TE + opt TM (the validation solve)
scripts/run_comsol_batch.sh --study "opt TM" -p a=600[nm] -p t=400[nm] --outdir results/my_run
scripts/run_comsol_batch.sh --check-license              # licences the template needs; solves nothing
scripts/run_comsol_batch.sh --help
```
It runs one study of an `.mph` with `comsolbatch` (no GUI, no MPh, no LiveLink),
saves the solved model and writes `eig_<study>.csv` and `summary.json`
(Re/Im of every eigenfrequency, Q, below/above the light line, share of the
electric energy in the beam) to a new `results/batch_<time stamp>/` folder.
Defaults: `comsol/omc_unitcell_rib.mph`, study `std2` (`opt TE`), `kF = pi/a`.
Study tags of the rib template: `std1` mech sym, `std3` mech antisym, `std2`
opt TE, `std4` opt TM. Paths come from `COMSOL_ROOT` (or `--comsol-root`), so
the same command works on a machine with another install location.
For an optical study it refuses (exit 2) a rib template with the old air box
(no parameter `air_r`; see "Rib cell"), unless `--allow-old-air` is given, e.g.
to compare the Q of the old box and the half-disk deliberately. The air region
used is logged and saved as `air_region` in `summary.json`.

Which licences it takes (measured 2026-11-23, COMSOL 6.3, with `--watch-license`):

| Step | Time | Features checked out |
|---|---|---|
| solve (`-usebatchlic`, the default) | ~3 min | `COMSOLBATCH`, `WAVEOPTICSBATCH`, `MEMSBATCH`, `CADIMPORTBATCH`, `CLUSTERNODE` (+ `COMSOLUSER`, `CADIMPORTUSER`) |
| extract the eigenfrequencies | ~45 s, seats for ~30 s | `COMSOL`, `WAVEOPTICS`, `MEMS`, `CADIMPORT` (+ `COMSOLUSER`, `CADIMPORTUSER`) |

- The solve uses **no interactive seat**, so it runs while `WAVEOPTICS` is 3/3.
  It is COMSOL's `-usebatchlic` option that selects the batch pool, not batch
  mode itself (`--license interactive` leaves it out).
- The extract step **does need an interactive `WAVEOPTICS` and `MEMS` seat**,
  for a few seconds: it runs a small Java class, and COMSOL does not allow batch
  licences for class files. If no seat is free, the solved `.mph` is kept and
  the numbers can be fetched later with `--extract-only --outdir <folder>`.
- `MEMS[BATCH]` (2 seats) and `CADIMPORT[BATCH]` are taken even for an optical
  study, because the template lists them as its licences.
- All seats were back within 3 s of each process ending.

**2. Wait for a seat, then run:** `scripts/wait_for_comsol_seat.sh`
```bash
scripts/wait_for_comsol_seat.sh -r 3 -- python scripts/run_one.py --config configs/run_one_rib.yaml
scripts/wait_for_comsol_seat.sh --batch -r 2 -- scripts/run_comsol_batch.sh --study std2
scripts/wait_for_comsol_seat.sh -f WAVEOPTICS -f MEMS --all -r 3 -- \
    scripts/run_comsol_batch.sh --extract-only --outdir results/my_run --study std2
```
It polls `check_comsol_license.sh -p` (query only) until a watched feature has
a free seat, then runs the command and returns its exit status.
- `-f FEATURE` (repeatable; default `WAVEOPTICS` and `RF`, either of which can
  run ewfd), `--batch` (`WAVEOPTICSBATCH`, `RFBATCH`), `--any` (default) or
  `--all`.
- `-i` poll interval (60 s, minimum 30 s), `-t` timeout (exit 124).
- `-r N`: if the command fails with a licence error (someone took the seat
  between the check and the checkout), wait again and retry, at most N times.
- Ctrl-C stops the wait at once (exit 130).

**3. Several designs on one seat:** `scripts/run_batch_rib.py`
```bash
python scripts/run_batch_rib.py --config configs/run_batch_rib.yaml --dry-run   # list the queue
python scripts/run_batch_rib.py --config configs/run_batch_rib.yaml
scripts/wait_for_comsol_seat.sh -r 3 -- python scripts/run_batch_rib.py --config configs/run_batch_rib.yaml
```
Separate `run_one.py` calls each start COMSOL and give the seat back at exit, so
someone else can take it in between. `run_batch_rib.py` evaluates a queue (u
vectors in the YAML, a CSV of u vectors, and/or `run_one`-style configs) in one
process and one MPh session. `run_loop.py` already shares one session over its
iterations (`src/comsol_client.py` caches it); this is its counterpart for an
explicit list of designs.
- Every item is saved with `save_result(rec, db_path(celltype))`, exactly as
  `run_one.py` does, failed and infeasible ones included.
- COMSOL starts at the first item that needs a solver and is released
  (`comsol_client.release_client()`) when no queued item needs it any more, on
  an error, and on Ctrl-C. There is no keep-alive between or after items.
- A licence error stops the queue (exit 3) instead of failing every remaining
  item; `stop_on_license_error: false` carries on.

## Optical engine: COMSOL vs Tidy3D
The optical solver is chosen with the one config key `optical:`; mechanics is
always `mech: comsol | surrogate_stub`.

| `optical:` | hole | rib | `xsect` | notes |
|---|---|---|---|---|
| `surrogate`, `mpb` | yes | no | rect only | pre-screen |
| `comsol` | yes | yes | rect, wedge | same template as mechanics |
| `tidy3d` | yes | yes | rect, wedge | FDTD on the Tidy3D cloud, **costs FlexCredits** |
| `none` | no | yes | rect, wedge | mechanics only |

`optical: tidy3d` (`src/optical_tidy3d.py`) drives `diamond_bandStructure_sim.py`
(needs `pip install tidy3d` and an API key; that file stays a stand-alone CLI,
with the new option `--xsect`). Examples: `configs/run_one_tidy3d.yaml`,
`configs/run_one_rib_tidy3d.yaml`, `configs/run_opt_comsol_rib_tidy3d.yaml`.
It is available in `run_one.py`, `run_loop.py`, `run_batch_rib.py` and the rib
loop of `run_opt_comsol.py`; the hole loop of `run_opt_comsol.py` is COMSOL-only.

**Refractive index.** Diamond is n = 2.386 (telecom) for every optical engine,
from `configs/materials.yaml` (also `diamond_telecom` of the Tidy3D module and
`P.nbeam` of `../omc-comsol-chang`). The COMSOL templates still store
relpermittivity 5.76 (n = 2.40): the Python drivers overwrite it before every
optical solve (`comsol_client.set_diamond_index`), so the binary hole templates
need no rebuild. Mechanics does not use the index.

What still reads the number stored in the `.mph`: `comsol batch`
(`scripts/comsol_batch.py` can only set global parameters and refuses a
template without the parameter `n_diamond` unless `--keep-template-index` is
given), and any session that opens the template directly (the COMSOL GUI,
`mphopen` in MATLAB). To make the file itself carry the index:
- Rib templates: rebuild with `build_omc_unitcell_rib.m` (it defines
  `n_diamond = 2.386` and `relpermittivity = n_diamond^2`). Delete the old
  `.mph` first; the builder does not overwrite.
- Hole templates (`omc_unitcell.mph`, `omc_unitcell_iso.mph`; made by hand,
  no builder): patch them in place, in MATLAB with LiveLink:
  ```matlab
  cd python-scripts/comsol
  patch_hole_templates_index('-selftest')   % optional, no COMSOL needed
  report = patch_hole_templates_index();    % both hole templates, n = 2.386
  ```
  It adds the same parameter and sets the diamond material's relative
  permittivity to `n_diamond^2`, nothing else: no solve, and the geometry,
  mesh, studies, solver sequences and datasets are kept. Each file is backed
  up to `<name>.bak.mph` first, a file that is already patched is skipped, and
  the saved file is reopened and checked. A solution stored in the template
  is not recomputed.
  By hand in the GUI instead: Global Definitions > Parameters 1: add
  `n_diamond`, expression `2.386`; Component 1 > Materials > Diamond (100) >
  Material Contents > Relative permittivity: `n_diamond^2`; File > Save.
Once a template has `n_diamond`, the Python drivers and `comsol_batch.py` set
that parameter (from `configs/materials.yaml`) instead of the material
property.

**Cross-section (`xsect: rect | wedge`, default `rect`).**
- `rect`: the slab only (vertical walls, flat underside).
- `wedge`: the slab plus the triangular underside of an angled etch, as in
  `comsol/omc_unitcell_iso.mph`: under a strip of width b, an isosceles triangle
  on the bottom face of the slab with base angles 30 deg, apex b/(2*sqrt(3))
  below the slab.
  - Hole cell: one prism along the beam (b = w), the hole cut through both.
  - Rib cell: the "roof" of the FILLETED outline. Under a point p the wedge is
    `d(p)*tan(30 deg)` deep, with `d(p)` the in-plane distance from p to the
    filleted outline (neighbouring cells included;
    `geometry.rib_roof_distance`). Straight parts of the spine (b = s) and of
    the rib (b = t) keep the triangle above; the tip is hipped; under each tip
    fillet (r2) the wedge is a cone with its apex below the arc centre; around
    each root fillet (r1), and over the whole rib/spine junction, it is a cone
    `d = rho - r1` about the arc centre, so the wedge is deepest under the
    junction centre (`wedge_depths` in the record: `spine`, `rib`, `junction`,
    `max`).
  - COMSOL builds the rib roof exactly from planes and cones. Tidy3D builds it
    as `wedge_layers` (default 12) stacked slabs of exactly offset outlines:
    exact in plan view, a staircase in depth with an error of at most
    (deepest point)/(2*wedge_layers), about 6 nm for the default cell (set
    `tidy3d: {params: {wedge_layers: N}}` for more).

| cell | `rect` template | `wedge` template |
|---|---|---|
| hole | `comsol/omc_unitcell.mph` | `comsol/omc_unitcell_iso.mph` |
| rib | `comsol/omc_unitcell_rib.mph` | `comsol/omc_unitcell_rib_wedge.mph` |

The same template serves optics and mechanics, so `xsect` changes both. The
rib wedge template is not shipped; build it once in MATLAB with LiveLink:
```matlab
cd python-scripts/comsol
model = build_omc_unitcell_rib('', 'wedge');   % in memory first: check it
mphgeom(model)                                 % or mphlaunch(model)
build_omc_unitcell_rib('wedge');     % writes comsol/omc_unitcell_rib_wedge.mph
```
The wedge part of the builder has not been built yet. Check on the first
build: one beam domain whose footprint is the filleted outline; depth below
the slab `s/2*tan(30deg)` under the spine centreline between ribs,
`t/2*tan(30deg)` under the rib centreline mid-arm, the parameter `wd` under
the junction centre (the deepest point) and zero along the outline; the rib
prism lies under the rib (not rotated); root cones widen downwards, tip cones
point downwards; the builder's own checks (periodic faces, selections) and
the mesh pass. The list is repeated in the builder at the wedge block.
A missing template is an error that names this step. `surrogate` and `mpb`
model the rectangle only and refuse `xsect: wedge` before any solve. For the
hole cell, `wedge` adds the feasibility rule `t >= w/(2*sqrt(3))`
(`WEDGE_HOLE_BLIND`): the template's hole cutter is 3*t long, so in a thinner
slab the hole would not go through the wedge. The rib wedge adds no rule.
`xsect` is saved in every record and is part of the record id.
`run_opt_comsol.py` (hole loop) keeps its old behaviour when the key is absent:
the wedge template and all optical modes.

**Mid-thickness mirror symmetry (optics only).** A rectangular cell has a
mirror plane at mid-thickness, a wedge cell has not:
- `rect`: only the thickness-even (TE-like) modes are scored for `opt TE`
  (the odd ones for `opt TM`). Tidy3D uses a symmetry plane, symmetry tuple
  `(0, -1, +1)` for TE. The COMSOL templates model the full thickness, so the
  drivers read the parity of each solved mode from its field: the share of the
  in-plane field `|Ey|^2 + |Ez|^2` in `|E|^2`, integrated over the layer
  `|x| < thickness/8` around the mid-plane, is near 1 for an even mode and near
  0 for an odd one (threshold 0.5, `comsol_optics.py`). Twice as many
  eigenvalues are requested, and the value of every mode is saved in
  `optical_comsol_detail`.
- `wedge`: all modes are scored. Tidy3D uses `(0, -1, 0)` and spreads its
  dipoles over the slab thickness; COMSOL applies no filter.
- Mechanics always uses the full thickness and all modes.

**Cost safety (Tidy3D).** Nothing is uploaded unless the config has
`tidy3d: {allow_cloud: true, max_credits: <FlexCredits>}`.
- `mode: dry_run` builds and validates the simulations locally (plan.json,
  grid_report.txt), no upload. `mode: estimate` uploads drafts, records the
  maximum cost and deletes them. `mode: run` submits.
- `max_credits` is the budget of the whole process (a loop or queue stops
  cleanly, with everything saved, before the batch that would exceed it);
  `max_credits_total` is an optional budget over all runs
  (`results/tidy3d/credit_ledger.jsonl`). Both count Tidy3D's maximum estimates.
- Downloads go to `results/tidy3d/<celltype>/<record id>/` and are re-analyzed
  instead of paid again.
- `keep_hdf5: all | best | none` (default `all`) controls the downloaded `.hdf5`
  files, a few MB per cell. `best`: the loops (`run_opt_comsol.py` rib loop,
  `run_loop.py`, `run_batch_rib.py`) keep them for the best-scoring cell only
  and delete the others once their record is saved; `run_one.py` keeps its one
  cell. `none`: deleted as soon as a cell is analyzed. `params.json`,
  `plan.json`, `bands.csv` and `evaluation.json` stay for every cell, each
  deletion is logged (stderr and `hdf5_note.json` in the record folder), and
  only folders downloaded under `best`/`none` are ever pruned. A repeated cell
  (common on resume) is re-scored from the saved `evaluation.json`, so it is
  not uploaded or charged again (`result_source: cached_analysis` in the
  record). What is lost without the files: the resonance filters cannot be
  re-applied and `plot_tidy3d_bands.py` plots the saved resonances only.
- `k_mode: x_only` = 1 simulation per cell (gap at X; default for the hole);
  `k_mode: sweep` = `n_k` simulations (default for the rib, `n_k = rib.n_k_opt`).

**Progress (Tidy3D).** Everything the engine prints goes to **stderr**, never
to stdout (the JSON of `run_one.py` stays clean): Tidy3D's web progress (upload,
task status and progress bars, cost, download) and the engine's own steps
(`[2/4] optical: building ...`, `uploading`, `estimated maximum cost`,
`running`, `loading`, `analyzing resonances`), also in `mode: dry_run` and
`mode: estimate`. `quiet: true` (`run_one.py`, `run_batch_rib.py`) silences
both; errors remain. Loops without a `quiet` key show Tidy3D's progress.

**Scoring.** Tidy3D resonances go through the gap functions used with COMSOL:
`guided_interval_gap` (rib sweep) and `optical_gap_tracked` (hole, and any
x_only run). The module's own "two lowest guided resonances at X" gap is stored
in `optical_tidy3d_detail.x_point`. Records carry `optical_engine`,
`optical_index`, `xsect` and `optical_z_symmetry`; the id depends on engine,
index and cross-section, so such results never overwrite each other (nor the
records made before 2026-11-23 with n = 2.40).

**Looking at downloaded Tidy3D results.** `scripts/plot_tidy3d_bands.py`
re-analyzes the `.hdf5` files of a record offline (no cloud call; it blocks the
network while it runs) and plots the band diagram with the gap as the pipeline
scores it. It needs tidy3d, so run it with the `tidy` conda env:
```bash
python scripts/plot_tidy3d_bands.py cand_f85eb9955630       # a record id
python scripts/plot_tidy3d_bands.py --all                   # every record + index CSV
python scripts/plot_tidy3d_bands.py cand_f85eb9955630 --min-guided-q 1000 --light-line-margin 0.03
```
Inputs: record ids, `<celltype>/<record id>`, run folders or `.hdf5` files.
Output goes to `<record folder>/analysis/` (or `--outdir`):
`band_diagram.png`, `resonances.csv` (k, f, Q, amplitude, error, guided) and
`summary.json` (scored gap, weight, gap at X, counts, and whether it matches
the run's `evaluation.json`). `--minQ`, `--minamp`, `--maxerr` and the
frequency window can be overridden; `--min-guided-q` / `--light-line-margin`
additionally report the gap without low-Q or light-line "guided" peaks, which
harmonic inversion produces and which can bound or split a gap.

**How comparable are the two engines?**
- `rect`: same cell, same index, same thickness-even family. Left: fillets are
  8-point polygons in Tidy3D, PML vs scattering boundary, and FDTD resolution.
- `wedge`: same cell and index, all parities in both. The rib wedge is exact
  in COMSOL and a 12-layer staircase in Tidy3D. An FDTD run finds only the
  modes its Ey dipoles excite, so a mode COMSOL returns can be missing and the
  Tidy3D gap can come out larger. Confirm finalists in COMSOL.
- Neither the rib wedge template nor any Tidy3D simulation has been run yet.

## Layout
```
configs/   bounds.yaml (hole + rib sections)  targets.yaml  materials.yaml
           run_one.yaml  run_loop.yaml  run_opt_comsol.yaml  (per-script --config)
           run_one_rib.yaml  run_loop_rib.yaml  run_opt_comsol_rib.yaml
           run_batch_rib.yaml  (queue of designs for run_batch_rib.py)
           run_one_tidy3d.yaml  run_one_rib_tidy3d.yaml
           run_opt_comsol_rib_tidy3d.yaml  (optical: tidy3d examples)
           plot_bands.yaml  plot_geometry.yaml  plot_mode.yaml
src/       geometry (hole + rib, fabricability)  bandgap  comsol_client
           optical_{surrogate,mpb,comsol,tidy3d}  comsol_optics  acoustic_comsol
           objective  objective_rib  optimizer  database
scripts/   run_one.py  run_loop.py  run_opt_comsol.py  run_opt_comsol_rib.py
           run_batch_rib.py  (queue in one COMSOL session)
           check_comsol_license.sh  wait_for_comsol_seat.sh
           run_comsol_batch.sh + comsol_batch.py  (comsol batch, batch licences)
           plot_bands.py  plot_geometry.py  plot_mode.py
           plot_tidy3d_bands.py  (offline band diagrams from Tidy3D .hdf5)
comsol/    hole template recipe + MATLAB LiveLink driver;
           patch_hole_templates_index.m (index parameter for the hole templates)
           build_omc_unitcell_rib.m + omc_mechanical_livelink_rib.m (rib)
diamond_bandStructure_sim.py   stand-alone Tidy3D band-structure tool
tests/     test_pipeline.py  (numpy-only, hole + rib)
           test_optical_tidy3d.py  test_xsect_index.py  (no cloud, no COMSOL)
           test_batch_queue.py  (run_batch_rib.py + session release, no COMSOL)
results/   runs.sqlite3 (hole; + .jsonl fallback)  runs_rib.sqlite3 (rib)
docs/      strategy.md
```
