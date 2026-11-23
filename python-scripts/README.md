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
- **Frame:** beam along z, thickness along x, modelled as the y >= 0 half with an air box.
- **Material:** anisotropic cubic diamond (C11/C12/C44), unrotated, so the beam
  runs along <100>. The hole templates are isotropic.
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
  `surrogate` and `mpb` are not available for the rib.
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

## Layout
```
configs/   bounds.yaml (hole + rib sections)  targets.yaml  materials.yaml
           run_one.yaml  run_loop.yaml  run_opt_comsol.yaml  (per-script --config)
           run_one_rib.yaml  run_loop_rib.yaml  run_opt_comsol_rib.yaml
           run_batch_rib.yaml  (queue of designs for run_batch_rib.py)
           plot_bands.yaml  plot_geometry.yaml  plot_mode.yaml
src/       geometry (hole + rib, fabricability)  bandgap  comsol_client
           optical_{surrogate,mpb,comsol}  comsol_optics  acoustic_comsol
           objective  objective_rib  optimizer  database
scripts/   run_one.py  run_loop.py  run_opt_comsol.py  run_opt_comsol_rib.py
           run_batch_rib.py  (queue in one COMSOL session)
           check_comsol_license.sh  wait_for_comsol_seat.sh
           run_comsol_batch.sh + comsol_batch.py  (comsol batch, batch licences)
           plot_bands.py  plot_geometry.py  plot_mode.py
comsol/    hole template recipe + MATLAB LiveLink driver;
           build_omc_unitcell_rib.m + omc_mechanical_livelink_rib.m (rib)
tests/     test_pipeline.py  (numpy-only, hole + rib)
           test_batch_queue.py  (run_batch_rib.py + session release, no COMSOL)
results/   runs.sqlite3 (hole; + .jsonl fallback)  runs_rib.sqlite3 (rib)
docs/      strategy.md
```
