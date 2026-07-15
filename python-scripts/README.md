# 1D Diamond OMC Unit-Cell Bandgap Optimization

Closed-loop optimization of a 1D diamond optomechanical-crystal **unit cell**
(rectangular cross-section nanobeam + elliptical hole) for a simultaneous
optical bandgap near **1550 nm** and a mechanical bandgap near **8 GHz**.

See `docs/strategy.md` for the full strategy, field review, and physics notes.

## Quick start (numpy-only smoke test, runs anywhere)
```bash
pip install -r requirements.txt           # numpy/pyyaml/pandas/matplotlib (+optuna)
python tests/test_pipeline.py             # all tests pass with numpy only
python scripts/run_one.py  --u 0.5 0.6 0.5 0.6 0.5 --optical surrogate --mech surrogate_stub
python scripts/run_loop.py --n-init 40 --n-iter 80 --optical surrogate --mech surrogate_stub
```
The surrogate is a coarse numpy TMM model for loop debugging only -- not physically accurate.

`--u` takes 5 normalized values in `[0,1]` -- `ua uw uhx uhy ut` -- mapping to
`(a, w, hx, hy, t)` via `src/geometry.py:u_to_geometry` (physical ranges in
`configs/bounds.yaml`). `t` (beam thickness) is a free design variable
optimized by `run_loop.py` (`src/optimizer.py:N_DIM = 5`). A legacy 4-value
`--u` is still accepted by `run_one.py` and pins `t` to its minimum.

`run_one.py` prints per-stage progress (feasibility / optical / mechanical / score,
with elapsed time) to **stderr** as it runs, so `stdout` stays pure JSON even when
redirected to a file. Pass `--quiet` to suppress the progress lines.

## Visualization (runs anywhere)
```bash
python scripts/plot_geometry.py --u 0.5 0.6 0.5 0.6 --periods 5   # unit cell: top view + cross-section
python scripts/plot_bands.py    --u 0.5 0.6 0.5 0.6 --kind optical-surrogate   # band diagram + gap
python scripts/plot_mode.py     --demo                            # mode-field rendering (synthetic)
```
For real band/mode figures, the MPB and COMSOL backends expose `save_bands(...)`
and `export_mode_grid(...)`; feed the resulting `.npz` to `plot_bands.py --npz`
/ `plot_mode.py --npz`.

## Real physics (on the Mac)
COMSOL runs on your Mac, not in any sandbox. **Setup + diagnostic:** see
`docs/comsol_setup.md`, then `python scripts/check_comsol.py`.

Optical (fast, ideal geometry) -- MPB via conda:
```bash
conda create -n omc -c conda-forge pymeep pymeep-extras && conda activate omc
python scripts/run_one.py --u 0.5 0.6 0.5 0.6 0.5 --optical mpb --mech surrogate_stub
```
Mechanical + fabrication-aware optical -- COMSOL 6.2:
1. Build `comsol/omc_unitcell.mph` once (see `comsol/README_template.md`).
2. Python driver: `pip install MPh`, start a COMSOL server, then
   `python scripts/run_loop.py --optical comsol --mech comsol`.
   Or MATLAB LiveLink: `comsol/omc_mechanical_livelink.m`.

## Layout
```
configs/   bounds.yaml  targets.yaml  materials.yaml
src/       geometry  bandgap  optical_{surrogate,mpb,comsol}  acoustic_comsol
           objective  optimizer  database
scripts/   run_one.py  run_loop.py
comsol/    template recipe + MATLAB LiveLink driver
tests/     test_pipeline.py  (numpy-only)
results/   runs.sqlite (+ .jsonl fallback)
docs/      strategy.md
```
