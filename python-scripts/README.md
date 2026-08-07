# 1D Diamond OMC Unit-Cell Bandgap Optimization

Closed-loop optimization of a 1D diamond optomechanical-crystal **unit cell**
(rectangular cross-section nanobeam + elliptical hole) for a simultaneous
optical bandgap near **1550 nm** and a mechanical bandgap near **8 GHz**.

See `docs/strategy.md` for the full strategy, field review, and physics notes.

## Quick start (numpy-only smoke test, runs anywhere)
```bash
pip install -r requirements.txt           # numpy/pyyaml/pandas/matplotlib (+optuna)
python tests/test_pipeline.py             # all tests pass with numpy only
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

## Layout
```
configs/   bounds.yaml  targets.yaml  materials.yaml
           run_one.yaml  run_loop.yaml  run_opt_comsol.yaml  (per-script --config)
           plot_bands.yaml  plot_geometry.yaml  plot_mode.yaml
src/       geometry  bandgap  optical_{surrogate,mpb,comsol}  acoustic_comsol
           objective  optimizer  database
scripts/   run_one.py  run_loop.py  run_opt_comsol.py
           plot_bands.py  plot_geometry.py  plot_mode.py
comsol/    template recipe + MATLAB LiveLink driver
tests/     test_pipeline.py  (numpy-only)
results/   runs.sqlite (+ .jsonl fallback)
docs/      strategy.md
```
