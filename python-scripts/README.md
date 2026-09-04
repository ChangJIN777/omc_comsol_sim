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
           objective  optimizer  database  cli_progress  comsol_client
scripts/   see "Scripts reference" below
comsol/    template recipe + MATLAB LiveLink driver
tests/     test_pipeline.py  (numpy-only)
results/   runs.sqlite3 (+ .jsonl fallback)
docs/      strategy.md  comsol_setup.md  run_opt_comsol.md
```

## Scripts reference

Every script in `scripts/`. **Mac/COMSOL** = requires COMSOL + MPh (runs only on
the Mac with LiveLink). **anywhere** = pure numpy/stdlib, runs on any machine.
Scripts that overwrite result files in place are flagged ⚠.

### Optimization drivers
| Script | Purpose | Runs |
|---|---|---|
| `run_one.py` | Evaluate ONE candidate `--u` and save the record; per-stage progress to stderr, pure JSON to stdout. Flags: `--u`, `--optical {surrogate,mpb,comsol}`, `--mech {comsol,surrogate_stub}`, `--quiet`. | anywhere (surrogate); Mac/COMSOL if `--optical mpb`/`--mech comsol` |
| `run_loop.py` | Closed ask→evaluate→save→tell loop via `src/`. Flags: `--n-init`, `--n-iter`, `--optical`, `--mech`, `--backend {auto,random,optuna,botorch}`, `--quiet`. | anywhere (surrogate); Mac/COMSOL with `--mech comsol` |
| `run_opt_comsol.py` | Production COMSOL co-optimization driver (optical TE gap near 1550 nm + breathing mechanical gap 5–10 GHz); Optuna TPE or Halton fallback; spawns `characterize_best.py` at the end. See `docs/run_opt_comsol.md`. | Mac/COMSOL |

### Monitoring
| Script | Purpose | Runs |
|---|---|---|
| `watch_progress.py` | Live polling dashboard (progress bar, rate/ETA, best-so-far, top-5) over a results JSON. Flags: `--file`, `--n-iter`, `--interval`, `--once`. | anywhere (read-only) |
| `replot.py` | Regenerate the progress figure from the saved JSON without re-running COMSOL. Flags: `--json`, `--fig`. | anywhere (read-only) |

### COMSOL setup / diagnostics
| Script | Purpose | Runs |
|---|---|---|
| `check_comsol.py` | Step-by-step COMSOL↔Python(MPh) connectivity diagnostic (arch, imports, `mph.start()` + license, trivial eval, optional template load). | Mac/COMSOL |
| `inspect_comsol.py` | Physical sanity check: Γ→X sweep + zone-edge mode shapes rendered as a 4-panel figure (optical bands+mode, mech bands colored by fy + breathing mode). Flags: `--u`, `--n-k`, `--study-*`, `--out`. | Mac/COMSOL |

### Evaluation / parameter sweeps
All COMSOL-driven, reuse `evaluate()` from `run_opt_comsol.py`.
| Script | Purpose | Runs |
|---|---|---|
| `eval_points.py` | Evaluate a hardcoded list of physics-seeded (large-hx / 1550 nm-tuned) u-vectors; appends to `results/opt_results.json`. | Mac/COMSOL |
| `eval_t_sweep.py` | Sweep beam thickness t (220→420 nm) at fixed a/w/hx/hy to map G_o vs t; caches done points → `results/t_sweep.json`. | Mac/COMSOL |
| `eval_zone_edge.py` | Evaluate a grid around the 1550 nm zone-edge sweet spot; dedups, appends to `results/opt_results.json`. | Mac/COMSOL |
| `eval_large_hx.py` | Evaluate large-hx geometries at a≈840–900 nm (strong optical modulation); dedups, appends to `results/opt_results.json`. | Mac/COMSOL |

### Plotting / visualization
| Script | Purpose | Runs |
|---|---|---|
| `plot_geometry.py` | Draw the unit cell: top view + y–x cross-section. Flags: `--u` (4 values), `--periods`, `--out`. | anywhere |
| `plot_bands.py` | Γ→X band diagram + gap shading; from a saved `.npz` (`--kind {optical,mechanical}`) or the numpy surrogate (`--u --kind optical-surrogate`). | anywhere |
| `plot_mode.py` | Render a solved mode on a 2D cut plane from a grid `.npz`; `--demo` synthesizes a breathing mode with no COMSOL. | anywhere |
| `plot_gap_vs_tw.py` | Scatter G_o and G_m vs (t, w) from existing records (marks fy-verified and ≥target points). Flags: `--files`. | anywhere (read-only) |
| `plot_tpe_landscape.py` | Corner plot of the TPE "good region" 5-D KDE density + mode-finding of density islands; writes a figure and a new `_islands.json`. Flags: `FILE...`, `--good-frac`, `--n-starts`. | anywhere |

### Reprocessing / maintenance
| Script | Purpose | Runs |
|---|---|---|
| `reanalyze_opt_gap.py` | Re-run optical-gap detection (zone-edge + light-line filter) on stored `_freqs_o`, recompute score, print diagnostics. Read-only — emits `<name>_reanalyzed.json`. Flag: `--json`. | anywhere |
| `rescore_results.py` ⚠ | Re-score existing records with the current `score_result()` (optional require/g-min overrides); stashes original as `score_prev`. **Mutates `FILE...` in place.** | anywhere |
| `migrate_t_bounds.py` ⚠ | One-time re-map of stored `u[4]` + rescore after changing `t_max` in bounds.yaml. **Mutates `--in` in place by default.** Flags: `--t-max-old`, `--t-max-new`, `--in`, `--out`. | anywhere |

### Characterization
| Script | Purpose | Runs |
|---|---|---|
| `characterize_best.py` | End-of-run characterization of the best stored result (band sweeps, mode fields, figures); launched as a subprocess by `run_opt_comsol.py`, also usable standalone. Flags: `--mph`, `--json`, `--out-dir`, `--no-cache`. | Mac/COMSOL |
