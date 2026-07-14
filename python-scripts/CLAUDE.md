# Project: 1D diamond OMC mirror-cell bandgap optimization

## Goal
Optimize a 1D diamond photonic/phononic crystal **unit cell** (rectangular
cross-section nanobeam + one elliptical hole) to obtain simultaneously:
- optical normalized bandgap (TE-like / y-even) **>= 20%**, mid-gap near **1550 nm**
- mechanical normalized bandgap **>= 20%** near **8 GHz**
  - A *symmetry-restricted* (quasi) gap for one mode family (even-even
    breathing-like) is acceptable -- a complete phononic gap is NOT required.

Material: single-crystal diamond (n ~ 2.40; cubic elastic constants in
`configs/materials.yaml`). For now: unit cells only (no cavity taper, no Q,
no g0). Find cells with both gaps first; seed cavity design later.

## Backends
- **Optical, fast pre-screen:** MPB (`src/optical_mpb.py`). Ideal geometry.
- **Optical, fabrication-aware:** COMSOL Wave Optics (`src/optical_comsol.py`).
- **Optical, dependency-free smoke test:** numpy TMM surrogate
  (`src/optical_surrogate.py`). NOT physically accurate; loop debugging only.
- **Mechanical:** COMSOL Solid Mechanics eigenfrequency + Floquet
  (`src/acoustic_comsol.py`), driven from Python via `MPh`, on the user's Mac.

## Important files
- `src/geometry.py`   normalized u in [0,1]^4 -> physical geometry + feasibility.
- `src/bandgap.py`    extract complete or symmetry-restricted gaps from bands.
- `src/optical_*.py`  optical band-structure backends.
- `src/acoustic_comsol.py`  mechanical band-structure backend (COMSOL/MPh).
- `src/objective.py`  geometry -> (G_o, G_m, score) with penalties + caching.
- `src/optimizer.py`  ask/tell: numpy fallback -> Optuna -> BoTorch.
- `src/database.py`   SQLite (+ JSONL fallback). Never lose a simulation.
- `scripts/run_one.py`, `scripts/run_loop.py`  CLI entry points.
- `comsol/`           .mph template + MATLAB LiveLink alternative + recipe.

## Rules
- Preserve the result-record schema (back-compatible additions only).
- Save EVERY candidate, including infeasible/failed (status field).
- Optimizer variables are normalized u in [0,1]; map in `geometry.py`.
- Enforce minimum feature size BEFORE launching any solver.
- Do NOT change solver physics (materials, BCs, parities) without explicit ask.
- Prefer small, testable changes. Run `tests/test_pipeline.py` (numpy-only)
  before long simulations.
- Optical pre-screen with MPB/surrogate; confirm finalists in COMSOL.
