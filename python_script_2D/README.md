# python_script_2D -- 2D diamond phononic-crystal ("boomerang" cell) optimization

Python/MPh optimization pipeline for the 2D hexagonal-lattice "boomerang"
phononic-crystal unit cell, adapted from two sources:

- **Physics/geometry reference**: the MATLAB COMSOL pipeline in
  `../omc-comsol-chang/` -- specifically `test_Boomerang.m` (front-panel
  script) calling `solveBands.m` -> `runBands_2D.m` (2D Brillouin-zone sweep,
  Floquet boundary conditions, z-parity symmetry BCs) -> `buildBoomerangUnitCell.m`
  (the boomerang-hole geometry itself).
- **Software architecture reference**: `../python-scripts/` -- the 1D
  diamond nanobeam optimization pipeline this package's ask/tell loop,
  SQLite result store, and config-driven CLI scripts are patterned after.

## What's different from `../python-scripts/` (read this first)

| | `../python-scripts/` (1D nanobeam) | `python_script_2D` (this package) |
|---|---|---|
| Periodicity | 1D, along the beam axis (`kz`) | **True 2D**, hexagonal Brillouin zone (`kx,ky`), Gamma->M->K->Gamma |
| Symmetry family | in-plane y-even/y-odd (breathing mode, `fy`) | **z-parity** (evenz/oddz about the slab midplane) |
| Design variables | `a, w, hx, hy, t` (5) | `a, w, r, r1, r2, th` (6) -- boomerang leg geometry |
| Mechanical backend | implemented (`acoustic_comsol.py`) | implemented (`acoustic_comsol_2d.py`) |
| Optical backend | implemented (MPB / COMSOL / surrogate) | **not yet implemented** -- see `docs/optical_2d_plan.md` |

If you're coming from `../python-scripts/`, the module names and CLI
structure will look familiar (`run_one_2d.py`/`run_loop_2d.py`,
`objective2d.py`, `optimizer.py`, `database.py`) -- several modules
(`bandgap.py`, `optimizer.py`, `database.py`, `comsol_client.py`) are ported
with little or no change, since that logic is genuinely backend/geometry
agnostic. Everything geometry- and physics-specific (`geometry2d.py`,
`acoustic_comsol_2d.py`, `objective2d.py`) was rewritten for the 2D case.

## Quickstart

```bash
pip install -r requirements.txt        # numpy, pyyaml, optuna

# 1. Numpy-only smoke test (no COMSOL/MPh needed):
python tests/test_pipeline_2d.py

# 2. Build the COMSOL template ONCE (see comsol/README_template_2d.md):
#    save it as comsol/omc2d_boomerang.mph

# 3. On the Mac, with a COMSOL server running (`comsol mphserver`):
pip install MPh
python scripts/run_one_2d.py           # evaluate configs/run_one_2d.yaml's `u`
python scripts/run_loop_2d.py          # closed-loop optimization
```

Both CLI scripts read all settings from YAML (`configs/run_one_2d.yaml` /
`configs/run_loop_2d.yaml`); pass `--config path/to/file.yaml` to override.

## Layout

```
comsol/    README_template_2d.md  -- .mph template build recipe (READ THE CAVEATS)
           omc2d_boomerang.mph    -- you build this once, in the COMSOL GUI
configs/   bounds_2d.yaml, materials.yaml, targets_2d.yaml,
           run_one_2d.yaml, run_loop_2d.yaml
src/       geometry2d.py, bandgap.py, acoustic_comsol_2d.py,
           optical_comsol_2d.py (stub), objective2d.py,
           optimizer.py, database.py, comsol_client.py
scripts/   run_one_2d.py, run_loop_2d.py
tests/     test_pipeline_2d.py     -- numpy-only, no COMSOL/MPh required
docs/      optical_2d_plan.md      -- future optical-bandgap implementation plan
```

## Status / known limitations

- **Mechanical**: implemented, but the `.mph` template (`comsol/
  omc2d_boomerang.mph`) has not been built yet -- do that first, per
  `comsol/README_template_2d.md`.
- **Feasibility check** (`src/geometry2d.check_feasibility`) uses a
  conservative heuristic for the boomerang hole's minimum remaining solid
  material -- not re-derived exactly from the Boolean geometry. See its
  docstring; validate visually (`mphgeom`) near the bounds' extremes.
- **Optical**: not implemented. `src/optical_comsol_2d.py` raises
  `NotImplementedError`; `require_opt: false` is the default everywhere.
  See `docs/optical_2d_plan.md` for the design writeup.
- No plotting/characterization scripts yet (the 1D project has
  `plot_bands.py`, `characterize_best.py`, etc. in `../python-scripts/
  scripts/`) -- worth porting once the mechanical loop has produced results
  worth visualizing.
