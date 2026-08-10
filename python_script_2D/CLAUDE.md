# Project: 2D diamond phononic-crystal "boomerang" cell optimization

## Goal

Optimize a 2D diamond phononic-crystal unit cell -- a hexagonal lattice with
a "boomerang"-shaped air hole (3 fillet-rounded rectangular legs at 120 deg) --
for a mechanical bandgap. This is the Python/MPh counterpart of the MATLAB
pipeline in `omc-comsol-chang/` (see `test_Boomerang.m` ->
`buildBoomerangUnitCell.m` for geometry, `runBands_2D.m` for physics/BCs and
the Brillouin-zone sweep), reworked to run as a closed optimization loop the
way `../python-scripts/` does for the 1D nanobeam project.

Two things are genuinely different from `../python-scripts/` (read
`docs/optical_2d_plan.md` and `comsol/README_template_2d.md` before touching
solver physics):
- **True 2D Bloch periodicity**: the wavevector is `(kx, ky)`, swept around
  the hexagonal Brillouin zone Gamma -> M -> K -> Gamma, not a single 1D `kz`.
- **Z-parity family split, not a breathing-mode fraction**: the slab has a
  z=0 mirror plane; mechanical modes split into z-even ("mech evenz") and
  z-odd ("mech oddz") families. A gap can be required across BOTH
  (`gap_mode: "complete"`) or within just one (`gap_mode: "symmetry"`) --
  see `configs/targets_2d.yaml`.

Optical (photonic) bandgap optimization is **planned but not implemented**
-- see `docs/optical_2d_plan.md`. `src/optical_comsol_2d.py` is a stub that
raises `NotImplementedError`; keep `require_opt: false` /
`optical_backend: none` in run configs until it's built.

## Backend

- **Mechanical:** COMSOL Solid Mechanics eigenfrequency + Floquet, driven
  from Python via `MPh` (`src/acoustic_comsol_2d.py`), on the user's Mac.
  Requires a hand-built COMSOL template -- see
  `comsol/README_template_2d.md` for the exact node-by-node recipe.
- **Optical:** not yet implemented (`src/optical_comsol_2d.py`).

## Important files

- `src/geometry2d.py`   normalized u in [0,1]^6 -> physical geometry + feasibility.
- `src/bandgap.py`      extract complete or symmetry-restricted gaps from bands
                        (ported unchanged from `../python-scripts/src/bandgap.py`
                        -- fully backend/dimension-agnostic).
- `src/acoustic_comsol_2d.py`  mechanical band-structure backend (COMSOL/MPh),
                        including the hexagonal-BZ `kx,ky` sweep (`bz_path`).
- `src/optical_comsol_2d.py`   stub; see `docs/optical_2d_plan.md`.
- `src/objective2d.py`  geometry -> (G_m, [future] G_o, score), with
                        feasibility gating and the z-parity gap_mode logic.
- `src/optimizer.py`, `src/database.py`, `src/comsol_client.py`  ported
  near-unchanged from `../python-scripts/src/` (generic ask/tell optimizer,
  SQLite result store, shared MPh client singleton).
- `scripts/run_one_2d.py`, `scripts/run_loop_2d.py`  CLI entry points.
- `comsol/README_template_2d.md`  .mph template recipe (READ THE CAVEATS
  about the rhombic cell's periodic-BC face pairs before building it).

## Rules

(Same spirit as `../python-scripts/CLAUDE.md`; repeated here since this is a
separate package.)

- Preserve the result-record schema (back-compatible additions only).
- Save EVERY candidate, including infeasible/failed (status field).
- Optimizer variables are normalized u in [0,1]; map in `geometry2d.py`.
- Enforce minimum feature size BEFORE launching any solver
  (`geometry2d.check_feasibility` -- but see its docstring caveat: the
  boomerang hole's exact minimum-solid-material check is a conservative
  heuristic, not re-derived exactly from the Boolean geometry).
- Do NOT change solver physics (materials, BCs, z-parity convention) without
  explicit ask -- and cross-check against `omc-comsol-chang/
  buildBoomerangUnitCell.m` / `runBands_2D.m` first if in doubt.
- Prefer small, testable changes. Run `tests/test_pipeline_2d.py`
  (numpy-only, no COMSOL/MPh needed) before long simulations:
  `conda run -n omc python3 tests/test_pipeline_2d.py` (or your own
  environment with numpy + pyyaml installed).
- Optical is future work -- don't wire it into scored runs until
  `docs/optical_2d_plan.md`'s validation plan has been carried out.
