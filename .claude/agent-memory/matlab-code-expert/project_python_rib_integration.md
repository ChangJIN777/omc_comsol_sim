---
name: project-python-rib-integration
description: 2026-11-23 proposal to add celltype rib to python-scripts (run_one/run_loop → objective_rib.py); rib t = rib width, th = thickness
metadata:
  type: project
---

On 2026-11-23 a proposal was delivered to add `celltype: rib` to python-scripts. Rib runs go via `run_one.py` / `run_loop.py` → new `src/objective_rib.py`, not `run_opt_comsol.py`. In rib records and `RibGeometry`, `t` = rib width and `th` = thickness; in hole records `t` = thickness. Rib gaps use the new `bandgap.interval_gap` (a frequency-interval port of the MATLAB gapReportUnion). The objective sector is configurable (`rib.mech_gap_mode`: sym/antisym/complete), and optics defaults to `optical: none` for rib.

**Why:** the user named t/th deliberately, and each template study returns the modes nearest the 8 GHz shift, so band-index gaps are unreliable.

**How to apply:** check whether it was implemented before re-proposing. Never translate rib `t` into thickness. Still open: search bounds, fixed th/r1/r2, sector, the optics approach, porting rib into run_opt_comsol.py, and shared vs separate DB.

Related: [[project-python-mph-templates]], [[reference-python-envs]].

**User decisions 2026-11-23:**
- Bounds a [300,1000], w [50,3000], t [50,1000], s [50,1000] nm, plus a separate fabricability check before any solver.
- th = 250 nm and r1 = r2 = 40 nm, fixed.
- Score the complete mechanical gap.
- Add rib optics, with light-line tracking and pseudo-mode filtering moved to a shared module.
- Add rib to run_opt_comsol.py.
- Use a separate sqlite DB for rib (no hole schema change).
A revised proposal was requested.
- As of 2026-11-23, local In_dev (merge b963194 with the rib template + driver) had diverged from origin/In_dev (9c9054e). Don't push without the user. Reference rib bands are in python-scripts/comsol/rib_bands/rib_bands_20261123_133400.mat. The localhost:2036 COMSOL server was taken over mid-run by another client once, so always ask before MPh/LiveLink runs.
- 2026-11-23: the revised proposal was approved and handed to an implementer. It adds a new src/comsol_optics.py (shared optics helpers moved verbatim from run_opt_comsol.py), src/objective_rib.py, scripts/run_opt_comsol_rib.py (reached via a celltype dispatch) and results/runs_rib.sqlite3. Before any rib optimization with optics, the optical validation solve (opt TE/TM at k = π/a) is required.
- 2026-11-23: IMPLEMENTED (uncommitted). No-COMSOL checks passed: test_pipeline, the verbatim-move check, the hole run_one JSON identical to HEAD, and a dry rib loop. The gated MPh steps B0–B4 haven't run yet. Known gap: the progress figure draws the rib only for 'success' records.
