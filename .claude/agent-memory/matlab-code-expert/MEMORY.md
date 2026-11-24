# Agent Memory Index

- [User role](user_role.md) — diamond OMC photonics researcher, high MATLAB/COMSOL fluency, wants decidable plans before implementation.
- [Preserve legacy geometry behavior](feedback_preserve_legacy_geometry_behavior.md) — reproduce builder quirks bit-for-bit, warn about defects instead of fixing them.
- [Scope discipline](feedback_scope_discipline.md) — one file only, no unrequested tests or feature gates; how to verify with headless MATLAB.
- [geom.runAll never finalizes](project_geom_runall_not_finalized.md) — latent defect: builders ending in runAll return empty xEnd/yEnd/bndSel indices.
- [Rib celltype](project_visible_rib_template.md) — rib = visible design merged with LN rib; template provenance, d/fillet/L-precedence rules.
- [Band-tracking proposal](project_band_tracking_proposal.md) — 2026-11-22 proposal for tracked-band fitness in optimize_rib_unitCell; decisions open, check if implemented.
- [Bayesopt unit-cell proposal](project_bayesopt_unitcell_proposal.md) — 2026-11-22 shared driver + hole/rib adapters; decisions open, check if implemented.
- [findGaps touching bands](project_findgaps_touching_bands.md) — strict test at findGaps.m:31 reports spurious gaps; surrogates on ds.sym must overlap bands.
- [Python .mph templates](project_python_mph_templates.md) — omc_unitcell*.mph: 4 studies, half model, isotropic; dset1 is optical in _iso; rib livelink proposal 2026-11-23.
- [Python rib integration](project_python_rib_integration.md) — 2026-11-23 proposal: celltype rib via objective_rib.py; rib t = rib width, th = thickness.
- [Python envs](reference_python_envs.md) — WSL python3 lacks numpy; conda env comsol_sim has MPh/optuna, no pytest; mph.start() = separate COMSOL session.
- [Fab check separate](feedback_fab_check_separate.md) — fabricability = own pre-solver gate with reason codes + YAML thresholds, apart from validity; cite sources.
- [Rib air half-disk](project_rib_air_halfdisk.md) — 2026-11-23 proposal: air box -> half-disk (air_r) in rib builder; check if implemented / Q improved.
