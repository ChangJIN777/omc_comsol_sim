---
name: result-db-quirks
description: Non-obvious facts about the hole/rib result databases (runs.sqlite3, runs_rib.sqlite3) that affect any analysis of saved records
metadata:
  type: project
---

Found while building scripts/plot_mech_bands.py (2026-11-24); verified on the real databases.

- Hole loop ids (run_opt_comsol._make_id) depend on geometry (+xsect, index) only, and save_result does INSERT OR REPLACE. So runs.sqlite3 keeps only the LAST evaluation of a cell. 18 records of the fy campaign (opt_results_t5d_iso_fy.json) were overwritten by later --no-fy / optical-only re-evaluations, so the out_json and SQLite values differ for those ids.
- 27 hole records saved 2026-07-01 ~17:20-18:15 carry the LEGACY breathing gap: fy < 0.5 set to NaN point by point, then largest_gap(nan_safe), so both edges had to be breathing. The current breathing_gap reproduces every other record exactly.
- Hole records from src/objective.py (run_one/run_loop, ids cand_*) store no band arrays. Only run_opt_comsol.py records (opt_*) have _k_m/_freqs_m/_fy_arr, from the study_mech study, which is not recorded.
- Rib records: _freqs_m = {sym, antisym}. All 29 successful records match interval_gap exactly. Records with score_version None or 2 have no pen_window_m.

**Why:** these make recorded-vs-recomputed comparisons look broken when they are not.
**How to apply:** when analysing or rescoring saved records, expect these cases and do not "fix" them by rewriting the databases (results/ is read-only, see [[python-scripts-safety-rules]]).
