# `run_opt_comsol.py` — COMSOL co-optimization driver

Co-optimizes the optical **and** mechanical bandgaps of a 1D diamond OMC unit
cell by driving COMSOL directly, iteration by iteration, with a Bayesian
optimizer (Optuna TPE, or a Halton random-search fallback).

- **Optical target:** TE bandgap ≥ 20%, center near 1550 nm.
- **Mechanical target:** breathing-mode (`fy > 0.5`) gap ≥ 20%, center in 5–10 GHz.

This is the **production driver** for the isotropic-template campaign. It is a
self-contained program: unlike `run_loop.py` (which composes `src/objective.py`,
`src/optimizer.py`, and the `src/*_comsol.py` backends), `run_opt_comsol.py`
carries its own inline candidate evaluator, scoring function, optimizer wrapper,
COMSOL band sweeps, gap detectors, and plotting. The two are **parallel
pipelines targeting the same physics problem**, not layered components — see
[Relationship to `run_loop.py`](#relationship-to-run_looppy).

---

## Prerequisites

- **COMSOL 6.x + MPh** reachable from Python (`comsol_client.get_model`).
  Runs on the user's Mac where COMSOL LiveLink is available.
- **The `.mph` template:** `comsol/omc_unitcell_iso.mph` (hardcoded as
  `_TEMPLATE`). It must expose studies named `mech sym` and `opt TE` (the
  `--study-mech` / `--study-opt` defaults) and geometry parameters
  `a, w, t, hx, hy, kF`. The script exits early if the template is missing.
- **Optuna** (optional). If importable, TPE is used; otherwise the script falls
  back to a built-in Halton sampler. No code change needed either way.
- **matplotlib** (progress figures), **numpy**. `shapely`/`scipy` are used by
  some characterization paths.

---

## Quick start

```bash
# From the python-scripts/ directory:
python scripts/run_opt_comsol.py --n-iter 20
python scripts/run_opt_comsol.py --n-init 8 --n-iter 30 --n-k 9
python scripts/run_opt_comsol.py --n-iter 5 --n-k 5 --skip-characterization   # quick smoke test
```

By default it writes to the `_iso`-suffixed result files (see
[Outputs](#outputs)) and **auto-resumes** from any existing
`results/opt_results_t5d_iso.json`.

---

## How it works

Each candidate is a **5-D normalized vector** `u ∈ [0,1]⁵` mapping to physical
`(a, w, hx, hy, t)` via `src/geometry.py:u_to_geometry`. The main loop runs
**`n_iter` total iterations** (see the note on `n_init` below):

```
for it in range(n_iter):
    u        = opt_ask(optimizer)          # Optuna TPE or Halton
    g        = u_to_geometry(u)
    ok       = check_feasibility(g)        # stage: feasibility
    if not ok: save infeasible record; continue
    rec      = evaluate(g):                # runs COMSOL:
                 _mech_sweep + breathing_gap    # stage: mechanical
                 _opt_sweep  + optical_gap       # stage: optical
                 score_result(...)               # stage: score
    save_result(rec, path=--db)            # → SQLite (+ JSONL fallback)
    opt_tell(optimizer, u, score)          # feed the optimizer
    _save_json / _plot_progress / _save_progress   # every iteration
# end of run:
characterize_best.py (fresh subprocess)    # unless --skip-characterization
```

**`n_init` is a subset of `n_iter`, not additive.** With `--n-iter 30
--n-init 8` you get **30** total recorded evaluations, the first ~8 being
random/Halton startup asks before TPE's model takes over. (This differs from
`run_loop.py`, whose total is `n_init + n_iter`.)

**All trials are recorded** — startup/init asks, TPE asks, and
infeasible/failed candidates alike. There is no separate unsaved
initialization phase. Infeasible candidates get `status="infeasible"`,
`score=-10.0`; caught COMSOL exceptions get `status="failed"`, `score=-5.0`.

**Resume / warm-start.** On startup, prior results in `--out-json` are replayed
into the optimizer (re-scored with the *current* constants and fed as
`add_trial`) to seed its model — this replay is **not** re-saved. Pass
`--fresh-start` to discard prior results (the old JSON is backed up to
`<file>.bak` first).

**Progress bar.** Per-stage progress prints to **stderr**
(`iter 3/30 [2/4] mechanical ... running`), separate from the stdout summary
and the on-disk progress file. Suppress with `--quiet`.

---

## Scoring

The optimizer **maximizes** `score_result()` (higher = better). In the default
co-optimized mode (both gaps required):

```
score = G_o + G_m
        − 20.0 · max(0, 0.15 − G_o)²                    # optical gap shortfall
        −  5.0 · max(0, 0.20 − G_m)²                    # mechanical gap shortfall
        −  9.0 · ((f_o_c − f_1550) / OPT_F_TOLERANCE)²  # optical center off 1550 nm
        − 2000 · max(0, (f_o_c − c/2a) / (c/2a))²       # light-cone penalty
```

- **Base** is the *sum* of the two normalized gaps (not `min`).
- **Gap shortfalls** are one-sided quadratics that only bite below the minimum
  (`G_MIN_OPT=0.15`, `G_MIN_MECH=0.20`). Optical is weighted 4× heavier.
- **Optical center penalty** normalizes the detuning from the 1550 nm target by
  an **absolute tolerance** `OPT_F_TOLERANCE` (a keyword arg of `score_result`,
  default **10 THz**), not by the target frequency. So a fixed absolute detuning
  costs the same regardless of target; a smaller tolerance stiffens centering.
  (Previously the denominator was `OPT_F_TARGET` ≈ 193.4 THz, making this penalty
  ~370× weaker — see the note below.)
- **Light-cone penalty** discourages an optical gap sitting *above* the
  zone-edge light line `c/(2a)` (not a true guided-mode gap). Weight is large
  (2000) because it's quadratic in a small fractional excess.
- **Mechanical frequency penalty is disabled** (`LAMBDA_MECH_F=0`).
- **Optical-only mode** (`--no-require-mech`) drops `G_m` and its penalty and
  instead adds thinness + narrowness rewards (favoring smaller `t`/`w`).
  Mechanical data is still computed and recorded, just not optimized for.
- Non-finite gaps return the sentinel `-5.0`.

Constants live at the top of the file (~lines 43–74) and have been hand-tuned
(dated 2026-07 comments). **Scores are not comparable to `run_loop.py`**, which
uses a different scoring function (`min`-based, weights from `targets.yaml`, no
light-cone term).

> **Scoring change (2026-07-16):** the optical center penalty now divides by an
> absolute `OPT_F_TOLERANCE` (default 10 THz) instead of `OPT_F_TARGET`
> (≈193.4 THz). This makes the penalty roughly **370× stronger** for the same
> detuning, so scores from this version are **not comparable** to runs made
> before the change. When resuming an older `--out-json`, prior results are
> re-scored on load with the new tolerance, so the warm-start stays internally
> consistent — but absolute score values will shift versus the old files.

---

## CLI reference

| Flag | Default | Meaning |
|---|---|---|
| `--n-init` | `6` | Random/Halton startup trials before TPE engages (subset of `--n-iter`). |
| `--n-iter` | `14` | **Total** iterations, including the init trials. |
| `--n-k` | `7` | k-points for the optical zone-edge sweep and the mechanical sweep. |
| `--n-bands-mech` | `8` | Mechanical bands per k-point (optimization loop). |
| `--n-bands-opt` | `6` | Optical bands per k-point (optimization loop). |
| `--study-mech` | `mech sym` | COMSOL study name for the mechanical sweep. |
| `--study-opt` | `opt TE` | COMSOL study name for the optical sweep. |
| `--no-fy` | off | Skip per-mode `fy` field evaluation (saves ~30–50% COMSOL time when `mech sym` BCs already isolate breathing modes). |
| `--out-fig` | `results/figures/opt_progress_iso.png` | Progress figure (4-panel). |
| `--out-json` | `results/opt_results_t5d_iso.json` | Results JSON (resume source of truth; band arrays stripped). |
| `--db` | `results/runs_iso.sqlite3` | SQLite DB for candidate records. Overrides `OMC_DB_PATH` for this run; JSONL fallback uses the same stem. |
| `--progress-file` | `results/progress_iso.txt` | Human-readable status snapshot; `tail -f` it. |
| `--explore-every` | `0` | Force a random exploration sample every N iterations (0 = off). |
| `--resume` | off | **Legacy no-op** — prior results always load automatically if `--out-json` exists. |
| `--fresh-start` | off | Discard prior results (backs up the old JSON to `.bak` first). |
| `--skip-characterization` | off | Skip the end-of-run characterization subprocess. |
| `--no-require-mech` | off | Drop mechanical gap from the score (still recorded). Enables thinness/narrowness rewards. |
| `--no-require-opt` | off | Symmetric to `--no-require-mech`, for the optical gap. |
| `--g-min-opt` | `0.15` | Minimum optical gap for scoring. Ignored under `--no-require-opt`. |
| `--g-min-mech` | `0.20` | Minimum mechanical gap for scoring. Ignored under `--no-require-mech`. |
| `--quiet` | off | Suppress per-stage progress on stderr. |

---

## Outputs

| Data | Default path | Flag | When | Notes |
|---|---|---|---|---|
| Results JSON (resume source of truth) | `results/opt_results_t5d_iso.json` | `--out-json` | every iteration | `_`-prefixed band arrays **stripped**. |
| SQLite DB (full record incl. band arrays) | `results/runs_iso.sqlite3` | `--db` | every iteration | `record` column holds the full dict; `--db` overrides `OMC_DB_PATH`. |
| JSONL fallback | `results/runs_iso.jsonl` | (`--db` stem) | on SQLite lock error | safety net for network/FUSE mounts. |
| Progress figure | `results/figures/opt_progress_iso.png` | `--out-fig` | every iteration | score history + Pareto + gap-vs-geometry. |
| Progress text | `results/progress_iso.txt` | `--progress-file` | every iteration | `tail -f` friendly dashboard. |
| JSON backup | `results/opt_results_t5d_iso.json.bak` | `--fresh-start` | once, at start | only with `--fresh-start`. |
| Characterization log / figure / cache | `results/figures/{char_log.txt, best_characterization.png, best_char_data.npz}` | via `--out-fig` dir | once, at end | produced by the `characterize_best.py` subprocess. |
| Mode field exports + manifest | `results/mode_exports/*.txt` | — | once, at end | raw COMSOL field profiles. |

> **Note:** the JSON and the SQLite record are *not* identical — the raw
> band-structure arrays (`_k_m`, `_freqs_m`, `_fy_arr`, `_k_o`, `_freqs_o`)
> survive only in the SQLite `record` column, never in the JSON.

---

## Monitoring a run

```bash
# Live text dashboard (updates each iteration):
tail -f results/progress_iso.txt

# Rich polling dashboard in a second terminal (point it at the same JSON):
python scripts/watch_progress.py --file results/opt_results_t5d_iso.json --n-iter 30

# Regenerate the progress figure without touching the running job:
python scripts/replot.py           # defaults to the _iso files
```

---

## Companion scripts

- **`characterize_best.py`** — launched automatically at end of run (unless
  `--skip-characterization`) in a fresh subprocess for a clean COMSOL state;
  produces the best-geometry figure, cache, and mode exports.
- **`watch_progress.py`** — standalone polling dashboard over the results JSON.
- **`replot.py`** — regenerates the progress figure from the JSON.
- **`migrate_t_bounds.py`** — one-time re-mapping of stored `u[4]` (thickness)
  and re-scoring after changing `t_max` in `configs/bounds.yaml`.

---

## Changing the target SQLite database

Use `--db` (recommended, per-campaign):

```bash
python scripts/run_opt_comsol.py --n-iter 30 \
    --db results/runs_mycampaign.sqlite3 \
    --out-json results/opt_results_mycampaign.json
```

`--db` passes an explicit path to `save_result`, so it **takes precedence over
the `OMC_DB_PATH` environment variable** for this script. (The other scripts —
`run_loop.py`, `run_one.py` — have no `--db` flag and still honor
`OMC_DB_PATH`.) A brand-new filename just starts a fresh DB; the `candidates`
table is auto-created.

---

## Relationship to `run_loop.py`

| | `run_loop.py` | `run_opt_comsol.py` |
|---|---|---|
| Size / role | ~65-line modular loop | ~2750-line COMSOL production driver |
| Objective | `src/objective.py` (`min`-based, `targets.yaml`) | inline `score_result` (sum-based + light-cone + geometry rewards) |
| Optimizer | `src/optimizer.py` ask/tell | inline Optuna TPE / Halton, warm-start from JSON |
| Solver access | via `optical_*`/`acoustic_comsol` backends | drives COMSOL directly via `comsol_client` |
| Persistence | `database.py` only | JSON (source of truth) **+** SQLite dual-write + figures + progress + characterization |
| `n_init` meaning | added to total (`n_init + n_iter`) | subset of total (`n_iter`) |

Because the scoring functions and weight sources differ, **a numeric score
means different things under each script** — don't compare across them.

---

## Gotchas

- **Study names must match the template.** The defaults (`mech sym`, `opt TE`)
  must exist in `omc_unitcell_iso.mph`, or `model.solve()` won't resolve.
- **Windows encoding.** Progress/log file writes are UTF-8 (dashboard glyphs
  like `█ ░ │`); the characterization subprocess runs with
  `PYTHONIOENCODING=utf-8` so its redirected stdout doesn't crash on Windows.
- **Three separate scoring definitions** exist in the repo
  (`run_opt_comsol.py`, `src/objective.py`, `migrate_t_bounds.py`), each with
  its own constants. Keep this in mind when rescoring old records.
- Not listed in `CLAUDE.md`'s "Important files," though it is the larger,
  actively-maintained driver.
