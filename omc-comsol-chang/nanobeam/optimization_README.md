# Nanobeam cavity taper optimization

Optimizes the defect taper of a 1-D diamond optomechanical-crystal nanobeam cavity over
two parameters, where every objective evaluation is a full COMSOL mechanical + optical +
PML solve. Because a single evaluation costs minutes, the run is **persisted to SQLite as
it goes** and can be stopped and resumed without losing a solve.

---

## Contents

- [Files](#files)
- [Quick start](#quick-start)
- [Requirements](#requirements)
- [The design space](#the-design-space)
- [The objective](#the-objective)
  - [Fitness function](#fitness-function)
  - [Why the objective is `-log10(fitness)`](#why-the-objective-is--log10fitness)
  - [The penalty ladder](#the-penalty-ladder)
  - [Alternative fitness functions](#alternative-fitness-functions)
- [Optimizers](#optimizers)
- [Stop and resume](#stop-and-resume)
  - [How it works](#how-it-works)
  - [SQLite backends](#sqlite-backends)
  - [Inspecting a run](#inspecting-a-run)
  - [Schema](#schema)
- [Outputs](#outputs)
- [Configuration reference](#configuration-reference)
- [Known gotchas](#known-gotchas)
- [Testing](#testing)
- [Related files](#related-files)

---

## Files

| File | Role |
|---|---|
| `optimize_oblong_maxdef.m` | The optimizer. Nelder–Mead today, `bayesopt` reserved. |
| `OptimStore.m` | SQLite persistence + stop/resume. Reusable, not specific to this script. |
| `test_OptimStore_resume.m` | COMSOL-free proof that resume works. Runs in ~1 s. |
| `bayesopt_oblong_maxdef_PROPOSAL.md` | Design for the Bayesian backend; decisions still open. |
| `sweep_oblong_maxdef.m` | Grid-sweep sibling over the same two parameters. |

---

## Quick start

```matlab
% 1. Start COMSOL with LiveLink for MATLAB.
% 2. Sanity-check without running anything:
checkcode('optimize_oblong_maxdef.m')
OptimStore.selfTest          % verifies the SQLite backend end-to-end
test_OptimStore_resume       % verifies stop/resume, no COMSOL needed

% 3. Edit the TUNABLE KNOBS block — at minimum OPT.runId, the bounds,
%    OPT.QmechMin, and MaxFunEvals in the optimset call.

% 4. Run:
[best_oblong, best_maxdef, fitnessBest, dsBest] = optimize_oblong_maxdef();
```

**Time a single evaluation before committing to a study.** `MaxFunEvals` defaults to
200; at even 10 minutes per solve that is well over a day. Run
`test_nanobeamRectFEM_withPML.m` once at the same `P.mMesh` / `P.oMesh` to get the
per-evaluation cost, and multiply.

To stop a run, Ctrl-C is safe — the run is marked `interrupted` and every completed
evaluation is already in the database. To continue, re-run with the same `OPT.runId`.

---

## Requirements

- **COMSOL with LiveLink for MATLAB, already running.** `RunNanobeamFEM` drives the
  COMSOL Java API directly.
- **Base MATLAB only** for the Nelder–Mead path — `fminsearch` and `optimset` are not
  toolbox functions.
- **A SQLite backend** for persistence. One of Database Toolbox, MATLAB's Python
  interface, or a `sqlite3` executable — see [SQLite backends](#sqlite-backends). Set
  `OPT.useStore = 0` to run without persistence.
- Must be run from `nanobeam/` so `RunNanobeamFEM`, `CreateNanobeamGeom`,
  `LoadMaterialParams` and `CreateFileBase` are on the path.
- `OPT.optimizer = 'bayesopt'` additionally needs the **Statistics and Machine Learning
  Toolbox**, which is *licensed but not installed* on this machine. See the proposal.

---

## The design space

Two continuous variables:

| `x` | Symbol | Meaning |
|---|---|---|
| `x(1)` | `defectAspectRatio` | hole aspect ratio at the cavity centre |
| `x(2)` | `maxdef` | fraction the lattice constant shrinks at the centre |

The script does **not** search `oblong` directly. It searches the *defect aspect ratio*
and inverts the taper relation for each candidate:

```
defectAspectRatio = (hy/hx) · (1-maxdef)^(2·oblong)

  =>  r      = defectAspectRatio / (P.hy/P.hx)
      oblong = log(r) / (2·log(1-maxdef))
```

This is the more physical coordinate: it fixes the *shape* of the hole at the cavity
centre rather than the exponent that gets you there. One consequence worth keeping in
mind — because `oblong` depends on **both** variables, a rectangular box in
`(defectAspectRatio, maxdef)` is a curved region in `(oblong, maxdef)`.

`hx_0` / `hy_0` are defined once at the top of the knobs block and `P.hx` / `P.hy` are
assigned from them, because the start point and the objective's mapping both depend on
`hy/hx`. Writing the ratio out twice lets them drift, and `x0` then stops corresponding
to `oblong_0`.

---

## The objective

### Fitness function

Active fitness is `fitness_coopProduct` — a cooperativity-like figure of merit with a
**binary** mechanical-Q gate:

```
Q_mech <  OPT.QmechMin  ->  fitness = 0                       (rejected)
Q_mech >= OPT.QmechMin  ->  fitness = g_OM² · Q_opt · L(λ)
Q_mech unavailable      ->  fitness = NaN                     (not evaluable)
```

where the wavelength factor is a **Gaussian** in the detuning:

```
L(λ) = exp( -((λ_opt - OPT.targetLambda) / OPT.lambdaTol)² )
```

so `lambdaTol` is the 1/e half-width — a `lambdaTol` miss costs a factor 1/e, a
`2·lambdaTol` miss costs 1/e⁴. Gaussian rather than exponential because it is flat
inside the tolerance (a 20 nm FEM discretization error should not cost 18 % of the
score) and falls off a cliff outside it (1750 nm is a different device, not a slightly
worse one). It is also smooth at the target, whereas `exp(-|Δλ|/tol)` has a kink there
— exactly the ridge geometry that makes a simplex collapse prematurely.

`Q_mech` is a pass/fail gate, not a factor. That makes the landscape piecewise-flat:
every design below the threshold scores exactly 0 and is indistinguishable from every
other. If a run stalls immediately, either lower `OPT.QmechMin` or switch to
`fitness_QxQxLambda`, which varies smoothly with `Q_mech` and so always supplies a
gradient.

### Why the objective is `-log10(fitness)`

`optimset`'s `TolFun` is an **absolute** tolerance on the spread of objective values
across the simplex. Raw fitness here runs to ~10¹⁶–10¹⁸ (`g_OM²·Q_opt`), so minimizing
`-fitness` makes any sane `TolFun` unreachable and the run always burns through
`MaxFunEvals` regardless of convergence.

Minimizing `-log10(fitness)` puts the objective on an O(10) scale, which turns `TolFun`
into an effectively *relative* test on the fitness. `TolFun = 1e-3` is about a 0.23 %
change. Tightening it much further is pointless: adaptive meshing (`P.mAdjMesh`,
`P.oAdjMesh`) makes the objective slightly discontinuous as the geometry changes, so the
numerical noise floor sits above a 1e-6 tolerance.

> **Nelder–Mead has no noise model.** If a run converges onto a point whose fitness does
> not reproduce when re-run, that is a mesh artifact, not an optimum.

### The penalty ladder

Infeasible outcomes are ordered by badness so the simplex prefers drifting back through
solvable territory rather than through the dead zone:

| Outcome | `J` | Meaning |
|---|---|---|
| feasible | `-log10(fitness)` ≈ −10 … −18 | usable design |
| gated | `OPT.Jpenalty` = 1000 | solved cleanly, but `Q_mech` below threshold |
| unusable | `2 · OPT.Jpenalty` = 2000 | no localized mode, or the solve failed |
| out of bounds | `OPT.Jpenalty · (2 + d)` | never solved; `d` = normalized distance outside the box |

The out-of-bounds penalty is a **sloped barrier**, not a flat plateau: a constant
penalty gives the simplex no direction to move once a vertex leaves the box, so it can
collapse against a bound and report false convergence. Scaling by distance — normalized
per parameter so the two axes are comparable — restores a downhill direction back into
the feasible set. No FEM call is made for an out-of-bounds point.

### Alternative fitness functions

Three ready-made objectives sit at the bottom of the file. Switch by editing the single
`OPT.fitnessFcn` line:

| Function | Formula | When to use |
|---|---|---|
| `fitness_coopProduct` **(active)** | `g_OM²·Q_opt·L_gauss(λ)`, binary `Q_mech` gate | default |
| `fitness_QxQxLambda` | `Q_opt·Q_mech·exp(-\|Δλ\|/tol)` | when the binary gate flattens the landscape; `Q_mech` is smooth here |
| `fitness_maxGOM` | `max\|g_OM\|` | quick "where is coupling strongest" scan; ignores both Q's |
| `fitness_maxLSiV` | `max λ_SiV` | SiV strain coupling; **requires `P.calcS = 1`** |

Anything returning NaN is treated as unusable; anything returning 0 is treated as a
rejected-but-evaluated design.

---

## Optimizers

Selected by `OPT.optimizer`:

| Value | Status | Needs |
|---|---|---|
| `'neldermead'` | **live** — `fminsearch` simplex | base MATLAB |
| `'bayesopt'` | **reserved** — raises a specific error | Statistics and ML Toolbox |

The `'bayesopt'` branch deliberately reports *why* it cannot run:

```
bayesopt is not available in this MATLAB install.
  licence permits it (1=yes) : 1
  installed on disk (1=yes)  : 0
Licence 1 with install 0 means: install the Statistics and Machine Learning
Toolbox from the Add-On Explorer. No re-licensing is needed.
```

A bare `Unrecognized function 'bayesopt'` reads like a licence problem and sends you
debugging the wrong thing. See `bayesopt_oblong_maxdef_PROPOSAL.md` for the design and
the open decisions — notably that the penalty ladder above should **not** carry over to
a GP surrogate.

---

## Stop and resume

### How it works

`fminsearch` is **deterministic** given `x0` and `options`. A resumed session therefore
retraces its earlier path exactly before extending it — so the objective checks the
database before spending a solve:

```matlab
hit = store.lookup(OPT.runId, [defectAspectRatio_i, maxdef_i]);
if ~isempty(hit)
    J = hit.J;          % no COMSOL call at all
    return;
end
```

Points are matched on `OptimStore.paramKeyOf` — `sprintf('%.10g|%.10g', ...)` — rounded
so a value that has round-tripped through SQLite's `REAL` storage still matches.

Replayed points are appended to the plotting history too, so the convergence figure
shows the whole search, not just the part solved in the current session.

Practical shape of a resumed run:

```
Run store    : .\test\1D_OMC_hole\optim_runs.sqlite3  [backend: python]
Run id       : oblongMaxdef_trial1  (13 evaluation(s) already stored)
Resuming     : best stored fitness 9.9756e+15 at eval 11 (dAR=1.1804, maxdef=0.2071)
  eval  1: oblong=1.9600  maxdef=0.2200  [CACHED  fitness=8.4e+15  ok]
  eval  2: ...                           [CACHED  ...]
  ...
  eval 14: oblong=1.8873  maxdef=0.2093  ->  .\test\...\eval_0014_ob1.8873_md0.2093\
```

The incumbent is seeded from the store *before* the search starts, so a resumed run can
never report a worse design than one already paid for — even if every evaluation in the
new session fails.

A run is marked `interrupted` automatically by an `onCleanup` guard if the function exits
without reaching the end (an error, or Ctrl-C), and `complete` on the happy path.

### SQLite backends

Auto-detected in this order. Database Toolbox is **not** required.

| Backend | Requires | On this machine |
|---|---|---|
| `dbtoolbox` | Database Toolbox `sqlite()` | licensed, not installed |
| `python` | MATLAB Python interface (`pyenv`) | **active** — Python 3.13, SQLite 3.45.3 |
| `cli` | `sqlite3` on PATH (Anaconda ships one) | available, tested |

An append-only `<db>.jsonl` mirror is written **before** each SQLite insert, so a failed
database write can never discard a solve you have already paid hours for.

### Inspecting a run

```matlab
store = OptimStore(fullfile('.','test','1D_OMC_hole','optim_runs.sqlite3'));

store.listRuns()                          % every study in the database
T = store.loadEvals('oblongMaxdef_trial1');   % table, one row per evaluation
[x, fitness, evalIdx] = store.best('oblongMaxdef_trial1');
n = store.evalCount('oblongMaxdef_trial1');
```

`T` carries `eval, ts, dAR, maxdef, oblong, wM, gOM, LSiV, Qmech, Qopt, lambdaNm,
fitness, J, status`, so the whole search history is available for plotting or filtering
without re-parsing the CSV.

Because it is plain SQLite, `sqlite3` on the command line and any SQL client work too.

### Schema

```sql
runs  (run_id PK, created_ts, updated_ts, script, optimizer, status, config_json)
evals (run_id, eval_idx, ts, param_key, dar, maxdef, oblong,
       wm_hz, gom_hz, lsiv_hz, q_mech, q_opt, lambda_nm,
       fitness, objective_j, status, ev_loc,
       PRIMARY KEY (run_id, eval_idx))
INDEX idx_evals_lookup ON evals(run_id, param_key)
```

`config_json` snapshots the whole `OPT` struct (function handles stringified) so a run's
settings can be recovered months later. `ev_loc` points at the per-evaluation output
folder, so a replayed point still tells you where its mode profiles were written.

`status` values: `ok`, `Qmech_gated`, `no_mech_modes`, `no_opt_modes`,
`no_mech_modes+no_opt_modes`, `failed`.

> `Qmech_gated` is distinct from `failed` on purpose. Conflating "solved cleanly but the
> mechanical Q was too low" with "the solve blew up" makes the log impossible to
> interpret afterwards.

---

## Outputs

Everything lands under `OPT.rootLoc`, except the database, which lives at `OPT.dbPath`:

| Path | Contents |
|---|---|
| `<rootLoc>/eval_NNNN_obX_mdY/` | per-evaluation FEM output — `.mat`, mode-profile `.png`/`.fig` |
| `<rootLoc>/BEST/` | final re-run of the best design, with geometry and mode plots on |
| `<rootLoc>/optimize_log.csv` | one row per evaluation, appended immediately; survives resume |
| `<rootLoc>/optimize_convergence.png` / `.fig` | fitness trace + simplex trajectory |
| `<dbPath>` | SQLite run database (shared across studies) |
| `<dbPath minus extension>.jsonl` | append-only mirror of every evaluation |

`<rootLoc>` is `./test/1D_OMC_hole/optimize_<OPT.runId>/` — keyed on the run id, **not**
date-stamped, so a study resumed weeks later keeps all of its output in one place.

CSV columns: `eval, oblong, maxdef, wM_GHz, gOM_kHz, LSiV_MHz, Q_mech, Q_opt,
lambda_opt_nm, fitness, status`. The file is reopened and appended per evaluation so it
survives a mid-run crash. It is opened `'a'`, never `'w'`: because `<rootLoc>` is stable
across sessions, truncating would erase every earlier session's history. The header is
written only when the file is created, and each resumed session inserts a seam marker:

```
eval,oblong,maxdef,wM_GHz,...,fitness,status
1,1.9600,0.2200,7.1,...,4.5e16,ok
2,1.9418,0.2213,7.0,...,3.9e16,Qmech_gated
# resumed 2026-09-05 17:01:16
3,1.8873,0.2093,7.2,...,6.1e16,ok
```

Read it back with `readtable(logFile, 'CommentStyle', '#')`, which skips the seam lines
and returns the evaluations as one table. The CSV duplicates the database and is kept
for quick eyeballing; `store.loadEvals` is the better programmatic route.

The convergence figure has two panels:

1. **Fitness vs. evaluation** — every evaluated point, a running-best trace, and
   gated/failed points marked along the bottom (they cannot go on a log axis, and
   hiding them would make a mostly-rejected run look like gaps).
2. **Simplex trajectory** — the path through the `(defectAspectRatio, maxdef)` box,
   coloured by evaluation order, with the start point, the best point, and the bounds.

---

## Configuration reference

### Search and objective

| Knob | Default | Meaning |
|---|---|---|
| `hx_0`, `hy_0` | 196 nm, 578 nm | nominal mirror hole dims; `P.hx`/`P.hy` derive from these |
| `oblong_0`, `maxdef_0` | 1.96, 0.22 | start point, expressed in the physical parameters |
| `OPT.defectAspectRatio_min/max` | ±50 % of start | box bounds on `x(1)` |
| `OPT.maxdef_min/max` | ±50 % of start | box bounds on `x(2)` |
| `OPT.targetLambda` | 1550 nm | optical resonance target |
| `OPT.lambdaTol` | 100 nm | **Gaussian 1/e half-width**, not a decay length |
| `OPT.QmechMin` | 1e7 | binary gate; fitness = 0 below this |
| `OPT.targetFreq` | 7 GHz | mechanical solver target (`P.freq`) |
| `OPT.Jpenalty` | 1e3 | penalty scale (see the ladder above) |
| `OPT.fitnessFcn` | `fitness_coopProduct` | swap for an alternative |

### Optimizer and persistence

| Knob | Default | Meaning |
|---|---|---|
| `OPT.optimizer` | `'neldermead'` | `'neldermead'` or `'bayesopt'` (reserved) |
| `OPT.useStore` | `1` | 0 disables SQLite persistence entirely |
| `OPT.runId` | `'oblongMaxdef_trial1'` | names the study — **keep it to resume, change it to start fresh** |
| `OPT.dbPath` | `./test/1D_OMC_hole/optim_runs.sqlite3` | shared database across studies |
| `OPT.rootLoc` | `./test/1D_OMC_hole/optimize_<runId>/` | output folder, derived from `OPT.runId` |
| `OPT.resume` | `1` | 0 re-solves everything (database still written) |
| `optimset(...)` | `TolX 1e-4`, `TolFun 1e-3`, `MaxIter 100`, `MaxFunEvals 200` | simplex stopping rules |

`OPT.runId` is deliberately **not** date-stamped: a date-stamped id would silently start
a new study tomorrow instead of continuing today's.

### Physics (the expensive settings)

`P.solveMech = 1`, `P.solveOpt = 1`, `P.calcG = 1`, `P.solveMechPML = 1`, `P.calcS = 0`.

`P.solveMechPML = 1` with `P.PMLLen = 10e-6` is what makes `Q_mech` available — and it
is a substantial share of each evaluation's cost. It is load-bearing: `Q_mech` is the
acceptance gate. Turning the PML off makes `Q_mech` NaN, which the objective treats as
*unusable* rather than *gated*, so every point would be penalized.

`P.plotMech`/`P.plotOpt`/`P.saveplots` are on during the search, so each evaluation
writes mode profiles into its own folder. Set them to 0 for a faster, quieter run.

---

## Known gotchas

**Renaming `OPT.runId` mid-study orphans the old output.** Output paths are derived
from it (`optimize_<runId>/`), and lookups are keyed on it, so changing it starts a
genuinely new study: a fresh folder, a fresh CSV, and no replay of the stored points.
That is the intended way to start over — just be aware the old run is still in the
database under its old id, not deleted.

**Deleting `<rootLoc>` but keeping the database leaves dangling `ev_loc` pointers.**
Replayed points still resolve their `J` correctly (that comes from the database), but
the folder their mode profiles were written to is gone.

**This script and `sweep_oblong_maxdef.m` use different device geometry** — `w = 750 nm`,
`hy = 578 nm` here versus `w = 800 nm`, `hy = 651 nm` in the sweep. Their results are
**not** directly comparable, and sweep points cannot be reused to seed this optimizer
without re-solving.

**The optics-unavailable branch of the fitness returns `g_OM²` alone**, which is not on
the same scale as the full expression. With `P.solveOpt = 1` this should never trigger,
but a run that mixes the two branches cannot be ranked meaningfully.

**A cache hit does not re-solve, so it does not re-create its output folder.** If you
delete `<rootLoc>` but keep the database, replayed points will have `ev_loc` pointing at
folders that no longer exist.

---

## Testing

Neither test needs COMSOL.

```matlab
OptimStore.selfTest
% Exercises the active backend: schema, insert, lookup, best, NULL round-trip,
% reopen, JSONL mirror.

test_OptimStore_resume
% Runs three optimizer sessions against an analytic fitness and asserts that
% session 2 replays session 1's path *identically* (to 1e-12) from the database
% before extending it, and that a fresh runId shares nothing.
```

Expected output of the resume test:

```
session 1 : 13 solved,  0 cached, best fitness 9.97556e+15
session 2 : 18 solved, 13 cached, best fitness 9.99996e+15
            replayed 13 point(s) identically, then extended
session 3 :  8 solved,  0 cached (fresh runId)
ALL RESUME ASSERTIONS PASSED
```

The identical-path assertion is the important one: if the replay diverged, the
fast-forward would be serving stored answers for a trajectory the optimizer is no longer
on.

Static analysis:

- `checkcode('OptimStore.m')` — clean.
- `checkcode('optimize_oblong_maxdef.m')` — clean apart from three "function might be
  unused" warnings for the intentionally-kept alternative fitness functions, and one
  "value might be unused" for a defensive `status_i` initialization.

---

## Related files

| File | Relationship |
|---|---|
| `sweep_oblong_maxdef.m` | grid sweep over the same two parameters (different geometry) |
| `RunNanobeamFEM.m` | the per-evaluation pipeline: geometry → mesh → solve → post-process |
| `CreateNanobeamGeom.m` | builds the hole arrays; line 77 is the lithography-gap check |
| `SolveNanobeamFEM.m` | source of `ds.mfem.QAll`, `ds.ofem.Q`, `ds.cpl.gMax` |
| `../bayesopt_boomerang.m`, `../bayesopt_README.md` | MATLAB `bayesopt` reference implementation |
| `../optimize_hole_unitCell.m` | Nelder–Mead optimizer for hole *unit cells* |
| `../../python-scripts/src/database.py` | the SQLite pattern `OptimStore` follows |
| `../CLAUDE.md` | repo conventions — `P` struct, path handling, file naming |
