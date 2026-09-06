# Nanobeam cavity taper optimization

Optimizes the defect taper of a 1-D diamond optomechanical-crystal nanobeam cavity over
two parameters, where every objective evaluation is a full COMSOL mechanical + optical +
PML solve. Because a single evaluation costs minutes, the run is **persisted to SQLite as
it goes** and can be stopped and resumed without losing a solve.

Two optimizers are selectable — **Bayesian** (`bayesopt`, the default) and **Nelder–Mead**
(`fminsearch`). Both write the same physics to the same database, so a study started with
one can be resumed or seeded with the other.

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
  - [Why Bayesian by default](#why-bayesian-by-default)
  - [The gate becomes a constraint, not a penalty](#the-gate-becomes-a-constraint-not-a-penalty)
  - [Budget semantics](#budget-semantics--read-this-before-resuming)
  - [Resume, and cross-optimizer seeding](#resume-and-cross-optimizer-seeding)
- [Fabrication pre-filter](#fabrication-pre-filter)
  - [How each optimizer uses it](#how-each-optimizer-uses-it)
  - [Startup scan](#startup-scan)
- [Stop and resume](#stop-and-resume)
  - [How it works](#how-it-works)
  - [SQLite backends](#sqlite-backends)
  - [Inspecting a run](#inspecting-a-run)
  - [Schema](#schema)
- [Outputs](#outputs)
- [Coupling breakdown](#coupling-breakdown)
- [Configuration reference](#configuration-reference)
  - [Search and objective](#search-and-objective)
  - [Optimizer and persistence](#optimizer-and-persistence)
  - [Physics (the expensive settings)](#physics-the-expensive-settings)
- [Known gotchas](#known-gotchas)
- [Testing](#testing)
- [Related files](#related-files)

---

## Files

| File | Role |
|---|---|
| `optimize_oblong_maxdef.m` | The optimizer. Bayesian by default, Nelder–Mead switchable. |
| `OptimStore.m` | SQLite persistence + stop/resume. Reusable, not specific to this script. |
| `isFabricable.m` | Lithography feasibility of a design, without any solve. |
| `seedBayesFromStore.m` | Turns stored evaluations into `bayesopt` seed data. |
| `test_OptimStore_resume.m` | COMSOL-free proof that resume works. Runs in ~1 s. |
| `test_bayesopt_wiring.m` | COMSOL-free proof of the Bayesian plumbing. |
| `test_isFabricable.m` | Cross-validates the pre-filter against `CreateNanobeamGeom`. |
| `bayesopt_oblong_maxdef_PROPOSAL.md` | The design this backend was built from. |
| `sweep_oblong_maxdef.m` | Grid-sweep sibling over the same two parameters. |

---

## Quick start

```matlab
% 1. Start COMSOL with LiveLink for MATLAB.
% 2. Sanity-check without running anything (none of these need COMSOL):
checkcode('optimize_oblong_maxdef.m')
OptimStore.selfTest          % SQLite backend, end to end
test_OptimStore_resume       % stop/resume
test_isFabricable            % pre-filter vs. the real geometry builder
test_bayesopt_wiring         % Bayesian plumbing (needs Statistics Toolbox)

% 3. Edit the TUNABLE KNOBS block — at minimum OPT.runId, the bounds,
%    OPT.QmechMin, and the evaluation budget:
%      bayesopt    -> OPT.bo.maxEvals   (default 60)
%      neldermead  -> MaxFunEvals in the optimset call (default 200)

% 4. Run:
[best_oblong, best_maxdef, fitnessBest, dsBest] = optimize_oblong_maxdef();
```

**Time a single evaluation before committing to a study.** The default budget is 60
COMSOL solves under `bayesopt`, 200 under Nelder–Mead; at even 10 minutes each that is
half a day to well over a day. Run `test_nanobeamRectFEM_withPML.m` once at the same
`P.mMesh` / `P.oMesh` to get the per-evaluation cost, and multiply.

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
- **Statistics and Machine Learning Toolbox** for `OPT.optimizer = 'bayesopt'` (the
  default). Verify with `exist('bayesopt')` — 2 or 5 means available, 0 means not. Set
  `OPT.optimizer = 'neldermead'` to run with base MATLAB only.

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

Both optimizers minimize `-log10(fitness)` rather than `-fitness`. The transform is
**monotonic**, so it cannot change which design wins — `argmax f = argmax log f`. What it
changes is how well each optimizer's machinery works on the way there, and it helps the
two for partly different reasons.

#### The shared reason: a product spanning many decades

```
fitness = g_OM² · Q_opt · L(λ)      →      log f = 2·log g_OM + log Q_opt + log L
```

Every factor is a gain or a ratio, so the physics combines **multiplicatively**. Across a
representative box the fitness spans about **4.6 decades** (2×10¹³ → 7×10¹⁷). In linear
space that landscape is one enormous spike on a floor that is numerically flat: a mediocre
design at 2×10¹³ and a bad one at 2×10¹¹ differ by 0.004 % of the peak. Log compresses the
same landscape to a range of 4.6, where every region has structure an optimizer can act
on — and makes progress additive, so `-16.65` can be read off as "2 from `g_OM²`, 5.7 from
`Q_opt`, −0.05 from detuning" to see which factor is limiting.

It also fixes *what counts as progress*. Doubling `g_OM` is the same engineering
achievement whether `Q_opt` is 10⁵ or 10⁶, but in linear space those two cases differ by a
factor of 10 in objective change, so the optimizer weights them completely differently.

#### For Nelder–Mead: it makes `TolFun` reachable, and scale-free

`fminsearch` terminates only when **both** conditions hold:

```
max |f(vertex) − f(best)| ≤ TolFun     AND     simplex size ≤ TolX
```

`TolFun` is an **absolute** tolerance on the spread of function values. Converging on a
peak of 4.5×10¹⁶:

| simplex within | spread with `J = -fitness` | spread with `J = -log10(f)` |
|---|---|---|
| 50 % of the optimum | 2.3e16 | 0.176 |
| 10 % | 4.5e15 | 0.041 |
| 1 % | 4.5e14 | 0.0043 |
| 0.2 % | 9.0e13 | 0.00087 |

A simplex clustered within **0.2 %** of the optimum — converged by any reasonable standard
— still has an absolute spread of 9×10¹³. For `TolFun` to fire you would need it around
1e13, and that value is only correct for *this* fitness magnitude: switch to
`fitness_maxGOM` (which returns ~10⁵ instead of ~10¹⁶) and it is wrong by eleven orders of
magnitude. Without the transform the run simply never satisfies `TolFun` and always burns
through `MaxFunEvals` regardless of convergence.

In log space `TolFun` becomes a *relative* test, because `Δlog₁₀f = Δf / (f·ln10)`:

| `TolFun` | change in fitness |
|---|---|
| 1e-2 | 2.33 % |
| 1e-3 | 0.23 % |
| 1e-4 | 0.023 % |

The same number now means the same thing whichever fitness function is plugged in.
Tightening much past `1e-3` is pointless anyway: adaptive meshing (`P.mAdjMesh`,
`P.oAdjMesh`) puts the numerical noise floor above a `1e-6` tolerance.

There is a second, subtler benefit. Nelder–Mead's reflect / expand / contract decisions are
pure **comparisons** between vertex values, so the flat-floor problem above hits it
directly: in the low-fitness region every vertex looks identically worthless and the
simplex has nothing to descend.

#### For bayesopt: it makes the GP's assumptions approximately true

A GP with a stationary kernel assumes the function varies on comparable scales everywhere,
and `IsObjectiveDeterministic = false` fits a **single** noise variance for the whole
domain. A 4.6-decade objective violates the first assumption; relative mesh noise violates
the second. One transform fixes both.

Mesh-adaptation noise is **relative** — a geometry re-solves to within a few percent, not
to within a few Hz. With 3 % mesh noise, here is what each optimizer actually sees:

| fitness | σ with `-fitness` | σ with `-log10(f)` |
|---|---|---|
| 1e11 | 3.0e9 | **0.0130** |
| 1e14 | 3.0e12 | **0.0128** |
| 1e17 | 3.0e15 | **0.0125** |

In linear space the noise is **heteroscedastic**, tracking the value across six orders of
magnitude — something a single fitted σ cannot represent. The GP would badly overestimate
uncertainty near the peak and underestimate it in the tails, and expected-improvement would
chase noise. In log space the noise is **constant**, which is exactly the homoscedastic
model the GP assumes. `OPT.bo.deterministic = false` and this transform are really one
decision seen from two ends.

The same argument applies to the kernel: a stationary kernel is a poor fit to a function
varying over 4.6 decades, and a good one for a function varying over 4.6 units.

#### One wrinkle the transform forces

`log(0)` is undefined, so a gated design (`fitness = 0`) cannot be expressed in the
transformed objective at all. That is the transform *exposing* something true rather than
creating a problem: "rejected" was never a fitness value, it is a feasibility statement.
Hence the separate handling — a penalty tier under Nelder–Mead
([the penalty ladder](#the-penalty-ladder)) and a coupled constraint under `bayesopt`
([the gate becomes a constraint](#the-gate-becomes-a-constraint-not-a-penalty)).

**Why base 10 rather than `ln`?** Legibility only — `-16.65` reads directly as
"fitness ≈ 10^16.65". Any base is equivalent up to a constant factor that just rescales
`TolFun`.

> **Nelder–Mead still has no noise model.** The transform makes the noise *uniform across
> the domain*; only `bayesopt` actually models it. If a Nelder–Mead run converges onto a
> point whose fitness does not reproduce when re-run, that is a mesh artifact, not an
> optimum.

### The penalty ladder

Infeasible outcomes are ordered by badness so the simplex prefers drifting back through
solvable territory rather than through the dead zone:

| Outcome | `J` | Meaning |
|---|---|---|
| feasible | `-log10(fitness)` ≈ −10 … −18 | usable design |
| gated | `OPT.Jpenalty` = 1000 | solved cleanly, but `Q_mech` below threshold |
| unusable | `2 · OPT.Jpenalty` = 2000 | no localized mode, or the solve failed |
| unfabricable | `OPT.Jpenalty · (2 + v)` | rejected by the pre-filter; `v = -margin/50 nm`, capped at 10 |
| out of bounds | `OPT.Jpenalty · (2 + d)` | `d` = normalized distance outside the box |

The last two are **sloped barriers**, not flat plateaus: a constant penalty gives the
simplex no direction to move once a vertex leaves the feasible region, so it can collapse
against a bound and report false convergence. Scaling by distance — normalized so the two
axes are comparable — restores a downhill direction back in. Neither costs a FEM call.

This ladder is **Nelder–Mead only**. Under `bayesopt` the same outcomes become coupled
constraints and an `XConstraintFcn` instead — see [Optimizers](#optimizers).

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

Selected by `OPT.optimizer`. Both are live:

| Value | Method | Needs |
|---|---|---|
| `'bayesopt'` **(default)** | GP surrogate + expected-improvement-plus | Statistics and ML Toolbox |
| `'neldermead'` | `fminsearch` simplex | base MATLAB |

### Why Bayesian by default

The binding constraint is the **evaluation budget**, not the convergence rate — each
evaluation is a multi-minute COMSOL solve. Two properties of this objective also defeat a
simplex, and a GP handles both:

1. **The `Q_mech` gate makes the landscape piecewise-flat.** Every design below
   `OPT.QmechMin` scores exactly 0, so a simplex inside that region has no direction to
   move. A GP models the plateau *and its boundary*, and proposes points at the edge.
2. **The objective is mildly noisy** (adaptive meshing). Nelder–Mead has no noise model
   and will contract onto a mesh artifact; `IsObjectiveDeterministic = false` fits a
   noise term instead.

Against that: the search is only 2-D, where BO's advantage is smallest. The case rests on
the budget and the plateau, not on dimensionality — which is why the switch exists rather
than a replacement.

### The gate becomes a constraint, not a penalty

This is the one part that is **not** a mechanical port of the Nelder–Mead path.

Feeding the penalty ladder (1000, 2000) to a GP that also sees feasible values around −16
wrecks the kernel's length scale: the surrogate spends its capacity modelling the penalty
cliff instead of the physics. Under `bayesopt` the infeasible cases become **coupled
constraints** (`≤ 0` means satisfied):

| | Nelder–Mead | bayesopt |
|---|---|---|
| feasible | `-log10(fitness)` | `-log10(fitnessRaw)`, both constraints `< 0` |
| `Q_mech` gated | `J = 1000` | real objective + `con(1) = 1 - Qmech/QmechMin > 0` |
| solve failed | `J = 2000` | `objective = NaN`, `con(2) = +1` |
| out of bounds | sloped barrier | never proposed — the box *is* the variable range |

So the GP learns **where the gate boundary is**, which is the actual design question,
rather than learning that a wall exists.

That requires the **ungated** fitness, which is why `OPT.fitnessFcnRaw` exists (the same
`fitness_coopProduct` with `QmechMin = -Inf`) and why the store keeps both `fitness` and
`fitness_raw`. A gated `0` is the same number everywhere below the threshold and tells a
surrogate nothing.

### Budget semantics — read this before resuming

`OPT.bo.maxEvals` is the **total** study size, and points seeded from the store **count
against it** — `bayesopt`'s `MaxObjectiveEvaluations` includes `InitialX`. Resuming a
study that already holds `maxEvals` points therefore does *nothing*. Raise `maxEvals` to
extend it, exactly as `resume()` would. The script warns rather than silently no-op'ing:

```
Warning: The store already holds 60 point(s) for run "oblongMaxdef_trial1", which
meets or exceeds OPT.bo.maxEvals = 60. bayesopt counts seeded points against its
budget, so NO new evaluations will be made. Raise OPT.bo.maxEvals to extend the study.
```

### Resume, and cross-optimizer seeding

`seedBayesFromStore` reads the run's stored evaluations and returns `InitialX` /
`InitialObjective` / `InitialConstraintViolations`. This is better than a `.mat`
checkpoint of the study object: it survives a MATLAB version change, it is queryable with
plain SQL, and **a Nelder–Mead run can seed a Bayesian run** — the comparison worth doing.

Rows written before `fitness_raw` existed cannot seed the GP (their ungated score is
unknown, and inventing one would poison the surrogate). They stay in the database and
still serve the Nelder–Mead cache.

---

## Fabrication pre-filter

`CreateNanobeamGeom` throws on three lithography guards (`CreateNanobeamGeom.m:72-84`).
Reaching that throw means the FEM pipeline has already started building geometry, so an
unfabricable candidate costs a real fraction of a solve. `isFabricable` reproduces those
guards as pure arithmetic — microseconds instead of minutes:

```matlab
[ok, margin, why] = isFabricable(defectAspectRatio, maxdef, P);
```

`margin` is the tightest rule's slack **in metres** (negative = violated, and by how
much), which is what lets the Nelder–Mead path build a *graded* penalty rather than a
flat one. `why` names the binding rule.

Checked at both the nominal mirror cell and the cavity centre, where the taper is most
extreme:

| Rule | Default | Override |
|---|---|---|
| beam edge to hole edge | 148 nm | `P.minSidewallGap` |
| gap between adjacent holes | 50 nm | `P.minHoleGap` |
| smallest hole dimension | 50 nm | `P.minFeature` |

Defaults match `CreateNanobeamGeom`'s hard-coded values exactly, so the filter can never
disagree with the code that actually throws.

The binding rule is normally `a_def - hx_def`: into the taper the lattice contracts
(`a_def = a(1-maxdef)`) while the hole *widens* (`hx_def = hx(1-maxdef)^(1-oblong)` grows
for `oblong > 1`), so the solid bridge between holes is what runs out first.

### How each optimizer uses it

| Optimizer | Mechanism | Cost of an unfabricable candidate |
|---|---|---|
| `bayesopt` | `XConstraintFcn` | never proposed — rejected before the objective is called |
| `neldermead` | graded penalty in `evalPoint` | one arithmetic evaluation |

`XConstraintFcn` is the better of the two: it prunes *candidates*, so an unfabricable
design never becomes an evaluation at all. This is distinct from the coupled constraints
of the `Q_mech` gate, which need an evaluation before the GP can learn from them.
`isFabricable` is vectorized because `bayesopt` passes a table of many candidate rows.

### Startup scan

Every run scans a 41×41 grid of the box before starting and reports what fraction is
reachable:

```
Fabricable   : 100.0% of the search box (tightest margin 2.2 nm)
```

A fully infeasible box is an **error** — the study is not started, rather than burning a
day discovering it. Below 25 % is a warning, because seed points are drawn from the whole
box and most would be rejected.

> **At the current bounds this filter prunes nothing** — the whole box is fabricable, with
> 2.2 nm to spare at the tightest corner. It is insurance rather than a saving today: it
> becomes load-bearing the moment the bounds are widened, and the startup line tells you
> immediately if they were widened too far. The 2.2 nm margin is also worth knowing on its
> own — the box is closer to the lithography limit than it looks.

---

## Stop and resume

### How it works

Both optimizers resume from the same table, by different routes:

| Optimizer | Resume mechanism |
|---|---|
| `neldermead` | memoized objective — the deterministic path is replayed from the store for free |
| `bayesopt` | `seedBayesFromStore` → `InitialX` / `InitialObjective` / `InitialConstraintViolations` |

**Nelder–Mead.** `fminsearch` is **deterministic** given `x0` and `options`. A resumed
session therefore retraces its earlier path exactly before extending it — so the objective
checks the database before spending a solve:

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

**bayesopt.** A GP has no path to retrace, so the stored evaluations are handed to it as
seed data instead — see [Resume, and cross-optimizer seeding](#resume-and-cross-optimizer-seeding).
Note the budget caveat in the same section: seeded points count against `OPT.bo.maxEvals`.

Either way the incumbent is seeded from the store *before* the search starts, so a resumed
run can never report a worse design than one already paid for — even if every evaluation
in the new session fails.

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

`T` carries `eval, ts, dAR, maxdef, oblong, wM, gOM, gMB, gPE, LSiV, Qmech, Qopt,
lambdaNm, fitness, fitnessRaw, J, status`, so the whole search history is available for plotting or filtering
without re-parsing the CSV.

Because it is plain SQLite, `sqlite3` on the command line and any SQL client work too.

### Schema

```sql
runs  (run_id PK, created_ts, updated_ts, script, optimizer, status, config_json)
evals (run_id, eval_idx, ts, param_key, dar, maxdef, oblong,
       wm_hz, gom_hz, gmb_hz, gpe_hz, lsiv_hz, q_mech, q_opt,
       lambda_nm, fitness, fitness_raw, objective_j, status, ev_loc,
       PRIMARY KEY (run_id, eval_idx))
INDEX idx_evals_lookup ON evals(run_id, param_key)
```

A database created before `gmb_hz` / `gpe_hz` / `fitness_raw` existed is migrated in place
on open —
`OptimStore.migrateSchema` attempts the `ALTER TABLE ADD COLUMN` and treats a
duplicate-column error as success. No rows are rewritten; historical evaluations keep
their data and read back `NaN` for the new columns.

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
| `<rootLoc>/optimize_convergence.png` / `.fig` | fitness trace + sampled designs |
| `<rootLoc>/bayesopt_results.mat` | the `BayesianOptimization` object (`bayesopt` only) |
| `<dbPath>` | SQLite run database (shared across studies) |
| `<dbPath minus extension>.jsonl` | append-only mirror of every evaluation |

`<rootLoc>` is `./test/1D_OMC_hole/optimize_<OPT.runId>/` — keyed on the run id, **not**
date-stamped, so a study resumed weeks later keeps all of its output in one place.

CSV columns: `eval, oblong, maxdef, wM_GHz, gOM_kHz, gMB_kHz, gPE_kHz, LSiV_MHz,
Q_mech, Q_opt, lambda_opt_nm, fitness, status`. `gMB_kHz` and `gPE_kHz` are the
moving-boundary and photoelastic parts of the coupling — see
[Coupling breakdown](#coupling-breakdown). The file is reopened and appended per evaluation so it
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

Unfabricable candidates are **not** logged — they are not evaluations. They are counted
and reported at the end of the run:

```
Total evaluations      = 47
Rejected unfabricable  = 0 (no solve attempted)
Served from store      = 13 (COMSOL solves skipped)
```

The convergence figure has two panels:

1. **Fitness vs. evaluation** — every evaluated point, a running-best trace, and
   gated/failed points marked along the bottom (they cannot go on a log axis, and
   hiding them would make a mostly-rejected run look like gaps).
2. **Where the search went** — the `(defectAspectRatio, maxdef)` box, coloured by
   evaluation order, with the start point, the best point, and the bounds. Under
   Nelder–Mead the points are joined into the simplex path; under `bayesopt` they are
   not, because a GP jumps wherever the acquisition function points and a connecting
   line would imply a trajectory that does not exist.

---

## Coupling breakdown

`CalcGOM` splits the optomechanical coupling into two physical mechanisms, and both are
carried through to the CSV and the database alongside the net. **`README_calcGOM.md` is the
reference for the physics** — every equation, its COMSOL expression, and the sign
conventions; what follows is only what the optimizer does with the results.

| Quantity | Source | Meaning |
|---|---|---|
| `gMB` | `cpl.gMBmax` | **moving boundary** — the dielectric interface shifting |
| `gPE` | `cpl.gPEmax` | **photoelastic** — strain changing the refractive index |
| `gOM` | `cpl.gMax` | net, `= gMB + gPE` |

`gMB` and `gPE` are stored **signed**; `gOM` is a magnitude. The two mechanisms carry
opposite signs (`CalcGOM.m:220` and `:278`) and routinely partially cancel, so the net
alone hides something important: a design can have a small `gOM` while both mechanisms
are individually large. Such a design sits in a **cancellation notch** and is fragile —
a small fabrication error moves it a long way. Taking `abs()` of each part would erase
exactly that signal, which is why `local_getGOM` returns them signed.

Upstream, `CalcGOM` also builds `ds.cpl.breakdown` — a table with one row per evaluated
(optical, mechanical) mode pair, sorted by `|gOM|`:

```matlab
>> ds.cpl.breakdown(1:3, {'oSol','mSol','wM_Hz','gMB_Hz','gPE_Hz','gOM_Hz','cancelRatio'})
```

| Column | Meaning |
|---|---|
| `oSol`, `mSol` | absolute COMSOL solution numbers for the pair |
| `wM_Hz`, `lambda_nm` | mechanical frequency, optical wavelength |
| `gMB_Hz`, `gPE_Hz`, `gOM_Hz` | the three contributions, signed |
| `fracMB`, `fracPE` | each part as a fraction of the net; they sum to 1 |
| `cancelRatio` | `(\|gMB\|+\|gPE\|)/\|gOM\|` — 1 = reinforcing, >1 = opposing |
| `gPE_p11_Hz`, `gPE_p12_Hz`, `gPE_p44_Hz` | photoelastic part split by tensor component |

The table exists because `cpl.gMB`/`gPE`/`gOM` are indexed by *absolute* solution number
and are therefore mostly zeros — reading `cpl.gOM(3,7)` cannot distinguish "no coupling"
from "never evaluated". `cpl.breakdown` lists only the pairs actually computed. Full
column reference: `README_calcGOM.md` §6.5.

A `cancelRatio` above 1.5 is also called out on the console during the run.

**Why this matters for the optimizer.** The active fitness is `g_OM² · Q_opt · L(λ)`, so
nothing in the objective penalizes a cancellation notch: the search will happily settle in
one if the optics improve enough to compensate. The stored `gMB_kHz` / `gPE_kHz` columns
are how you detect that after a run —

```matlab
T = store.loadEvals('oblongMaxdef_trial1');
R = (abs(T.gMB) + abs(T.gPE)) ./ abs(T.gOM);   % cancellation ratio per evaluation
[~, k] = max(T.fitness);
fprintf('best design sits at cancelRatio = %.1f\n', R(k));
```

A large `R` at the optimum means the design is fragile to fabrication error even though
its nominal `g_OM` looks fine. If that turns out to matter, `R` is a natural constraint to
add — see `README_calcGOM.md` §6.4 for the physics.

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
| `OPT.fitnessFcn` | `fitness_coopProduct` | gated score — Nelder–Mead |
| `OPT.fitnessFcnRaw` | same, `QmechMin = -Inf` | ungated score — the GP objective |

### Optimizer and persistence

| Knob | Default | Meaning |
|---|---|---|
| `OPT.optimizer` | `'bayesopt'` | `'bayesopt'` or `'neldermead'` |
| `OPT.bo.maxEvals` | `60` | **total** study size; seeded points count against it |
| `OPT.bo.numSeedPoints` | `8` | random seed points, reduced by what the store supplies |
| `OPT.bo.acquisition` | `'expected-improvement-plus'` | acquisition function |
| `OPT.bo.deterministic` | `false` | fits a noise term — correct with adaptive meshing; see [why the objective is `-log10(fitness)`](#why-the-objective-is--log10fitness) |
| `OPT.bo.plotFcn` | `{@plotMinObjective, @plotObjectiveModel}` | live `bayesopt` plots |
| `P.minSidewallGap` | 148 nm | lithography: beam edge to hole edge |
| `P.minHoleGap` | 50 nm | lithography: gap between adjacent holes |
| `P.minFeature` | 50 nm | lithography: smallest hole dimension |
| `OPT.useStore` | `1` | 0 disables SQLite persistence entirely |
| `OPT.runId` | `'oblongMaxdef_trial1'` | names the study — **keep it to resume, change it to start fresh** |
| `OPT.dbPath` | `./test/1D_OMC_hole/optim_runs.sqlite3` | shared database across studies |
| `OPT.rootLoc` | `./test/1D_OMC_hole/optimize_<runId>/` | output folder, derived from `OPT.runId` |
| `OPT.resume` | `1` | 0 re-solves everything (database still written) |
| `optimset(...)` | `TolX 1e-4`, `TolFun 1e-3`, `MaxIter 100`, `MaxFunEvals 200` | **Nelder–Mead only** — simplex stopping rules |

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
A cache hit does not re-solve, so it does not re-create its output folder. Replayed
points still resolve correctly (that comes from the database), but the folder their mode
profiles were written to is gone.

**This script and `sweep_oblong_maxdef.m` use different device geometry** — `w = 750 nm`,
`hy = 578 nm` here versus `w = 800 nm`, `hy = 651 nm` in the sweep. Their results are
**not** directly comparable, and sweep points cannot be reused to seed this optimizer
without re-solving.

**The optics-unavailable branch of the fitness returns `g_OM²` alone**, which is not on
the same scale as the full expression. With `P.solveOpt = 1` this should never trigger,
but a run that mixes the two branches cannot be ranked meaningfully.


---

## Testing

None of these need COMSOL. `test_bayesopt_wiring` needs the Statistics and Machine
Learning Toolbox; the rest run on base MATLAB.

```matlab
OptimStore.selfTest
% Exercises the active backend: schema, insert, lookup, best, NULL round-trip,
% reopen, JSONL mirror.

test_OptimStore_resume
% Runs three optimizer sessions against an analytic fitness and asserts that
% session 2 replays session 1's path *identically* (to 1e-12) from the database
% before extending it, and that a fresh runId shares nothing.

test_isFabricable
% Cross-validates the pre-filter against CreateNanobeamGeom over two grids,
% asserting ZERO false positives (accepting a design the builder rejects,
% which would buy a failed solve) and zero false negatives. The wide-region
% scan is asserted to reject something, so the comparison cannot go vacuous.

test_bayesopt_wiring    % needs the Statistics and ML Toolbox
% Drives bayesopt against an analytic fitness and a synthetic Q_mech field,
% asserting: both fitness flavours persist; gated rows keep a finite ungated
% fitness; the Q_mech constraint is positive exactly on the gated designs;
% seedBayesFromStore returns data bayesopt accepts; and a Nelder-Mead cache
% read can re-derive its own J from rows a Bayesian run wrote.
```

Expected output of the resume test:

```
session 1 : 13 solved,  0 cached, best fitness 9.97556e+15
session 2 : 18 solved, 13 cached, best fitness 9.99996e+15
            replayed 13 point(s) identically, then extended
session 3 :  8 solved,  0 cached (fresh runId)
ALL RESUME ASSERTIONS PASSED
```

of the pre-filter test:

```
nominal design: ok=1, margin=24.0 nm
search box   21x21 : filter 100.0% ok, builder 100.0% ok, agreement 100.0%
wide region  25x25 : filter  86.9% ok, builder  86.9% ok, agreement 100.0%
ALL ISFABRICABLE ASSERTIONS PASSED
```

and of the bayesopt wiring test:

```
session 1 : 14 evaluation(s), 14 objective row(s) in the study
           2 of 14 design(s) gated by Q_mech
           seed: 14 point(s), 14 objective(s), 14x2 constraint(s)
session 2 : 6 new evaluation(s), study holds 20 point(s)
           NM re-derives J = -13.4 from a BO-written row
ALL BAYESOPT WIRING ASSERTIONS PASSED
```

The identical-path assertion is the important one: if the replay diverged, the
fast-forward would be serving stored answers for a trajectory the optimizer is no longer
on.

Static analysis — `checkcode` on every file in this directory:

| File | Result |
|---|---|
| `isFabricable.m` | clean |
| `seedBayesFromStore.m` | clean |
| `test_isFabricable.m`, `test_bayesopt_wiring.m`, `test_OptimStore_resume.m` | clean |
| `OptimStore.m` | one "value might be unused" — a defensive init in `selfTest` |
| `optimize_oblong_maxdef.m` | three "function might be unused" for the intentionally-kept alternative fitness functions, one defensive `status_i` init |

---

## Related files

| File | Relationship |
|---|---|
| `sweep_oblong_maxdef.m` | grid sweep over the same two parameters (different geometry) |
| `RunNanobeamFEM.m` | the per-evaluation pipeline: geometry → mesh → solve → post-process |
| `CreateNanobeamGeom.m` | builds the hole arrays; lines 72-84 are the lithography guards `isFabricable` mirrors |
| `SolveNanobeamFEM.m` | source of `ds.mfem.QAll`, `ds.ofem.Q` |
| `CalcGOM.m` | computes `ds.cpl` — `gMax`, `gMBmax`, `gPEmax`, `breakdown` |
| `README_calcGOM.md` | **the reference for the coupling physics** — equations, signs, limitations |
| `../bayesopt_boomerang.m`, `../bayesopt_README.md` | MATLAB `bayesopt` reference implementation |
| `../optimize_hole_unitCell.m` | Nelder–Mead optimizer for hole *unit cells* |
| `../../python-scripts/src/database.py` | the SQLite pattern `OptimStore` follows |
| `../CLAUDE.md` | repo conventions — `P` struct, path handling, file naming |
