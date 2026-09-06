# Proposal: Bayesian optimization for `optimize_oblong_maxdef.m`

Port the Bayesian-optimization pattern already used elsewhere in this repo onto the
nanobeam cavity taper search, **keeping Nelder–Mead selectable** so the two can be
compared on the same objective, the same store, and the same fitness function.

> **Status: phases 0 and 3 are implemented.** The Statistics and Machine Learning Toolbox
> was installed on 2026-09-06, and `OPT.optimizer = 'bayesopt'` is now the default. This
> document is kept as the design rationale; `optimization_README.md` is the user guide.

| Piece | State |
|---|---|
| SQLite stop/resume infrastructure | **Built and tested** (`OptimStore.m`, `test_OptimStore_resume.m`) |
| `OPT.optimizer` switch | **Live both ways**; `'bayesopt'` is the default |
| Bayesian backend | **Built and tested** (`test_bayesopt_wiring.m`) — coupled constraints, store seeding |
| Feasibility pre-filter | **Still proposed** (phase 1) — not written |

---

## 1. What already exists (the recall)

There is no Bayesian optimizer for hole-shaped unit cells in the MATLAB tree. What
exists is three separate things that together are the template:

**`omc-comsol-chang/bayesopt_boomerang.m`** (2746 lines, + a 751-line README) — the
only MATLAB `bayesopt` implementation here. It optimizes **boomerang** unit cells for
a complete mechanical bandgap over three *integer* variables (`a`, `r`, `w`), and it
is where the mechanics to copy live:

- `optimizableVariable` array + `bayesopt(...)` with `expected-improvement-plus`
- `XConstraintFcn` — prunes infeasible candidates **before** any COMSOL solve
- `OutputFcn` → `checkpointState`, saving the `BayesianOptimization` object every
  iteration so `resume(results, 'MaxObjectiveEvaluations', n)` can extend it
- Three solver backends selected by string (`'comsol' | 'surrogate' | 'stub'`) so the
  whole loop is testable in seconds without COMSOL
- A provenance guard so a synthetic dry-run result can never be mistaken for a real one

> Its README flags the file as **not re-verified** since the three-variable and
> constraint changes — treat it as a design reference, not as known-good code.

**`omc-comsol-chang/optimize_hole_unitCell.m`** (243 lines) — *is* the hole-cell
optimizer, but it is **Nelder–Mead**, not Bayesian: a 1-D area-preserving search over
the aspect ratio `r = hx/hy` at fixed ellipse area `A0`.

**`python-scripts/`** — the closest thing to "Bayesian optimization for hole-shaped
unit cells", and the source of the persistence design:

- `src/optimizer.py` — Optuna **TPE** (a Bayesian method) over 5 normalized variables
  `[a, w, hx, hy, t]` of a rectangular nanobeam with one **elliptical hole** per period
- `src/database.py` — every candidate written to SQLite immediately, with an
  append-only JSONL fallback for filesystems that cannot do SQLite locking
- `_replay()` — rebuilds the sampler from completed rows; that *is* the resume path
- Backend tiering `auto → optuna → random`, so a missing dependency degrades instead
  of failing

The proposal below is: **MATLAB `bayesopt` mechanics from the boomerang script, resume
semantics from `database.py`, applied to the oblong/maxdef taper search.**

---

## 2. Why Bayesian optimization fits this problem

The binding constraint is the **evaluation budget**, not the convergence rate. Each
evaluation is a COMSOL mechanical + optical + PML solve — minutes to tens of minutes.
The current `MaxFunEvals = 200` is a multi-day run.

Three properties of *this* objective make Nelder–Mead a poor match, and all three are
things a GP surrogate handles natively:

1. **The binary `Q_mech` gate makes the landscape piecewise-flat.** Every design below
   `OPT.QmechMin = 1e7` scores exactly 0. A simplex sitting inside that region has no
   direction to move in. A GP models the plateau and the boundary and proposes points
   *at the edge* of it, where the information is.
2. **Infeasible points are currently paid for in full.** A candidate that violates the
   lithography gap check throws from `CreateNanobeamGeom.m:77` — but only after the
   geometry build starts. `XConstraintFcn` rejects it with zero solver cost.
3. **The objective is mildly noisy** (adaptive meshing, `P.mAdjMesh`/`P.oAdjMesh`).
   Nelder–Mead has no noise model and will contract onto a mesh artifact and call it
   convergence. A GP with a noise term will not.

Against that, one honest caveat: the search is only **2-D**. Bayesian optimization
earns its keep most clearly in 5–20 dimensions. In 2-D with a generous budget, a
coarse `sweep_oblong_maxdef` grid plus a short local polish is genuinely competitive.
The case for BO here rests on the *budget* and the *plateau*, not on dimensionality —
which is exactly why the switch should exist rather than a wholesale replacement.

---

## 3. Toolbox reality on this machine (checked, not assumed)

**Resolved.** As of 2026-09-06 the Statistics and Machine Learning Toolbox is installed
and `bayesopt` is callable:

```
exist('bayesopt')                                = 2     % callable
license('test','Statistics_Toolbox')             = 1     % licensed
isfolder(fullfile(matlabroot,'toolbox','stats')) = 1     % installed
```

It had previously been *licensed but not installed*, which surfaces as a bare
`Unrecognized function 'bayesopt'` and reads like a licence problem. The dispatch still
distinguishes the two cases explicitly, so the same confusion cannot recur on another
machine. MATLAB here is **R2025b**.

**Database Toolbox remains licensed but not installed**, which is why `OptimStore` does
not depend on it — see below.

---

## 4. What is already built

### `OptimStore.m` — SQLite persistence with stop/resume

One row per evaluation, written the moment it completes.

```
runs   (run_id PK, created_ts, updated_ts, script, optimizer, status, config_json)
evals  (run_id, eval_idx, ts, param_key, dar, maxdef, oblong,
        wm_hz, gom_hz, gmb_hz, gpe_hz, lsiv_hz, q_mech, q_opt,
        lambda_nm, fitness, fitness_raw, objective_j, status, ev_loc,
        PRIMARY KEY (run_id, eval_idx))
INDEX  idx_evals_lookup ON evals(run_id, param_key)
```

**Backends, auto-detected in order** — because Database Toolbox is not installed:

| Backend | Requires | Status here |
|---|---|---|
| `dbtoolbox` | Database Toolbox `sqlite()` | not installed |
| `python` | MATLAB Python interface (`pyenv`) | **active** — Python 3.13, SQLite 3.45.3 |
| `cli` | `sqlite3` on PATH | available (Anaconda), **tested** |

An append-only `<db>.jsonl` mirror is written *before* the SQLite insert, so a failed
database write can never discard a solve you have already paid hours for. That is
`database.py`'s rule, carried over.

### How resume works, and why it is sound

`fminsearch` is **deterministic** given `x0` and `options`. So a resumed session
retraces its earlier path exactly before extending it. `objective()` therefore checks
the store first:

```matlab
hit = store.lookup(OPT.runId, [defectAspectRatio_i, maxdef_i]);
if ~isempty(hit)
    J = hit.J;              % no COMSOL call at all
    return;
end
```

Matching is on `OptimStore.paramKeyOf` — `sprintf('%.10g|%.10g', ...)` — rounded so a
value that has round-tripped through SQLite's `REAL` storage still matches.

**Verified** by `test_OptimStore_resume.m` (analytic fitness, no COMSOL, ~1 s):

```
session 1 : 13 solved,  0 cached, best fitness 9.97556e+15
session 2 : 18 solved, 13 cached, best fitness 9.99996e+15
            replayed 13 point(s) identically, then extended
session 3 :  8 solved,  0 cached (fresh runId)
```

The test asserts the replayed path is identical to the original to 1e-12 — if it were
not, the fast-forward would be serving answers for a trajectory the optimizer is no
longer on.

### Usage

```matlab
% Start a study
OPT.runId  = 'oblongMaxdef_trial1';
OPT.resume = 1;
optimize_oblong_maxdef

% ... interrupted (Ctrl-C, LiveLink crash, reboot). The run is marked
% 'interrupted' automatically by an onCleanup guard.

% Continue it: same runId, bigger budget. Prior points replay for free.
optimize_oblong_maxdef

% Inspect from the command line
store = OptimStore(fullfile('.','test','1D_OMC_hole','optim_runs.sqlite3'));
store.listRuns()
T = store.loadEvals('oblongMaxdef_trial1');
[x, fit] = store.best('oblongMaxdef_trial1');
```

---

## 5. Proposed Bayesian backend

### 5.1 Variables

```matlab
optVars = [ ...
    optimizableVariable('dAR',    [OPT.defectAspectRatio_min, OPT.defectAspectRatio_max], 'Type','real'), ...
    optimizableVariable('maxdef', [OPT.maxdef_min,            OPT.maxdef_max],            'Type','real')];
```

Both continuous, both already normalized to sensible ranges. No transform — neither
spans more than a factor of 3.

### 5.2 Replace the penalty ladder with coupled constraints

This is the one change that is **not** cosmetic, and the main reason the BO branch is
a proposal rather than a patch.

The current objective encodes three infeasible outcomes as large positive numbers:

```
gated (Q_mech < QmechMin) -> 1 * Jpenalty  = 1000
unusable (no mode/error)  -> 2 * Jpenalty  = 2000
out of bounds             -> Jpenalty*(2+d)
```

Feeding 1000 and 2000 into a GP that also sees feasible values around −16 destroys the
kernel's length scale — the surrogate spends its capacity modelling the penalty cliff
instead of the physics. `bayesopt` has first-class machinery for this:

```matlab
function [objective, coupledConstraints] = bayesObjective(x, ...)
    % coupledConstraints < 0 means satisfied.
    %   c(1): Q_mech gate      -> 1 - Qmech/QmechMin
    %   c(2): usable solve     -> -1 on success, +1 on failure
    % objective is NaN when the point could not be evaluated; bayesopt models
    % those errors with a separate classifier rather than treating NaN as a value.
```

So:

| outcome | today (`fminsearch`) | proposed (`bayesopt`) |
|---|---|---|
| feasible | `-log10(fitness)` | `-log10(fitness)`, both constraints < 0 |
| `Q_mech` gated | `J = 1000` | real objective + `c(1) > 0` |
| solve failed | `J = 2000` | `objective = NaN`, `c(2) > 0` |
| out of bounds | sloped barrier | never proposed — the box is the variable range |

The GP then learns *where the gate boundary is* — which is the actual design question
— instead of learning that a wall exists.

### 5.3 Pre-solve feasibility filter (worth doing for **both** optimizers)

`CreateNanobeamGeom.m:77` rejects a geometry when

```
a*(1-maxdef) - hx*(1-maxdef)^(1-oblong) < 50 nm
```

That check is pure arithmetic on `(a, hx, maxdef, oblong)` — no COMSOL. Lifting it
into a shared `isFabricable(dAR, maxdef, P)` helper lets:

- `bayesopt` use it as `XConstraintFcn` (candidate never proposed), and
- `fminsearch` use it as an instant penalty (no wasted FEM attempt).

Cheapest win in this whole document, and it is optimizer-independent.

### 5.4 Determinism setting

`bayesopt_boomerang.m` sets `'IsObjectiveDeterministic', true`. **Recommend `false`
here.** Adaptive meshing makes repeated evaluations at the same geometry differ
slightly; declaring determinism forces the GP to interpolate that jitter exactly and
produces a badly conditioned surrogate. `false` fits a noise term, which is what is
actually happening.

### 5.5 Resume for the BO branch

`OptimStore` already holds everything needed. Seeding is a table read, not a
`.mat` reload:

```matlab
T = store.loadEvals(OPT.runId);
ok = isfinite(T.J) & T.J < OPT.Jpenalty;
results = bayesopt(objFcn, optVars, ...
    'InitialX',         table(T.dAR(ok), T.maxdef(ok), 'VariableNames', {'dAR','maxdef'}), ...
    'InitialObjective', T.J(ok), ...
    'NumSeedPoints',    max(0, cfg.numSeedPoints - nnz(ok)), ...);
```

This is strictly better than the boomerang script's `.mat` checkpoint: it survives a
MATLAB version change, it is queryable with plain SQL, and **it lets a Nelder–Mead run
seed a Bayesian run** — which is the comparison worth doing.

### 5.6 Dry-run backend

Borrow `cfg.solverBackend = 'comsol' | 'surrogate' | 'stub'` from the boomerang script
before the first real study. The BO loop has considerably more machinery than the
simplex path (constraints, seeding, acquisition), and validating it against hours of
COMSOL is the wrong order of operations.

---

## 6. Phasing

| Phase | Work | Cost | Depends on |
|---|---|---|---|
| **0 — done** | `OptimStore`, resume, `OPT.optimizer` switch | — | — |
| **1** | `isFabricable` pre-filter, shared by both optimizers | small | — |
| **2** | Surrogate dry-run backend | small | — |
| ~~**3**~~ | ~~`bayesopt` branch~~ — **done**: `optimizableVariable` pair, coupled constraints, `seedBayesFromStore` | — | — |
| **4** | Head-to-head on equal budget, seeded from the same store | 1 study | 1–2 |

Phases 1 and 2 pay off on both optimizer paths and need no toolbox, so they remain
worth doing. Phase 1 is now the cheapest remaining win: lifting the lithography check
into `isFabricable` lets it become `bayesopt`'s `XConstraintFcn`, which prunes
unfabricable candidates *before* any COMSOL solve rather than paying for a failed one.

---

## 7. Decisions needed

1. ~~**Install the Statistics and Machine Learning Toolbox?**~~ **Done** (2026-09-06).
2. ~~**Coupled constraints, or keep the penalty ladder?**~~ **Resolved: constraints.**
   The ladder actively degrades a GP. Note the consequence — the two optimizers score
   infeasible points differently, so a feasible-point comparison stays apples-to-apples
   but "how many evaluations were wasted" does not.
3. *Still open.* **Is `Q_mech ≥ 1e7` a hard requirement or a preference?** If it is really a
   preference, a smooth `Q_mech` factor (`fitness_QxQxLambda`, already in the file)
   suits both optimizers better than the binary gate.
4. **Seed from `sweep_oblong_maxdef`?** A 5×5 grid is 25 evaluations of coverage that
   could initialize the GP instead of random seed points — but the sweep script uses a
   different device geometry (`w = 800 nm`, `hy = 651 nm` vs `750`/`578` here), so its
   points are not directly reusable without a re-run.

---

## 8. Related files

| File | Role |
|---|---|
| `nanobeam/optimize_oblong_maxdef.m` | the script being extended |
| `nanobeam/OptimStore.m` | SQLite persistence + resume |
| `nanobeam/test_OptimStore_resume.m` | COMSOL-free proof that resume works |
| `nanobeam/sweep_oblong_maxdef.m` | grid-sweep sibling (different geometry) |
| `bayesopt_boomerang.m`, `bayesopt_README.md` | MATLAB `bayesopt` reference |
| `optimize_hole_unitCell.m` | hole-cell Nelder–Mead optimizer |
| `python-scripts/src/optimizer.py` | Optuna TPE over elliptical-hole cells |
| `python-scripts/src/database.py` | the SQLite pattern `OptimStore` follows |
