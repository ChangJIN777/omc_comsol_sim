# `bayesopt_boomerang.m`

Bayesian optimization of the boomerang unit cell for a complete **mechanical** bandgap centred on a target frequency.

Searches four geometry parameters — lattice constant `a`, hole arm length `r`, hole arm width `w`, and slab thickness `th` — using MATLAB's `bayesopt`, where each objective evaluation is a full COMSOL eigenfrequency band structure solve.

> **Status: written but never executed.** No MATLAB or COMSOL was available in the environment where this was authored, so it has not been run and `checkcode` has not been passed over it. It has been reviewed by hand and structurally checked (all nine local functions balance), but treat the first run as a shakedown. Run `checkcode('bayesopt_boomerang.m')` before committing to a long study.

---

## Contents

- [Why bayesopt](#why-bayesopt)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Read this before your first real run](#read-this-before-your-first-real-run)
- [Unit cell geometry](#unit-cell-geometry)
- [The objective](#the-objective)
- [Constraints](#constraints)
- [Configuration reference](#configuration-reference)
- [Outputs](#outputs)
- [Resuming an interrupted study](#resuming-an-interrupted-study)
- [How caching works](#how-caching-works)
- [`solveBands` quirks this script works around](#solvebands-quirks-this-script-works-around)
- [Local function reference](#local-function-reference)
- [Tuning guidance](#tuning-guidance)
- [Related files](#related-files)

---

## Why bayesopt

Each objective evaluation is a multi-minute COMSOL solve, so the **evaluation budget** — not the convergence rate — is the binding constraint. `bayesopt` fits a Gaussian-process surrogate to the objective and picks each next design point where it is most informative.

Compared with the `fminsearch` approach in `optimization_client.m`:

| | `fminsearch` | `bayesopt` |
|---|---|---|
| Spends solves informatively | No — fixed simplex moves | Yes, via GP + acquisition function |
| Integer design variables | No | Yes, natively |
| Skips infeasible points | No — must evaluate, then penalise | Yes, `XConstraintFcn` prunes before solving |
| Resumable after a crash | No | Yes, via `resume` |
| Handles flat/zero regions | Poorly — a local simplex method stalls | Yes, the GP models them |

That last row matters here: the objective is exactly zero everywhere no complete gap exists, and a simplex crawling across such a plateau has no gradient information to work with.

---

## Requirements

- **Statistics and Machine Learning Toolbox** — for `bayesopt` and `optimizableVariable`. There is no base-MATLAB fallback; without it the script cannot run.
- **COMSOL with LiveLink for MATLAB, already running** before you start the script. `solveBands` → `runBands_2D` calls the COMSOL Java API directly.
- MATLAB R2017b or newer (`isfile`).
- Must be run from the `omc-comsol-chang/` directory, so `solveBands`, `runBands_2D`, `buildBoomerangUnitCell`, `findGaps`, and `LoadMaterialParams` are all on the path.

---

## Quick start

```matlab
% 1. Start COMSOL with LiveLink for MATLAB.
% 2. Sanity-check the script parses cleanly:
checkcode('bayesopt_boomerang.m')

% 3. Edit the cfg block near the top — at minimum cfg.targetFreq,
%    cfg.bounds.*, and cfg.maxEvaluations. See the caveat below.

% 4. Run:
run('bayesopt_boomerang.m')
```

The script prints a per-evaluation summary, keeps two live `bayesopt` plots open, and leaves `results` (a `BayesianOptimization` object) in the base workspace.

**Before spending a long run, time a single solve** with `test_Boomerang.m` at the same `meshSize` and `nbands`. Default settings cost `cfg.maxEvaluations × 2` COMSOL band-structure runs — 40 evaluations is **80 runs**, since each design point solves two symmetry sectors (even and odd about z). Multiply your single-solve time by that.

---

## Read this before your first real run

**The default bounds do not contain your current design point.** `cfg.bounds` was taken from the grid ranges in `sweep_boomerang_code.m`, which predate the geometry now in `test_Boomerang.m` (`a=480, w=140, r=177, th=220` nm):

| variable | `cfg.bounds` default | `test_Boomerang.m` | |
|---|---|---|---|
| `a` | 600 – 1000 nm | 480 | **below the range** |
| `w` | 50 – 120 nm | 140 | **above the range** |
| `th` | 275 – 350 nm | 220 | **below the range** |
| `r` | 150 – 250 nm | 177 | inside ✓ |

So as shipped, the optimizer cannot reach — or even model — the region you are presumably most interested in. Decide deliberately which design space you want searched. To centre it on the current geometry:

```matlab
cfg.bounds.a  = [400  600];
cfg.bounds.r  = [140  210];
cfg.bounds.w  = [100  180];
cfg.bounds.th = [180  260];
```

Whatever you choose, keep `r < sqrt(3)*a/4` reachable across the box or a large fraction of candidates will be pruned as infeasible (see [Constraints](#constraints)). That ceiling is 173 nm at `a`=400 nm and 260 nm at `a`=600 nm.

---

## Unit cell geometry

Built by `buildBoomerangUnitCell.m` — the builder `runBands_2D` dispatches for `celltype = 'boomerang'`. (Note `buildBoomerangUnitCell_2D.m` exists but is *not* on this path.)

**The cell** is the primitive cell of a hexagonal (triangular) lattice: a rhombus of side `a` with a 60° interior angle, spanned by lattice vectors `a1 = (a, 0)` and `a2 = (a/2, a·√3/2)`. The cell centre sits at `(3a/4, a·√3/4)`. **`a` alone sets the cell footprint.**

**The hole** is a three-pointed star: three rectangular arms, each `w` wide and `r` long, radiating from the cell centre at 120° spacing (θ = 90°, 210°, 330°). This tri-arm shape is what gives the cell its name.

```
              arm 1 (90 deg)
                   |
                   |
                   O   <- cell centre, (3a/4, a*sqrt(3)/4)
                  / \
                 /   \
          arm 2 /     \ arm 3
       (210 deg)       (330 deg)
```

| `P` field | What it actually is |
|---|---|
| `a` | Lattice constant; side of the rhombic cell |
| `w` | Hole **arm width** — the narrowest etched feature |
| `r` | Hole **arm length**, cell centre to arm tip |
| `th` | **Full** slab thickness in z |
| `r1` | Fillet radius at the three **inner** corners, where the arms meet near the centre |
| `r2` | Fillet radius at the three **outer** arm tips |

⚠️ Comments in several other scripts in this directory label `w` as "unit cell width (along x)" and `r` as "unit cell height (along y)". **Both are wrong** — they are hole dimensions, not cell dimensions. Inherited `r1`/`r2` comments describing a "cross" or "inner block" belong to the `'cross'`/`'hollow'`/`'solid'` cell types and do not apply here. `test_Boomerang.m` and `sweep_boomerang_code.m` have since been corrected.

On `th` and symmetry: because `P.mbevenz` is nonzero, the builder subtracts the lower half of the slab and meshes only `z ∈ [0, th/2]` against a mirror plane at `z = 0`. `th` remains the full physical thickness; the mesh covers half of it.

---

## The objective

`bayesopt` minimizes, so the script returns a negated fitness:

```
objective = -max over all complete gaps of
              (gapSize / midGap) * exp( -((f_target - midGap) / sigma)^2 )
```

Two factors: the **fractional gap width** `gapSize/midGap`, which is the thing you actually want maximised, and a **Gaussian penalty** pulling the mid-gap frequency toward `cfg.targetFreq` with width `cfg.sigma`.

Every complete gap returned by `findGaps` is scored and the best is kept. This deliberately replaces the gap-*selection* step in `optimization_client.m`, which had two defects:

```matlab
gap_ind = find(abs(target_freq-mFreqs) < mGaps);          % can match SEVERAL gaps
if isempty(gap_ind)
    gap_ind = find(min(abs(target_freq-mFreqs) - mGaps/2));  % a no-op, not an argmin
end
```

The first line can match more than one gap, leaving `mGap`/`mFreq` as vectors and the fitness a vector. The second looks like an argmin but is not: `find(min(v))` evaluates the scalar minimum and returns `1` if it is nonzero, empty if it is exactly zero — so it always picks gap 1 regardless of which gap is nearest the target. Scoring all gaps and taking the max removes the need to select at all.

### Return values and what they mean

| value | meaning |
|---|---|
| negative | A complete gap was found; more negative is better |
| `0` | Valid solve, but **no usable complete gap** — the worst real value, since fitness is non-negative by construction |
| `NaN` | A genuine solver **failure** (exception from `solveBands`). Returned deliberately so `bayesopt` records an error point rather than learning the region is merely mediocre |

The `status` field distinguishes `target-in-gap` (the target frequency falls inside the winning gap) from `nearest-gap` (it does not), mirroring the message `optimization_client.m` printed.

---

## Constraints

Both are deterministic and enforced through `XConstraintFcn`, so infeasible candidates are pruned **before** any COMSOL time is spent — as opposed to being evaluated and then penalised, which also distorts the surrogate.

1. **Minimum feature size:** `w >= cfg.minFeatureNm` (default 50 nm). For the tri-arm hole the narrowest feature is simply the arm width, since that is an etched slit.

2. **Arms stay inside the cell:** `r < sqrt(3)*a/4`. The cell centre sits `√3·a/4` from the nearest cell edge and an arm reaches `r`, so violating this punches the hole through the cell boundary and produces a shape the periodic BCs no longer describe. The default bounds clear this by only ~10 nm at `a`=600, `r`=250, so it matters as soon as ranges widen.

`boomerangFabConstraint` is **vectorized over a multi-row table**, as `bayesopt` requires, and compares in integer nm so a design sitting exactly on a limit is not decided by floating-point round-off.

### On `a - r`

The `min(a - r, ...)` fabrication proxy used in `boomerang_optimize_sweep_diamond.m` is **deliberately not reproduced here** — it has no geometric meaning for this shape. Sampling the tri-arm outline against its periodic images gives a true inter-hole diamond ligament of ~161 nm for the current design, where `a - r` reports 303 nm; it overstates by roughly 2×. (`√3a/4 - r` ≈ 31 nm is the other tempting expression, but that is only the gap to the cell *boundary* — material continues past it into the neighbouring cell, so it is not a feature either.) Over the default bounds the true ligament never drops below ~210 nm, so an `a - r` term would neither bind nor measure anything real.

---

## Configuration reference

All settings live in the `cfg` struct at the top of the script.

### Objective

| field | default | meaning |
|---|---|---|
| `cfg.targetFreq` | `13e9` | Target mechanical mid-gap frequency [Hz] |
| `cfg.sigma` | `5e9` | Width of the Gaussian frequency penalty [Hz] |

### Fabrication

| field | default | meaning |
|---|---|---|
| `cfg.minFeatureNm` | `50` | Minimum arm width [nm]. Integer nm so the comparison is exact |

### Search space (integer nm)

| field | default |
|---|---|
| `cfg.bounds.a` | `[600 1000]` |
| `cfg.bounds.r` | `[150 250]` |
| `cfg.bounds.w` | `[50 120]` |
| `cfg.bounds.th` | `[275 350]` |

See [the caveat above](#read-this-before-your-first-real-run) — these do not contain the current design point.

### Fixed geometry

| field | default | meaning |
|---|---|---|
| `cfg.r1` | `10e-9` | Fillet radius, inner corners (arms meet at centre) |
| `cfg.r2` | `10e-9` | Fillet radius, outer arm tips |

### Solver fidelity

Kept identical to `sweep_boomerang_code.m` so results are comparable with the existing sweeps.

| field | default | meaning |
|---|---|---|
| `cfg.kpts` | `9` | k-points **excluding** Γ. The solve covers `3*kpts` points on the Γ–M–K–Γ circuit |
| `cfg.nbands` | `15` | Bands per k-point. Must be high enough to bracket the gap you care about |
| `cfg.meshSize` | `4` | COMSOL `autoMeshSize` level — **higher means coarser.** `test_Boomerang.m` uses 3, which is *finer*, so frequencies are not directly comparable between the two |
| `cfg.maxDof` | `3e6` | DOF cap. `runBands_2D` coarsens the mesh in a loop until the estimate falls below this, so it bounds solve time rather than erroring |

### Optimizer budget

| field | default | meaning |
|---|---|---|
| `cfg.maxEvaluations` | `40` | Total design points `bayesopt` may evaluate (× 2 sectors = COMSOL runs) |
| `cfg.numSeedPoints` | `8` | Random points evaluated before the GP takes over |

Other `bayesopt` options are set inline: `'expected-improvement-plus'` acquisition, and `'IsObjectiveDeterministic', true` (a COMSOL solve is repeatable, so the GP need not model observation noise).

### I/O and plotting

| field | default | meaning |
|---|---|---|
| `cfg.datLoc` | `.\test\boomerang_bayesopt\<mmddyyyy>\` | Output folder. **The trailing separator is required** |
| `cfg.plotgeom` | `0` | 1 to plot geometry every evaluation (slow, noisy) |
| `cfg.savebndplot` | `1` | 1 to save a band diagram per evaluation |
| `cfg.closeSolveFigures` | `1` | Close figures the solve opened, keeping the `bayesopt` live plots alive |

That last option exists because a plain `close all` would also destroy the `bayesopt` plot windows. Instead the figure list is differenced against a snapshot taken before each solve, so only figures the solve created get closed.

---

## Outputs

Everything lands under `cfg.datLoc`.

| file | contents |
|---|---|
| `<fileBase>_bds.mat` | Per-evaluation band structure `ds` struct, written by `solveBands`. Also the cache — see below |
| `<fileBase>_fullBands.png` / `.fig` | Band diagram per evaluation, when `cfg.savebndplot = 1` |
| `bayesopt_boomerang_log.txt` | One tab-separated row per evaluation |
| `bayesopt_boomerang_state.mat` | Checkpoint of the `BayesianOptimization` object, rewritten every iteration |
| `bayesopt_boomerang_results.mat` | Final `results`, `cfg`, and `xBest` |

`<fileBase>` is e.g. `boomerang_a_700nm_r_262nmw_99nm_th_300nm_r1_10nm_r2_10nm_`. The script sets `P.fileBase` explicitly, reproducing byte-for-byte the string `solveBands` builds for `celltype = 'boomerang'` — including its missing `_` after the `r` field. That is intentional: it lets the script predict the cache path, and keeps filenames interchangeable with data already produced by the sweep scripts.

### Log columns

`bayesopt_boomerang_log.txt` is written and closed once per evaluation rather than buffered, so a run killed mid-study still leaves a complete record.

| column | notes |
|---|---|
| `a_nm`, `r_nm`, `w_nm`, `th_nm` | Design point |
| `midGap_GHz`, `gapSize_GHz` | Winning gap |
| `gapRatio` | `gapSize/midGap` |
| `penalty` | Gaussian frequency penalty, 0–1 |
| `fitness` | `gapRatio * penalty` |
| `objective` | `-fitness` — what `bayesopt` sees |
| `cached` | 1 if the band structure was reloaded from an existing `.mat` |
| `evalTime_min` | Wall-clock minutes for this evaluation |
| `status` | `target-in-gap`, `nearest-gap`, `no-complete-gap`, or `solver-error` |

Load it with `readtable('bayesopt_boomerang_log.txt')`.

---

## Resuming an interrupted study

The `BayesianOptimization` object is checkpointed after **every** iteration, because losing a partly finished study to a dropped LiveLink connection is expensive at these evaluation costs.

```matlab
S = load('.\test\boomerang_bayesopt\<date>\bayesopt_boomerang_state.mat');
results = resume(S.results, 'MaxObjectiveEvaluations', 20);
```

To extend a study that finished cleanly, `resume` the in-memory object directly:

```matlab
results = resume(results, 'MaxObjectiveEvaluations', 20);
```

Note that `resume` continues with the surrogate already fitted — it does not re-spend seed points.

---

## How caching works

`solveBands` **skips the solve entirely** when `[datLoc fileBase '_bds.mat']` already exists. This script leans on that deliberately:

- Re-running an interrupted study re-uses everything already solved.
- Existing sweep data under the same `datLoc` is picked up for free, since the filenames match.

**This is why the design variables are declared as integer nanometres rather than continuous metres.** `solveBands` names its files by rounding every dimension to whole nm (`'%.0f'`). With continuous variables, many distinct GP query points would alias onto one cached result, silently feeding the surrogate stale data for design points it never actually evaluated. On an integer-nm grid the rounding is exact, so a repeated point is a genuine cache hit.

The flip side: **editing physics settings without changing geometry will silently reuse old results.** Clear the dated folder — or point `cfg.datLoc` somewhere new — whenever you change `kpts`, `nbands`, `meshSize`, or the symmetry settings.

---

## `solveBands` quirks this script works around

Recorded here because they are easy to trip over again.

1. **Two different return shapes.** A fresh solve returns `ds` with a `.full` field, but a cache hit prints `Data folder exists in working directory` and returns a stub with only empty `.sym`/`.asym`. Reading `ds.full` off that stub throws `Unrecognized field name` — which over an optimization run, where revisiting a design point is routine, would abort the whole study. `solveBoomerangBands` detects this and reloads the cached `.mat` instead.

2. **`P.datLoc` needs its trailing separator.** `solveBands` normalises a separator into a local variable but then writes the `.mat` using the raw `P.datLoc`, so a missing separator drops output into the parent folder under a mangled name.

3. **`P.savedat` must stay `1`.** It is what writes the `.mat`, and therefore what makes caching work at all.

4. **`TwoSymPlanes` is set to `0`, deliberately.** On the 2D path `runBands_2D` applies a Floquet periodic BC in y and reads `P.mbeveny` without ever using it, so only `mbevenz` changes the problem posed. With `TwoSymPlanes = 1`, `solveBands` runs four sectors that collapse to two distinct solves run twice each — double the cost, and `full.F` then holds every band twice, making `findGaps` report every complete gap twice over. Two sectors (even/odd about z) is already the complete set. Note `test_Boomerang.m` currently sets `1`.

5. **`P.bandStruct_2D` must stay `1`.** `runBands` (the `0` path) has no `'boomerang'` branch in its celltype dispatch — it handles only `boomerang_strip_v2` and `boomerang_lower`.

---

## Local function reference

| function | role |
|---|---|
| `boomerangFabConstraint` | Vectorized deterministic feasibility test on a candidate table |
| `boomerangObjective` | Solve one design point and score it; returns `-fitness` |
| `boomerangParams` | Build the `P` struct for one design point, ported from `sweep_boomerang_thickness` |
| `solveBoomerangBands` | Call `solveBands` and always return a usable gap struct, handling the cache-hit stub |
| `closeNewFigures` | Close only the figures a solve opened, preserving `bayesopt` plots |
| `gapFitness` | Score every complete gap against the target, keep the best |
| `initIterationLog` | Create the log with its header row |
| `logIteration` | Append one evaluation, flushed immediately |
| `checkpointState` | `OutputFcn` saving the optimizer state every iteration |

---

## Tuning guidance

**`cfg.sigma`** controls how hard the frequency target is enforced. Too small and almost every design scores ~0, leaving the GP nothing to learn from; too large and the optimizer happily returns a wide gap at the wrong frequency. If early evaluations come back with `penalty` near zero across the board, widen it.

**`cfg.numSeedPoints`** should stay a reasonable fraction of `cfg.maxEvaluations` — with 4 variables, 8 seeds out of 40 is a sensible split. Too few and the GP extrapolates from almost nothing; too many and you spend the budget on random search.

**`cfg.nbands`** must be large enough to bracket the gap of interest. If gaps go missing at higher frequency, raise it — but note that COMSOL solve cost grows with band count.

**`cfg.maxEvaluations`** is a hard stop, not a convergence criterion. Prefer a modest first run followed by `resume`, since that lets you inspect the log before committing more solve time.

---

## Related files

| file | relationship |
|---|---|
| `test_Boomerang.m` | Single-point solve; the authoritative unit cell geometry write-up |
| `sweep_boomerang_code.m` | Grid sweeps of `a` / `th` / `(w,r)`; source of the default bounds and solver fidelity |
| `optimization_client.m` | The original `fminsearch` optimizer; source of the objective. Note it calls `solveBands(vars,P)`, a stale two-argument signature — the current one is `solveBands(P)` |
| `boomerang_optimize_sweep_diamond.m` | Combined optical + mechanical fitness; source of the `a - r` fab proxy not reproduced here |
| `solveBands.m` | Symmetry-sector driver, calls `findGaps` |
| `runBands_2D.m` | 2D band structure solve; builds the COMSOL model |
| `buildBoomerangUnitCell.m` | The geometry builder actually used for `celltype = 'boomerang'` |
| `findGaps.m` | Locates complete gaps in an assembled band matrix |
