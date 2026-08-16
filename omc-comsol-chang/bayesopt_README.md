# `bayesopt_boomerang.m`

Bayesian optimization of the boomerang unit cell for a complete **mechanical** bandgap centred on a target frequency.

Searches **three** geometry parameters — lattice constant `a`, hole arm length `r`, hole arm width `w` — using MATLAB's `bayesopt`, where each objective evaluation is a full COMSOL eigenfrequency band structure solve. Slab thickness `th` is **fixed** at `cfg.th`, and an optional **filling-factor constraint** restricts the search to designs whose air/dielectric area ratio sits in a band around `cfg.fillingFactor`.

> **Status: NOT re-verified since the three-variable / filling-factor change.** The end-to-end dry run described below was executed against the *four-variable* version of this file. Fixing `th`, adding the filling-factor constraint, replacing the containment test and adding the feasibility scan all postdate it, and **none of that has been run** — not even `checkcode`. Treat the next dry run as the verification, not as a formality.
>
> *Previously verified, and still expected to hold where untouched:* `checkcode` passed with no messages (MATLAB R2026a); a full dry run ran end-to-end with the surrogate backend, producing all four figures, the log, the checkpoint and the results `.mat`; cache isolation was demonstrated and the provenance guard tested by planting synthetic data at a real filename. COMSOL LiveLink has never been driven from this script, so `solveBandsViaComsol` is exercised only on its cache path.

---

## Contents

- [Why bayesopt](#why-bayesopt)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Dry-run mode](#dry-run-mode)
  - [The three backends](#the-three-backends)
  - [How cache isolation is enforced](#how-cache-isolation-is-enforced)
  - [How a synthetic result is labelled](#how-a-synthetic-result-is-labelled)
  - [What the surrogate actually computes](#what-the-surrogate-actually-computes)
- [Read this before your first real run](#read-this-before-your-first-real-run)
- [Unit cell geometry](#unit-cell-geometry)
- [The objective](#the-objective)
  - [The fitness function](#the-fitness-function)
  - [Symbols](#symbols)
  - [Properties worth knowing](#properties-worth-knowing)
- [Constraints](#constraints)
- [Configuration reference](#configuration-reference)
  - [Filling-factor constraint](#filling-factor-constraint)
- [Outputs](#outputs)
- [Cross-platform paths](#cross-platform-paths)
- [Visualization](#visualization)
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

- **Statistics and Machine Learning Toolbox** — for `bayesopt` and `optimizableVariable`. There is no base-MATLAB fallback; without it the script cannot run. **This is required for a [dry run](#dry-run-mode) too** — the dry run replaces the *solver*, not the optimizer.
  > If you get `Unrecognized function or variable 'optimizableVariable'` (or the same for `bayesopt`), the toolbox is not necessarily missing — check whether the **licence** simply failed to check out, which is intermittent on a shared or network licence:
  > ```matlab
  > license('test','Statistics_Toolbox')             % 1 = licensed
  > [ok,msg] = license('checkout','Statistics_Toolbox')
  > isfolder(fullfile(matlabroot,'toolbox','stats'))  % 1 = installed on disk
  > ```
  > All three returning 1 while the script still fails means a transient checkout failure: just re-run.
- **COMSOL with LiveLink for MATLAB, already running** before you start the script. `solveBands` → `runBands_2D` calls the COMSOL Java API directly. **Not** needed for a dry run.
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

Before that, though, do a **dry run** — the whole loop in seconds instead of hours:

```matlab
% Edit one line near the top of the cfg block:
cfg.solverBackend = 'surrogate';
run('bayesopt_boomerang.m')
% ... then set it back to 'comsol' for the real study.
```

---

## Dry-run mode

A mode that drives the **entire optimization loop with no COMSOL solves at all**, so the plumbing around the solve — the objective, `XConstraintFcn`, the iteration log, the checkpoint, the report and all four summary figures — can be validated in seconds rather than over a multi-hour study.

It is selected by an explicit **string**, `cfg.solverBackend`, not a boolean `dryRun` flag. That follows `mech_backend` in `python-scripts/src/objective.py` (`'comsol'` vs `'surrogate_stub'`): three tiers do not fit in a boolean, a fourth would not either, and a name self-documents everywhere it is printed, logged or saved instead of becoming an anonymous `1`.

An unrecognised name is an **error**, never a silent fall-through — the MATLAB counterpart of that file's `raise ValueError(f"unknown mech_backend {mech_backend}")`. There is no safe default: falling back to `'comsol'` would spend hours of solve time on a typo, and falling back to `'surrogate'` would hand back a folder of fabricated results that look exactly like a study. The check runs **twice** — once at configuration time, before any path is derived from the name, and again in the backend dispatcher, for a `cfg` that was hand-edited or reloaded from an old results `.mat` afterwards.

```matlab
cfg.solverBackend = 'surrogate';    % or 'comsol' (default), or 'stub'
cfg.dryRunDelay     = 0;            % artificial seconds per evaluation
cfg.dryRunFailEvery = 0;            % >0 injects deterministic pseudo-failures
cfg.dryRunSaveBands = 1;            % 0 = write no .mat at all
```

`cfg.isDryRun` is derived once, as `~strcmp(cfg.solverBackend,'comsol')`, and is the single flag everything downstream tests — so there is exactly one place in the file where the real backend is distinguished from the cheap ones.

### The three backends

| `cfg.solverBackend` | Cost / evaluation | Returns | What it exercises |
|---|---|---|---|
| `'comsol'` *(default)* | minutes | real band structure | the actual study. **The only backend whose numbers mean anything.** |
| `'surrogate'` | microseconds | analytic fake `ds` | the **success** path: `gapFitness` scores a correctly-shaped gap, the log fills, the checkpoint writes, all four figures draw from data |
| `'stub'` | microseconds | nothing (`objective = NaN`) | the **failure/skip** path: `solver-error` status, the red ✗ markers on the convergence panel, and `visualizeBestDesign`'s `band-data-unavailable` degradation |

Two cheap tiers are not redundant, and the distinction is the same one `objective.py` draws between its analytic `optical_surrogate` and its `surrogate_stub`. A dry run that only ever *succeeded* would leave the failure paths untested — and those are the branches that only execute on a bad day, which makes them the ones most likely to be broken. `'stub'` is this file's shorter name for `objective.py`'s `'surrogate_stub'`.

`cfg.dryRunFailEvery = N` lets the surrogate tier reach the failure path too: any design whose `a+r+w+th` (in nm) is divisible by `N` reports a simulated solver failure. It is keyed on the **design**, never on an evaluation counter or a random draw, because `bayesopt` is told `'IsObjectiveDeterministic', true` and a counter would make the second visit to a design succeed where the first failed — inventing observation noise the real solver does not have. Off by default.

`cfg.dryRunDelay` inserts an artificial `pause` per evaluation, on the cheap backends only. Default `0`; set `0.5`–`2` s when the thing being debugged is the live `bayesopt` plots themselves, which are hard to watch fill in at thousands of evaluations a minute.

### How cache isolation is enforced

⚠️ **This is the hazard that matters.** `solveBands` caches by **filename** — it skips the solve entirely whenever `[datLoc fileBase '_bds.mat']` exists — and `boomerangParams` deliberately reproduces the sweep scripts' filename byte-for-byte so real sweep data gets reused (see [How caching works](#how-caching-works)). A dry run that wrote those same filenames would make **every later real run at those geometries silently load fabricated band data**: invisible, permanent corruption of the one thing this directory exists to produce.

Four independent guards prevent it. Any one would be sufficient; they are layered because the first two are conventions about filenames, and a convention is only as strong as the last person who moved a file.

| # | Guard | Mechanism |
|---|---|---|
| 1 | **Separate folder** | Dry-run output goes to `test/boomerang_bayesopt/DRYRUN_<backend>_<date>/`. A real run never looks inside it, so the two namespaces do not intersect at all |
| 2 | **Filename prefix** | Every synthetic `_bds.mat` is named `DRYRUN_boomerang_a_...`. `cfg.dryRunPrefname` is also set as `P.prefname`, and the prefix is applied by hand exactly as `solveBands` would (`[prefname '_' fBase]`) — necessary because `solveBands` only applies `P.prefname` when `P.fileBase` is unset, and this script sets it |
| 3 | **In-file provenance marker** | Every synthetic `ds` carries `ds.isSynthetic = true` plus `syntheticBackend`, `syntheticWarning` and `syntheticCreated`. `assertBandDataProvenance` checks **every** cache load, on both the objective path and the figure path. Synthetic data reaching a real run is an **error**, not a warning — stopping a study that can be resumed from its checkpoint is a far smaller loss than finishing one whose numbers cannot be trusted. The reverse (real data reaching a dry run) only warns: it wastes the dry run but corrupts nothing |
| 4 | **Write nothing** | `cfg.dryRunSaveBands = 0` writes no `.mat` at all. The default is `1` because saving is what lets the dry run also exercise the cache-hit branch and `loadCachedBandData`, i.e. the code paths that read band data back off disk |

The same separation covers the log, the checkpoint, the results `.mat` and every figure, all of which are named through **`cfg.filePrefix`** (`bayesopt_boomerang` → `bayesopt_boomerang_DRYRUN`). For `'comsol'` that prefix reproduces exactly the literal that used to be hard-coded at each use site, so a real study writes byte-identical filenames to before this option existed.

### How a synthetic result is labelled

A result must be unmistakably synthetic from every direction someone could approach it from:

- **Console** — a boxed banner at the start *and* at the end of the run (the start banner has scrolled away by the time the best design prints), plus a per-evaluation `[DRY RUN] backend = surrogate : SYNTHETIC result, not physics`.
- **Log** — a trailing `backend` column reading `surrogate-SYNTHETIC` / `stub-SYNTHETIC` (a real row is the bare `comsol`).
- **`bayesopt` record** — `userData.backend`, so `results.UserDataTrace{k}.backend` answers "was this number real?" for any evaluation in a saved study.
- **Folder** — a `DRYRUN_README.txt` dropped beside the data, in plain text, readable without MATLAB.
- **Figures** — a red `*** DRY RUN: SYNTHETIC ... DATA - NOT PHYSICS - DO NOT USE FOR DESIGN ***` line in every figure headline, in the band panel's title, and at the top of the summary text panel; plus a `[DRY RUN - SYNTHETIC]` tag in each figure window name. Colour is load-bearing here, not decoration: a PNG pasted into an email arrives with no folder name and no console output attached.
- **Saved data** — `ds.isSynthetic` and `ds.syntheticWarning` inside every `_bds.mat`, and `cfg` (hence `cfg.solverBackend`) inside the results `.mat`.

### What the surrogate actually computes

**Nothing physical.** `surrogateBoomerangBands` solves no eigenproblem and is a fit to nothing; it evaluates a handful of closed-form expressions chosen to *look* like a band structure. The same disclaimer, for the same reasons, as `python-scripts/src/optical_surrogate.py` ("a cheap, dependency-free estimator … use it only as a cheap pre-screen and to exercise the loop") — except weaker still, since that one at least models a Bragg stack.

What it does get right is what makes it useful as a test:

1. **Shape.** Returns the same `ds` a real two-sector solve does: `ds.sym` / `ds.asym` each with `F` `[nk × nbands]` in Hz and `k_norm` `[nk × 1]` running 0→3, and `ds.full` with `F`, `midGap` and `gapSize`. `nk = 3*cfg.kpts + 1` (28 by default), matching the Γ–M–K–Γ circuit `runBands_2D` walks plus the wrapped final Γ.
2. **Determinism.** No `rand`, no clock, no counters. Required: `'IsObjectiveDeterministic', true` is passed to `bayesopt`, and a surrogate that answered differently on a repeated point would corrupt the GP fit — exactly the class of bug a dry run exists to find, not to introduce.
3. **Self-consistency.** `midGap`/`gapSize` are **not** hand-written. The synthesized band matrix is passed through the real `findGaps`, the way `solveBands` does it, which keeps the shaded gaps consistent with the plotted bands and exercises `findGaps` too.
4. **A non-trivial optimum.** A flat or monotone surface would make the convergence trace and the four design slices meaningless as tests — a plot of nothing renders perfectly well.

The spectrum is built as a ladder of `2*cfg.nbands` pass bands with centres going as `m^0.55` (so high orders crowd, as a folded spectrum does), split odd/even between the two symmetry sectors, and given a periodic `k` dependence that returns to its starting value at `k = 3`. A gap is opened between rungs 3 and 4 whose width is driven by a `gapStrength ∈ [0,1]`, itself a product of five Gaussian factors:

| factor | optimum | role |
|---|---|---|
| air filling fraction `fill` | 0.0823 | gap width |
| slab aspect ratio `th/a` | 0.390 | gap width |
| arm aspect ratio `w/r` | 0.425 | gap width |
| hole reach `r/a` | 0.250 | gap width; breaks the `(a,r,w)` ridge the other ratios leave |
| absolute thickness `th` | 312 nm | gap width; the only **non**-dimensionless factor, and the only thing that breaks the `(a,th)` scaling degeneracy |

plus one scaling law — frequencies go as **1/a** — so the Gaussian target-frequency penalty in `gapFitness` has real work to do and *fights* the gap-width factors instead of agreeing with them.

Every optimum is placed so the peak lands at the **centre of the default bounds**: `a = 800`, `r = 200`, `w = 85`, `th ≈ 312` nm, fitness ≈ 0.24, mid-gap 13.0 GHz (= `cfg.targetFreq`, by construction). Centring matters — an optimum near a bound makes the best design *on* that bound almost as good as the peak, and the corresponding design slice then looks flat and reports `[at bound]`. Measured on a grid scan of the feasible box, the peak beats the best design available on each bound face by **2.5× (`a`), 1.9× (`r`), 2.4× (`w`), 1.5× (`th`)**, and about **64%** of the box has no complete gap at all — so the gapless branch of `gapFitness` gets exercised too.

> If you change `cfg.bounds`, `cfg.targetFreq` or `cfg.sigma`, re-check that the surrogate optimum is still interior. It is tuned against the shipped values, and the absolute-thickness factor in particular does not follow `cfg.th`, which is now a single fixed value rather than a searched range — so if you move `cfg.th` away from ~312 nm the surrogate's thickness factor is simply off-peak everywhere, uniformly. A surrogate whose optimum has drifted onto a bound tests considerably less than it appears to.
>
> A measured example, since this is easy to hit in practice: retargeting to `cfg.targetFreq = 7e9`, `cfg.sigma = 3e9` with `cfg.nbands = 10` and `cfg.bounds.r = [100 250]` still leaves the surrogate working — 20 % of a 3969-design scan score nonzero, best fitness 0.0515 — but the optimum moves to `a = 1000`, `r = 250`, i.e. **pinned on both upper bounds**, versus an interior `a = 800, r = 200, w = 85, th = 312` at fitness 0.2427 under the shipped 13 GHz / 5 GHz. The cause is structural: surrogate frequencies scale as 1/a and the ladder is centred so that 13 GHz falls at a = 800, so a lower target can only be reached with a bigger cell and `a` runs to the ceiling. The loop is still fully exercised and the run is still valid as a *plumbing* test — but the `a` and `r` slices will look bound-limited, and that says nothing whatever about where the real COMSOL optimum lies.

---

## Read this before your first real run

**The default bounds do not contain your current design point.** `cfg.bounds` was taken from the grid ranges in `sweep_boomerang_code.m`, which predate the geometry now in `test_Boomerang.m` (`a=480, w=140, r=177, th=220` nm):

| variable | `cfg.bounds` default | `test_Boomerang.m` | |
|---|---|---|---|
| `a` | 600 – 1000 nm | 730 | inside ✓ |
| `w` | 50 – 120 nm | 125 | **above the range** |
| `r` | 150 – 250 nm | 300 | **above the range** |
| `th` | *fixed at `cfg.th` = 300 nm* | 220 | not searched |

So as shipped, the optimizer cannot reach — or even model — the region you are presumably most interested in. Decide deliberately which design space you want searched. To centre it on the current geometry:

```matlab
cfg.bounds.a  = [400  600];
cfg.bounds.r  = [140  210];
cfg.bounds.w  = [100  180];
cfg.th        = 220e-9;      % fixed, not a range
```

Whatever you choose, check the containment ceiling stays reachable across the box or candidates get pruned (see [Constraints](#constraints)). With the default `holeCentreFrac` = 0.4 that ceiling is ~379 nm at `a` = 730 nm — looser than the old centred-hole `sqrt(3)*a/4` = 316 nm, so the containment test binds less than it used to. **The filling-factor band is now the constraint that prunes**, by a wide margin; size the bounds against [its feasibility table](#filling-factor-constraint), not against containment.

---

## Unit cell geometry

Built by `buildBoomerangUnitCell.m` — the builder `runBands_2D` dispatches for `celltype = 'boomerang'`. (Note `buildBoomerangUnitCell_2D.m` exists but is *not* on this path.)

**The cell** is the primitive cell of a hexagonal (triangular) lattice: a rhombus of side `a` with a 60° interior angle, spanned by lattice vectors `a1 = (a, 0)` and `a2 = (a/2, a·√3/2)`. The cell centre sits at `(3a/4, a·√3/4)`. **`a` alone sets the cell footprint.**

**The hole** is a three-pointed star: three rectangular arms, each `w` wide and `r` long, radiating from a common centre at 120° spacing (θ = 90°, 210°, 330°). This tri-arm shape is what gives the cell its name.

⚠️ **The hole is no longer centred on the cell.** `resolveHoleCentreFrac.m` places it at `f·(a1 + a2)` with `f` = `P.holeCentreFrac`, **default 0.4** — the centroid is `f` = 0.5. Sliding down the diagonal trades clearance from the two down-pointing arms, which have room to spare, to the single up-pointing arm, which is what binds: max `r` before the hole crosses the cell boundary goes from **316 nm to ~379 nm** at `a` = 730 nm. This is a translation of the motif inside a periodic cell, so the crystal it tiles into — and therefore the band structure — is unchanged. Set `P.holeCentreFrac = 0.5` to restore the centred geometry.

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

### The fitness function

For a design $\mathbf{x} = (a,\, r,\, w,\, t_h)$, `solveBands` returns a set of complete mechanical bandgaps. Index them by $k$, and let gap $k$ have lower and upper edges $f_k^-$ and $f_k^+$. `findGaps` reports each gap as a mid-gap frequency and a width:

$$\Omega_k \;=\; \frac{f_k^+ + f_k^-}{2}, \qquad \Delta_k \;=\; f_k^+ - f_k^-$$

Each gap is scored by its **fractional gap width** times a **Gaussian frequency penalty**:

$$\eta_k \;=\; \frac{\Delta_k}{\Omega_k} \;=\; \frac{2\,(f_k^+ - f_k^-)}{f_k^+ + f_k^-}, \qquad p_k \;=\; \exp\!\left[-\left(\frac{f_\mathrm{target} - \Omega_k}{\sigma}\right)^{\!2}\right]$$

The **fitness** is the best score over all admissible gaps, and since `bayesopt` minimizes, the **objective** is its negation:

$$\boxed{\;F(\mathbf{x}) \;=\; \max_{k \in \mathcal{G}(\mathbf{x})} \; \frac{\Delta_k}{\Omega_k} \cdot \exp\!\left[-\left(\frac{f_\mathrm{target} - \Omega_k}{\sigma}\right)^{\!2}\right] \;}$$

$$\text{objective}(\mathbf{x}) \;=\; -F(\mathbf{x})$$

where $\mathcal{G}(\mathbf{x}) = \\{\, k : \Omega_k > 0 \ \wedge\ \Delta_k > 0 \,\\}$ is the set of admissible complete gaps, and $F(\mathbf{x}) = 0$ by definition when $\mathcal{G}(\mathbf{x}) = \varnothing$.

In plain text, for anyone reading this file unrendered:

```
eta_k     = gapSize(k) / midGap(k)                                  % fractional gap width
p_k       = exp( -((f_target - midGap(k)) / sigma)^2 )              % frequency penalty
F         = max_k ( eta_k * p_k )                                   % fitness,  F >= 0
objective = -F                                                      % bayesopt minimizes
```

### Symbols

| symbol | code | meaning |
|---|---|---|
| $a,\ r,\ w,\ t_h$ | `x.a`, `x.r`, `x.w`, `x.th` | Design variables [integer nm] |
| $f_k^\pm$ | — | Lower / upper edge of gap $k$ [Hz] |
| $\Omega_k$ | `gapData.midGap(k)` | Mid-gap frequency of gap $k$ [Hz] |
| $\Delta_k$ | `gapData.gapSize(k)` | Width of gap $k$ [Hz] |
| $\eta_k$ | `gapRatios(k)` | Fractional gap width, dimensionless |
| $p_k$ | `penalties(k)` | Gaussian frequency penalty, $p_k \in (0,\,1]$ |
| $f_\mathrm{target}$ | `cfg.targetFreq` | Target mid-gap frequency [Hz], default 13 GHz |
| $\sigma$ | `cfg.sigma` | Penalty width [Hz], default 5 GHz |
| $F$ | `fitness` | Fitness, maximised |
| $-F$ | `objective` | What `bayesopt` minimizes |

### Properties worth knowing

- **$F \geq 0$ always**, since $\eta_k > 0$ and $p_k > 0$ on $\mathcal{G}$. So `objective` $\leq 0$, and `objective = 0` is the worst attainable value.
- **$p_k \leq 1$, with equality iff $\Omega_k = f_\mathrm{target}$.** The penalty can only ever discount the fractional gap width, never inflate it — so $F \leq \max_k \eta_k$, and $F = \eta_k$ exactly when gap $k$ is centred on the target.
- **$\eta_k$ is the standard fractional bandwidth** $2(f^+ - f^-)/(f^+ + f^-)$, so it is scale-free: two geometrically similar cells differing only in overall scale have the same $\eta$ but different $\Omega$. The penalty term is what breaks that degeneracy and pins down an absolute frequency.
- **$\sigma$ sets how far off-target a gap may sit before it stops counting.** At $|f_\mathrm{target} - \Omega_k| = \sigma$ the score is already down by $e^{-1} \approx 0.37$; at $2\sigma$ by $e^{-4} \approx 0.018$. With the defaults ($f_\mathrm{target} = 13$ GHz, $\sigma = 5$ GHz), a gap centred at 23 GHz retains under 2 % of its fractional width.
- **This is the same functional form as `optimization_client.m`** (`mech_cont = mGap/mFreq*freq_penalty; F = -1*mech_cont`). What changed is only *which* gap enters the formula — see below.

Two factors, then: the **fractional gap width** $\eta_k$, which is what you actually want maximised, and a **Gaussian penalty** $p_k$ pulling the mid-gap frequency toward `cfg.targetFreq` with width `cfg.sigma`.

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

2. **Arms stay inside the cell** — now **measured, not analytic.** This used to be `r < sqrt(3)*a/4`, which is the containment limit for a hole sitting on the cell *centroid*. That stopped being the geometry the builder produces when `resolveHoleCentreFrac.m` introduced a default **0.4** diagonal offset (the centroid is `0.5`), so the test is now `calcFillingFactor`'s `armsOverhang` flag — it compares hole area before and after the clip to the cell, which is exact for any hole position and stays correct if the offset is retuned.

   The offset is not cosmetic: it raises the largest usable `r` from ~316 nm to ~379 nm at `a` = 730 nm, against a physical merge limit of ~400 nm. See [`resolveHoleCentreFrac.m`](resolveHoleCentreFrac.m).

3. **Filling factor in a band:** `|ff - cfg.fillingFactor| <= cfg.fillingFactorTol`, where `ff = area(air) / area(dielectric)` from [`calcFillingFactor.m`](calcFillingFactor.m) — the same polygons `buildBoomerangUnitCell` hands COMSOL, so the constraint and the geometry cannot disagree. Set `cfg.fillingFactor = []` to disable.

   **Why a band and not an equality.** `a`, `r`, `w` live on an integer-nm grid, so the set where the ratio hits a target *exactly* has measure zero — every candidate would be pruned. The tolerance turns the equality into a thin feasible shell: roughly two-dimensional inside a three-dimensional box, which is the point. Combined with fixed `th`, the optimizer searches ~2 free directions instead of 4.

   **Why it is measured rather than computed from a formula.** The three arms overlap at the hole centre, so the naive `3*w*r` overstates the air area by about **6%** at the current geometry (99 000 vs 93 600 nm² at `w`=110, `r`=300 nm) — and the error varies nonlinearly with `w/r`, so it cannot be absorbed into a constant.

   **Cost.** `calcFillingFactor` is called with `'FeatureSizes',false`, which skips its `minSolidFeature` measurement. That measurement is O(n²) over a few thousand densified boundary points × six neighbours — fine once per design, ruinous in a constraint function `bayesopt` evaluates on thousands of candidates per iteration. The areas, the filling factor and `armsOverhang` all survive the fast path. Rows already failing test (1) are skipped before any geometry work.

`boomerangFabConstraint` is **vectorized over a multi-row table**, as `bayesopt` requires, and compares in integer nm so a design sitting exactly on a limit is not decided by floating-point round-off.

### On `a - r`

The `min(a - r, ...)` fabrication proxy used in `boomerang_optimize_sweep_diamond.m` is **deliberately not reproduced here** — it has no geometric meaning for this shape. Sampling the tri-arm outline against its periodic images gives a true inter-hole diamond ligament of ~161 nm for the current design, where `a - r` reports 303 nm; it overstates by roughly 2×. (`√3a/4 - r` ≈ 31 nm is the other tempting expression, but that is only the gap to the cell *boundary* — material continues past it into the neighbouring cell, so it is not a feature either.) Over the default bounds the true ligament never drops below ~210 nm, so an `a - r` term would neither bind nor measure anything real.

---

## Configuration reference

All settings live in the `cfg` struct at the top of the script.

### Solver backend / dry run

The first block in `cfg`, because it is the most consequential switch in the file. See [Dry-run mode](#dry-run-mode) for the full story.

| field | default | meaning |
|---|---|---|
| `cfg.solverBackend` | `'comsol'` | `'comsol'` \| `'surrogate'` \| `'stub'`. Anything else is an **error**, checked at configuration time and again in the dispatcher |
| `cfg.isDryRun` | *derived* | `~strcmp(cfg.solverBackend,'comsol')`. Single source of truth; do not set by hand |
| `cfg.dryRunDelay` | `0` | Artificial seconds per evaluation, cheap backends only. `0.5`–`2` to watch the live plots at human speed |
| `cfg.dryRunFailEvery` | `0` | `0` = off. `N` > 0 makes any design with `mod(a+r+w+th, N) == 0` report a simulated solver failure. Keyed on the design, not a counter, so `IsObjectiveDeterministic` still holds |
| `cfg.dryRunSaveBands` | `1` | `1` writes the synthetic `ds` to a `_bds.mat` under the dry-run folder, which is what exercises the cache-hit branch and `loadCachedBandData`. `0` touches nothing but the log, and the visualization then degrades to its `band-data-unavailable` path |
| `cfg.dryRunPrefname` | `'DRYRUN'` | Filename prefix stamped onto synthetic band data. Applied by hand rather than left to `solveBands`' `P.prefname` — see guard 2 above |
| `cfg.filePrefix` | *derived* | `bayesopt_boomerang`, or `bayesopt_boomerang_DRYRUN` in a dry run. Used for the log, checkpoint, results `.mat` and all four figures |

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
| `cfg.bounds.r` | `[100 250]` |
| `cfg.bounds.w` | `[50 120]` |

`cfg.bounds.th` is **gone** — `th` is no longer searched. See [the caveat above](#read-this-before-your-first-real-run) — these do not contain the current design point.

### Fixed geometry

| field | default | meaning |
|---|---|---|
| `cfg.th` | `300e-9` | **Full slab thickness [m]. Fixed, not searched.** Set it to whatever thickness the process actually delivers — there is no value in optimizing a dimension you cannot choose |
| `cfg.r1` | `10e-9` | Fillet radius, inner corners (arms meet at centre) |
| `cfg.r2` | `10e-9` | Fillet radius, outer arm tips |

### Filling-factor constraint

| field | default | meaning |
|---|---|---|
| `cfg.fillingFactor` | `0.25` | Target `area(air)/area(dielectric)`. `[]` disables the constraint entirely |
| `cfg.fillingFactorTol` | `0.01` | Accept `\|ff − target\| <= tol`. Too tight and everything is pruned; see the scan below |
| `cfg.fillingFactorMinFeasible` | `0.02` | Feasible fraction below which the startup scan warns. Zero feasible is always an error |

**A startup feasibility scan runs before any solve.** It scores a 12³ grid over the bounds through the *same* constraint function `bayesopt` will use, and reports what fraction survives. Zero → **error**, study not started. Below `cfg.fillingFactorMinFeasible` → **warning**.

That check exists because the failure it catches is silent and expensive: `bayesopt` never reports "your constraint pruned everything" — it keeps drawing candidates, finds almost none admissible, and either stalls or spends the whole budget in a sliver of the box.

**Choose the target with the bounds in mind.** Measured over the shipped bounds (`a`∈[600,1000], `r`∈[150,250], `w`∈[50,120] nm, `holeCentreFrac` = 0.4):

| statistic | filling factor |
|---|---|
| range | 0.025 – 0.368 |
| 10th / 25th percentile | 0.051 / 0.067 |
| median | 0.092 |
| 75th / 90th percentile | 0.128 / 0.174 |

and the resulting feasible fractions:

| target | ±0.005 | ±0.01 | ±0.02 |
|---|---|---|---|
| 0.15 | 3.3% | 6.6% | 13.4% |
| 0.20 | 1.3% | 2.5% | 5.3% |
| 0.25 | 0.5% | **0.9%** | 1.9% |

⚠️ The shipped default of **0.25 ± 0.01 leaves under 1% of the box** — it matches the geometry in `test_Boomerang.m` (`a`=730, `w`=125, `r`=300 → ff = 0.254), but that design sits *outside* these bounds, and the box as shipped is mostly low-fill. Expect the tight-feasibility warning. Either move `cfg.bounds` toward the shell you care about, widen `cfg.fillingFactorTol`, or drop the target toward the box median.

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
| `cfg.datLoc` | `test/boomerang_bayesopt/<mmddyyyy>/` (host separators) | Output folder, built with `fullfile` — see [Cross-platform paths](#cross-platform-paths). **The trailing separator is required.** A dry run gets `DRYRUN_<backend>_<mmddyyyy>/` instead |
| `cfg.plotgeom` | `0` | 1 to plot geometry every evaluation (slow, noisy) |
| `cfg.savebndplot` | `1` | 1 to save a band diagram per evaluation |
| `cfg.closeSolveFigures` | `1` | Close figures the solve opened, keeping the `bayesopt` live plots alive |

That last option exists because a plain `close all` would also destroy the `bayesopt` plot windows. Instead the figure list is differenced against a snapshot taken before each solve, so only figures the solve created get closed.

### Figures

| field | default | meaning |
|---|---|---|
| `cfg.makeSummaryFigures` | `1` | Master switch for the whole post-run [visualization](#visualization) |
| `cfg.figResolution` | `150` | PNG export resolution [dpi], matching the Python scripts |
| `cfg.figNPeriods` | `3` | Lattice periods tiled in the geometry panel. `1` draws the bare cell; ≥2 makes the inter-hole ligament visible |
| `cfg.closeSummaryFigures` | `0` | `0` leaves the saved figures on screen; `1` closes **only** the handles it created |

---

## Outputs

Everything lands under `cfg.datLoc`. `<prefix>` below is `cfg.filePrefix`: `bayesopt_boomerang` for a real study, `bayesopt_boomerang_DRYRUN` for a [dry run](#dry-run-mode).

| file | contents |
|---|---|
| `<fileBase>_bds.mat` | Per-evaluation band structure `ds` struct, written by `solveBands`. Also the cache — see below |
| `<fileBase>_fullBands.png` / `.fig` | Band diagram per evaluation, when `cfg.savebndplot = 1`. Real runs only — the cheap backends never call `solveBands` |
| `<prefix>_log.txt` | One tab-separated row per evaluation |
| `<prefix>_state.mat` | Checkpoint of the `BayesianOptimization` object, rewritten every iteration |
| `<prefix>_results.mat` | Final `results`, `cfg`, and `xBest` |
| `<prefix>_summary.png` / `.fig` | Composite 3×4 figure — see [Visualization](#visualization) |
| `<prefix>_geometry.png` / `.fig` | Standalone best-cell geometry, tiled over `cfg.figNPeriods` |
| `<prefix>_bestbands.png` / `.fig` | Standalone band structure of the best design, gaps shaded |
| `<prefix>_progress.png` / `.fig` | Standalone convergence trace + per-variable design slices |
| `DRYRUN_README.txt` | **Dry runs only.** Plain-text warning that everything in the folder is synthetic |

`<fileBase>` is e.g. `boomerang_a_700nm_r_262nmw_99nm_th_300nm_r1_10nm_r2_10nm_`. The script sets `P.fileBase` explicitly, reproducing byte-for-byte the string `solveBands` builds for `celltype = 'boomerang'` — including its missing `_` after the `r` field. That is intentional: it lets the script predict the cache path, and keeps filenames interchangeable with data already produced by the sweep scripts.

In a dry run `<fileBase>` additionally gains a `DRYRUN_` prefix and the whole folder changes, so nothing here can collide with a real study — see [How cache isolation is enforced](#how-cache-isolation-is-enforced).

### Log columns

`<prefix>_log.txt` is written and closed once per evaluation rather than buffered, so a run killed mid-study still leaves a complete record.

| column | notes |
|---|---|
| `a_nm`, `r_nm`, `w_nm` | Design point |
| `th_nm` | **Constant** — `cfg.th`, no longer searched. The column is kept so the 14-column format and any existing `readtable` call still work |
| `midGap_GHz`, `gapSize_GHz` | Winning gap |
| `gapRatio` | `gapSize/midGap` |
| `penalty` | Gaussian frequency penalty, 0–1 |
| `fitness` | `gapRatio * penalty` |
| `objective` | `-fitness` — what `bayesopt` sees |
| `cached` | 1 if the band structure was reloaded from an existing `.mat` |
| `evalTime_min` | Wall-clock minutes for this evaluation |
| `status` | `target-in-gap`, `nearest-gap`, `no-complete-gap`, or `solver-error` |
| `backend` | `comsol` for a real solve; `surrogate-SYNTHETIC` / `stub-SYNTHETIC` for a [dry run](#dry-run-mode) |

Load it with `readtable('<prefix>_log.txt')` — 14 columns.

The `backend` column is appended last, so a log written before the column existed still lines up for its first 13 fields. A dry run writes to a different folder anyway, so the column is not what keeps the two apart — but a row that has been copied, pasted or `readtable`'d out of its folder has lost every other clue, and a table of gap ratios with no provenance is precisely the artifact that gets believed.

---

## Cross-platform paths

Every path in this script is built with `fullfile`, never by pasting a separator in by hand, so the same file produces a proper nested directory on Windows, macOS and Linux alike:

```matlab
cfg.datLoc  = [fullfile('.','test','boomerang_bayesopt',datedFolder), filesep];
cfg.logPath = fullfile(cfg.datLoc,[cfg.filePrefix,'_log.txt']);
```

This matters because the rest of this directory hard-codes `'.\test\...'`, and that form is only correct on Windows. On macOS and Linux a backslash is an ordinary filename character rather than a separator, so those scripts create a **single file literally named** `test\boomerang_sweep\08132026\...` instead of a folder tree — and that name cannot even be represented on Windows, so the artifacts stop being portable in both directions. `fullfile` emits `filesep` for the host and collapses duplicate separators.

**On Windows the result is byte-identical to the old hard-coded literal** (`.\test\boomerang_bayesopt\<date>\`), which is asserted in the verification, so nothing about an existing Windows workflow changes.

Two deliberate exceptions:

- **`cfg.datLoc` keeps a trailing separator** (a `filesep`, not a literal `'\'`), because `solveBands` needs one — see [quirk 2](#solvebands-quirks-this-script-works-around).
- **Extensions are still appended by concatenation** (`[pathNoExt '.png']`). That is a suffix, not a path join, so `fullfile` would be wrong there.

> **The sibling scripts have not been converted.** `test_Boomerang.m`, `sweep_boomerang_code.m` and the other `test_*`/`sweep_*` scripts still hard-code `'.\test\...'`, so their output is still Windows-only. Only `bayesopt_boomerang.m` was in scope here.

Note that `runBands_2D.m` independently normalises `P.datLoc` with `strrep(...,'\',filesep)` before creating directories, while `solveBands.m` does not — which is why, before this change, a macOS run scattered its band diagrams into a real directory tree and its `.mat` files into a backslash-named file beside it.

---

## Visualization

After the study finishes, the script draws a **composite 3×4 figure** plus three standalone figures, modelled on the multi-panel characterization figure in `python-scripts/scripts/run_opt_comsol.py`.

The composite layout:

| row | contents |
|---|---|
| 1 | Best-cell **geometry** (2 tiles, tiled over `cfg.figNPeriods` periods) \| **band structure** of the best design (2 tiles, complete gaps shaded) |
| 2 | **Convergence** trace (2 tiles) \| **summary** text panel (2 tiles, monospace) |
| 3 | **Design slices** — fitness vs `a`, `r`, `w`, one tile each, with the search bounds drawn. The fourth tile is now **empty**: `drawDesignSpaceRow` reads variable names off `XTrace` rather than hard-coding them, so dropping `th` needs no change there, but the layout still reserves four tiles |

Notes on reading these:

- Everything is plotted as **fitness** (= `-objective`) so up is better, matching the sign convention you actually care about. The axis labels say so.
- The **design slices** exist to show whether the optimizer pinned a variable against a bound. If the best point sits on a bound, the bound — not the physics — is setting your design, and the box should be widened.
- Errored evaluations (`NaN`) are marked with a red ✗ on the convergence panel rather than silently dropped.
- The band panel marks the target frequency and the specific mid-gap the fitness scored. With several complete gaps shaded it is otherwise ambiguous which one produced the number.
- The y-axis is extended to keep the target line on screen even when the solved bands stop below it, with an explicit note in that case — otherwise the marker vanishes in exactly the situation where it matters most.
- The geometry panel **omits the `r1`/`r2` fillets** (noted in grey inside the axes). It is drawn from `polyshape` in pure MATLAB, so it needs no COMSOL and costs nothing.

Band data is reloaded from the best design's cached `_bds.mat`. The whole section is **best-effort**: it runs *after* `bayesopt_boomerang_results.mat` is written, and it is wrapped at three levels (the block, each figure, each panel), so a plotting failure can only ever cost you a panel — never a finished multi-hour study. Missing cache, a cache-hit stub with no `.full`, a gapless result, and an all-`NaN` trace were all tested and each degrades to a warning plus a labelled empty panel.

It never calls a bare `close all`, which would destroy the live `bayesopt` plot windows; only figure handles it created are closed, and only when `cfg.closeSummaryFigures = 1`.

In a [dry run](#dry-run-mode) every figure gains a red banner line, and the summary text panel drops its blank separator lines to make room — the panel's 21 lines fill its tile exactly, and a 22nd is clipped by the tile below rather than scaled, landing on the design-slice titles. If you add a line to that panel, re-render and look at the PNG.

---

## Resuming an interrupted study

The `BayesianOptimization` object is checkpointed after **every** iteration, because losing a partly finished study to a dropped LiveLink connection is expensive at these evaluation costs.

```matlab
S = load(fullfile('test','boomerang_bayesopt','<date>','bayesopt_boomerang_state.mat'));
results = resume(S.results, 'MaxObjectiveEvaluations', 20);
```

To extend a study that finished cleanly, `resume` the in-memory object directly:

```matlab
results = resume(results, 'MaxObjectiveEvaluations', 20);
```

Note that `resume` continues with the surrogate already fitted — it does not re-spend seed points.

⚠️ A checkpoint written **before** the dry-run option existed cannot be resumed: `resume` re-resolves the stored objective handle against the current file, so the old `cfg` (which has no `solverBackend` field) meets the new dispatcher. You get an explicit error saying so. Re-run the script instead — it rebuilds `cfg` and the handle, and the existing `_bds.mat` files under `cfg.datLoc` are reused, so no solve time is lost.

---

## How caching works

`solveBands` **skips the solve entirely** when `[datLoc fileBase '_bds.mat']` already exists. This script leans on that deliberately:

- Re-running an interrupted study re-uses everything already solved.
- Existing sweep data under the same `datLoc` is picked up for free, since the filenames match.

**This is why the design variables are declared as integer nanometres rather than continuous metres.** `solveBands` names its files by rounding every dimension to whole nm (`'%.0f'`). With continuous variables, many distinct GP query points would alias onto one cached result, silently feeding the surrogate stale data for design points it never actually evaluated. On an integer-nm grid the rounding is exact, so a repeated point is a genuine cache hit.

The flip side: **editing physics settings without changing geometry will silently reuse old results.** Clear the dated folder — or point `cfg.datLoc` somewhere new — whenever you change `kpts`, `nbands`, `meshSize`, or the symmetry settings.

This filename-keyed cache is also the reason a dry run is dangerous if left unguarded, and why four separate guards keep synthetic band data out of it — see [How cache isolation is enforced](#how-cache-isolation-is-enforced). Note the contrast with the Python pipeline, which sidesteps the problem by folding the backend into the cache key itself (`cid = _hash_u(u, (optical_backend, mech_backend))` in `objective.py`); here the key is a filename that must stay byte-compatible with the sweep scripts, so isolation is enforced around the cache rather than inside it.

---

## `solveBands` quirks this script works around

Recorded here because they are easy to trip over again.

1. **Two different return shapes.** A fresh solve returns `ds` with a `.full` field, but a cache hit prints `Data folder exists in working directory` and returns a stub with only empty `.sym`/`.asym`. Reading `ds.full` off that stub throws `Unrecognized field name` — which over an optimization run, where revisiting a design point is routine, would abort the whole study. `solveBoomerangBands` detects this and reloads the cached `.mat` instead.

2. **`P.datLoc` needs its trailing separator.** `solveBands` normalises a separator into a local variable but then writes the `.mat` using the raw `P.datLoc`, so a missing separator drops output into the parent folder under a mangled name. This is why `cfg.datLoc` keeps a trailing `filesep` even though every other path here is built with `fullfile`.

6. **`solveBands` decides its own cache hit with `strcmp(P.datLoc(end),'\')`** — a hard-coded backslash. On macOS and Linux that test fails for a correctly-formed path, so it appends a literal backslash to the path it then uses for its `dir()` probe and `mkdir`, producing a stray directory named `\` and a probe that can never match the file it just wrote. This script therefore answers the cache question itself before calling `solveBands`, which makes caching behave identically on all three platforms; on Windows the only difference is that one no-op call into `solveBands` is skipped.

3. **`P.savedat` must stay `1`.** It is what writes the `.mat`, and therefore what makes caching work at all.

4. **`TwoSymPlanes` is set to `0`, deliberately.** On the 2D path `runBands_2D` applies a Floquet periodic BC in y and reads `P.mbeveny` without ever using it, so only `mbevenz` changes the problem posed. With `TwoSymPlanes = 1`, `solveBands` runs four sectors that collapse to two distinct solves run twice each — double the cost, and `full.F` then holds every band twice, making `findGaps` report every complete gap twice over. Two sectors (even/odd about z) is already the complete set. Note `test_Boomerang.m` currently sets `1`.

5. **`P.bandStruct_2D` must stay `1`.** `runBands` (the `0` path) has no `'boomerang'` branch in its celltype dispatch — it handles only `boomerang_strip_v2` and `boomerang_lower`.

---

## Local function reference

| function | role |
|---|---|
| `boomerangFabConstraint` | Vectorized deterministic feasibility test on a candidate table |
| `boomerangObjective` | Evaluate one design point and score it; returns `-fitness` |
| `boomerangParams` | Build the `P` struct for one design point, ported from `sweep_boomerang_thickness`; applies the `DRYRUN` filename prefix |
| `solveBoomerangBands` | **Dispatcher** on `cfg.solverBackend` — the single seam between the loop and whatever stands in for the physics. Errors on an unknown backend |
| `solveBandsViaComsol` | The real path: call `solveBands` and always return a usable gap struct, handling the cache-hit stub |
| `closeNewFigures` | Close only the figures a solve opened, preserving `bayesopt` plots |
| `gapFitness` | Score every complete gap against the target, keep the best |
| `initIterationLog` | Create the log with its header row |
| `logIteration` | Append one evaluation, flushed immediately |
| `logBackendLabel` | Value for the log's `backend` column |
| `checkpointState` | `OutputFcn` saving the optimizer state every iteration |

[Dry-run backends](#dry-run-mode) — inert unless `cfg.solverBackend` is a cheap tier, except for the provenance guard, which a real run calls precisely so it can refuse synthetic data:

| function | role |
|---|---|
| `solveBandsViaSurrogate` | Cheap tier 1: analytic fake bands, walking the same cache/save/return path as the real solve |
| `solveBandsViaStub` | Cheap tier 2: no bands at all, every evaluation reported as a failure |
| `surrogateBoomerangBands` | The analytic band-structure generator. **Not physics** — see [what it computes](#what-the-surrogate-actually-computes) |
| `simulatedSolverFailure` | Deterministic pseudo-failure predicate, keyed on the design |
| `assertKnownBackend` | Reject an unrecognised `cfg.solverBackend`, loudly |
| `formatBackendForError` | Describe whatever was found in `cfg.solverBackend`, without assuming it is printable |
| `isSyntheticBandData` | True if a `ds` carries the dry-run provenance marker |
| `assertBandDataProvenance` | Refuse to mix synthetic and real band data — error one way, warning the other |
| `printDryRunBanner` | Boxed console notice, printed at the start *and* the end of the run |
| `writeDryRunMarker` | Write `DRYRUN_README.txt` into the output folder |
| `dryRunBannerText` | One-line figure banner, or `{}` outside a dry run |
| `withDryRunHeadline` | Prefix a figure headline with the banner and grow its reserved strip |
| `dryRunHeadlineColor` | Red in a dry run, `[]` (no opinion) otherwise |

Visualization (all best-effort, none on the optimization path):

| function | role |
|---|---|
| `visualizeBestDesign` | Entry point; builds and saves all four figures |
| `buildCompositeFigure` | Assembles the 3×4 composite |
| `drawBoomerangCell` | Draws the rhombic cell minus the tri-arm hole, tiled, via `polyshape` |
| `plotBestBands` | Band structure with shaded complete gaps and reference lines |
| `plotBandSector` | Draws one symmetry sector's bands |
| `plotConvergenceTrace` | Best-so-far fitness plus per-evaluation points |
| `plotDesignSlice` | Fitness vs one design variable, with bounds |
| `drawDesignSpaceRow` | Lays out the four slices |
| `writeSummaryPanel` | Monospace text summary |
| `loadCachedBandData` | Reloads the best design's `_bds.mat` |
| `newSummaryFigure` | Creates a themed figure of a given size |
| `summaryFigureName` | Figure window title, tagged when the data behind it is synthetic |
| `addFigureHeadline` | Figure-level suptitle that cannot collide with axes titles; optional text colour |
| `saveFigurePair` | Writes `.png` + `.fig` |
| `tryPanel` | Wraps one panel so its failure cannot take down the figure |

---

## Tuning guidance

**`cfg.sigma`** controls how hard the frequency target is enforced. Too small and almost every design scores ~0, leaving the GP nothing to learn from; too large and the optimizer happily returns a wide gap at the wrong frequency. If early evaluations come back with `penalty` near zero across the board, widen it.

**`cfg.numSeedPoints`** should stay a reasonable fraction of `cfg.maxEvaluations` — with 3 variables, 8 seeds out of 40 is a sensible split, and one variable fewer than before means the same budget goes further. Raise it if the filling-factor scan warns that the feasible shell is thin, since seeds are drawn from the whole box and most will be rejected. Too few and the GP extrapolates from almost nothing; too many and you spend the budget on random search.

**`cfg.nbands`** must be large enough to bracket the gap of interest. If gaps go missing at higher frequency, raise it — but note that COMSOL solve cost grows with band count.

**`cfg.maxEvaluations`** is a hard stop, not a convergence criterion. Prefer a modest first run followed by `resume`, since that lets you inspect the log before committing more solve time.

**`cfg.fillingFactorTol`** is the knob that decides whether there is anything to search. Start from the [feasibility table](#filling-factor-constraint), and let the startup scan tell you rather than guessing — it is cheap and it runs before any solve. If the scan warns, raising `cfg.numSeedPoints` helps: seeds are drawn from the whole box, so a thin shell means most are rejected and the GP starts with very little.

**Do a dry run first.** `cfg.solverBackend = 'surrogate'` walks the whole loop in seconds, so every mistake that is not about physics — a bad bound, a constraint that prunes everything, a `datLoc` that cannot be written, a budget split that spends everything on seed points — surfaces before you spend a solve on it. Then `'stub'` for the failure paths, then `'comsol'`.

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
| `findGaps.m` | Locates complete gaps in an assembled band matrix — used by the surrogate too, so the synthetic gaps are found the same way the real ones are |
| `calcFillingFactor.m` | Measures air/dielectric area from the reconstructed polygons; supplies both the filling-factor constraint and the `armsOverhang` containment test. Called with `'FeatureSizes',false` on the constraint path |
| `resolveHoleCentreFrac.m` | Single source for the hole's position in the cell (`f` = 0.4 by default). Shared with `buildBoomerangUnitCell` so the constraint measures the cell COMSOL builds |
| `plotBoomerangCell.m` | Draws cell + hole + landmarks from the same `calcFillingFactor` result; useful for eyeballing a candidate before committing solve time |
| `python-scripts/src/objective.py` | Source of the string-valued backend selection, the fail-loud-on-unknown-backend rule, and the backend-in-the-cache-key idea |
| `python-scripts/src/optical_surrogate.py` | Source of the "analytic surrogate, loop debugging only" pattern the surrogate backend follows |
