---
name: feedback-scope-discipline
description: Keep scope to what was asked — no bonus test scripts or speculative feature gates; plus how to verify MATLAB changes headless without COMSOL
metadata:
  type: feedback
---

Write only the artifact that was asked for. Do not touch existing files, do not
add a test script "for convenience", and do not add optional features behind a
flag on the theory that they might be wanted later. If something seems essential
but was not requested, ask rather than assume.

**Why:** this repo has many near-duplicate builders (`buildHoleStrip_2D/_3D`,
`buildSnowflakeStrip_2D/_3D`, `buildBoomerangUnitCellStrip`/`_v2`) and several
`.asv` files, so unrequested additions compound an existing sprawl problem.
Downstream solvers also index COMSOL domains by hard-coded number, so a
speculative extra geometry feature can silently misassign materials.

**How to apply:**
- When asked for a builder targeting the mechanical path, *omit* the optical air
  region entirely rather than gating it behind `P.add_airDisk` — an extra domain
  renumbers the diamond that `runBands_2D` pins as `b_domind = 1`. Note the
  omission in the header instead.
- Verify with `checkcode` before reporting. MATLAB R2025b is installed at
  `C:\Program Files\MATLAB\R2025b\bin\matlab`; `matlab -batch "checkcode(...)"`
  runs headless without COMSOL, and validation/warning paths can be exercised by
  calling a builder with `model = []` (it throws at the first `model.*` call,
  after all the pure-MATLAB checks have run).
- To exercise a builder's *local* functions (param validation, warn-only
  geometry checks) with many inputs, copy the local function's text out of the
  `.m` file into a temp file on the path and call it directly. Better coverage
  than the `model = []` trick, which stops at the first `model.*` call, and it
  runs the real code rather than a re-typed copy.
- Do not attempt a COMSOL run unless asked; LiveLink must already be running.
- Verified 2026-11-21: headless `matlab.exe -batch` has NO LiveLink on the path
  (`mphstart`/`mphmeasure` absent, `com.comsol.model.util.ModelUtil` unresolvable),
  so COMSOL-dependent tests skip there. `omc-comsol-chang/RibUnitCellTest.m` uses
  an assumption for this; COMSOL-tagged tests in it have never actually executed.
- The MATLAB MCP tools (`mcp__matlab__*`) are not always exposed to subagents;
  fall back to the `matlab -batch` CLI above rather than reporting unverified.
- Nuance: "don't touch other files" is about *unrequested* edits. When a task
  names the files to change (e.g. add a `P` field read by `runBands.m` and
  documented in a `test_*.m`), edit them — but keep new behaviour opt-in with a
  default that reproduces the old path exactly, since many sweep/optimize
  scripts call the same solver.
- Proposal-only mode (since 2026-11-22): headless `matlab.exe -batch "<multi-line
  script>"` from WSL can check numerics inline without writing files, but
  function files can't be exercised that way. State clearly which parts were
  executed and which weren't.
- Inline `matlab -batch` can also run real repo functions (`findGaps`,
  `validateRibParams`, `surrogateCrossBands`) and a headless `bayesopt` with
  `PlotFcn = []`, `Verbose = 0`, without writing files.
- Update 2026-11-23: when the user has a COMSOL 6.3 server running, headless
  `matlab.exe -batch` CAN reach it via `addpath('C:\Program Files\COMSOL\COMSOL63\Multiphysics\mli'); mphstart(2036)`.
  The 2026-11-21 "no LiveLink" note only means the mli path isn't added by default.
  That server is the user's own live session (their `model` is on it), so create
  temp models with `ModelUtil.uniquetag`, remove them afterwards, never touch
  `model`, avoid `mphopen` on repo files (it leaves .lock files), and report
  any solves you ran.
- Inline checks of classifier/matcher logic with synthetic inputs catch real flaws
  (2026-11-23: a 0.45/0.55 hybrid was classified 'lost' instead of 'ambiguous'). In
  proposal mode, run those bodies inline before delivering, and say which final
  version was or wasn't re-run.

Related: [[user-role]], [[feedback-preserve-legacy-geometry-behavior]].
