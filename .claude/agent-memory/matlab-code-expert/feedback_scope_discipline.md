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
- The MATLAB MCP tools (`mcp__matlab__*`) are not always exposed to subagents;
  fall back to the `matlab -batch` CLI above rather than reporting unverified.
- Nuance: "don't touch other files" is about *unrequested* edits. When a task
  names the files to change (e.g. add a `P` field read by `runBands.m` and
  documented in a `test_*.m`), edit them — but keep new behaviour opt-in with a
  default that reproduces the old path exactly, since many sweep/optimize
  scripts call the same solver.

Related: [[user-role]], [[feedback-preserve-legacy-geometry-behavior]].
