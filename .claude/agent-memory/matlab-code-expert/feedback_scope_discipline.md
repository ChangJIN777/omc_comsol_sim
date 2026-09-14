---
name: feedback-scope-discipline
description: Deliver exactly the one file asked for — no edits to existing files, no test scripts, no speculative feature gates; ask instead of assuming
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
- Do not attempt a COMSOL run unless asked; LiveLink must already be running.

Related: [[user-role]], [[feedback-preserve-legacy-geometry-behavior]].
