---
name: feedback-house-style
description: How to write code in this repo - hard errors over silent degradation, worked arithmetic in comments, warn-only audits, never edit shared pipeline files unasked
metadata:
  type: feedback
---

Rules for writing in `omc-comsol-chang/nanobeam/`, confirmed by how the user asks
for work and by what they accepted without pushback.

**1. A geometric or physical precondition gets a `error`, not a `warning`.**
**Why:** in this pipeline a wrong geometry does not crash - it produces a
perfectly plausible-looking Q or g_OM. An air cylinder that swallows the phononic
shield still solves. A boundary condition applied to an empty selection still
solves. Silence is the failure mode.
**How to apply:** guard with a hard error carrying the offending numbers and the
name of the `P` field to change. Reserve `warning` for things that are
approximations rather than mistakes, and make those audits **warn-only, never
geometry-altering** (see `auditCrossShield`, `assertNonEmptySel`).

**2. Comments carry worked numbers, not symbols.**
**Why:** the user checks the arithmetic. "the cylinder clears the shield" is not
reviewable; "x-gap = 8951.8 - 7401.8 = 1550.0 nm, and x is the ONLY separation
because y and z overlap" is.
**How to apply:** compute the values headless first, then write them into the
comment block with units in nm. Explain *why* a bound is load bearing and which
term would break if a parameter moved.

**3. Do not edit the shared pipeline files unless explicitly authorised.**
`SetupNanobeamFEM.m`, `SolveNanobeamFEM.m`, `CalcGOM.m`, `RunNanobeamFEM.m`,
`BuildNanobeamFEM.m` and anything in `omc-comsol-master/` are consumed by many
test scripts.
**Why:** a change there silently moves results for every existing caller. The
repo already documents known bugs in `BuildNanobeamFEM.m` that were left in place
for exactly this reason.
**How to apply:** treat those files as a contract to satisfy, not to change. If a
task seems to need one changed, stop and report instead. New behaviour goes in a
new builder plus `P`-field gates.

**4. New behaviour must be bit-identical when its gate is off.**
**Why:** existing mechanical-only cross-shield runs must not move when optics is
added.
**How to apply:** one boolean gate derived once (`useOpt = isfield(P,'solveOpt')
&& P.solveOpt`), every addition inside it, and verify by dry-running both
branches and diffing the console output.

**5. `checkcode` must be clean, measured against a `git show HEAD:` baseline.**
Both edited files should report nothing new. Pre-existing findings in other files
stay.

Also: SI units in `P`, `*1e9` only for display; `isfield` guards on every optional
field; local helper functions at the bottom of the file; path separators always
from `fullfile`/`filesep` (see CLAUDE.md for why).
