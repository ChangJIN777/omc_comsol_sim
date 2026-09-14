---
name: feedback-preserve-legacy-geometry-behavior
description: When porting/extending a COMSOL geometry builder, reproduce existing quirks bit-for-bit and surface problems as warnings instead of fixing them silently
metadata:
  type: feedback
---

When extending or porting an existing geometry builder, reproduce the original
behaviour exactly — even behaviour I have demonstrated is buggy — and report the
defect as a build-time **warning** plus a doc-block note, rather than correcting
it.

**Why:** simulation results and fabricated devices already exist for the current
geometry. A "fix" that changes the mesh or the filleting silently invalidates
comparisons against that body of work. The user's words: reproduce it
"bit-for-bit so results stay comparable with already-simulated/fabricated
designs." Correctness of the *reported* diagnosis still matters — they want to
know, they just do not want the geometry moving underneath them.

**How to apply:**
- Change only what is structurally required by the new use (e.g. per-instance
  coordinate offsets and unique feature tags when looping a single-cell motif
  into an N-cell strip). Leave radii, tolerances, conditions, and orderings alone.
- Where I found a real defect, add: (a) a pure-MATLAB pre-build check that WARNS
  and names the offending entities, (b) a KNOWN ISSUE block in the header saying
  it is reproduced deliberately and what the correct value would be.
- Warnings, never errors, for reproduced legacy defects. Hard errors are still
  right for genuinely invalid input (e.g. a dimension that disconnects the solid).
- Concrete precedent: `buildCrossStrip.m` reproduces `buildCrossUnitCell`'s
  `addFillet` annuli `[w/2, 1.5w]` and `[h/2 +/- 5nm]` unchanged, even though the
  second selects nothing for typical parameters and the first bleeds into
  neighbouring cells.

Related: [[user-role]], [[feedback-scope-discipline]].
