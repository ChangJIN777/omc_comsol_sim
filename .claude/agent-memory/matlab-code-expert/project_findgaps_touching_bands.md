---
name: project-findgaps-touching-bands
description: findGaps reports a spurious gap when two bands merely touch; matters for any surrogate or fitness scored on ds.sym
metadata:
  type: project
---

`findGaps` reports a spurious gap when two bands merely touch, because of the strict comparison at `findGaps.m:31`. `surrogateCrossBands` shows this in `ds.sym` for gapless designs (about 0.893 GHz). It's harmless for `bayesopt_cross`, which reads `ds.full`.

**Why:** found while checking surrogates for the hole/rib bayesopt proposal on 2026-11-22.

**How to apply:** any surrogate or fitness that scores `ds.sym` must make band groups overlap, not just touch, for gapless designs. Consider fixing the comparison in `findGaps` if real COMSOL bands can also touch.

Related: [[project-bayesopt-unitcell-proposal]].
