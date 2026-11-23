---
name: project-band-tracking-proposal
description: Band-tracking proposal for optimize_rib_unitCell.m (2026-11-22). Check whether it was implemented before re-proposing; key design decisions are still open
metadata:
  type: project
---

On 2026-11-22 a proposal was delivered for tracking one mechanical band during `optimize_rib_unitCell.m`. It is opt-in via `OPT.track` / `P.exportModes`, and proposes new helpers `ribModeGrid`, `exportBandModes`, `bandOverlap`, `modeFingerprints`, `connectBandsAlongK`, `trackedBandGap`, `bandTrackCandidate`, `trackBand` and `surrogateRibBandsTracked`, plus a `BandTrackingTest` class. Matching combines mode-shape overlap (MAC) on a normalized grid, parity veto, fingerprints, a golden mode, and global assignment with `matchpairs`.

**Why:** the nearest-gap fitness picks gaps by frequency alone, so it can silently hop to a different band. `runBands` returns frequencies only, in COMSOL's order, with no sort and no eigenvectors.

**How to apply:** before proposing band tracking again, grep for those helper names to see whether it was implemented. Decisions the user still has to make: which band edge (lower or upper at X), anchor k-point (X or Γ), the fitness (`'margin'` vs gap size), whether bridge solves and `refineK` are acceptable, and `rxtal` (the parity veto is valid only at 0° or 45°).

Related: [[project-visible-rib-template]], [[feedback-scope-discipline]].

**Update 2026-11-23, user decisions:**
- Track both gap-edge bands, identified at X.
- Make it a standalone tool; it is NOT part of the optimizer for now.
- Resolve ambiguous matches by bridging (midpoint solves).
- rxtal = 0 in all rib scripts.
A revised proposal was requested on this basis.

**Update 2026-11-23 (later):** a revised proposal was delivered and approved for implementation: the standalone `trackRibBands(designs, opts)` plus `matchBandModes`, `plotRibBandTrack`, a redesigned `surrogateRibBandsTracked` and a rewritten `BandTrackingTest`. `trackBand`, `trackedBandGap` and all optimizer edits were dropped. Defaults: follow the continuous eigen-branch through avoided crossings (character via `macGold_*`), bridge only on 'ambiguous', refineK off. Still open: branch vs character following, bridging on 'lost', a default targetFreq, re-solving when a band leaves the eigenvalue window.
