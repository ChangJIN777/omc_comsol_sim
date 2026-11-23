---
name: project-bayesopt-unitcell-proposal
description: 2026-11-22 bayesopt proposal for hole + rib unit cells (shared driver + adapters). Check whether it was implemented before re-proposing; decisions open
metadata:
  type: project
---

On 2026-11-22 a proposal was delivered for `bayesoptUnitCell.m` (shared driver) + `unitCellAdapter_hole.m` / `unitCellAdapter_rib.m` + entry scripts `bayesopt_hole.m` / `bayesopt_rib.m`, templated on `bayesopt_cross.m`. Feasibility functions `holeFabOK` / `ribFabOK` share crossFabOK's metres-in/broadcast contract. Both cells score `ds.sym`.

**Why:** the user wanted Bayesian optimization for both unit cells. The request came via another Claude session ("matlab-code-expert-2 session routing").

**How to apply:** before re-proposing, check whether these files exist. Open decisions: rxtal (rib 0° vs 45°; hole uses 45°), the hole targetFreq (5/7/3 GHz) and baseline (618/600/380 vs 528.6/800/500), hole free vs fixed area, whether to search rib `a`/`d`, sigma/maxFreq, gapSelect ('bestScore' vs 'closest'), gap ratio vs absolute Hz, and the hole sidewall floor (148 vs 150 nm).

Related: [[project-band-tracking-proposal]], [[project-findgaps-touching-bands]].
