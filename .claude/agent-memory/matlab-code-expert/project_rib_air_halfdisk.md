---
name: project-rib-air-halfdisk
description: 2026-11-23 proposal to replace the rib builder's air box with a parametric half-disk (air_gap/air_xc/air_r); check if implemented and whether Q rose
metadata:
  type: project
---

On 2026-11-23 a proposal was delivered to replace the air box in `python-scripts/comsol/build_omc_unitcell_rib.m` with a half-disk:
- parameters `air_gap = 1.55 um`, `air_xc` (0 for rect, −wd/2 for wedge) and `air_r = sqrt((w/2)^2+(th/2-air_xc)^2)+air_gap`;
- built as a Cylinder minus a Block;
- `selScat` becomes a DifferenceSelection, cross-checked against probe boxes;
- stale templates are detected by the absence of `air_r`.

**Why:** guided modes showed Q of 200–340 with the box boundary only about 2 decay lengths from the slab.

**How to apply:** before re-proposing, check whether it was implemented, and whether Q rose in the re-run optical validation solve.

Related: [[project-python-mph-templates]], [[project-python-rib-integration]].
