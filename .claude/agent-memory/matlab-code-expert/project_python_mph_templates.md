---
name: project-python-mph-templates
description: python-scripts/comsol/omc_unitcell*.mph don't match README_template.md — 4 studies, half model, isotropic, and dset1 is optical in the _iso file
metadata:
  type: project
---

`python-scripts/comsol/omc_unitcell.mph` and `omc_unitcell_iso.mph` don't match `README_template.md`:
- They're y ≥ 0 half models with an air box.
- They have four studies: `std1` "mech sym", `std3` "mech antisym", `std2` "opt TE", `std4` "opt TM". There's no "Study 1"/"Study 2".
- Both use isotropic diamond (E = 1050 GPa, ν = 0.1), not the cubic C11/C12/C44 the README describes.
- `omc_unitcell.mph` has no solver sequence or dataset.
- In `omc_unitcell_iso.mph`, `dset1` is "opt TE/Solution 1".

**Why:** a driver that runs `std1` and then reads `mphglobal(...,'dataset','dset1')` can silently get optical data from the `_iso` file. Found 2026-11-23 while proposing the rib counterpart of `omc_mechanical_livelink.m`.

**How to apply:** read `smodel.json` / `dmodel.xml` inside the `.mph` zip rather than trusting the README. In build scripts, pin datasets by calling `createAutoSequences('sol')` in study order. A rib template + driver proposal (`build_omc_unitcell_rib.m`, `omc_mechanical_livelink_rib.m`) was delivered on 2026-11-23; check whether it was implemented before re-proposing.
