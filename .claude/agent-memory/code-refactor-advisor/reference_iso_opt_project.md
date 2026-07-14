---
name: reference-iso-opt-project
description: Sibling Python OMC-optimization project on Desktop that uses a template-driven COMSOL geometry pattern worth borrowing
metadata:
  type: reference
---

A related project lives at `/Users/changjin/Desktop/Research/OMC iso design optimization Claude/` (Python + COMSOL via the `mph`/MPh library, not MATLAB LiveLink). It targets the SAME physics as this repo (1D diamond OMC unit cell, elliptical hole, optical gap near 1550 nm + mechanical breathing gap near 8 GHz) but is architected declaratively:

- Geometry topology is baked into ONE pre-built parametric `.mph` template (`comsol/omc_unitcell_iso.mph`) with named COMSOL parameters `a,w,t,hx,hy,kF`. Python never constructs geometry — it only does `model.parameter(name, val)` then solves.
- `src/geometry.py`: a `Geometry` dataclass + `u_to_geometry`/`phys_to_u` (normalized [0,1] <-> physical meters) + `check_feasibility` (min-feature-size gate).
- `configs/*.yaml`: `bounds.yaml` (variable ranges + feasibility), `targets.yaml` (objectives/weights), `materials.yaml`.
- A new design point = a different params dict / u-vector. Zero code edits.

**Why:** Useful reference when discussing how to make this repo's geometry generation fully parameter-driven (see the user's Task-3 goal: sweeps with no per-design-point code edits).
**How to apply:** When proposing sweep/optimization ergonomics for omc_comsol_sim, contrast against this template-driven + single-config pattern. Note it is a separate external dir (may move/change) — verify before relying on specific file contents.
