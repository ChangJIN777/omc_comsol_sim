---
name: optical-engine-comparability
description: User decisions of 2026-11-23 on diamond index, cross-section (xsect) and thickness-parity rule, and what is still unverified
metadata:
  type: project
---

User decisions (2026-11-23) for python-scripts and omc-comsol-chang:
- Diamond index is n = 2.386 everywhere (Python and MATLAB), even in scripts at visible wavelengths.
- `xsect: rect | wedge` for both cell types; the mid-thickness mirror rule applies to OPTICS only (rect: thickness-even family, wedge: all modes); mechanics is never filtered.

**Why:** makes COMSOL and Tidy3D optical gaps comparable; the hole templates differ (omc_unitcell.mph is the rectangle, omc_unitcell_iso.mph has the wedge).

**How to apply:** still unverified as of that date, so do not present as validated: the rib wedge template (builder written blind, never built), the COMSOL field-parity filter (never run on a live model), and every Tidy3D simulation (tidy3d not installed). The rib wedge is, by user decision, the roof of the FILLETED outline (depth = tan30 * distance to outline; deepest under the junction): exact planes/cones in the COMSOL builder, a 12-layer stack in Tidy3D. A rib wedge .mph built before that (straight prisms, no parameter wdj) is stale and refused by template_path. Related: [[tidy3d-no-cloud-spend]].

Rib optical scoring, user decision 2026-11-23 ("score version 2"): no window on the gap centre; the optical gap enters the rib score as G_o * exp(-((f_c - f_target)/optical_tol)^2) with `rib.optical_tol_thz` = 2 THz, and the old quadratic optical centre penalty is dropped. I flagged that 2 THz (~1 % of 193 THz) is about the solver's own frequency accuracy and makes a needle-like landscape; the user asked for exactly this, so do not widen it unasked. Hole scoring is unchanged.
