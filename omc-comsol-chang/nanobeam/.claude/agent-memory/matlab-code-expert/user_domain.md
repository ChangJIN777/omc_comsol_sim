---
name: user-domain
description: User works on diamond optomechanical-crystal nanobeam cavity FEM in COMSOL LiveLink for MATLAB; expert level, wants physics reasoning not hand-holding
metadata:
  type: user
---

Works on **diamond optomechanical crystal (OMC) nanobeam cavities** simulated with
COMSOL LiveLink for MATLAB: optical modes near 1550 nm, mechanical breathing modes
near 7 GHz, optomechanical coupling g_OM, radiative mechanical Q via phononic
shields plus a mechanical PML, and SiV strain coupling.

How to be useful to them:

- They reason at the level of "is this mode actually inside the bandgap", "is the
  PML stretch keyed to the right wave speed", "does this selection resolve to the
  faces I think it does". Answer at that level. Do not explain what a PML is.
- They want the **numbers**, computed, not estimated in prose. Run
  `CreateNanobeamGeom` / `LoadMaterialParams` headless via `matlab -batch` (both are
  pure MATLAB, no COMSOL needed) and put real values in comments rather than
  approximations.
- They are comfortable being told plainly that something is infeasible or that an
  approximation is uncontrolled. Say so with the arithmetic attached.
- Never attempt a COMSOL solve on their behalf - it needs a live LiveLink session
  and takes a very long time. Geometry arithmetic run headless is welcome.

See [[feedback-house-style]] and [[project-crossshield-optics]].
