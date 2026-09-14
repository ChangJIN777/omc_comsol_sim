---
name: user-role
description: User is a photonics/optomechanics researcher building COMSOL LiveLink MATLAB geometry for diamond optomechanical crystals; deep domain and MATLAB fluency
metadata:
  type: user
---

Works on diamond optomechanical crystal (OMC) cavities — phononic/photonic band
structure and nanobeam FEM — driven from MATLAB via COMSOL LiveLink.

Fluency level: high. Reads COMSOL's Java API idioms directly, reasons about
Floquet/Bloch BCs, symmetry sectors, and lattice geometry without explanation.
Do not explain what a fillet or a periodic boundary condition is. Do explain
*specific* non-obvious findings in their own code (e.g. a DiskSelection whose
radius does not match the vertex radius it was meant to catch).

Reviews plans before implementation on non-trivial geometry work, and gives
precise, itemised decisions in response. Write plans that are decidable —
enumerate the open questions explicitly rather than picking defaults silently.

Cares about *comparability with already-simulated and fabricated designs*, which
repeatedly outranks code correctness when the two conflict. See
[[feedback-preserve-legacy-geometry-behavior]].

Runs optimization loops (`bayesopt_cross.m`, `*_optimize_sweep_*.m`) that call
geometry builders once per evaluation, so build-time cost and unguarded
`figure`/`disp` calls inside builders are real problems, not style nits.
