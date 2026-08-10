---
name: user-role
description: User is a researcher building a COMSOL-driven phononic/photonic crystal optimizer; expert in the MATLAB pipeline, treats it as physics ground truth
metadata:
  type: user
---

Changjin (changjin@g.harvard.edu) works on diamond optomechanical-crystal
simulation. Two parallel pipelines exist in this repo: a mature MATLAB +
COMSOL LiveLink pipeline (`omc-comsol-chang/`) and newer Python/MPh
optimization loops (`python-scripts/` 1D nanobeam, `python_script_2D/` 2D
boomerang cell).

How to collaborate:
- Treat the MATLAB implementation as the physics/geometry/BC ground truth.
  When Python and MATLAB disagree, MATLAB is almost always right and the
  Python side is the port that drifted.
- The user already knows COMSOL's Java API, Voigt-notation stiffness
  tensors, Floquet/Bloch BCs, and eigenfrequency sweeps. Skip tutorials;
  go straight to specifics.
- The user often arrives with a numbered list of findings they have already
  established and wants them *verified and corrected*, not re-derived. Do the
  independent check anyway (compute the tensor, count the solves), and say
  plainly when one of their claims is wrong or incomplete.

See [[feedback-audit-style]].
