---
name: feedback-audit-style
description: For audit/review tasks here - analysis only, no file edits, file:line anchors, and explicitly classify each finding as crash vs silent-wrong-answer
metadata:
  type: feedback
---

When asked to audit code in this repo: do NOT modify files. Return findings
directly in the final message (not as a written .md report). Anchor every
finding to `absolute/path.py:LINE`, and for each one state explicitly whether
the failure mode is a **crash** or a **silent wrong answer**.

**Why:** the user is auditing a numerical-simulation pipeline where a wrong
answer that still "runs" (e.g. an optimizer scoring every candidate on a
frozen geometry) is far more dangerous than an exception. They asked for that
classification by name.

**How to apply:** applies to any review of the COMSOL/Python coupling. Also:
quantify waste in concrete units the user can act on (number of eigenvalue
solves per candidate, GB of .mph, GHz of frequency error), not adjectives.
Verify claims numerically where cheap - e.g. recompute a rotated stiffness
tensor rather than eyeballing it.

See [[user-role]].
