---
name: python-scripts-safety-rules
description: Standing rules when working on python-scripts/ (Tidy3D cloud, COMSOL server, result files, line endings, test envs)
metadata:
  type: feedback
---

When working on `python-scripts/`: never submit to the Tidy3D cloud (no `--run`, no `--estimate-cost`, no `web.Batch`/`web.run`/upload), never start COMSOL/MPh or touch the user's LiveLink server, never install packages, never commit/push unless asked, and never write to `results/` (sqlite databases, `results/tidy3d/*`, out_json files).

**Why:** Tidy3D runs cost FlexCredits and even a cost estimate uploads a draft; the user often has a live optimization campaign writing to `results/` and a COMSOL server in use while the agent works.

**How to apply:** Test with fakes or with the real tidy3d package offline (sockets blocked); write test output to a temp dir and delete it (use literal paths, `rm -rf $VAR/*` is blocked). Use the Windows conda pythons from WSL: env `comsol_sim` (no tidy3d) and env `tidy` (tidy3d). Keep .py/.yaml/.md CRLF and README ASCII. Fakes of tidy3d must mimic the real module layout (`td.web` is not an attribute; a MagicMock once hid that bug). When a relayed numeric example contradicts the stated formula, implement the formula and flag the mismatch in the report.
