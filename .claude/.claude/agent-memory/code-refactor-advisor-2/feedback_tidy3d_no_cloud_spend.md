---
name: tidy3d-no-cloud-spend
description: Never submit to the Tidy3D cloud or spend FlexCredits when working on python-scripts optical engine; test offline only
metadata:
  type: feedback
---

Never upload to or run on the Tidy3D cloud from an agent session: no `--run`, no `--estimate-cost` (it uploads drafts), no `web.Batch` / `web.run`.

**Why:** Tidy3D charges FlexCredits per simulation; the user set this as a strict constraint when the `optical: tidy3d` engine was added (2026-11-23).

**How to apply:** Test `python-scripts/src/optical_tidy3d.py` and `diamond_bandStructure_sim.py` only offline: `tests/test_optical_tidy3d.py` mocks `tidy3d` in `sys.modules`. The `comsol_sim` conda env has no tidy3d; the `tidy` env (/mnt/c/Users/alber/anaconda3/envs/tidy/python.exe) has tidy3d 2.12.2 and may be used for OFFLINE checks only (import, inspect, build_simulations, geometry.inside; block sockets). Do not install packages without asking. A bare MagicMock hid a real bug once (`td.web` does not exist until `import tidy3d.web`): test fakes must mimic the real package layout, and real-API assumptions should be checked in the `tidy` env. Related: [[optical-engine-comparability]].
