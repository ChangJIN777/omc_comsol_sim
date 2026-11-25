---
name: tidy3d-conda-env
description: Where a real tidy3d install lives for reading its source (the comsol_sim env has neither tidy3d nor rich)
metadata:
  type: reference
---

The Windows conda env `comsol_sim` (used for the tests) has no tidy3d and no rich.
A separate env `tidy` has tidy3d 2.12.2 + rich 13.9.4:
`/mnt/c/Users/alber/anaconda3/envs/tidy/Lib/site-packages/tidy3d`.
Read its source there to check Tidy3D APIs (e.g. `log.py`, `web/api/container.py`) instead of guessing.
Never import it to call `web.*`: any upload/estimate/run reaches the paid cloud.
