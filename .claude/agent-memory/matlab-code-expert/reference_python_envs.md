---
name: reference-python-envs
description: Python environments for python-scripts: WSL python3 lacks numpy; the Windows conda env comsol_sim has MPh/numpy/optuna but no pytest
metadata:
  type: reference
---

- WSL `python3` has no numpy.
- The Windows conda env `C:\Users\alber\anaconda3\envs\comsol_sim` has MPh 1.3.1, numpy, scipy and optuna, but no pytest. Run `tests/test_pipeline.py` directly with python.
- `mph.start()` launches its own COMSOL session, separate from the live LiveLink server on port 2036, and needs a free licence seat. Ask the user before starting one.

Related: [[project-python-rib-integration]].
- WSL has no `unzip`. Read `.mph` files with `python3 -I -c "import zipfile..."`.
- Windows MATLAB can `checkcode` a WSL `/tmp` file via the `wslpath -w` UNC path.
- The conda env `tidy` has tidy3d 2.12.2; `comsol_sim` does not.
