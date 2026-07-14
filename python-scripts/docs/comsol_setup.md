# Running the COMSOL loop on your MacBook Pro M2 — step by step

## First, the honest constraint
My sandbox is an **isolated Linux container** with no network route to your
laptop and no COMSOL install or license. It **cannot** reach or run COMSOL on
your Mac — that boundary can't be opened. So the simulation loop runs **on your
Mac**, natively, right next to COMSOL. I write and validate the Python; you run
it locally. This is also the normal, fast way to do it: the optimizer is tiny;
COMSOL does the heavy lifting in-process.

There are two ways for Python (or MATLAB) to drive COMSOL. Pick **A** if you want
the pure-Python loop (recommended, matches your preference); **B** is a rock-solid
fallback since you already have MATLAB.

---

## Path A — Python drives COMSOL via MPh (recommended)

`MPh` is a Python package that talks to COMSOL through its Java API (via JPype).
COMSOL 6.2 supports Apple Silicon natively, and MPh works with it. The **one**
thing that bites people: the **Python architecture must match the COMSOL build**,
because JPype loads COMSOL's Java VM into your Python process.

### Step 1 — find out which COMSOL build you have
COMSOL 6.2 ships as either a **native Apple Silicon** build or an **Intel**
build (run under Rosetta). Check:
```bash
# About box: COMSOL > Help > About COMSOL Multiphysics  (look for "Apple Silicon")
# or inspect the binary:
file /Applications/COMSOL62/Multiphysics/bin/comsol
```
`arm64` in the output → native Apple Silicon. `x86_64` → Intel build.

### Step 2 — install a matching Python
Use a clean environment. Miniforge gives you an arm64 Python by default.
```bash
# install miniforge (arm64) if you don't have conda:
#   https://github.com/conda-forge/miniforge  (Apple Silicon installer)

# --- if COMSOL is native Apple Silicon (arm64) ---
conda create -n omc python=3.12
conda activate omc

# --- if COMSOL is the Intel build ---
# run an x86_64 Python under Rosetta instead:
#   CONDA_SUBDIR=osx-64 conda create -n omc-x86 python=3.12
#   conda activate omc-x86
#   conda config --env --set subdir osx-64
```

### Step 3 — install the Python packages
```bash
pip install MPh numpy pyyaml pandas matplotlib optuna
# MPh pulls in JPype1 automatically.
```

### Step 4 — run the diagnostic (this is the "make sure it works" step)
```bash
cd "<this repo>"
python scripts/check_comsol.py
```
It checks, in order: Python arch, MPh/JPype import, `mph.start()` (COMSOL
discovery + license), a trivial parameter/evaluate, and your template if present.
Each line prints `[ OK ]` / `[warn]` / `[FAIL]` with a fix hint. Get all the
critical steps green before running simulations.

Common fixes if `mph.start()` fails:
- Close any open COMSOL GUI/sessions (license contention).
- Non-standard install path → `export COMSOL_HOME=/Applications/COMSOL62/Multiphysics`.
- Architecture mismatch (step 1 ≠ step 2) → recreate the env with the right subdir.

### Step 5 — build the unit-cell template once
Follow `comsol/README_template.md` in the COMSOL GUI and save as
`comsol/omc_unitcell.mph` (parameters `a,w,t,hx,hy,kF`; mechanical
Eigenfrequency study "Study 1"; optional optical study "Study 2"). Re-run
`check_comsol.py` — it will now also report the template's studies.

### Step 6 — run it
```bash
# mechanical only, single cell:
python scripts/run_one.py --u 0.5 0.6 0.5 0.6 --optical surrogate --mech comsol

# full loop (MPB optical pre-screen + COMSOL mechanical):
python scripts/run_loop.py --n-init 40 --n-iter 120 --optical mpb --mech comsol
```

---

## Path B — MATLAB LiveLink (fallback, very robust on your setup)

You have MATLAB, and COMSOL 6.2 includes a LiveLink update for Apple Silicon.
1. Launch the linked session: from COMSOL, *"COMSOL Multiphysics with MATLAB"*,
   or start `comsol mphserver matlab`.
2. In MATLAB, with `comsol/omc_unitcell.mph` built:
   ```matlab
   bands = omc_mechanical_livelink(500e-9,700e-9,220e-9,150e-9,250e-9,15,12);
   ```
   (`comsol/omc_mechanical_livelink.m`). It sweeps Γ→X, returns
   `[n_k × n_bands]` eigenfrequencies, and prints the gap near 8 GHz.
3. To keep the Python optimizer as the driver, have it shell out to MATLAB
   (`matlab -batch ...`) per candidate and read back a results file — I can wire
   this if you choose Path B.

---

## Which to use?
Start with **Path A** (matches your Python preference, one process, fast). If the
JPype/arch setup gives trouble, **Path B** will almost certainly work since your
MATLAB+COMSOL link is already a supported, tested combination. Either way the
geometry, scoring, database, and plotting code are shared — only the solver call
differs.

## Sources
- [COMSOL 6.2 macOS Apple Silicon Native Support (KB 1307)](https://www.comsol.com/support/knowledgebase/1307)
- [LiveLink for MATLAB on macOS Apple Silicon, COMSOL 6.2 (KB 1308)](https://www.comsol.com/support/knowledgebase/1308)
- [MPh — Pythonic scripting interface for COMSOL](https://github.com/MPh-py/MPh)
- [MPh installation docs](https://mph.readthedocs.io/en/stable/installation.html)
