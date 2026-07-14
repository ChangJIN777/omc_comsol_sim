# Preliminary strategy: 1D diamond OMC unit-cell bandgap optimization

**Scope of this document.** A concrete starting strategy for optimizing a single
diamond optomechanical-crystal (OMC) *unit cell* — rectangular cross-section
nanobeam with one elliptical hole per period — toward a simultaneous optical
bandgap near 1550 nm and a mechanical bandgap near 8 GHz. It covers what the
field has done, which tools to use (and why), the diamond-specific physics, the
objective/loop design, and a staged plan. A runnable repo implementing this plan
accompanies the document.

---

## 1. What the field has done (1D nanobeam OMCs)

The canonical 1D OMC is a suspended nanobeam with a periodic row of holes, the
architecture Painter's group established in silicon: the same 1D photonic-crystal
mirror that opens an *optical* bandgap also opens a *mechanical* bandgap, so a
central defect can co-localize a ~200 THz optical mode and a few-GHz acoustic
mode with strong radiation-pressure coupling. The mirror cell is the object you
optimize first; the cavity defect (quadratic taper of period/width/hole size)
comes later.

Two facts from the literature shape the plan:

First, **diamond OMCs already exist and sit right at your targets.** Burek,
Lončar and co-workers (Harvard) demonstrated diamond OMCs with optical modes at
λ ≈ 1529 nm (≈196 THz) and GHz mechanical modes in a 1D nanobeam — though their
beams used a *triangular* cross-section from angled etching. Your rectangular
cross-section is the more standard thin-film geometry and is well suited to a
clean width/thickness parameterization.

Second, and most important for your expectations: **a plain rectangular beam with
elliptical holes generally does not have a *complete* phononic bandgap.** A
complete 1D phononic gap usually requires extra geometric degrees of freedom —
"cross" or "stub" shaped cells, or strong width modulation (Gomis-Bresco /
Navarro-Urrios et al., *Nat. Commun.* 2014). What the simple ellipse-in-beam cell
*does* reliably provide is a **symmetry-restricted ("quasi") bandgap**: a gap for
one mechanical mode family of a given symmetry (the even–even "breathing"
family), even while modes of other symmetry classes cross the gap. This is
exactly the relaxation you said you're comfortable with, and it is how the
original silicon OMC mirror cells are actually designed — the gap is taken within
a chosen parity/symmetry sector, not over all modes.

**Consequence for the design.** Define the mechanical gap on the breathing
(even-y, even-z) family. If later you want a *complete* phononic gap, the cheapest
geometric upgrade is to add a "stub"/cross feature or width modulation to the
cell — keep that as a fallback parameter set, not the starting point.

---

## 2. Tooling — can this be done in pure Python?

**Optical: yes, comfortably.** The community-standard tool for photonic *band
structures* is **MPB (MIT Photonic Bands)** — a fast plane-wave eigensolver with
a Python API, used routinely by the OMC groups for exactly this mirror-cell
calculation (sweep Γ→X, read the TE-like gap). For a single unit cell it runs in
seconds. The trade-off: MPB assumes *ideal* geometry (sharp corners, vertical
sidewalls) and a piecewise-constant dielectric, so it does not represent
fabrication corrections (fillets, sidewall angle, surface layers).

Other optical options and why they rank lower for you:

- **FDTD (Meep / Lumerical FDTD / Tidy3D)** — better for cavity Q and spectra,
  overkill and slower for a unit-cell band structure.
- **Guided-mode expansion (legume)** — very fast, but it is built for photonic-
  crystal *slabs*; a nanobeam is a 1D-periodic waveguide, and GME does not
  cleanly handle arbitrary cross-sections or the fabrication detail you care
  about. Not recommended as the primary tool here.
- **COMSOL Wave Optics (FEM eigenfrequency with Floquet BC)** — slower than MPB
  but represents the *real* fabricated geometry, and uses the *same* mesh/geometry
  as the mechanical model.

**Mechanical: pure Python is possible but not turnkey.** There is no
MPB-equivalent mature, pip-installable phononic band-structure package. The
honest options are (a) write a Bloch-Floquet 3D-elasticity eigensolver yourself
in FEniCSx/scikit-fem (works, but you build and validate the complex
phase-shifted periodic BC, anisotropic diamond stiffness, and SLEPc eigensolve),
or (b) use **COMSOL Solid Mechanics, Eigenfrequency study with Floquet
periodicity** — the de-facto standard in the OMC literature. Since you already
have COMSOL 6.2 with the MATLAB LiveLink, (b) is clearly the right call: it is
faster to get correct, handles anisotropy and fabrication geometry natively, and
shares one geometry with the optical model.

### Recommendation

> **Optical:** MPB as the fast pre-screen on ideal geometry; **COMSOL Wave
> Optics** for fabrication-aware confirmation of finalists.
> **Mechanical:** **COMSOL Solid Mechanics + Floquet**, automated from Python
> (`MPh`) or MATLAB LiveLink.
> **One COMSOL geometry serves both physics** — so fabrication corrections
> (fillets, sidewall draft, etched surface) apply identically to optics and
> mechanics. This single-geometry coupling is the strongest reason to run optics
> in COMSOL too, and it directly answers your concern that the geometry is too
> involved for slab modal-expansion methods.

The optimizer itself stays in Python (Optuna → BoTorch); it calls whichever
backend you select. The repo lets you switch backends with a flag without
touching the loop.

A note on this sandbox: I could not run MPB/COMSOL here (no conda, PyPI blocked,
and COMSOL lives on your Mac), so the *running* validation I did was the numpy
TMM surrogate plus the full loop plumbing. The MPB and COMSOL drivers are written
and ready to run on your machine.

---

## 3. Diamond physics and target scaling

Material constants are in `configs/materials.yaml`: optical index n ≈ 2.40 at
1550 nm; density 3515 kg/m³; cubic stiffness C11 = 1076, C12 = 125, C44 = 577
GPa. Diamond is elastically anisotropic — fix the crystal orientation relative to
the beam axis early (start with beam axis along ⟨100⟩) and keep it constant, since
orientation shifts the mechanical bands by a few percent.

**The two targets are mutually consistent at one lattice constant — which is the
whole point of an OMC.** Optical mid-gap normalized frequency for a 1D nanobeam
mirror is typically a/λ ≈ 0.3, so

  a ≈ 0.3 × 1550 nm ≈ 0.46 µm.

The mechanical breathing-mode frequency scales as roughly f ≈ v / (2 w_eff) with
a diamond shear velocity v ≈ 1.3 × 10⁴ m/s, so for a beam width w ≈ 0.7 µm

  f ≈ 1.3e4 / (2 × 0.7e-6) ≈ 9 GHz,

i.e. ~8 GHz is reachable in the same size cell that puts the optical gap at
1550 nm. That is why the starting bounds in `configs/bounds.yaml` are
a ∈ [0.40, 0.60] µm, w ∈ [0.45, 0.90] µm, hole semi-axes hx ∈ [0.08, 0.24] µm,
hy ∈ [0.10, 0.35] µm, thickness fixed at 0.22 µm for the first campaign.

---

## 4. Objective and the closed loop

Design vector (normalized, internal to the optimizer): u = (u_a, u_w, u_hx, u_hy)
∈ [0,1]⁴, mapped to physical dimensions in `geometry.py`. Thickness is fixed
first; promote to a 5th variable once the 4-D search is mapped.

Normalized gaps:

  G_o = (ω_o2 − ω_o1) / ((ω_o2 + ω_o1)/2),   G_m likewise (breathing family).

Score (maximize):

  S = min(G_o, G_m)
      − λ_o · max(0, 0.20 − G_o)²
      − λ_m · max(0, 0.20 − G_m)²
      − λ_of · ((f_o,c − f_o,target)/f_o,target)²
      − λ_mf · ((f_m,c − f_m,target)/f_m,target)².

The first term rewards making *both* gaps large; the squared hinge terms push
each gap past 20%; the frequency terms keep the gaps at 1550 nm / 8 GHz rather
than at some large gap in the wrong spectral range. Weights live in
`configs/targets.yaml`.

Feasibility is checked *before* any solve: minimum feature size, minimum
solid bridge between holes (a − 2·hx), minimum sidewall (w/2 − hy). Infeasible
candidates return a penalized record — they never crash or launch a simulation.

Loop (in `scripts/run_loop.py`):

  load cached results → init optimizer → for each iteration:
  ask(u) → evaluate (optical + mechanical) → save EVERY result → tell(score).

Every candidate (success, infeasible, or solver-failed) is persisted to SQLite
(with a JSONL fallback), so the loop is resumable and no simulation is ever lost.

---

## 5. Staged plan

1. **Loop plumbing (done, validated here):** geometry → surrogate optical gap →
   score → save → propose. Confirms feasibility gating, scoring, persistence,
   and resume work end-to-end with numpy only.
2. **Real optical pre-screen (MPB, on the Mac):** swap `--optical mpb`. Run a
   ~50–200 point Sobol/Halton sweep to map where a ≥20% TE-like gap sits near
   1550 nm. Cheap; identifies the feasible optical region and broken geometries.
3. **Mechanical in COMSOL:** build `comsol/omc_unitcell.mph` once (recipe in
   `comsol/README_template.md`), define the breathing-family gap via mode parity,
   then `--mech comsol`. Now S = min(G_o, G_m) is meaningful.
4. **Bayesian optimization:** switch the optimizer from Optuna (TPE) to BoTorch
   (constrained, optionally multi-objective max (G_o, G_m)) — same ask/tell
   interface, no loop changes.
5. **Fabrication-aware + robust:** turn on COMSOL optics on the same geometry,
   add fillets/sidewall angle, and optimize robustness
   S_robust = E_δ[S(x+δ)] − β·Std_δ[S(x+δ)] so finalists tolerate fabrication
   spread. Optionally add the stub/cross feature if a *complete* phononic gap is
   wanted.
6. **Hand off to cavity design:** use the best mirror cells as seeds for the
   quadratic taper / defect optimization (Q, g0) — out of scope for unit cells.

---

## 6. Open choices for you

- **Crystal orientation** (beam axis ⟨100⟩ vs ⟨110⟩) — affects mechanical bands
  by a few %. Pick one and fix it.
- **Complete vs symmetry gap** — start with the symmetry (breathing) gap as
  agreed; decide later whether a complete gap (stub/cross cell) is worth the
  extra parameters.
- **Optical backend for production** — MPB-only (fast, ideal) vs COMSOL optics
  (fabrication-faithful, one geometry). Given your fabrication-correction goal,
  COMSOL-for-both is recommended for finalists, MPB for bulk screening.

---

## Sources
- [Diamond optomechanical crystals — Burek et al., *Optica* 3, 1404 (2016)](https://opg.optica.org/optica/fulltext.cfm?uri=optica-3-12-1404&id=354754) ([arXiv:1512.04166](https://arxiv.org/abs/1512.04166))
- [Harvard Lab for Nanoscale Optics — Optomechanical Crystals](https://nano-optics.seas.harvard.edu/optomechanical-crystal)
- [Si₃N₄ nanobeam optomechanical crystals — arXiv:1411.5996](https://arxiv.org/pdf/1411.5996)
- [A one-dimensional optomechanical crystal with a complete phononic band gap — *Nat. Commun.* 5, 4452 (2014)](https://www.nature.com/articles/ncomms5452)
- [Engineering multiple GHz mechanical modes in OMC cavities — arXiv:2208.00890](https://arxiv.org/pdf/2208.00890)
- [Phononic bandgap nano-acoustic cavity with ultralong phonon lifetime — arXiv:1901.04129](https://arxiv.org/pdf/1901.04129)
- [MPB (MIT Photonic Bands) — NanoComp/mpb](https://github.com/NanoComp/mpb)
- [Spin-embedded diamond optomechanical resonator, Q_m > 1e6 — arXiv:2508.05906](https://arxiv.org/pdf/2508.05906)
