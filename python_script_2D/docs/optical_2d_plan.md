# Plan: 2D optical (photonic) bandgap -- future work

Status: **not implemented**. `src/optical_comsol_2d.py` is a documented stub
that raises `NotImplementedError`. This file is the fuller design writeup;
the module docstring carries a condensed version next to the code itself.

## Why this isn't just a copy of the 1D optical backend

`python-scripts/src/optical_comsol.py` (and `optical_mpb.py`,
`optical_surrogate.py`) handle a 1D nanobeam: one Bloch direction (kz along
the beam axis), one relevant guided-mode family ("TE-like" = Ey-dominant,
y-even), and no light-line bookkeeping because MPB/the surrogate already
restrict to guided (below-index) solutions by construction.

This project's photonic crystal is 2D-periodic in-plane, exactly like the
mechanical problem (same hexagonal Brillouin zone, same `(kx,ky)` sweep via
`src/acoustic_comsol_2d.py:bz_path`). Two things make the optical case harder
than the mechanical one, and harder than the 1D project's optical case:

1. **Light-line filtering at every k-point.** A finite-thickness slab
   surrounded by air only *guides* photonic-crystal-slab modes below the
   light line (`omega < c*|k_parallel|`); above it, "bands" computed by an
   eigenfrequency solver with an artificial outer boundary are numerical
   artifacts of that boundary, not real guided modes. The mechanical problem
   has no analog of this -- elastic waves in vacuum-clad diamond don't leak
   the way photons do, so `acoustic_comsol_2d.py` never needs to check
   anything like this.
2. **Two lattice directions instead of one.** The light line itself becomes
   a light CONE `omega < c*sqrt(kx^2+ky^2)`, evaluated at every swept
   `(kx,ky)` along the BZ path, not a single per-k-point scalar comparison.

## Planned physical setup

Extend `comsol/omc2d_boomerang.mph` (do not fork a separate geometry -- the
whole point, per `python-scripts/comsol/README_template.md`'s stated
rationale, is that ONE geometry serves both physics so fabrication details
like fillets apply identically to both):

- Physics: **Electromagnetic Waves, Frequency Domain** (or Mode Analysis),
  added to the same extruded boomerang-cell domain.
- Same Floquet periodicity face pairs as the mechanical studies
  (`kFloquet = (kx, ky, 0)`), so both physics share the exact same BZ sweep.
- Air cladding above/below the slab, with either a scattering boundary
  condition or a PML, so leaky (above-light-line) solutions are damped
  rather than spuriously reflected/confined by an artificial hard boundary --
  this matters more here than for a simple guided-mode 1D waveguide because
  every k-point needs a correct answer about which modes are real.
- z-mirror boundary condition (PEC/PMC on the z=0 plane, the photonics
  analog of the mechanical Symmetry/Antisymmetry pair) to classify modes into
  z-even ("TE-like", by the convention already declared in
  `configs/targets_2d.yaml`) vs z-odd families -- two studies, e.g.
  `"opt evenz"` / `"opt oddz"`, mirroring `STUDY_EVENZ`/`STUDY_ODDZ` in
  `acoustic_comsol_2d.py`.
- Eigenfrequency study per parity, searching `n_bands` modes near the target
  (`configs/targets_2d.yaml: optical.target_wavelength_nm`).

## Planned software interface

```python
# src/optical_comsol_2d.py, once implemented
def run_optical_comsol_2d(g, *, n_per_segment=15, n_bands=8,
                          unitcell="hexagonal", template=_TEMPLATE):
    """Mirrors run_mechanical_comsol_2d's signature and return shape:
    dict(k_norm, kx, ky, freqs_evenz [n_k,n_bands], freqs_oddz [n_k,n_bands],
         a). PLUS a light-line mask per (k, band) -- e.g. an additional
    `guided_evenz`/`guided_oddz` boolean array of the same shape, computed as
    freq < c * sqrt(kx**2 + ky**2) / (2*pi) at each swept point -- so
    src/objective2d.py can NaN-out (or otherwise exclude) non-guided entries
    before calling bandgap.gap_near_frequency, the same way the 1D project's
    bandgap extraction only ever sees physically meaningful bands.
    """
```

Then in `src/objective2d.py`:
- Add an `_optical_gap` helper analogous to `_mechanical_gap`, applying the
  light-line mask (set masked entries to NaN) before calling
  `bandgap.gap_near_frequency(..., nan_safe=True)`.
- Wire `optical_backend="comsol"` to actually call `run_optical_comsol_2d`
  instead of raising -- the `try/except` around it in `evaluate_candidate`
  already exists and does not need to change.
- Flip `require_opt: true` in `configs/targets_2d.yaml`-consuming run configs
  once validated on a few hand-picked geometries.

## Validation plan before trusting the optimizer with it

1. Build the template additions above; solve one geometry by hand at a few
   k-points; sanity-check the light line visually against a plotted band
   diagram (frequency vs `sqrt(kx^2+ky^2)`) before wiring it into the loop.
2. Add `tests/test_pipeline_2d.py` cases for the light-line mask logic using
   synthetic `(kx,ky,freq)` triples (no COMSOL needed) once that function
   exists, the same way `test_bz_path_hexagonal_shape_and_closure` tests
   `bz_path` without a COMSOL server.
3. Cross-check a handful of candidates' optical gap against the MATLAB
   pipeline (`omc-comsol-chang/solveOpticalBands.m`, if/when a 2D-hexagonal
   optical study exists there) before trusting the Python driver's numbers.
