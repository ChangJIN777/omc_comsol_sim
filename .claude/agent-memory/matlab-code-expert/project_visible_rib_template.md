---
name: project-visible-rib-template
description: 'rib' celltype (formerly visible_rib, merged with the LN rib 2026-11-21) - template provenance, user decisions, LN compatibility rules
metadata:
  type: project
---

`comsol_templates/visible_rib_omc.m` (COMSOL 6.3) is a GUI-session history.
Its final rib cells came in via `geom.insertFile('.../omc rib cavities/model.mph')`,
which does not export their settings. The visible set's fillet radii remain
unrecoverable (TODO, 40 nm borrowed from the IR set).
s = t = 100 nm are corroborated by the diamond rib sweep logs in
`test/diamond_ribUnitCell*/`.

User decisions, 2026-11-21. The device was first `'visible_rib'`, then the user had it
renamed to `'rib'` and merged with the LN rib cell: `buildRibUnitCell_LN.m` keeps its
name but holds the merged builder (helpers ribDefaultParams / ribLength /
validateRibParams / addRibCell / addRibPrisms; nanobeam createRibCavityGeom /
buildRibCavity; reference doc rib_README.md). Merged-device rules: P.d in radians,
fillets default 0 when d ~= 0 and are rejected with d ~= 0, z-symmetry rejected
with d ~= 0, explicit P.L > P.w > L(a) rule (P.ribDerived marks written-back values),
legacy solveBands_noSym / solveOpticalBands 'rib' file names kept.
Other decisions:
- th defaults to 110 nm = the template's EFFECTIVE thickness (z-offset
  defect: work plane at -125 nm, cut at -60 nm); airZOffset 60 nm reproduces
  its air-cylinder offset. Nominal design 120 nm.
- L(a) = L_ref + dLda*(a - a_ref): 400 nm @ 240 nm, slope -0.4998 (IR set).
- taperConvention default 'centerspacing' (matches CreateNanobeamGeom exactly,
  verified). 'cellwidth' is the only one that reproduces template positions.
- mEndBC default 'lowreflecting' (visible_rib only); rxtal default 0;
  rib-at-centre maps to holeatctr = 0 parity (oevenx forced 1).

**Why:** results get compared against the GUI model, so the 110 nm / offset
choice and convention switch matter.
**How to apply:** when touching visible_rib, keep defaults in
`ribDefaultParams.m` (the single source). nanobeam tests reach the
root helpers via `addpath(root,'-end')`. Ask for the source model.mph before
trusting fillet radii.

Related: [[feedback-preserve-legacy-geometry-behavior]], [[user-role]].
