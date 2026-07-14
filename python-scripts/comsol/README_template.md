# Building the COMSOL unit-cell template (`omc_unitcell.mph`)

Build this ONCE in the COMSOL 6.2 GUI on your Mac, then the Python/MATLAB
drivers just set parameters and solve. One geometry serves BOTH the optical and
mechanical studies, so fabrication corrections (fillets, sidewall angle) are
applied identically to both physics -- this is the main reason to do optics in
COMSOL too rather than only MPB.

## 1. Global parameters  (Model > Parameters)
| name | expr | meaning |
|------|------|---------|
| a  | 500[nm]  | lattice constant (period along z) |
| w  | 700[nm]  | beam width (y) |
| t  | 220[nm]  | beam thickness (x) |
| hx | 150[nm]  | hole semi-axis along z |
| hy | 250[nm]  | hole semi-axis along y |
| kF | 0[1/m]   | Floquet wavevector along z (swept by driver) |

## 2. Geometry  (one period; meshable as a single unit cell)
1. Block: width = w (y), depth/length = a (z), height = t (x), centered at origin.
2. Elliptic cylinder (or revolved ellipse) through the thickness x, semi-axes
   hx (z) and hy (y), centered at origin; subtract (Boolean Difference) from the
   block to make the air hole.
3. (Fabrication option) add fillets to hole/edges; add a slight sidewall angle
   by drafting the block faces. Keep these parametric.

## 3. Material
Diamond: density 3515[kg/m^3]; for accuracy use the anisotropic elasticity
matrix (cubic): C11=1076[GPa], C12=125[GPa], C44=577[GPa]. Set the material
coordinate system so the beam axis z maps to your chosen crystal axis
(<100> to start). For optics set n = 2.40 (or relative permittivity 5.76).

## 4. Mechanical study  ("Study 1")
- Physics: **Solid Mechanics**.
- The two faces normal to z (z = +-a/2): **Periodic Condition > Floquet
  periodicity**, k-vector = (0,0,kF). The x and z... (free surfaces elsewhere:
  the beam is suspended, so y- and x-normal outer faces are FREE/traction-free).
- Study: **Eigenfrequency**. Search for `n_bands` (e.g. 12) eigenfrequencies
  "around" 8[GHz] (use the shift/search-around-value option).
- (Optional) define integration coupling `intop1` over the domain to evaluate
  `ux^2, uy^2, uz^2` so the driver can classify the breathing (even-even) family.

## 5. Optical study  ("Study 2")  -- optional, for fabrication-aware optics
- Physics: **Electromagnetic Waves, Frequency Domain** (or Mode Analysis).
- Floquet periodicity on the z-faces, k = (0,0,kF). Outer x/y faces:
  scattering / PML padding region of air so guided modes decay.
- Study: **Mode Analysis / Eigenfrequency** searching `n_bands` modes near
  193.4[THz] (1550 nm). Filter TE-like (Ey-dominant, y-even) modes.

## 6. Save as `comsol/omc_unitcell.mph`
The drivers (`src/acoustic_comsol.py`, `src/optical_comsol.py`,
`comsol/omc_mechanical_livelink.m`) reference parameter names a,w,t,hx,hy,kF and
study names "Study 1" / "Study 2". Keep these names or update the drivers.

## Sweep convention
The driver sets `kF = pi*s/a` for `s` in (0, 1] to walk Gamma->X, where
`s = k_z*a/pi` is the normalized wavevector and the zone edge X is at `s = 1`.
