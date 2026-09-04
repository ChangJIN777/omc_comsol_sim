# `CalcGOM.m` — Optomechanical Coupling Rate

Reference documentation for `nanobeam/CalcGOM.m`. This file states **explicitly** every
equation the code evaluates, and maps each equation to the MATLAB/COMSOL expression that
implements it.

---

## 1. Purpose and Signature

```matlab
[ds, model] = CalcGOM(ds, model, oModes, mModes)
```

Computes the single-photon optomechanical coupling rate $g_0$ between each optical mode in
`oModes` and each mechanical mode in `mModes` of a diamond nanobeam photonic-crystal
cavity, as the sum of two physically distinct contributions:

$$
g_0 \;=\; g_{\mathrm{MB}} \;+\; g_{\mathrm{PE}}
$$

| Term | Name | Physical origin | Integral type |
|---|---|---|---|
| $g_{\mathrm{MB}}$ | Moving boundary (a.k.a. `MB`, `Bnd`) | Mechanical motion displaces the dielectric interface | **Surface** integral over beam boundaries |
| $g_{\mathrm{PE}}$ | Photoelastic (a.k.a. `PE`, `Str`) | Strain modifies the refractive index in the bulk | **Volume** integral over the beam |

### Inputs

| Argument | Description |
|---|---|
| `ds` | Data struct; must contain `ds.ofem` (optical FEM results), `ds.mfem` (mechanical FEM results), and `ds.P` (parameter struct). Errors out if `ofem` or `mfem` is missing. |
| `model` | Live COMSOL model handle (LiveLink Java object) that already holds **solved** optical (`emw`) and mechanical (`solid`) eigenmode solutions and the datasets `odset` / `mdset`. |
| `oModes` | Vector of optical solution numbers (`solnum`) to evaluate. Indexed over *all* optical solutions, including low-Q modes. |
| `mModes` | Vector of mechanical solution numbers to evaluate. Indexed over *all* mechanical solutions, including unlocalized modes. Typically `ds.mfem.locInd`. |

Called from `RunNanobeamFEM.m:143` as `CalcGOM(ds, model, 1, ds.mfem.locInd)` when `P.calcG` is set,
and from the `OptimizeNanobeam*.m` sweep drivers.

### Output

Returns `ds` with a new/updated substructure `ds.cpl` (see §8) and `model` with the added
result datasets `odset_bnd`, `mdset_bnd`, `jdset_bnd`, `odset_vol`, `mdset_vol`, `jdset_vol`
and the geometry selection `beamBndsAll`.

---

## 2. Notation

| Symbol | Code | Meaning |
|---|---|---|
| $\mathbf{E}$, $\mathbf{D}$ | `emw.Ex…`, `emw.Dx…` | Complex optical eigenmode fields |
| $\mathbf{u}=(u,v,w)$ | `u`, `v`, `w` | Complex mechanical eigenmode displacement |
| $S_{ij}$ | `solid.eXX…` | Mechanical strain tensor |
| $n$ | `n = ofem.n{end}` | Refractive index of the beam material (diamond); air is $n=1$ |
| $\varepsilon_0$ | `epsilon0_const` | Vacuum permittivity |
| $\varepsilon_r^{\mathrm{av}}$ | `emw.epsrAv` | Local relative permittivity (COMSOL average of the diagonal) |
| $\hat{\mathbf{n}}$ | `nx,ny,nz` (optical frame), `nX,nY,nZ` (mechanical/material frame) | Outward boundary normal |
| $\nu_O$ | `wO(oi)` | Optical frequency, **cyclic** — see §7 |
| $\nu_M$ | `wM = mfem.freqs(mi)` | Mechanical frequency [Hz] |
| $x_{\mathrm{zpf}}$ | `cpl.xzpf(mi)` | Zero-point displacement amplitude |
| $\rho$ | `mfem.rho` | Mass density |
| $p_{11}, p_{12}, p_{44}$ | `P.p11, P.p12, P.p44` | Cubic photoelastic (Pockels) coefficients |
| $\theta$ | `P.rxtal` [deg] | In-plane crystal rotation about the $\langle 100\rangle$ surface normal |

Material constants pre-computed at `CalcGOM.m:44-45`:

$$
\Delta n^2 \;\equiv\; n^2 - 1 \quad(\texttt{Dn2}),
\qquad
\Delta\!\left(\tfrac{1}{n^2}\right) \;\equiv\; \frac{1}{n^2} - 1 \quad(\texttt{Dn2\_1})
$$

These are the diamond-minus-air jumps in $\varepsilon_r$ and $\varepsilon_r^{-1}$ respectively:
$\Delta\varepsilon = \varepsilon_0\,\Delta n^2$ and
$\Delta(\varepsilon^{-1}) = \Delta(1/n^2)/\varepsilon_0$.

---

## 3. Symmetry Bookkeeping

Only a $1/2^N$ wedge of the physical structure is meshed, with symmetry planes at $x=0$,
$y=0$, $z=0$ selected by `P.mevenx`, `P.meveny`, `P.mevenz` ($+1$ even, $-1$ odd, $0$ fixed BC)
and their optical counterparts `P.oevenx/oeveny/oevenz`.

### 3.1 Volume/mass multiplier `symFac` (`CalcGOM.m:48`)

$$
\texttt{symFac} \;=\; 2^{\,|m_x| \,+\, |m_y| \,+\, |m_z|\cdot\delta_{\mathrm{rect}}}
$$

where $m_\alpha = $ `P.mevenα` and $\delta_{\mathrm{rect}} = 1$ only if `P.xsect == 'rect'`.
This converts a wedge integral into a full-structure integral. It is used for the optical
energy denominator $L_V$ (§5.1) and, in `SolveNanobeamFEM.m:439`, for the effective mass.

### 3.2 Octant sign sum `sgnCpl` (`CalcGOM.m:52-61`)

The coupling integrands are *linear* in $\mathbf{u}$, so mirrored octants can add or cancel.
The code enumerates all 8 octants and sums the product of mirror signs:

```matlab
XsymVec = [1 P.mevenx];   YsymVec = [1 P.meveny];   ZsymVec = [1 P.mevenz];
[XV,YV,ZV] = meshgrid(XsymVec,YsymVec,ZsymVec);
sgnAll = prod(transpose([XV(:),YV(:),ZV(:)]));
sgnCpl = sum(sgnAll);
```

$$
\texttt{sgnCpl} \;=\; \sum_{o=1}^{8} s_x^{(o)} s_y^{(o)} s_z^{(o)},
\qquad s_\alpha^{(o)} \in \{1,\; m_\alpha\}
$$

Resulting values:

| Mechanical symmetry | `sgnCpl` | Interpretation |
|---|---|---|
| All three even ($+1,+1,+1$) | **8** | All octants add — full-structure coupling |
| Any one odd ($-1$) | **0** | Mirror octants cancel exactly ⇒ $g_0 = 0$ by symmetry |
| One plane a fixed BC ($0$) | **4** (or 2, 1) | Only the meshed half contributes |

`sgnCpl` therefore plays the same role for the (linear-in-$\mathbf{u}$) numerators that
`symFac` plays for the (quadratic-in-field) denominators.

When `sgnCpl == 0` the loop still records the localized-mode entry so the sweep drivers
have a well-defined `cpl.gMax` (`CalcGOM.m:310-320`).

---

## 4. Integration Regions and Joined Datasets

### 4.1 Domains (`CalcGOM.m:67-70`)

- `bdomM = mfem.dia_domind` — beam (diamond) domains in the mechanical model
- `bdomO = ofem.dia_domind` — beam domains in the optical model
- `adomO = ofem.air_domind` — surrounding air cylinder

### 4.2 Boundaries (`CalcGOM.m:77-98`)

An `AdjacentSelection` named `beamBndsAll` extracts every 2-D boundary adjacent to the
3-D `beamSel` selection. Symmetry-plane faces are then **removed** whenever *both* the
mechanical and the optical simulation exploit that plane:

```matlab
if (abs(P.mevenx) && abs(P.oevenx)), bndsM = setdiff(bndsM,P.bndSel.cylXsym); end
if (abs(P.meveny) && abs(P.oeveny)), bndsM = setdiff(bndsM,P.bndSel.cylYsym); end
if (abs(P.mevenz) && abs(P.oevenz)), bndsM = setdiff(bndsM,P.bndSel.cylZsym); end
```

A symmetry plane is not a physical dielectric interface, so it must not contribute to the
moving-boundary integral.

### 4.3 Join datasets (`CalcGOM.m:106-156`)

A single integrand needs *both* the optical and the mechanical solution simultaneously.
COMSOL's `Join` dataset provides this: `data1(...)` evaluates on the first (optical)
solution, `data2(...)` on the second (mechanical) solution.

| Dataset | Type | Selection | Used for |
|---|---|---|---|
| `odset_bnd` / `mdset_bnd` | duplicate of `odset` / `mdset`, dim 2 | `bndsM` | operands of the boundary join |
| `jdset_bnd` | `Join`, method `explicit` | — | moving-boundary surface integral |
| `odset_vol` / `mdset_vol` | duplicate of `odset` / `mdset`, dim 3 | `bdomO` / `bdomM` | operands of the volume join |
| `jdset_vol` | `Join`, method `explicit` | — | photoelastic volume integral |

Inside the mode loops the specific eigenmode pair is selected by
`jdset.set('solnum', oi)` and `jdset.set('solnum2', mi)`.

---

## 5. Normalizations

### 5.1 Optical energy denominator (`CalcGOM.m:161-162`)

```matlab
epsE2Str = 'abs(emw.normE)^2*emw.epsrAv*epsilon0_const';
LV = symFac*mphint2(model,epsE2Str,'volume','dataset','odset', ...
                    'selection',[adomO, bdomO],'solnum',oModes);
```

$$
\boxed{\;
L_V \;=\; \texttt{symFac} \int_{V_{\mathrm{beam}} \cup V_{\mathrm{air}}}
      \varepsilon_0\,\varepsilon_r(\mathbf{r})\,\bigl|\mathbf{E}(\mathbf{r})\bigr|^2 \; dV
\;}
$$

This is (twice) the electric energy of the optical mode and serves as the common
denominator of both coupling terms. Note the PML domains are excluded.

### 5.2 Optical effective volume (`CalcGOM.m:163-164`)

$$
L_V^{\max} \;=\; \max_{V}\;\varepsilon_0 \varepsilon_r |\mathbf{E}|^2,
\qquad
\boxed{\;
V_{\mathrm{eff}} \;=\; \frac{L_V / L_V^{\max}}{\bigl(\lambda / 2n\bigr)^{3}}
\;}
$$

i.e. the mode volume $\int \varepsilon|E|^2 dV / \max(\varepsilon|E|^2)$ expressed in units
of a cubic half-wavelength in the material. Stored as `cpl.Veff`.

### 5.3 Zero-point displacement (`SolveNanobeamFEM.m:437-440`, read at `CalcGOM.m:167`)

$$
m_{\mathrm{eff}} \;=\; \texttt{symFac}\;\rho\;
   \frac{\displaystyle\int_{V_{\mathrm{beam}}} |\mathbf{u}|^2\,dV}
        {\displaystyle\max_{V_{\mathrm{beam}}} |\mathbf{u}|^2}
\qquad\Longrightarrow\qquad
\boxed{\;
x_{\mathrm{zpf}} \;=\; \sqrt{\frac{\hbar}{2\,m_{\mathrm{eff}}\;\bigl(2\pi \nu_M\bigr)}}
\;}
$$

This is the standard $\sqrt{\hbar/2m_{\mathrm{eff}}\omega_M}$ for a mode normalized so that
its maximum displacement is the generalized coordinate.

### 5.4 Displacement field normalization (`CalcGOM.m:168`)

```matlab
maxDisp = mphmax(model,'abs(solid.disp)','volume','dataset','mdset','selection','all','solnum',mModes);
```

$$
u_{\max} \;=\; \max_V \bigl|\mathbf{u}\bigr|
$$

Every coupling expression below carries the factor $x_{\mathrm{zpf}}/u_{\max}$, which rescales
the arbitrarily-normalized COMSOL eigenvector so that its peak displacement equals the
zero-point amplitude — consistent with the $m_{\mathrm{eff}}$ definition in §5.3.

---

## 6. The Two Contributions

Both terms follow the perturbation-theory result of Johnson *et al.* for the frequency shift
of a dielectric cavity mode, evaluated at the zero-point displacement:

$$
g_0 \;=\; \frac{\partial \nu_O}{\partial x}\, x_{\mathrm{zpf}}
$$

### 6.1 Moving-boundary term $g_{\mathrm{MB}}$

**Integrand** (`CalcGOM.m:172-177`):

```matlab
mDispExpr = 'data2(u*nX + v*nY + w*nZ)';
oEtExpr   = '(Δn²)*epsilon0_const*data1( abs(emw.normE)^2 - abs(nx*emw.Ex+ny*emw.Ey+nz*emw.Ez)^2 )';
oDnExpr   = '(Δ(1/n²))/epsilon0_const*data1( abs(nx*emw.Dx)^2 + abs(ny*emw.Dy)^2 + abs(nz*emw.Dz)^2 )';
MB = [mDispExpr,'*(',oEtExpr,'-',oDnExpr,')'];
```

In equation form, with $\mathbf{E}_\parallel$ the field component tangential to the surface
and $D_\perp = \hat{\mathbf{n}}\cdot\mathbf{D}$ the normal displacement field:

$$
\boxed{\;
L_{\mathrm{MB}} \;=\; \oint_{\partial V_{\mathrm{beam}}}
\bigl(\mathbf{u}\cdot\hat{\mathbf{n}}\bigr)
\Bigl[\;
\underbrace{\varepsilon_0\,\Delta n^2 \Bigl(|\mathbf{E}|^2 - |\hat{\mathbf{n}}\cdot\mathbf{E}|^2\Bigr)}_{\Delta\varepsilon\,\left|\mathbf{E}_\parallel\right|^2}
\;-\;
\underbrace{\frac{\Delta(1/n^2)}{\varepsilon_0}\Bigl(|n_x D_x|^2 + |n_y D_y|^2 + |n_z D_z|^2\Bigr)}_{\Delta(\varepsilon^{-1})\,\left|D_\perp\right|^2\ \text{(see note)}}
\;\Bigr]\, dA
\;}
$$

evaluated as `mphint2(model, MB, 'surface', 'dataset','jdset_bnd', 'selection','all')`.

**Coupling rate** (`CalcGOM.m:220`):

$$
\boxed{\;
g_{\mathrm{MB}}(o,m) \;=\; -\,\texttt{sgnCpl}\;\cdot\;
\frac{x_{\mathrm{zpf}}(m)}{u_{\max}(m)}\;\cdot\;\frac{\nu_O(o)}{2}\;\cdot\;
\frac{L_{\mathrm{MB}}}{L_V(o)}
\;}
$$

```matlab
cpl.gMB(oi,mi) = -sgnCpl*cpl.xzpf(mi)*0.5*wO(oi)*(LMB/LV(oIdx)/maxDisp(mIdx));
```

> **Note on $|D_\perp|^2$.** The canonical Johnson form uses
> $|\hat{\mathbf{n}}\cdot\mathbf{D}|^2 = |n_x D_x + n_y D_y + n_z D_z|^2$. The code
> instead sums the squares of the three products,
> $|n_x D_x|^2 + |n_y D_y|^2 + |n_z D_z|^2$. These agree only where the normal is
> axis-aligned (the flat top/bottom and sidewall faces), and differ on curved hole
> sidewalls. See §9.

### 6.2 Photoelastic term $g_{\mathrm{PE}}$

**Step 1 — strain in Voigt (engineering) notation** (`CalcGOM.m:180-186`):

$$
S = \bigl(S_1,\dots,S_6\bigr)^{\!\top} =
\bigl(S_{XX},\; S_{YY},\; S_{ZZ},\; 2S_{YZ},\; 2S_{XZ},\; 2S_{XY}\bigr)^{\!\top}
$$

```matlab
S{1}='data2(solid.eXX)';  S{2}='data2(solid.eYY)';  S{3}='data2(solid.eZZ)';
S{4}='data2(2*solid.eYZ)';S{5}='data2(2*solid.eXZ)';S{6}='data2(2*solid.eXY)';
```

**Step 2 — cubic photoelastic tensor** (`CalcGOM.m:230-235`):

$$
p \;=\;
\begin{pmatrix}
p_{11} & p_{12} & p_{12} & 0 & 0 & 0\\
p_{12} & p_{11} & p_{12} & 0 & 0 & 0\\
p_{12} & p_{12} & p_{11} & 0 & 0 & 0\\
0 & 0 & 0 & p_{44} & 0 & 0\\
0 & 0 & 0 & 0 & p_{44} & 0\\
0 & 0 & 0 & 0 & 0 & p_{44}
\end{pmatrix}
$$

Diamond defaults from `LoadMaterialParams.m`: $p_{11} = -0.25$, $p_{12} = 0.043$,
$p_{44} = -0.172$ (an alternate set $-0.094 / 0.017 / -0.051$ exists for the second material branch).

**Step 3 — crystal rotation** (`CalcGOM.m:236`, implemented in `RotateXtalTensor.m`):

For a counter-clockwise in-plane rotation by $\theta = $ `P.rxtal` about the
$\langle 100\rangle$ surface normal, with $c=\cos\theta$, $s=\sin\theta$:

$$
K \;=\;
\begin{pmatrix}
c^2 & s^2 & 0 & 0 & 0 & 2sc\\
s^2 & c^2 & 0 & 0 & 0 & -2sc\\
0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & c & -s & 0\\
0 & 0 & 0 & s & c & 0\\
-sc & sc & 0 & 0 & 0 & c^2-s^2
\end{pmatrix}
\qquad\Longrightarrow\qquad
\boxed{\;p^R \;=\; K\,p\,K^{\!\top}\;}
$$

$K$ is the stress-type Bond matrix $M_\sigma$; since the strain-type matrix satisfies
$M_\varepsilon^{-1} = M_\sigma^{\!\top}$, the mixed stress/strain index structure of $p$
transforms as $p^R = M_\sigma\, p\, M_\varepsilon^{-1} = K p K^{\!\top}$. The same routine is
used for the elastic stiffness tensor `P.D`.

**Step 4 — contract with strain** (`CalcGOM.m:243-269`):

$$
(pS)_i \;=\; \sum_{j=1}^{6} p^R_{ij}\, S_j , \qquad i = 1\ldots 6
$$

The nested loop builds these as COMSOL *text* expressions, emitting only the non-zero terms
of each row (and the literal `'0'` for an entirely zero row) to keep the expression strings short.

**Step 5 — the $\mathbf{E}\!\cdot\!(pS)\!\cdot\!\mathbf{E}$ integrand** (`CalcGOM.m:189-194`, `272-274`):

With $\varepsilon_0|E_\alpha|^2$ and $\varepsilon_0\,\mathrm{Re}(E_\alpha E_\beta^*)$ as the field
factors, the symmetric contraction splits into a diagonal and a shear part:

$$
\begin{aligned}
\mathcal{I}_{\mathrm{div}} &=\; \varepsilon_0\Bigl[(pS)_1|E_x|^2 + (pS)_2|E_y|^2 + (pS)_3|E_z|^2\Bigr] \\[4pt]
\mathcal{I}_{\mathrm{shear}} &=\; 2\varepsilon_0\Bigl[(pS)_4\,\mathrm{Re}\bigl(E_y E_z^*\bigr)
      + (pS)_5\,\mathrm{Re}\bigl(E_x E_z^*\bigr)
      + (pS)_6\,\mathrm{Re}\bigl(E_x E_y^*\bigr)\Bigr]
\end{aligned}
$$

$$
\boxed{\;
L_{\mathrm{PE}}^{D} = \int_{V_{\mathrm{beam}}} \mathcal{I}_{\mathrm{div}}\,dV ,
\qquad
L_{\mathrm{PE}}^{S} = \int_{V_{\mathrm{beam}}} \mathcal{I}_{\mathrm{shear}}\,dV
\;}
$$

Together these are exactly $\int \varepsilon_0\, \mathbf{E}^*\!\cdot\!(p^R\!:\!S)\cdot\mathbf{E}\, dV$
for the symmetric 2-tensor $p^R\!:\!S$ in stress-Voigt ordering
($1\!\to\!xx,\,2\!\to\!yy,\,3\!\to\!zz,\,4\!\to\!yz,\,5\!\to\!xz,\,6\!\to\!xy$).

**Step 6 — coupling rate** (`CalcGOM.m:278`):

The photoelastic permittivity perturbation is
$\delta\varepsilon_{ij} = -\varepsilon_0 n^4 (p\!:\!S)_{ij}$, hence

$$
\boxed{\;
g_{\mathrm{PE}}(o,m) \;=\; +\,\texttt{sgnCpl}\;\cdot\;
\frac{x_{\mathrm{zpf}}(m)}{u_{\max}(m)}\;\cdot\;\frac{\nu_O(o)}{2}\;\cdot\;
\frac{n^4\Bigl(L_{\mathrm{PE}}^{D} + L_{\mathrm{PE}}^{S}\Bigr)}{L_V(o)}
\;}
$$

```matlab
cpl.gPEc(oi,mi,pci) = sgnCpl*cpl.xzpf(mi)*0.5*wO(oi) ...
                      *(n^4*(LPEDc(oi,mi,pci)+LPESc(oi,mi,pci))/LV(oIdx)/maxDisp(mIdx));
```

Note the sign is **opposite** to $g_{\mathrm{MB}}$, which follows from the extra minus in
$\delta\varepsilon = -\varepsilon_0 n^4 p\!:\!S$.

**Step 7 — per-coefficient decomposition** (`CalcGOM.m:196-200`, `224-283`):

`pcompts` is a $3\times 3$ *diagonal* matrix, so the `pci` loop runs three times with only
one Pockels coefficient non-zero at a time:

```matlab
pcompts = [P.p11, 0,     0;
           0,     P.p12, 0;
           0,     0,     P.p44];
```

| `pci` | Active coefficient | Stored in |
|---|---|---|
| 1 | $p_{11}$ only | `cpl.gPEc(oi,mi,1)` |
| 2 | $p_{12}$ only | `cpl.gPEc(oi,mi,2)` |
| 3 | $p_{44}$ only | `cpl.gPEc(oi,mi,3)` |

Because the whole chain $p \to p^R \to (pS)_i \to \int \mathcal{I}\,dV \to g$ is **linear in
$p$**, summing the three partial results reproduces the full-tensor answer exactly:

$$
g_{\mathrm{PE}} \;=\; \sum_{\text{pci}=1}^{3} g_{\mathrm{PE}}^{(\text{pci})}
$$

So this loop costs 3× the integrals but yields a free, exact breakdown of which Pockels
coefficient dominates the coupling (`CalcGOM.m:281-283`).

### 6.3 Total

$$
\boxed{\;
g_{\mathrm{OM}}(o,m) \;=\; g_{\mathrm{MB}}(o,m) \;+\; g_{\mathrm{PE}}(o,m)
\;}
$$

Printed per mode pair as (`CalcGOM.m:287-290`):

```
wM = 5.83 GHz, g0 = -12.4 + 31.7 = 19.3 kHz
```

---

## 7. Sign and Unit Conventions — read before comparing to literature

1. **`wO` is a cyclic frequency, not an angular frequency.** `CalcGOM.m:160` sets
   `wO = c./lambdaAll`, i.e. $\nu_O = c/\lambda$, **not** $\omega_O = 2\pi c/\lambda$. The
   returned `cpl.gOM` is therefore $g_0/2\pi$ in **Hz**, which is the quantity conventionally
   quoted in kHz. The display line converts with `*1e-3` to kHz.

2. **`wM` is likewise a cyclic frequency** ($\nu_M$, Hz). The $2\pi$ *is* correctly applied in
   the $x_{\mathrm{zpf}}$ expression (`2.*pi.*mfem.freqs`).

3. **Only the real part is reported.** `cpl.gOM(oi,mi)` may be complex (lossy/PML eigenmodes);
   the max-tracking and printing use `real(...)`.

4. **Overall MB sign** depends on COMSOL's boundary normal orientation on interior faces
   (`nx` points from the "down" to the "up" domain). The leading $-$ in `gMB` and the sign of
   `sgnCpl` assume the outward-from-diamond convention. Validate against a known mode
   (e.g. a breathing mode, for which $g_{\mathrm{MB}}$ and $g_{\mathrm{PE}}$ should have
   opposite signs in diamond) before trusting the absolute sign.

5. **`symFac` vs `sgnCpl`** are *not* interchangeable: the former multiplies quadratic-in-field
   denominators, the latter linear-in-displacement numerators.

---

## 8. Output Reference — `ds.cpl`

| Field | Size | Description |
|---|---|---|
| `LMB` | `(oi,mi)` | Raw moving-boundary surface integral $L_{\mathrm{MB}}$ |
| `gMB` | `(oi,mi)` | Moving-boundary coupling [Hz] |
| `LPED` | `(oi,mi)` | $\sum_{\text{pci}} L_{\mathrm{PE}}^{D}$ — diagonal photoelastic integral |
| `LPES` | `(oi,mi)` | $\sum_{\text{pci}} L_{\mathrm{PE}}^{S}$ — shear photoelastic integral |
| `gPEc` | `(oi,mi,3)` | Photoelastic coupling split by $p_{11}$ / $p_{12}$ / $p_{44}$ [Hz] |
| `gPE` | `(oi,mi)` | Total photoelastic coupling [Hz] |
| `gOM` | `(oi,mi)` | $g_{\mathrm{MB}} + g_{\mathrm{PE}}$ [Hz] |
| `Veff` | `1×numel(oModes)` | Optical mode volume in units of $(\lambda/2n)^3$ |
| `xzpf` | `1×nMech` | Copy of `mfem.xzpf` [m] |
| `gMax` | scalar | $\max \bigl|\mathrm{Re}\,g_{\mathrm{OM}}\bigr|$ over all evaluated pairs |
| `gOMmax`, `gMBmax`, `gPEmax` | scalar | Components at the maximum |
| `oSol`, `mSol` | scalar | Solution numbers of the maximizing pair |
| `optWvl` | scalar | Optical wavelength at the maximum [m] |
| `mechFreq` | scalar | Mechanical frequency at the maximum [Hz] |
| `Q` | scalar | Optical Q at the maximum (`ofem.QAll(oi)`) |

`gMax` is initialized to `0` at `CalcGOM.m:37`, so repeated calls on the same `ds` re-run the
search from scratch while other `cpl` fields are preserved.

---

## 9. Assumptions, Limitations, and Observations

**Physics assumptions**

- First-order perturbation theory: valid for $|\mathbf{u}| \ll$ feature size, which zero-point
  motion ($\sim$fm) amply satisfies.
- Isotropic refractive index for the optical problem (`ofem.n = {1, nbeam}`); the photoelastic
  response *is* treated as anisotropic (cubic, rotated by `P.rxtal`).
- Cubic crystal symmetry for $p$; the rotation is restricted to in-plane about
  $\langle 100\rangle$.
- Air/vacuum outside ($n = 1$) is baked into `Dn2` and `Dn2_1`.

**Code-level observations worth checking if numbers look wrong**

1. **`oDnExpr` component-wise square** (`CalcGOM.m:175-176`) — as flagged in §6.1, this
   computes $\sum_\alpha |n_\alpha D_\alpha|^2$ rather than $|\hat{\mathbf{n}}\cdot\mathbf{D}|^2$.
   The two differ on non-axis-aligned faces (curved hole sidewalls, boomerang facets).
   The Johnson-formula form is `abs(nx*emw.Dx+ny*emw.Dy+nz*emw.Dz)^2`.

2. **`maxDisp` selection is `'all'`** (`CalcGOM.m:168`) whereas the `m_eff` normalization in
   `SolveNanobeamFEM.m:438` uses `'selection',bdom` (beam only). If a mechanical PML is
   active (`P.solveMechPML`) and carries larger displacement than the beam, the two
   normalizations become inconsistent.

3. **`symFac` uses only mechanical flags** (`CalcGOM.m:48`) but multiplies the *optical*
   integral $L_V$. This is correct as long as the optical and mechanical simulations use the
   same set of symmetry planes (`|mevenα| ≠ 0 ⟺ |oevenα| ≠ 0`), which is also the condition
   under which symmetry faces are excluded from `bndsM` (§4.2). Mixed configurations would
   break the normalization.

4. **`LPEDc` / `LPESc` are not preallocated** and grow inside the triple loop
   (`CalcGOM.m:276-277`). Harmless but flagged by `checkcode`.

5. **Dead code** at `CalcGOM.m:73-75`: `geomname` / `beam` are extracted from
   `mphmodel(model.geom)` but the subsequent calls use `P.geomname` instead (`beam.runCurrent`
   at line 83 is the one live use).

6. **`cpl` arrays are indexed by absolute solution number** (`oi`, `mi`) while `LV` and
   `maxDisp` are indexed by position within `oModes` / `mModes` (`oIdx`, `mIdx`). Sparse
   `mModes` therefore produces `gOM` arrays with zero-filled gaps — expected, but don't
   `sum`/`mean` over them blindly.

---

## 10. References

1. S. G. Johnson *et al.*, "Perturbation theory for Maxwell's equations with shifting material
   boundaries," *Phys. Rev. E* **65**, 066611 (2002) — the moving-boundary formula of §6.1.
2. J. Chan, "Laser cooling of an optomechanical crystal resonator to its quantum ground state
   of motion," PhD thesis, Caltech (2012) — the $pS$ / photoelastic formulation referenced
   directly in the code comment at `CalcGOM.m:240`.
3. M. Eichenfield, J. Chan, R. M. Camacho, K. J. Vahala, O. Painter, "Optomechanical crystals,"
   *Nature* **462**, 78 (2009).
4. T. C. T. Ting, *Anisotropic Elasticity: Theory and Applications* — source of the Bond
   rotation matrix $K$ used in `RotateXtalTensor.m`.

**Related files:** `SolveNanobeamFEM.m` (produces `ofem`/`mfem`, computes $m_{\mathrm{eff}}$ and
$x_{\mathrm{zpf}}$), `RotateXtalTensor.m` (tensor rotation), `LoadMaterialParams.m`
(Pockels coefficients), `CalcStrCplSiV.m` (the analogous strain–SiV coupling calculation),
`RunNanobeamFEM.m` (caller).
