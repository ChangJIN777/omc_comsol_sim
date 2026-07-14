"""FAST optical surrogate (numpy-only) for the 1D nanobeam photonic gap.

Effective-index 1D transfer-matrix (TMM) model: the perforated nanobeam is
reduced to a 2-layer Bragg stack per period -- a solid section (waveguide
effective index n_hi) and a holed section (lower effective index n_lo). This is
a coarse, dependency-free estimator of the TE-like photonic stop-band and its
center wavelength. It is NOT a substitute for the full vectorial solver
(MPB / COMSOL); use it only as a cheap pre-screen and to exercise the loop.

Returns the first stop-band: normalized gap (in frequency) and center
wavelength, so the optimizer can be debugged end-to-end without MPB/COMSOL.
"""
from __future__ import annotations

import numpy as np

C0 = 299_792_458.0  # m/s


def _effective_indices(g, n_core=2.40):
    """Crude effective indices for the solid and holed z-sections.

    n_hi: waveguide effective index of the solid beam (geometry-weighted).
    n_lo: index at the hole, reduced by the air fraction across the width.
    """
    # solid-section effective index: scale core index by transverse confinement.
    # Larger cross-section -> closer to n_core. Simple monotone heuristic.
    area = g.w * g.t
    conf = area / (area + (0.30e-6) ** 2)        # ->1 for large beams
    n_hi = 1.0 + (n_core - 1.0) * (0.55 + 0.45 * conf)

    # holed section: ellipse removes a fraction of the width at the midplane.
    frac_air = min(1.0, (2.0 * g.hy) / g.w)
    n_lo = np.sqrt(max(1.0, n_hi ** 2 * (1.0 - frac_air) + 1.0 * frac_air))
    return float(n_hi), float(n_lo)


def _bloch_trace(lam, n_hi, L_hi, n_lo, L_lo):
    """0.5*Tr of the period transfer matrix at free-space wavelength lam."""
    phi_hi = 2 * np.pi * n_hi * L_hi / lam
    phi_lo = 2 * np.pi * n_lo * L_lo / lam

    def M(phi, n):
        return np.array([[np.cos(phi), 1j * np.sin(phi) / n],
                         [1j * n * np.sin(phi), np.cos(phi)]])

    Mtot = M(phi_hi, n_hi) @ M(phi_lo, n_lo)
    return 0.5 * np.real(Mtot[0, 0] + Mtot[1, 1])


def optical_gap_surrogate(g, n_core=2.40, n_lam=4000,
                          lam_min=0.9e-6, lam_max=2.6e-6):
    """Compute the first photonic stop-band of the effective 1D Bragg stack.

    Returns dict: normalized_gap, center_wavelength_nm, center_freq_Hz,
    f_lower_Hz, f_upper_Hz, n_hi, n_lo.
    """
    n_hi, n_lo = _effective_indices(g, n_core)
    L_hi = max(g.a - 2 * g.hx, 1e-9)   # solid bridge length along z
    L_lo = 2 * g.hx                    # holed length along z

    lam = np.linspace(lam_min, lam_max, n_lam)
    tr = np.array([_bloch_trace(l, n_hi, L_hi, n_lo, L_lo) for l in lam])
    freq = C0 / lam                    # Hz (descending as lam ascends)

    forbidden = np.abs(tr) > 1.0       # evanescent -> inside a stop-band
    # find contiguous forbidden bands; pick the one with the highest center freq
    # (lowest order / first gap), which is the relevant TE-like mirror gap.
    gaps = []
    i = 0
    n = len(lam)
    while i < n:
        if forbidden[i]:
            j = i
            while j < n and forbidden[j]:
                j += 1
            f_hi = float(freq[i])       # at lam[i] (shorter lam side)
            f_lo = float(freq[j - 1])   # at lam[j-1] (longer lam side)
            fc = 0.5 * (f_hi + f_lo)
            ng = (f_hi - f_lo) / fc
            gaps.append((fc, ng, f_lo, f_hi))
            i = j
        else:
            i += 1

    if not gaps:
        return dict(normalized_gap=0.0, center_wavelength_nm=0.0,
                    center_freq_Hz=0.0, f_lower_Hz=0.0, f_upper_Hz=0.0,
                    n_hi=n_hi, n_lo=n_lo, found=False)

    # first (fundamental) gap = highest center frequency among low-order gaps;
    # restrict to telecom-ish range to avoid spurious high-order bands.
    gaps = [gp for gp in gaps if gp[0] < C0 / 1.0e-6]  # lam_center > 1 um
    fc, ng, f_lo, f_hi = max(gaps, key=lambda x: x[1])  # largest gap
    return dict(normalized_gap=float(ng),
                center_wavelength_nm=float(C0 / fc * 1e9),
                center_freq_Hz=float(fc),
                f_lower_Hz=float(f_lo), f_upper_Hz=float(f_hi),
                n_hi=n_hi, n_lo=n_lo, found=True)


def surrogate_bands(g, n_core=2.40, n_lam=6000, lam_min=0.9e-6, lam_max=2.6e-6):
    """Return (k_frac, bands_THz_list) for plotting the effective 1D band diagram.

    Inverts cos(Ka) = 0.5*Tr(omega) to get the Bloch wavevector for each allowed
    frequency, producing the folded dispersion k(omega) in [0, 0.5] (units pi/a).
    Returns k_frac per allowed point and the corresponding frequency (THz), split
    into contiguous propagating branches.
    """
    n_hi, n_lo = _effective_indices(g, n_core)
    L_hi = max(g.a - 2 * g.hx, 1e-9)
    L_lo = 2 * g.hx
    lam = np.linspace(lam_min, lam_max, n_lam)
    tr = np.array([_bloch_trace(l, n_hi, L_hi, n_lo, L_lo) for l in lam])
    freq_THz = (C0 / lam) / 1e12
    allowed = np.abs(tr) <= 1.0
    k_frac = np.full_like(tr, np.nan)
    k_frac[allowed] = np.arccos(tr[allowed]) / np.pi  # in units of pi/a
    return k_frac, freq_THz, allowed


if __name__ == "__main__":
    from geometry import u_to_geometry
    g = u_to_geometry([0.5, 0.6, 0.5, 0.6])
    print("geom (nm):", {k: round(v * 1e9, 1) for k, v in g.as_dict().items()})
    r = optical_gap_surrogate(g)
    for k, v in r.items():
        print(f"  {k}: {v}")
