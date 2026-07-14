#!/usr/bin/env python3
"""Inspect COMSOL band structure + mode shapes for physical sanity check.

Generates results/figures/comsol_inspect.png with four panels:
  [optical TE bands + light cone]  [optical TE mode |normE|² at zone edge]
  [mech bands (sym=red, antisym=blue)]  [mech breathing mode |u| at zone edge]

With symmetry BCs in the COMSOL template (quarter-domain simulations), each
study returns only one mode family — no post-hoc fy classification needed:
  Mech Sym    → breathing / u_y-odd-in-y family  [red]   (default: Study 1)
  Mech Antisym → other / u_y-even-in-y family    [blue]  (default: Study 3)
  Opt TE      → TE-like guided modes               (default: Study 2)
  Opt TM      → TM-like modes (reserved)          (default: Study 4)

For quarter-domain solves (y ≥ 0 only), mode profiles are symmetry-
reconstructed to show the full beam cross-section:
  Mech: u_y ODD in y, u_z EVEN in y → mirror with sign flip
  Opt : |E| EVEN in y (TE) → mirror without sign flip

Bands are scatter dots (no line-joining) to avoid stitching artefacts.

Usage:
    python scripts/inspect_comsol.py                      # full 15-pt sweep
    python scripts/inspect_comsol.py --n-k 1              # zone-edge quick test
    python scripts/inspect_comsol.py --u 0.4 0.5 0.4 0.6 --n-k 15
    python scripts/inspect_comsol.py \\
        --study-mech-sym "Study 1" --study-mech-antisym "Study 3" \\
        --study-opt-te  "Study 2" --study-opt-tm "Study 4"
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.interpolate import griddata

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from geometry import u_to_geometry, Geometry
from comsol_client import get_model
from bandgap import gap_near_frequency

C0 = 299_792_458.0
_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "omc_unitcell_iso.mph")


# ── COMSOL helpers ────────────────────────────────────────────────────────────

def _study_dataset(model, study):
    matches = [d for d in model.datasets() if d.startswith(study + '//')]
    return matches[0] if matches else None


def _find_inner_id(model, dset, target_hz):
    """(inner_id, freq_hz) for eigenmode closest to target_hz."""
    inner_ids = np.asarray(model.inner(dset)[0])
    freqs_hz = np.real(np.asarray(model.evaluate('freq', dataset=dset)))
    best = int(np.argmin(np.abs(freqs_hz - target_hz)))
    return int(inner_ids[best]), float(freqs_hz[best])


def _eval_list(model, dset, exprs, inner_id):
    """Evaluate expressions co-located at the same DOF grid.

    Returns array [n_nodes, n_exprs]. Coordinate columns are always valid;
    field columns are NaN outside the physics domain (solid mechanics only).
    """
    result = model.evaluate(exprs, dataset=dset, inner=[inner_id])
    return np.column_stack([np.real(np.asarray(r)).ravel() for r in result])


def _fy_from_nodes(model, dset, inner_id):
    """Fractional u_y energy for one eigenmode, from solid mesh nodes.

    Returns float in [0,1]: fraction of displacement energy in the y direction.
    Used to identify the breathing (pinch) mode, which has fy close to 1.
    """
    try:
        data = _eval_list(model, dset, ['real(u)', 'real(v)', 'real(w)'], inner_id)
        valid = ~np.isnan(data[:, 1])
        if valid.sum() == 0:
            return np.nan
        u2 = np.sum(data[valid, 0] ** 2)
        v2 = np.sum(data[valid, 1] ** 2)
        w2 = np.sum(data[valid, 2] ** 2)
        tot = u2 + v2 + w2
        return float(v2 / tot) if tot > 0 else np.nan
    except Exception:
        return np.nan


# ── Band sweeps ───────────────────────────────────────────────────────────────

def _solve_mech_sweep(model, g, study, n_k, n_bands):
    """Sweep Γ→X; return (k_frac, freqs [n_k,n_bands], fy [n_k,n_bands]).

    fy[i,j] = fractional u_y energy for each mode at each k-point, computed
    from mesh-node displacements.  Used for fy-coloured band plots.
    """
    k_frac = np.array([1.0]) if n_k == 1 else np.linspace(1e-3, 1.0, n_k)
    freqs  = np.full((n_k, n_bands), np.nan)
    fy_arr = np.full((n_k, n_bands), np.nan)

    for i, s in enumerate(k_frac):
        model.parameter("kF", f"{np.pi * s / g.a}[1/m]")
        model.solve(study)
        dset = _study_dataset(model, study)

        freqs_raw = np.real(np.asarray(model.evaluate('freq', dataset=dset)))
        inner_ids = np.asarray(model.inner(dset)[0])
        sort_idx  = np.argsort(freqs_raw)
        n_av = min(n_bands, len(freqs_raw))
        freqs[i, :n_av] = freqs_raw[sort_idx[:n_av]]
        sorted_ids = inner_ids[sort_idx[:n_av]]

        # fy per mode — cheap: just node retrieval, no re-solve
        for j, iid in enumerate(sorted_ids):
            fy_arr[i, j] = _fy_from_nodes(model, dset, int(iid))

        fy_str = ', '.join(f'{v:.2f}' for v in fy_arr[i, :4])
        print(f"  k={s:.2f}: [{', '.join(f'{f:.2f}' for f in freqs[i,:4]/1e9)}] GHz"
              f"  fy=[{fy_str}]")

    return k_frac, freqs, fy_arr


def _solve_opt_sweep(model, g, study, n_k, n_bands, min_THz=10.0):
    """Sweep Γ→X for optical; return (k_frac, freqs [n_k, n_bands]) in Hz.

    Filters spurious near-zero eigenvalues (< min_THz) that COMSOL EM solvers
    sometimes return before the first physical guided mode.
    """
    k_frac = np.array([1.0]) if n_k == 1 else np.linspace(1e-3, 1.0, n_k)
    freqs = np.full((n_k, n_bands), np.nan)
    for i, s in enumerate(k_frac):
        model.parameter("kF", f"{np.pi * s / g.a}[1/m]")
        model.solve(study)
        dset = _study_dataset(model, study)
        ev = model.evaluate("freq", dataset=dset) if dset else model.evaluate("freq")
        freqs_s = np.sort(np.real(np.asarray(ev)))
        freqs_s = freqs_s[freqs_s > min_THz * 1e12]   # drop spurious DC modes
        n_av = min(n_bands, len(freqs_s))
        freqs[i, :n_av] = freqs_s[:n_av]
        print(f"  k={s:.2f}: [{', '.join(f'{f:.1f}' for f in freqs[i, :4] * 1e-12)}] THz"
              f"  ({n_av} physical modes)")
    return k_frac, freqs


# ── Mode field extraction ─────────────────────────────────────────────────────

def _mech_mode(model, g: Geometry, study, target_GHz=8.0):
    """Zone-edge breathing mode on the y-z cut plane (x ≈ 0).

    Mode selection: iterates all eigenmodes in the study, computes the
    fractional u_y energy (fy) for each, and picks the mode with the highest
    fy (breathing/pinch character).  Within Mech Sym, this distinguishes the
    true breathing mode from acoustic flexural branches that also satisfy the
    Symmetry-at-y=0 BC but are dominated by u_z or u_x.

    Quarter-domain symmetry reconstruction (works for full-domain too):
      u_y ODD  in y → u_y(-y) = -u_y(y)   (pinch motion)
      u_z EVEN in y → u_z(-y) = +u_z(y)   (Floquet component)
    Takes y ≥ 0 nodes (= all nodes for quarter domain) and mirrors to y < 0.

    Uses all solid 3D nodes projected onto the y-z plane — far more data than
    the ~270-node x=0 boundary slice, giving smooth griddata without blurring.
    """
    model.parameter("kF", f"{np.pi / g.a}[1/m]")
    model.solve(study)
    dset = _study_dataset(model, study)

    # ── Pick highest-fy (breathing) mode ────────────────────────────────────
    inner_ids = np.asarray(model.inner(dset)[0])
    freqs_hz  = np.real(np.asarray(model.evaluate('freq', dataset=dset)))
    best_fy, best_iid, best_freq = -1.0, int(inner_ids[0]), float(freqs_hz[0])
    fy_report = []
    for iid, fhz in zip(inner_ids, freqs_hz):
        fy = _fy_from_nodes(model, dset, int(iid))
        fy_report.append(f"{fhz/1e9:.2f} GHz fy={fy:.2f}")
        if not np.isnan(fy) and fy > best_fy:
            best_fy, best_iid, best_freq = fy, int(iid), float(fhz)
    inner_id, freq_hz = best_iid, best_freq
    print(f"  Mech zone-edge mode selection:")
    for r in fy_report:
        print(f"    {r}")
    print(f"  → breathing mode: inner_id={inner_id}  {freq_hz/1e9:.3f} GHz  fy={best_fy:.2f}")

    data = _eval_list(model, dset,
                      ['x', 'y', 'z', 'real(u)', 'real(v)', 'real(w)'],
                      inner_id)
    x, y, z, ud, vd, wd = data.T

    # Use ALL solid nodes (not just an x=0 boundary slice).
    # The breathing mode is approximately uniform in x (thin beam), so projecting
    # all 3D nodes onto the y-z plane gives the correct cross-sectional profile
    # while providing ~6× more data points than the x=0 face alone.
    solid = ~np.isnan(ud)
    ys, zs = y[solid], z[solid]
    vs, ws = vd[solid], wd[solid]
    n_solid = solid.sum()

    # Remove nodes whose (y,z) sits inside the void ellipse
    outside_hole = (zs / g.hx)**2 + (ys / g.hy)**2 > 1.0
    yp, zp, vp, wp = ys[outside_hole], zs[outside_hole], vs[outside_hole], ws[outside_hole]

    # ── Symmetry reconstruction ──────────────────────────────────────────────
    # Take y ≥ 0 half (= all nodes for quarter-domain), mirror to y < 0.
    pos_half = yp >= 0
    yp_h, zp_h = yp[pos_half], zp[pos_half]
    vp_h, wp_h = vp[pos_half], wp[pos_half]

    yp_sym = np.concatenate([yp_h,  -yp_h])
    zp_sym = np.concatenate([zp_h,   zp_h])
    vp_sym = np.concatenate([vp_h,  -vp_h])    # u_y: odd  in y
    wp_sym = np.concatenate([wp_h,   wp_h])    # u_z: even in y
    disp_sym = np.sqrt(vp_sym**2 + wp_sym**2)

    # ── Griddata on fixed full-beam grid ─────────────────────────────────────
    n_grid = 150
    Gy, Gz = np.meshgrid(
        np.linspace(-g.w / 2, g.w / 2, n_grid),
        np.linspace(-g.a / 2, g.a / 2, n_grid),
    )
    interp_v = griddata((yp_sym, zp_sym), vp_sym, (Gy, Gz),
                        method='linear', fill_value=np.nan)
    interp_w = griddata((yp_sym, zp_sym), wp_sym, (Gy, Gz),
                        method='linear', fill_value=np.nan)
    mag = np.sqrt(np.nan_to_num(interp_v)**2 + np.nan_to_num(interp_w)**2).astype(float)

    # Mask void hole (Gz, Gy, g.hx, g.hy all in metres)
    hole = (Gz / g.hx)**2 + (Gy / g.hy)**2 <= 1.0
    mag[hole] = np.nan
    interp_v[hole] = np.nan
    interp_w[hole] = np.nan

    print(f"  Solid nodes: {n_solid}  outside hole: {outside_hole.sum()}  "
          f"after y-mirror: {len(yp_sym)}  hole masked: {hole.sum()}")

    return (Gy * 1e9, Gz * 1e9, mag, interp_v, interp_w,
            yp_sym * 1e9, zp_sym * 1e9, disp_sym, vp_sym, wp_sym, freq_hz / 1e9)


def _opt_mode(model, g: Geometry, study, target_THz=193.4):
    """Zone-edge TE mode on the y-z cut plane (x ≈ 0).

    Quarter-domain symmetry reconstruction:
      |E| EVEN in y for TE (PMC at x=0, PEC at y=0) → |E|(-y) = |E|(y).
    Takes y ≥ 0 nodes and mirrors to y < 0.
    """
    model.parameter("kF", f"{np.pi / g.a}[1/m]")
    model.solve(study)
    dset = _study_dataset(model, study)

    # Filter spurious near-zero modes before picking nearest physical mode
    inner_ids_all = np.asarray(model.inner(dset)[0])
    freqs_all = np.real(np.asarray(model.evaluate('freq', dataset=dset)))
    phys = freqs_all > 10e12   # > 10 THz: discard DC/spurious eigenvalues
    if phys.sum() == 0:
        print("  [warn] no physical optical modes found (all < 10 THz)")
        G = np.zeros((10, 10))
        return G, G, G, 0.0
    inner_ids_phys = inner_ids_all[phys]
    freqs_phys = freqs_all[phys]
    best = int(np.argmin(np.abs(freqs_phys - target_THz * 1e12)))
    inner_id = int(inner_ids_phys[best])
    freq_hz = float(freqs_phys[best])
    freq_THz = freq_hz / 1e12
    wl_nm = C0 / freq_hz * 1e9
    print(f"  Optical zone-edge: inner_id={inner_id}  {freq_THz:.3f} THz ({wl_nm:.0f} nm)"
          f"  [{phys.sum()} physical / {len(freqs_all)} total modes]")

    field_expr = None
    for expr in ['ewfd.normE', 'ewfd.Ey', 'ewfd.Ez', 'emw.normE']:
        try:
            data = _eval_list(model, dset, ['x', 'y', 'z', expr], inner_id)
            if not np.all(np.isnan(data[:, 3])):
                field_expr = expr
                print(f"  EM expression: '{expr}'")
                break
        except Exception:
            continue

    if field_expr is None:
        print("  [warn] no EM expression evaluated; panel will be blank")
        G = np.zeros((10, 10))
        return G, G, G, freq_THz

    x, y, z, field = data.T

    # Cut at x ≈ 0
    valid_x = ~np.isnan(x)
    x_thresh = max(np.percentile(np.abs(x[valid_x]), 8), 10e-9)
    cut = np.abs(x) <= x_thresh
    yp, zp = y[cut], z[cut]
    fp = np.abs(field[cut])

    # Symmetry reconstruction: |E| is even in y for TE modes
    pos_half = yp >= 0                          # all nodes for quarter domain
    yp_h, zp_h, fp_h = yp[pos_half], zp[pos_half], fp[pos_half]
    yp_sym = np.concatenate([yp_h,  -yp_h])
    zp_sym = np.concatenate([zp_h,   zp_h])
    fp_sym = np.concatenate([fp_h,   fp_h])    # |E|: even in y

    # Fixed grid spanning full beam width
    g1 = np.linspace(-g.w / 2, g.w / 2, 120)
    g2 = np.linspace(zp.min(), zp.max(), 120)
    Gy, Gz = np.meshgrid(g1, g2)
    interp_f = griddata((yp_sym, zp_sym), fp_sym, (Gy, Gz),
                        method='linear', fill_value=0.0)
    print(f"  Optical cut pts: {cut.sum()}  after y-mirror: {len(yp_sym)}")

    return Gy * 1e9, Gz * 1e9, interp_f**2, freq_THz


# ── Band scatter plots ────────────────────────────────────────────────────────

def _mech_scatter_fy(ax, k_frac, freqs_sym, fy_sym, freqs_anti, fy_anti,
                     gap_sym, title, target_line=None):
    """Mechanical band scatter colored by fy (fractional u_y energy).

    Color: blue (fy≈0, acoustic) → red (fy≈1, breathing/pinch).
    Sym study → filled circles.  Antisym study → hollow circles (edge only).
    """
    cmap = plt.cm.coolwarm
    norm = Normalize(vmin=0.0, vmax=1.0)

    def _plot(k_vals, freqs, fy_arr, filled):
        for b in range(freqs.shape[1]):
            for i in range(freqs.shape[0]):
                fv = freqs[i, b]
                if np.isnan(fv):
                    continue
                fy = fy_arr[i, b] if (fy_arr is not None and not np.isnan(fy_arr[i, b])) else 0.5
                color = cmap(norm(fy))
                if filled:
                    ax.scatter(k_vals[i], fv * 1e-9, c=[color], s=18,
                               zorder=3, linewidths=0)
                else:
                    ax.scatter(k_vals[i], fv * 1e-9, c='none', s=22,
                               zorder=3, edgecolors=[color], linewidths=1.2)

    _plot(k_frac, freqs_sym, fy_sym, filled=True)
    if freqs_anti is not None:
        _plot(k_frac, freqs_anti, fy_anti, filled=False)

    if gap_sym.found:
        lo, hi = gap_sym.f_lower * 1e-9, gap_sym.f_upper * 1e-9
        ax.axhspan(lo, hi, color='#ffd27f', alpha=0.55, zorder=0)
        ax.text(0.97, (lo + hi) / 2, f"{gap_sym.normalized_gap * 100:.1f}%",
                ha='right', va='center', fontsize=9, color='#8b6000',
                transform=ax.get_yaxis_transform())

    if target_line is not None:
        ax.axhline(target_line * 1e-9, ls='--', color='#555555', lw=1.2)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("fy  (u_y energy fraction)", fontsize=8)

    legend_handles = [
        Line2D([0], [0], marker='o', color='gray', markerfacecolor='gray',
               markersize=7, linestyle='None', label='Sym (filled)'),
        Line2D([0], [0], marker='o', color='gray', markerfacecolor='none',
               markersize=7, linestyle='None', markeredgewidth=1.2, label='Antisym (hollow)'),
    ]
    if target_line is not None:
        legend_handles.append(
            Line2D([0], [0], ls='--', color='#555555', lw=1.2, label='8 GHz target'))
    ax.legend(handles=legend_handles, fontsize=7.5, loc='upper left')

    ax.set_xlim(0, 1)
    ax.set_xlabel(r"$k_z a / \pi$  (Γ→X)", fontsize=9)
    ax.set_ylabel("Frequency [GHz]", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)


def _opt_scatter(ax, k_frac, freqs, g: Geometry, gap, title, target_line=None):
    """Optical band scatter with light-cone overlay.

    Guided modes (below the air light line) → solid blue.
    Leaky / radiation modes (above light cone) → faded grey.
    Y-axis restricted to the guided region.
    """
    f_THz = freqs * 1e-12
    f_light = C0 * k_frac / (2.0 * g.a) * 1e-12   # THz; f_light = C0·k/(2π)
    f_top = f_light[-1] * 1.08

    for b in range(f_THz.shape[1]):
        for i in range(f_THz.shape[0]):
            fv = f_THz[i, b]
            if np.isnan(fv):
                continue
            guided = fv < f_light[i]
            c = '#1f5fa6' if guided else '#999999'
            alpha = 1.0 if guided else 0.35
            ax.scatter(k_frac[i], fv, c=c, s=10, alpha=alpha, zorder=3, linewidths=0)

    ax.fill_between(k_frac, f_light, f_top, color='gray', alpha=0.13, zorder=0)
    ax.plot(k_frac, f_light, 'k-', lw=1.5, label='light line (air)', zorder=5)

    if gap.found:
        lo, hi = gap.f_lower * 1e-12, gap.f_upper * 1e-12
        if lo < f_top:
            ax.axhspan(lo, min(hi, f_top), color='#ffd27f', alpha=0.55, zorder=1)
            ax.text(0.97, (lo + min(hi, f_top)) / 2,
                    f"{gap.normalized_gap * 100:.1f}%",
                    ha='right', va='center', fontsize=9, color='#8b6000',
                    transform=ax.get_yaxis_transform())

    if target_line is not None:
        ax.axhline(target_line * 1e-12, ls='--', color='#cc3333', lw=1.2,
                   label='1550 nm target')

    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f5fa6',
               markersize=7, label='guided (below light cone)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#999999',
               markersize=7, alpha=0.5, label='leaky / radiation'),
        Line2D([0], [0], ls='-', color='k', lw=1.5, label='light line'),
    ]
    if target_line is not None:
        legend_handles.append(
            Line2D([0], [0], ls='--', color='#cc3333', lw=1.2, label='1550 nm'))
    ax.legend(handles=legend_handles, fontsize=7.5, loc='upper left')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, f_top)
    ax.set_xlabel(r"$k_z a / \pi$  (Γ→X)", fontsize=9)
    ax.set_ylabel("Frequency [THz]", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)


# ── Mode shape panel ──────────────────────────────────────────────────────────

def _add_mode_ax(ax, Gy_nm, Gz_nm, field_mag, vy=None, vz=None,
                 title="", cbar_label="",
                 beam_yw=None, beam_za=None,
                 hole_hy=None, hole_hz=None):
    """pcolormesh + optional quiver on the y-z cut plane."""
    pcm = ax.pcolormesh(Gz_nm, Gy_nm, field_mag, shading='auto', cmap='magma')
    plt.colorbar(pcm, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)

    if vy is not None and vz is not None:
        s = max(1, Gy_nm.shape[0] // 16)
        q1 = np.nan_to_num(vy[::s, ::s])
        q2 = np.nan_to_num(vz[::s, ::s])
        scale = max(np.max(np.abs([q1, q2])) * 14, 1e-30)
        ax.quiver(Gz_nm[::s, ::s], Gy_nm[::s, ::s],
                  q2 / scale, q1 / scale,
                  color='cyan', scale=8, width=0.005, alpha=0.7)

    if beam_yw is not None and beam_za is not None:
        rect = plt.Rectangle((-beam_za / 2, -beam_yw / 2), beam_za, beam_yw,
                              fill=False, edgecolor='white', lw=1.2, ls='--')
        ax.add_patch(rect)
    if hole_hy is not None and hole_hz is not None:
        ell = matplotlib.patches.Ellipse(
            (0, 0), 2 * hole_hz, 2 * hole_hy,
            fill=False, edgecolor='lime', lw=1.2)
        ax.add_patch(ell)

    ax.set_xlabel("z [nm]  (beam axis)", fontsize=9)
    ax.set_ylabel("y [nm]  (width)", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.set_aspect('equal', adjustable='datalim')


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Inspect COMSOL OMC unit-cell band structure and mode shapes.")
    ap.add_argument("--u", nargs=4, type=float, default=[0.5, 0.6, 0.5, 0.6],
                    metavar=("ua", "uw", "uhx", "uhy"))
    ap.add_argument("--n-k", type=int, default=15,
                    help="k-points per sweep (1 = zone-edge-only quick test)")
    ap.add_argument("--n-bands-mech", type=int, default=8,
                    help="Bands per mech study (quarter domain → fewer per family)")
    ap.add_argument("--n-bands-opt", type=int, default=4,
                    help="Bands per optical study (one polarisation per study)")
    ap.add_argument("--study-mech-sym", default="mech sym",
                    help="Mech study with Symmetry BCs → breathing family")
    ap.add_argument("--study-mech-antisym", default="mech antisym",
                    help="Mech study with Antisymmetry BCs → other families "
                         "(skipped if study name not found in template)")
    ap.add_argument("--study-opt-te", default="opt TE",
                    help="Optical TE study (PMC at x=0, PEC at y=0)")
    ap.add_argument("--study-opt-tm", default="opt TM",
                    help="Optical TM study (reserved; not plotted yet — "
                         "skipped if not found)")
    ap.add_argument("--out", default="results/figures/comsol_inspect.png")
    args = ap.parse_args()

    g = u_to_geometry(args.u)
    print("Geometry:", {k: f"{v * 1e9:.1f} nm" for k, v in g.as_dict().items()})

    template = os.path.abspath(_TEMPLATE)
    if not os.path.exists(template):
        print(f"[FAIL] Template not found: {template}")
        sys.exit(1)

    model = get_model(template)
    for name, val in g.as_dict().items():
        model.parameter(name, f"{val}[m]")

    avail = set(model.studies())
    print(f"Available studies: {sorted(avail)}")

    # ── 1. Mech Sym sweep (breathing family) ──────────────────────────────────
    print(f"\n=== Mech Sym sweep  ({args.n_k} k-pts, {args.n_bands_mech} bands)"
          f"  study: '{args.study_mech_sym}' ===")
    k_m, freqs_sym, fy_sym = _solve_mech_sweep(
        model, g, args.study_mech_sym, args.n_k, args.n_bands_mech)
    gap_m = gap_near_frequency(freqs_sym, 8e9)
    if gap_m.found:
        print(f"  -> Mech Sym gap: {gap_m.normalized_gap * 100:.1f}%  "
              f"center {gap_m.f_center / 1e9:.2f} GHz")
    else:
        print("  -> No gap near 8 GHz in Sym study")

    # ── 2. Mech Antisym sweep (other families) ────────────────────────────────
    freqs_anti = None
    fy_anti = None
    if args.study_mech_antisym in avail:
        print(f"\n=== Mech Antisym sweep  ({args.n_k} k-pts)"
              f"  study: '{args.study_mech_antisym}' ===")
        _, freqs_anti, fy_anti = _solve_mech_sweep(
            model, g, args.study_mech_antisym, args.n_k, args.n_bands_mech)
    else:
        print(f"\n[skip] Mech Antisym study '{args.study_mech_antisym}' not in template")

    # ── 3. Mech zone-edge breathing mode profile ──────────────────────────────
    print("\n=== Mech zone-edge breathing mode  (Mech Sym study) ===")
    try:
        (Gy_m, Gz_m, mag_m, vy_m, vz_m,
         yp_m, zp_m, disp_raw_m, vp_m, wp_m, freq_m) = _mech_mode(
            model, g, args.study_mech_sym, target_GHz=8.0)
        have_mech_mode = True
    except Exception as e:
        print(f"  [warn] mech mode extraction failed: {e}")
        have_mech_mode = False

    # ── 4. Optical TE sweep ───────────────────────────────────────────────────
    print(f"\n=== Opt TE sweep  ({args.n_k} k-pts, {args.n_bands_opt} bands)"
          f"  study: '{args.study_opt_te}' ===")
    k_o, freqs_o = _solve_opt_sweep(
        model, g, args.study_opt_te, args.n_k, args.n_bands_opt)
    gap_o = gap_near_frequency(freqs_o, C0 / 1550e-9)
    if gap_o.found:
        print(f"  -> TE gap: {gap_o.normalized_gap * 100:.1f}%  "
              f"center {C0 / gap_o.f_center * 1e9:.0f} nm")
    else:
        print("  -> No TE gap near 1550 nm")

    # ── 5. Optical TE zone-edge mode profile ──────────────────────────────────
    print("\n=== Optical TE zone-edge mode ===")
    try:
        Gy_o, Gz_o, normE2_o, freq_o = _opt_mode(
            model, g, args.study_opt_te, target_THz=193.4)
        have_opt_mode = True
    except Exception as e:
        print(f"  [warn] optical mode extraction failed: {e}")
        have_opt_mode = False

    # ── 6. Figure ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"COMSOL inspection — "
        f"a={g.a * 1e9:.0f}  w={g.w * 1e9:.0f}  t={g.t * 1e9:.0f}  "
        f"hx={g.hx * 1e9:.0f}  hy={g.hy * 1e9:.0f} nm",
        fontsize=11)

    _opt_scatter(axes[0, 0], k_o, freqs_o, g, gap_o,
                 f"Optical TE bands  ({args.study_opt_te})",
                 target_line=C0 / 1550e-9)

    _mech_scatter_fy(axes[1, 0], k_m, freqs_sym, fy_sym, freqs_anti, fy_anti, gap_m,
                     f"Mechanical bands  (color=fy, filled={args.study_mech_sym},"
                     f" hollow={args.study_mech_antisym})",
                     target_line=8e9)

    if have_opt_mode:
        wl = C0 / (freq_o * 1e12) * 1e9
        _add_mode_ax(axes[0, 1], Gy_o, Gz_o, normE2_o,
                     title=f"Optical TE |normE|²  k=X  ({wl:.0f} nm)\n"
                           f"y-reconstructed  (|E| even in y)",
                     cbar_label="|normE|²  [arb.]",
                     beam_yw=g.w * 1e9, beam_za=g.a * 1e9,
                     hole_hy=g.hy * 1e9, hole_hz=g.hx * 1e9)
    else:
        axes[0, 1].text(0.5, 0.5, "mode extraction failed",
                        ha='center', va='center',
                        transform=axes[0, 1].transAxes)
        axes[0, 1].set_title("Optical TE mode")

    if have_mech_mode:
        _add_mode_ax(axes[1, 1], Gy_m, Gz_m, mag_m,
                     vy=vy_m, vz=vz_m,
                     title=(f"Mech breathing mode  k=X  ({freq_m:.2f} GHz)\n"
                            f"y-reconstructed  (u_y odd, u_z even in y)"),
                     cbar_label="|u_y,z|  [arb.]",
                     beam_yw=g.w * 1e9, beam_za=g.a * 1e9,
                     hole_hy=g.hy * 1e9, hole_hz=g.hx * 1e9)
    else:
        axes[1, 1].text(0.5, 0.5, "mode extraction failed",
                        ha='center', va='center',
                        transform=axes[1, 1].transAxes)
        axes[1, 1].set_title("Mechanical mode")

    fig.tight_layout()
    out = args.out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\nSaved -> {out}")

    print("\n=== Summary ===")
    if gap_o.found:
        print(f"Optical TE gap : {gap_o.normalized_gap * 100:.1f}%  "
              f"center {C0 / gap_o.f_center * 1e9:.0f} nm")
    else:
        print("Optical TE gap : not found near 1550 nm")
    if gap_m.found:
        print(f"Mech breathing : {gap_m.normalized_gap * 100:.1f}%  "
              f"center {gap_m.f_center / 1e9:.2f} GHz")
    else:
        print("Mech breathing : not found near 8 GHz")


if __name__ == "__main__":
    main()
