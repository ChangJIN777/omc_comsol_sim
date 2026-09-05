#!/usr/bin/env python3
"""Plot a band diagram (Gamma -> X) and shade the bandgap.

Usage:
    python scripts/plot_bands.py                              # uses configs/plot_bands.yaml
    python scripts/plot_bands.py --config path/to/plot_bands.yaml

All settings live in the YAML config (default: configs/plot_bands.yaml).
Two input modes, selected by the `kind` key:

  (A) From a saved band array (any backend dumps freqs [n_k, n_bands]):
      set kind: optical (or mechanical) and npz: results/bands/cand_xxue.npz
      The .npz must contain arrays `k_frac` [n_k] and `freqs_hz` [n_k, n_bands].
      run_one.py --save-bands writes these for MPB / COMSOL runs.

  (B) Demo from the numpy optical surrogate (runs anywhere):
      set kind: optical-surrogate and u: [0.5, 0.6, 0.5, 0.6]

For mechanical bands, run COMSOL on the Mac and save an .npz (see
acoustic_comsol.save_bands), then use mode (A) with kind: mechanical.
"""
import argparse
import os
import sys

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from geometry import u_to_geometry            # noqa: E402
from bandgap import largest_gap, gap_near_frequency  # noqa: E402

C0 = 299_792_458.0

_DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                               "plot_bands.yaml")


def _plot_array(k_frac, freqs, unit_scale, unit_label, gap, title, out):
    fig, ax = plt.subplots(figsize=(6, 5))
    f = np.sort(np.asarray(freqs), axis=1) * unit_scale
    for b in range(f.shape[1]):
        ax.plot(k_frac, f[:, b], "-", color="#1f5fa6", lw=1.6)
    if gap.found:
        lo, hi = gap.f_lower*unit_scale, gap.f_upper*unit_scale
        ax.axhspan(lo, hi, color="#ffd27f", alpha=0.6, zorder=0)
        ax.text(0.02, 0.5*(lo+hi), f"gap {gap.normalized_gap*100:.1f}%",
                va="center", fontsize=10)
    ax.set_xlabel(r"Bloch wavevector  $k_z a/\pi$   ($\Gamma\rightarrow X$, X=1)")
    ax.set_ylabel(unit_label)
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"saved -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/plot_bands.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    npz = cfg.get("npz")
    u = cfg.get("u")
    kind = cfg["kind"]
    target_nm = cfg["target_nm"]
    target_GHz = cfg["target_GHz"]
    out = cfg["out"]

    if kind == "optical-surrogate":
        from optical_surrogate import surrogate_bands
        g = u_to_geometry(u or [0.5, 0.6, 0.5, 0.6])
        k_frac, freq_THz, allowed = surrogate_bands(g)
        # scatter the allowed (propagating) dispersion points
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(k_frac[allowed], freq_THz[allowed], ".", ms=2.5,
                color="#1f5fa6")
        # estimate gap from the surrogate
        from optical_surrogate import optical_gap_surrogate
        r = optical_gap_surrogate(g)
        if r["found"]:
            lo, hi = r["f_lower_Hz"]/1e12, r["f_upper_Hz"]/1e12
            ax.axhspan(lo, hi, color="#ffd27f", alpha=0.6, zorder=0)
            ax.text(0.02, 0.5*(lo+hi),
                    f"gap {r['normalized_gap']*100:.1f}%\n"
                    f"~{r['center_wavelength_nm']:.0f} nm", va="center")
        ax.axhline(C0/(target_nm*1e-9)/1e12, ls="--", color="grey",
                   lw=1, label=f"{target_nm:.0f} nm target")
        ax.set_xlabel(r"Bloch wavevector  $k_z a/\pi$   ($\Gamma\rightarrow X$, X=1)")
        ax.set_ylabel("Frequency [THz]")
        ax.set_title("Optical bands (numpy TMM surrogate — coarse)")
        ax.set_xlim(0, 1)
        if r["found"]:
            ax.set_ylim(max(80, lo - 60), hi + 60)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.savefig(out, dpi=150)
        print(f"surrogate gap {r['normalized_gap']*100:.1f}% near "
              f"{r['center_wavelength_nm']:.0f} nm  ->  {out}")
        return

    # mode (A): load saved bands
    d = np.load(npz)
    k_frac, freqs = d["k_frac"], d["freqs_hz"]
    if kind == "optical":
        gap = gap_near_frequency(freqs, C0/(target_nm*1e-9))
        _plot_array(k_frac, freqs, 1e-12, "Frequency [THz]", gap,
                    "Optical band structure", out)
    else:
        gap = gap_near_frequency(freqs, target_GHz*1e9)
        _plot_array(k_frac, freqs, 1e-9, "Frequency [GHz]", gap,
                    "Mechanical band structure", out)


if __name__ == "__main__":
    main()
