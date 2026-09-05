#!/usr/bin/env python3
"""Visualize a solved mode on a 2D cut plane (mechanical displacement or optical field).

Usage:
    python scripts/plot_mode.py                              # uses configs/plot_mode.yaml
    python scripts/plot_mode.py --config path/to/plot_mode.yaml

All settings live in the YAML config (default: configs/plot_mode.yaml).
Consumes a generic grid export (.npz) with:
    Y, Z      : 2D meshgrids of in-plane coords [m]   (top cut, x=0 plane)
    UX,UY,UZ  : 2D field components on that grid       (displacement or E)
    freq      : scalar mode frequency [Hz]
    kind      : 'mechanical' | 'optical' (optional, for labels)

The COMSOL driver writes this via acoustic_comsol.export_mode_grid(...); point
the `npz` key at that file. For a quick look at the rendering without COMSOL,
set `demo: true` (a synthetic breathing mode is generated and `npz` is ignored).
"""
import argparse
import os
import sys

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                               "plot_mode.yaml")


def _demo_field(n=120, a=500e-9, w=700e-9):
    z = np.linspace(-a/2, a/2, n)
    y = np.linspace(-w/2, w/2, n)
    Z, Y = np.meshgrid(z, y)
    # synthetic even-even "breathing" mode: width expands/contracts symmetrically
    UY = np.sin(np.pi * Y / w) * np.cos(np.pi * Z / a)
    UZ = 0.3 * np.sin(2 * np.pi * Z / a) * np.cos(np.pi * Y / w)
    UX = np.zeros_like(UY)
    return Y, Z, UX, UY, UZ, 8.0e9, "mechanical"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/plot_mode.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    npz = cfg.get("npz")
    demo = cfg["demo"]
    out = cfg["out"]

    if demo:
        Y, Z, UX, UY, UZ, freq, kind = _demo_field()
    else:
        d = np.load(npz)
        Y, Z = d["Y"], d["Z"]
        UX, UY, UZ = d["UX"], d["UY"], d["UZ"]
        freq = float(d["freq"]) if "freq" in d else 0.0
        kind = str(d["kind"]) if "kind" in d else "mode"

    mag = np.sqrt(np.abs(UX)**2 + np.abs(UY)**2 + np.abs(UZ)**2)
    nm = 1e9
    fig, ax = plt.subplots(figsize=(7, 4))
    pcm = ax.pcolormesh(Z*nm, Y*nm, mag, shading="auto", cmap="magma")
    # downsample arrows for the in-plane (y,z) displacement
    s = max(1, Y.shape[0] // 24)
    ax.quiver(Z[::s, ::s]*nm, Y[::s, ::s]*nm,
              np.real(UZ[::s, ::s]), np.real(UY[::s, ::s]),
              color="cyan", scale=20, width=0.003, alpha=0.8)
    fig.colorbar(pcm, ax=ax, label="|field|")
    funit = f"{freq/1e9:.2f} GHz" if kind == "mechanical" else f"{2.998e8/ (freq if freq else 1)*1e9:.0f} nm"
    ax.set_title(f"{kind} mode  ({funit})  — top cut (x=0)")
    ax.set_xlabel("z [nm] (beam axis)")
    ax.set_ylabel("y [nm] (width)")
    ax.set_aspect("equal")
    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
