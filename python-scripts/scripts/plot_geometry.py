#!/usr/bin/env python3
"""Visualize the simulated unit cell: top view (x-y... actually z-y) + cross-section.

Usage:
    python scripts/plot_geometry.py                              # uses configs/plot_geometry.yaml
    python scripts/plot_geometry.py --config path/to/plot_geometry.yaml

All settings live in the YAML config (default: configs/plot_geometry.yaml).
Draws (a) a top view of the beam along the periodic z-axis with the elliptical
hole(s) and the feasibility margins (bridge, sidewall), and (b) the rectangular
cross-section (y-x) with the hole footprint. Pure matplotlib; runs anywhere.
"""
import argparse
import os
import sys

import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from geometry import u_to_geometry, check_feasibility  # noqa: E402

_DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                               "plot_geometry.yaml")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/plot_geometry.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    u = [float(x) for x in cfg["u"]]
    periods = cfg["periods"]
    out = cfg["out"]

    g = u_to_geometry(u)
    ok, reasons = check_feasibility(g)
    nm = 1e9
    a, w, t, hx, hy = g.a*nm, g.w*nm, g.t*nm, g.hx*nm, g.hy*nm

    fig, (axt, axc) = plt.subplots(1, 2, figsize=(11, 4.2),
                                   gridspec_kw={"width_ratios": [3, 1]})

    # --- top view: beam along z, width y ---
    N = periods
    axt.add_patch(Rectangle((-N*a/2, -w/2), N*a, w, facecolor="#cfe8ff",
                            edgecolor="#1f5fa6", lw=1.5, zorder=1))
    for i in range(N):
        zc = -N*a/2 + (i + 0.5)*a
        axt.add_patch(Ellipse((zc, 0), 2*hx, 2*hy, facecolor="white",
                              edgecolor="#1f5fa6", lw=1.2, zorder=2))
    # annotate one period
    axt.annotate("", xy=(-N*a/2, -w/2-90), xytext=(-N*a/2+a, -w/2-90),
                 arrowprops=dict(arrowstyle="<->", color="k"))
    axt.text(-N*a/2+a/2, -w/2-150, f"a = {a:.0f} nm", ha="center")
    axt.annotate("", xy=(N*a/2+60, -w/2), xytext=(N*a/2+60, w/2),
                 arrowprops=dict(arrowstyle="<->", color="k"))
    axt.text(N*a/2+90, 0, f"w = {w:.0f} nm", rotation=90, va="center")
    bridge = a - 2*hx
    sidewall = w/2 - hy
    axt.set_title(f"Top view  (bridge={bridge:.0f} nm, sidewall={sidewall:.0f} nm)")
    axt.set_xlabel("z  [nm]  (periodic / beam axis)")
    axt.set_ylabel("y  [nm]  (beam width)")
    axt.set_aspect("equal")
    axt.set_xlim(-N*a/2-150, N*a/2+260)
    axt.set_ylim(-w/2-260, w/2+120)

    # --- cross-section: width y, thickness x ---
    axc.add_patch(Rectangle((-w/2, -t/2), w, t, facecolor="#cfe8ff",
                            edgecolor="#1f5fa6", lw=1.5))
    # hole is etched fully through the thickness -> through-slot of width 2*hy
    axc.add_patch(Rectangle((-hy, -t/2), 2*hy, t, facecolor="white",
                            edgecolor="#1f5fa6", lw=1.2))
    axc.set_title(f"Cross-section\nt = {t:.0f} nm")
    axc.set_xlabel("y  [nm]")
    axc.set_ylabel("x  [nm] (thickness)")
    axc.set_aspect("equal")
    axc.set_xlim(-w/2-80, w/2+80)
    axc.set_ylim(-t*2, t*2)

    status = "FEASIBLE" if ok else "INFEASIBLE: " + "; ".join(reasons)
    color = "#1a7f37" if ok else "#b00020"
    fig.suptitle(f"Diamond OMC unit cell  |  {status}", color=color, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"geometry (nm): a={a:.0f} w={w:.0f} t={t:.0f} hx={hx:.0f} hy={hy:.0f}")
    print(f"feasible={ok}  saved -> {out}")


if __name__ == "__main__":
    main()
