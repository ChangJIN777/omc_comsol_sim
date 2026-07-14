#!/usr/bin/env python3
"""Scatter-plot optical and mechanical bandgap vs. (t, w) from EXISTING
optimization records — no new COMSOL evaluations.

Every point comes from a real COMSOL evaluation already stored in one of
the campaign result files. Points are NOT on a regular grid (they're
optimizer-sampled across all 5 design variables), so no interpolation is
applied — each marker is one real design, colored by its own gap value.
a/hx/hy vary freely across points; this plot marginalizes over them, so
scatter (not a smooth heatmap) is the honest representation of what the
data actually supports.

Usage:
    python scripts/plot_gap_vs_tw.py [--files FILE [FILE ...]]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..")
RESULT_FILES = [
    "results/opt_results_t5d_iso_seeded.json",
]
OUT_DIR = os.path.join(ROOT, "results", "figures")


def load_records():
    seen = {}
    for fn in RESULT_FILES:
        path = os.path.join(ROOT, fn)
        if not os.path.exists(path):
            continue
        with open(path) as f:
            data = json.load(f)
        for r in data:
            if r.get("status") == "success" and r.get("params"):
                seen[r["id"]] = r
    return list(seen.values())


def _scatter_gap(ax, t_nm, w_nm, G, fy_mask, label, cmap, target=0.20):
    sc = None
    if (~fy_mask).any():
        sc = ax.scatter(t_nm[~fy_mask], w_nm[~fy_mask], c=G[~fy_mask] * 100,
                         cmap=cmap, vmin=0, vmax=max(20, np.nanmax(G) * 100),
                         s=45, marker="o", edgecolors="none", alpha=0.85,
                         label="no-fy (mech sym trusted)")
    if fy_mask.any():
        sc2 = ax.scatter(t_nm[fy_mask], w_nm[fy_mask], c=G[fy_mask] * 100,
                         cmap=cmap, vmin=0, vmax=max(20, np.nanmax(G) * 100),
                         s=60, marker="^", edgecolors="k", linewidths=0.5,
                         alpha=0.95, label="real fy (breathing-verified)")
        sc = sc if sc is not None else sc2
    hit = G >= target
    if hit.any():
        ax.scatter(t_nm[hit], w_nm[hit], facecolors="none",
                   edgecolors="red", linewidths=1.5, s=140,
                   label=f"≥{target*100:.0f}% target")
    cb = plt.colorbar(sc, ax=ax)
    cb.set_label(f"{label} [%]")
    ax.set_xlabel("t  (thickness) [nm]")
    ax.set_ylabel("w  (width) [nm]")
    ax.set_title(f"{label} vs. (t, w)\n(marginalized over a, hx, hy)")
    ax.legend(fontsize=7, loc="best")


def main():
    records = load_records()
    print(f"Loaded {len(records)} unique successful records")

    t_nm = np.array([r["params"]["t"] * 1e9 for r in records])
    w_nm = np.array([r["params"]["w"] * 1e9 for r in records])
    G_o = np.array([r.get("G_o", 0.0) for r in records])
    G_m = np.array([r.get("G_m", 0.0) for r in records])
    fy_mask = np.array([bool(r.get("fy_computed")) for r in records])

    print(f"  fy real: {fy_mask.sum()}   fy fallback (--no-fy): {(~fy_mask).sum()}")

    os.makedirs(OUT_DIR, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 6))
    _scatter_gap(ax, t_nm, w_nm, G_o, fy_mask, "Optical gap G_o", "viridis",
                target=0.15)
    fig.tight_layout()
    out1 = os.path.join(OUT_DIR, "gap_optical_vs_t_w.png")
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"Saved -> {out1}")

    fig, ax = plt.subplots(figsize=(7, 6))
    _scatter_gap(ax, t_nm, w_nm, G_m, fy_mask, "Mechanical gap G_m", "magma")
    fig.tight_layout()
    out2 = os.path.join(OUT_DIR, "gap_mechanical_vs_t_w.png")
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    print(f"Saved -> {out2}")


if __name__ == "__main__":
    main()
