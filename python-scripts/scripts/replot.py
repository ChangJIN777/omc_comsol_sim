#!/usr/bin/env python3
"""Regenerate opt_progress.png from the saved JSON without running COMSOL.

Useful when the optimization script is still running with an older version of
the plot function — just run this in a separate terminal to get the new plot.

Usage:
    python scripts/replot.py
    python scripts/replot.py --json results/opt_results_t5d.json --fig results/figures/opt_progress.png
"""
import argparse, json, os, sys, time
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--json", default="results/opt_results_t5d.json")
ap.add_argument("--fig",  default="results/figures/opt_progress.png")
args = ap.parse_args()

with open(args.json) as f:
    data = json.load(f)

all_scores = [r["score"] for r in data]
G_o = [r["G_o"] * 100 for r in data]
G_m = [r["G_m"] * 100 for r in data]

finite_scores = [s for s in all_scores if np.isfinite(s) and s > -9]
p5   = max(np.percentile(finite_scores, 5), -3.0) if finite_scores else -2.0
ymax = max(max(finite_scores) * 1.1, 0.05)        if finite_scores else 0.5

running_best, curr = [], -np.inf
for s in all_scores:
    curr = max(curr, s)
    running_best.append(curr)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

ax = axes[0]
iters = np.arange(len(all_scores))
clipped = np.clip(all_scores, p5 - 0.1, ymax + 0.1)
norm = plt.Normalize(vmin=p5, vmax=ymax)
ax.scatter(iters, clipped, c=clipped, cmap="RdYlGn", norm=norm,
           s=22, alpha=0.65, zorder=2, label="Individual eval")
ax.plot(iters, np.clip(running_best, p5 - 0.1, ymax + 0.1),
        lw=2.5, color="#1a7a3c", zorder=3, label="Running best")
ax.axhline(0.40, ls="--", color="#aa0000", lw=1.2, label="Target (0.40)")
ax.axhline(0.0,  ls=":",  color="gray",    lw=0.9, alpha=0.6)
ax.set_ylim(p5, ymax)
ax.set_xlabel("Iteration"); ax.set_ylabel("Score")
ax.set_title("Score history  (low tail clipped, best always visible)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[1]
cb_lo = max(np.percentile(finite_scores, 10), -2.0) if finite_scores else -2.0
cb_hi = max(finite_scores) if finite_scores else 0.5
sc = ax.scatter(G_o, G_m, c=np.clip(all_scores, cb_lo, cb_hi),
                cmap="RdYlGn", vmin=cb_lo, vmax=cb_hi,
                s=50, alpha=0.8, zorder=3, edgecolors="none")
plt.colorbar(sc, ax=ax, label=f"Score  (lower tail clipped at {cb_lo:.2f})")
ax.axvline(20, ls="--", color="#cc3333", lw=1.0, label="20% target")
ax.axhline(20, ls="--", color="#1f5fa6", lw=1.0)
ax.set_xlabel("Optical gap  G_o [%]"); ax.set_ylabel("Mechanical breathing gap  G_m [%]")
ax.set_title("Optical vs mechanical gap  (color = score)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.tight_layout()
os.makedirs(os.path.dirname(args.fig) or ".", exist_ok=True)
fig.savefig(args.fig, dpi=130)
plt.close(fig)
mtime = time.strftime("%H:%M:%S", time.localtime(os.path.getmtime(args.fig)))
print(f"Saved {args.fig}  ({os.path.getsize(args.fig):,} bytes, {mtime})")
print(f"  {len(data)} total results,  best score = {max(all_scores):.4f}")
