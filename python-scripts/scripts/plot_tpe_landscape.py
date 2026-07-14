#!/usr/bin/env python3
"""Visualize the TPE-implied "good region" density over all 5 design
variables, and numerically extract the most probable parameter setting(s).

Why this shape of plot: TPE (Tree-structured Parzen Estimator, what Optuna
uses here) works by splitting past trials into a "good" quantile and a
"bad" quantile, then fitting a Parzen-window / KDE density l(x) over the
good ones (and g(x) over the bad ones) and sampling where l(x)/g(x) is
high. The standard, most information-dense way domain scientists visualize
a multi-dimensional density like this is a CORNER PLOT (as used throughout
Bayesian/MCMC posterior analysis, e.g. the `corner` package): a lower-
triangle grid with 1-D marginal densities on the diagonal and 2-D joint
density contours off-diagonal. That is what this script produces, built
directly on `scipy.stats.gaussian_kde` (no extra dependency) rather than
Optuna's private/version-fragile internal `_ParzenEstimator` class.

"Good" trials = top `--good-frac` (default 25%, matching the common TPE
gamma convention) of successful records by score. The 5-D KDE is fit in
NORMALIZED u-space (all dims on [0,1], comparable scales — matches how TPE
itself operates) and only mapped to physical units for axis labels/output.

Mode-finding: the single best-scoring EVALUATED point is not necessarily
where the fitted density peaks (finite-sample noise, and TPE explores
around good regions rather than sitting exactly on them). This script
numerically finds the density's actual local maxima via multi-start
gradient-free optimization over the fitted 5-D KDE, clusters nearby optima
together (same "island"), and reports each surviving island's parameters —
i.e. a mode of the FITTED distribution, not just "the best row in the
table."

Usage:
    python scripts/plot_tpe_landscape.py FILE [FILE ...]
        [--good-frac 0.25] [--n-starts 200] [--out results/figures/tpe_landscape.png]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.stats import gaussian_kde
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from geometry import VARS, load_bounds, phys_to_u, u_to_geometry

UNITS = {"a": "nm", "w": "nm", "hx": "nm", "hy": "nm", "t": "nm"}
SCALE = 1e9  # meters -> nm for all 5 vars


def load_u_scores(files, bounds):
    seen = {}
    for fn in files:
        with open(fn) as f:
            data = json.load(f)
        for r in data:
            if r.get("status") == "success" and r.get("params"):
                seen[r["id"]] = r
    records = list(seen.values())
    u = np.array([phys_to_u(r["params"], bounds) for r in records])
    score = np.array([r.get("score", -99.0) for r in records])
    return u, score, records


def find_islands(kde, n_dim, n_starts=200, merge_dist=0.25, keep_rel_density=0.4,
                 seed=0):
    """Multi-start search for local maxima of `kde`, clustered into islands.

    merge_dist: two optima within this L2 distance in u-space are treated as
    the same island (keep the higher-density one).
    keep_rel_density: drop islands whose density is below this fraction of
    the best-found island's density — filters numerical noise / minor bumps.
    """
    rng = np.random.default_rng(seed)
    neg_log_dens = lambda x: -kde.logpdf(np.clip(x, 0, 1))[0]

    starts = rng.random((n_starts, n_dim))
    optima = []
    for x0 in starts:
        res = minimize(neg_log_dens, x0, method="L-BFGS-B",
                       bounds=[(0.0, 1.0)] * n_dim)
        if res.success:
            optima.append((res.x, -res.fun))

    if not optima:
        return []

    optima.sort(key=lambda t: -t[1])
    islands = []
    for x, logdens in optima:
        if all(np.linalg.norm(x - ix) > merge_dist for ix, _ in islands):
            islands.append((x, logdens))

    best_logdens = islands[0][1]
    islands = [(x, ld) for x, ld in islands
              if np.exp(ld - best_logdens) >= keep_rel_density]
    return islands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--good-frac", type=float, default=0.25,
                    help="Top fraction of successful trials (by score) "
                         "treated as TPE's 'good' set l(x). Default 0.25.")
    ap.add_argument("--n-starts", type=int, default=200,
                    help="Multi-start restarts for mode-finding.")
    ap.add_argument("--out", default="results/figures/tpe_landscape.png")
    ap.add_argument("--out-json", default=None,
                    help="Where to save found islands as JSON "
                         "(default: alongside --out, same basename .json)")
    args = ap.parse_args()

    bounds = load_bounds()
    u, score, records = load_u_scores(args.files, bounds)
    n = len(u)
    print(f"Loaded {n} unique successful records from {len(args.files)} file(s)")
    if n < 10:
        print("Too few points for a meaningful density estimate (<10). Aborting.")
        sys.exit(1)

    n_good = max(5, int(np.ceil(n * args.good_frac)))
    good_idx = np.argsort(-score)[:n_good]
    u_good = u[good_idx]
    print(f"Good set: top {n_good}/{n} by score "
          f"(score range {score[good_idx].min():.3f} to {score[good_idx].max():.3f})")

    d = len(VARS)
    kde_full = gaussian_kde(u_good.T)          # for mode-finding (joint 5-D)
    kde_1d = [gaussian_kde(u_good[:, i]) for i in range(d)]
    kde_2d = {}
    for i in range(d):
        for j in range(i):
            kde_2d[(i, j)] = gaussian_kde(np.vstack([u_good[:, i], u_good[:, j]]))

    # ── Mode-finding: islands of the fitted 5-D density ────────────────────────
    islands = find_islands(kde_full, d, n_starts=args.n_starts)
    print(f"\nFound {len(islands)} island(s) in the fitted 5-D density:")
    island_records = []
    v = bounds["variables"]
    for rank, (x, logdens) in enumerate(islands):
        g = u_to_geometry(list(x), bounds)
        phys = {k: getattr(g, k) * SCALE for k in VARS}
        rel = np.exp(logdens - islands[0][1])
        print(f"  #{rank+1}  rel.density={rel:.2f}  "
              + "  ".join(f"{k}={phys[k]:.0f}nm" for k in VARS))
        island_records.append({"rank": rank + 1, "u": list(map(float, x)),
                               "params_nm": {k: float(phys[k]) for k in VARS},
                               "relative_density": float(rel)})

    out_json = args.out_json or os.path.splitext(args.out)[0] + "_islands.json"
    os.makedirs(os.path.dirname(out_json) or ".", exist_ok=True)
    with open(out_json, "w") as f:
        json.dump({"good_frac": args.good_frac, "n_good": n_good, "n_total": n,
                  "islands": island_records}, f, indent=2)
    print(f"\nIslands saved -> {out_json}")

    # ── Corner plot ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(d, d, figsize=(2.6 * d, 2.6 * d))
    grid = np.linspace(0, 1, 200)

    def to_phys(i, ugrid):
        lo, hi = v[VARS[i]]["min"] * SCALE, v[VARS[i]]["max"] * SCALE
        return lo + ugrid * (hi - lo)

    island_colors = [plt.cm.tab10(k % 10) for k in range(len(islands))]

    for i in range(d):
        for j in range(d):
            ax = axes[i, j]
            if j > i:
                ax.set_axis_off()
                continue
            if i == j:
                dens = kde_1d[i](grid)
                ax.plot(to_phys(i, grid), dens, color="#2b6cb0", lw=1.6)
                ax.fill_between(to_phys(i, grid), dens, alpha=0.25, color="#2b6cb0")
                ax.scatter(to_phys(i, u[:, i]), np.zeros(n) - 0.02 * dens.max(),
                          c=score, cmap="RdYlGn", s=6, alpha=0.5, clip_on=False)
                for rank, (x, _) in enumerate(islands):
                    ax.axvline(to_phys(i, x[i]), color=island_colors[rank],
                              ls="--", lw=1.8)
                ax.set_yticks([])
                if i == d - 1:
                    ax.set_xlabel(f"{VARS[i]} [{UNITS[VARS[i]]}]")
                if i == 0:
                    ax.set_title(f"marginal density: {VARS[i]}", fontsize=9)
            else:
                # meshgrid: xi varies along COLUMNS, xj varies along ROWS.
                # kde_2d[(i,j)] was fit as vstack([u_good[:,i], u_good[:,j]]),
                # so xi.ravel() (cols) supplies the i-coordinate and xj.ravel()
                # (rows) supplies the j-coordinate -> in dens2d, i varies by
                # COLUMN and j varies by ROW. contourf(X,Y,Z) draws Z[row,col]
                # at (X[row,col],Y[row,col]), so X (the j-axis) must vary by
                # ROW -> use xj, and Y (the i-axis) must vary by COLUMN -> use
                # xi. Swapping these (as an earlier version did) silently
                # produces a warped/mislocated density blob that doesn't
                # overlay the actual scatter points.
                xi, xj = np.meshgrid(grid, grid)
                dens2d = kde_2d[(i, j)](np.vstack([xi.ravel(), xj.ravel()])).reshape(xi.shape)
                ax.contourf(to_phys(j, xj), to_phys(i, xi), dens2d, levels=12,
                           cmap="Blues", alpha=0.85)
                sc = ax.scatter(to_phys(j, u[:, j]), to_phys(i, u[:, i]),
                               c=score, cmap="RdYlGn", s=16, alpha=0.85,
                               edgecolors="k", linewidths=0.3)
                for rank, (x, _) in enumerate(islands):
                    ax.scatter(to_phys(j, x[j]), to_phys(i, x[i]), marker="*",
                              s=260, color=island_colors[rank], edgecolors="k",
                              linewidths=1.0, zorder=5)
                if i == d - 1:
                    ax.set_xlabel(f"{VARS[j]} [{UNITS[VARS[j]]}]")
                if j == 0:
                    ax.set_ylabel(f"{VARS[i]} [{UNITS[VARS[i]]}]")

    # ── Island legend in the unused upper-right corner ─────────────────────────
    if islands:
        legend_ax = fig.add_axes([0.62, 0.62, 0.35, 0.35])
        legend_ax.set_axis_off()
        handles = []
        for r, (c, (x, ld)) in enumerate(zip(island_colors, islands)):
            phys = u_to_geometry(list(x), bounds).as_dict()
            label = (f"#{r+1}  rel.dens={np.exp(ld - islands[0][1]):.2f}  "
                    + "  ".join(f"{k}={phys[k]*SCALE:.0f}" for k in VARS))
            handles.append(plt.Line2D([0], [0], marker="*", color="w",
                                      markerfacecolor=c, markeredgecolor="k",
                                      markersize=16, label=label))
        legend_ax.legend(handles=handles, loc="center left", fontsize=7.5,
                         title="Islands (rank, params in nm)", title_fontsize=8.5,
                         frameon=False)

    fig.suptitle(
        f"TPE 'good region' density (top {args.good_frac*100:.0f}% of {n} trials)  "
        f"—  ★ = found island mode(s), color-matched across panels", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=140)
    plt.close(fig)
    print(f"Corner plot saved -> {args.out}")


if __name__ == "__main__":
    main()
