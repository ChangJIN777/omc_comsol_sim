#!/usr/bin/env python3
"""Plot a 2D phononic band diagram over the hexagonal BZ path and shade the gap.

Usage:
    python scripts/plot_bands_2d.py                              # configs/plot_bands_2d.yaml
    python scripts/plot_bands_2d.py --config path/to/config.yaml

All settings live in the YAML config (default: configs/plot_bands_2d.yaml).
Two input modes, selected by the `kind` key:

  (A) kind: baseline-csv   -- RUNS TODAY, no COMSOL needed.
      Reads tests/data/baseline_oddz_reference.csv, the 27x10 odd-z band
      structure of the MATLAB reference design extracted from the solved
      pre-parameterization .mph. This is the only real band data the project
      has; use it to check the plotter and to see the target-frequency trap.

  (B) kind: npz            -- from a solved candidate.
      Reads an .npz written by acoustic_comsol_2d.save_bands, which
      objective2d.evaluate_candidate now calls for every successful mechanical
      solve (path recorded as `bands_npz` in the result record). Expected keys:
      `k_norm` [n_k], `freqs_evenz` and/or `freqs_oddz` [n_k, n_bands], and
      `a` (0-d) for the BZ annotation.

Deliberately NOT a copy of python-scripts/scripts/plot_bands.py. Four things
differ, and each one would be wrong if carried over:

  1. The x axis is a path through a 2D Brillouin zone, not a scalar k_frac in
     [0,1]. `k_norm` runs 0->3 over Gamma->M->K->Gamma (see
     acoustic_comsol_2d.bz_path, a port of runBands_2D.m:63-64), so the ticks
     are Gamma/M/K/Gamma at 0/1/2/3 with separators at the segment joins.
  2. There are TWO band families (z-even / z-odd about the slab midplane), not
     one. Both are drawn when present; either alone is fine.
  3. Gap shading follows targets_2d.yaml's `gap_mode`. In "complete" mode the
     truncation ceiling from objective2d._combined_bands is drawn as a
     horizontal line, because any apparent gap above it is an artifact of
     having solved a finite number of bands -- making that visible is the whole
     point of showing the line.
  4. There is no optical mode. src/optical_comsol_2d.py is a documented stub;
     carrying `optical`/`optical-surrogate` branches here would be dead code.

The figure also draws targets_2d.yaml's `target_frequency_GHz` and its
+/-rel_tol acceptance window, and -- when they differ -- BOTH the gap the
scorer would pick (inside the window) and the largest gap anywhere. On the
baseline data those are not the same gap, and seeing that is the fastest way to
understand why the reference design currently scores as a near-miss.
"""
import argparse
import os
import sys

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                          # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from bandgap import largest_gap, gap_near_frequency      # noqa: E402
from objective2d import _combined_bands                  # noqa: E402

_HERE = os.path.dirname(__file__)
_DEFAULT_CONFIG = os.path.join(_HERE, "..", "configs", "plot_bands_2d.yaml")
_DEFAULT_TARGETS = os.path.join(_HERE, "..", "configs", "targets_2d.yaml")
_DEFAULT_CSV = os.path.join(_HERE, "..", "tests", "data",
                            "baseline_oddz_reference.csv")

# High-symmetry points of the hexagonal (triangular-lattice) BZ, at the k_norm
# values runBands_2D.m's piecewise kx/ky expressions switch branches.
BZ_TICKS = [0.0, 1.0, 2.0, 3.0]
BZ_LABELS = [r"$\Gamma$", "M", "K", r"$\Gamma$"]

PARITY_STYLE = {
    "evenz": dict(color="#1f5fa6", label="z-even (Symmetry BC at $z{=}0$)"),
    "oddz": dict(color="#c1440e", label="z-odd (Antisymmetry BC at $z{=}0$)"),
}


# ── input ────────────────────────────────────────────────────────────────────

def load_baseline_csv(path=None):
    """Read tests/data/baseline_oddz_reference.csv -> (k_norm, {'oddz': F}, a).

    Same parsing as tests/test_pipeline_2d.py:load_baseline_oddz, kept separate
    so the plotter has no test-module dependency. `a` is the reference lattice
    constant, which the fixture header documents but does not tabulate.
    """
    path = path or _DEFAULT_CSV
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("k_norm"):
                continue
            if line.strip():
                rows.append([float(x) for x in line.split(",")])
    if not rows:
        raise ValueError(f"no data rows in {path}")
    arr = np.array(rows)
    return arr[:, 0], {"oddz": arr[:, 1:]}, 480e-9


def load_npz(path):
    """Read an acoustic_comsol_2d.save_bands .npz -> (k_norm, {parity: F}, a)."""
    d = np.load(path)
    if "k_norm" not in d:
        raise KeyError(f"{path} has no 'k_norm' array (keys: {list(d.keys())}). "
                       f"Was it written by acoustic_comsol_2d.save_bands?")
    k_norm = np.asarray(d["k_norm"], dtype=float)
    families = {p: np.asarray(d[f"freqs_{p}"], dtype=float)
                for p in ("evenz", "oddz") if f"freqs_{p}" in d}
    if not families:
        raise KeyError(f"{path} has neither 'freqs_evenz' nor 'freqs_oddz' "
                       f"(keys: {list(d.keys())})")
    a = float(d["a"]) if "a" in d else None
    return k_norm, families, a


def close_bz_loop(k_norm, families):
    """Append the Gamma point at k=3 if the data stops at 3 - 1/n.

    COMSOL's sweep is range(0, 1/kpts, 3 - 1/kpts), so it never samples k=3;
    runBands_2D.m:556 closes the path by copying row 0. Doing the same here
    makes the K->Gamma leg reach Gamma instead of stopping one step short. It
    cannot change any gap: duplicating an existing row leaves every
    max_k/min_k untouched.
    """
    if k_norm.size == 0 or float(k_norm[-1]) >= 3.0:
        return k_norm, families
    k2 = np.concatenate([k_norm, [3.0]])
    f2 = {p: np.vstack([F, F[0:1, :]]) for p, F in families.items()}
    return k2, f2


# ── gap selection ────────────────────────────────────────────────────────────

def select_gaps(families, f_target, rel_tol, gap_mode):
    """Return a dict describing what to shade.

    Keys: scored (Gap or None), scored_label, best (Gap or None), best_label,
    ceiling_hz (float or None), n_bands_usable (int or None).

    `scored` is the gap objective2d would actually score -- the one whose
    centre falls inside the +/-rel_tol window around f_target. `best` is the
    largest gap anywhere in the spectrum, drawn only when it is a different
    gap, so the figure shows the disagreement rather than hiding it.
    """
    if gap_mode == "complete":
        missing = {"evenz", "oddz"} - set(families)
        if missing:
            raise SystemExit(
                f"gap_mode 'complete' needs both parities; this data has only "
                f"{sorted(families)} (missing {sorted(missing)}). A complete "
                f"gap is defined as one no band of EITHER parity crosses, so "
                f"it cannot be evaluated from one family. Use gap_mode: "
                f"'symmetry' for single-parity data such as the baseline CSV.")
        bands, ceiling = _combined_bands(families["evenz"], families["oddz"])
        if bands.shape[1] < 2:
            return dict(scored=None, scored_label=None, best=None,
                        best_label=None, ceiling_hz=ceiling, n_bands_usable=0)
        scored = gap_near_frequency(bands, f_target, rel_tol=rel_tol)
        best = largest_gap(bands)
        return dict(scored=scored if scored.found else None,
                    scored_label="complete (both parities)",
                    best=best if best.found else None,
                    best_label="complete, largest anywhere",
                    ceiling_hz=float(ceiling),
                    n_bands_usable=int(bands.shape[1]))

    # symmetry mode: best single family near the target, and the family name.
    scored, scored_p = None, None
    for p, F in sorted(families.items()):
        gp = gap_near_frequency(F, f_target, rel_tol=rel_tol)
        if gp.found and (scored is None or gp.normalized_gap > scored.normalized_gap):
            scored, scored_p = gp, p
    best, best_p = None, None
    for p, F in sorted(families.items()):
        gp = largest_gap(F)
        if gp.found and (best is None or gp.normalized_gap > best.normalized_gap):
            best, best_p = gp, p
    return dict(scored=scored,
                scored_label=None if scored is None else f"{scored_p} symmetry gap",
                best=best,
                best_label=None if best is None else f"{best_p}, largest anywhere",
                ceiling_hz=None, n_bands_usable=None)


def _same_gap(g1, g2):
    if g1 is None or g2 is None:
        return False
    return (np.isclose(g1.f_lower, g2.f_lower)
            and np.isclose(g1.f_upper, g2.f_upper))


# ── rendering ────────────────────────────────────────────────────────────────

def render(k_norm, families, out, *, a=None, f_target_hz=8.0e9, rel_tol=0.5,
           gap_mode="symmetry", title=None, ylim_GHz=None, dpi=150,
           show_target_window=True):
    """Draw the band diagram and write it to `out`. Returns the gap info dict."""
    k_norm, families = close_bz_loop(np.asarray(k_norm, dtype=float), families)
    info = select_gaps(families, f_target_hz, rel_tol, gap_mode)

    fig, ax = plt.subplots(figsize=(8.6, 6.0))
    GHz = 1e-9

    # target acceptance window, behind everything -- this is what decides
    # whether a gap is scored at all (bandgap.gap_near_frequency's f_window).
    if show_target_window:
        lo = f_target_hz * (1 - rel_tol) * GHz
        hi = f_target_hz * (1 + rel_tol) * GHz
        ax.axhspan(lo, hi, color="#9aa5b1", alpha=0.13, zorder=0,
                   label=f"scoring window  {lo:.1f}-{hi:.1f} GHz "
                         f"(rel_tol={rel_tol:g})")
        ax.axhline(f_target_hz * GHz, ls="--", lw=1.2, color="#44515e", zorder=1,
                   label=f"target {f_target_hz*GHz:.2f} GHz")

    # the gap the scorer would use
    if info["scored"] is not None:
        g = info["scored"]
        ax.axhspan(g.f_lower * GHz, g.f_upper * GHz, color="#ffd27f", alpha=0.75,
                   zorder=1,
                   label=f"SCORED {info['scored_label']}: "
                         f"{g.normalized_gap*100:.1f}% @ {g.f_center*GHz:.2f} GHz")

    # the largest gap anywhere, when it is NOT the one being scored
    if info["best"] is not None and not _same_gap(info["best"], info["scored"]):
        g = info["best"]
        ax.axhspan(g.f_lower * GHz, g.f_upper * GHz, facecolor="none",
                   edgecolor="#2e8b57", hatch="///", lw=1.4, zorder=1,
                   label=f"NOT scored, {info['best_label']}: "
                         f"{g.normalized_gap*100:.1f}% @ {g.f_center*GHz:.2f} GHz")

    # bands
    for p in ("evenz", "oddz"):
        if p not in families:
            continue
        F = np.sort(np.asarray(families[p], dtype=float), axis=1) * GHz
        st = PARITY_STYLE[p]
        for b in range(F.shape[1]):
            ax.plot(k_norm, F[:, b], "-", lw=1.5, color=st["color"], zorder=3,
                    label=st["label"] if b == 0 else None)

    # truncation ceiling (complete mode only)
    if info["ceiling_hz"] is not None:
        ax.axhline(info["ceiling_hz"] * GHz, ls=":", lw=1.8, color="#8b0000",
                   zorder=4,
                   label=f"truncation ceiling {info['ceiling_hz']*GHz:.2f} GHz "
                         f"({info['n_bands_usable']} usable bands)")
        ax.text(0.015, info["ceiling_hz"] * GHz, " above: not solved in both "
                "families — any gap here is an artifact",
                transform=ax.get_yaxis_transform(), va="bottom", ha="left",
                fontsize=7.5, color="#8b0000", zorder=5)

    # BZ path axis
    for kk in (1.0, 2.0):
        ax.axvline(kk, color="0.55", lw=0.9, zorder=2)
    ax.set_xticks(BZ_TICKS)
    ax.set_xticklabels(BZ_LABELS)
    ax.set_xlim(0.0, 3.0)
    ax.set_xlabel(r"Bloch wavevector along $\Gamma\rightarrow$M$\rightarrow$"
                  r"K$\rightarrow\Gamma$   ($k$ = 0 … 3)")
    ax.set_ylabel("Frequency [GHz]")

    if ylim_GHz is not None:
        ax.set_ylim(float(ylim_GHz[0]), float(ylim_GHz[1]))
    else:
        fmax = max(np.nanmax(F) for F in families.values()) * GHz
        ax.set_ylim(0.0, 1.04 * fmax)

    if title is None:
        fams = " + ".join(sorted(families))
        title = f"2D boomerang phononic bands ({fams}, gap_mode='{gap_mode}')"
        if a is not None:
            title += f"   a = {a*1e9:.0f} nm"
    ax.set_title(title, fontsize=11)
    ax.grid(alpha=0.25, axis="y")
    # Legend below the axes: the entries are long (they carry the gap numbers)
    # and an in-axes legend hides the upper bands on real data.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2,
              fontsize=7.5, framealpha=0.95, borderaxespad=0.0)
    fig.tight_layout()

    d = os.path.dirname(os.path.abspath(out))
    os.makedirs(d, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return info


# ── CLI ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config "
                         "(default: configs/plot_bands_2d.yaml)")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh) or {}

    with open(cfg.get("targets", _DEFAULT_TARGETS)) as fh:
        targets = yaml.safe_load(fh)
    mech = targets["mechanical"]

    kind = cfg.get("kind", "baseline-csv")
    if kind == "baseline-csv":
        k_norm, families, a = load_baseline_csv(cfg.get("csv"))
    elif kind == "npz":
        npz = cfg.get("npz")
        if not npz:
            ap.error("kind: npz requires an `npz:` path in the config. "
                     "Successful mechanical solves record theirs as "
                     "`bands_npz` in the result record.")
        k_norm, families, a = load_npz(npz)
    else:
        ap.error(f"unknown kind {kind!r}; expected 'baseline-csv' or 'npz'. "
                 f"(There is no optical kind -- optical_comsol_2d is a stub.)")

    only = cfg.get("parities")
    if only:
        keep = {p: families[p] for p in only if p in families}
        if not keep:
            ap.error(f"config `parities` {list(only)} selects nothing; the data "
                     f"has {sorted(families)}")
        families = keep

    f_target = float(cfg.get("target_GHz", mech["target_frequency_GHz"])) * 1e9
    info = render(
        k_norm, families, cfg.get("out", "results/figures/bands_2d.png"),
        a=a,
        f_target_hz=f_target,
        rel_tol=float(cfg.get("rel_tol", 0.5)),
        gap_mode=cfg.get("gap_mode", mech.get("gap_mode", "symmetry")),
        title=cfg.get("title"),
        ylim_GHz=cfg.get("ylim_GHz"),
        dpi=int(cfg.get("dpi", 150)),
        show_target_window=bool(cfg.get("show_target_window", True)),
    )

    out = cfg.get("out", "results/figures/bands_2d.png")
    for key, lbl in (("scored", "scored gap"), ("best", "largest gap")):
        g = info[key]
        if g is None:
            print(f"{lbl:12}: none", file=sys.stderr)
        else:
            print(f"{lbl:12}: {g.normalized_gap*100:5.1f}%  "
                  f"{g.f_lower/1e9:.3f} -> {g.f_upper/1e9:.3f} GHz  "
                  f"(centre {g.f_center/1e9:.3f} GHz)", file=sys.stderr)
    if info["ceiling_hz"] is not None:
        print(f"{'ceiling':12}: {info['ceiling_hz']/1e9:.3f} GHz, "
              f"{info['n_bands_usable']} usable bands", file=sys.stderr)
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
