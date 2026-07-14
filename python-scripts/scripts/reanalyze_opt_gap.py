#!/usr/bin/env python3
"""Re-run optical gap detection (zone-edge-only) on stored band data.

This fixes the previous bug where largest_gap() used ALL k-points, causing
leaky radiation modes at k < π/a to fill the zone-edge guided-mode gap.

Usage:
    python scripts/reanalyze_opt_gap.py [--json results/opt_results.json]
"""
import argparse
import json
import os
import sys
import warnings
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from bandgap import largest_gap

C0 = 299_792_458.0
OPT_F_TARGET = C0 / 1550e-9  # 193.4 THz


def zone_edge_opt_gap(freqs_o_all, a):
    """Re-detect optical gap using zone-edge k-point only + light-line filter."""
    freqs_o_all = np.asarray(freqs_o_all, dtype=float)
    freqs_edge = freqs_o_all[-1:, :]        # last k-point = zone edge
    f_light = C0 / (2.0 * a)               # light line at k=π/a
    freqs_edge = np.where(freqs_edge < f_light, freqs_edge, np.nan)
    window = (OPT_F_TARGET * 0.5, OPT_F_TARGET * 1.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        gp = largest_gap(freqs_edge, f_window=window, nan_safe=True)
    return gp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="results/opt_results.json")
    args = ap.parse_args()

    with open(args.json) as fh:
        results = json.load(fh)

    print(f"Loaded {len(results)} results")
    print(f"\n{'a(nm)':>6} {'w(nm)':>6} {'hx(nm)':>6} {'hy(nm)':>6} "
          f"{'old G_o':>8} {'new G_o':>8} {'f_o(nm)':>8} "
          f"{'G_m':>8} {'f_m(GHz)':>9} {'score_new':>10}")
    print("-"*90)

    updated = []
    for r in results:
        if r.get("status") != "success":
            updated.append(r)
            continue
        if "_freqs_o" not in r:
            updated.append(r)
            continue

        a = r["params"]["a"]
        freqs_o_all = np.array(r["_freqs_o"])
        gp = zone_edge_opt_gap(freqs_o_all, a)

        G_o_new = gp.normalized_gap if gp.found else 0.0
        f_o_new = gp.f_center if gp.found else OPT_F_TARGET
        f_o_nm  = C0 / f_o_new * 1e9 if f_o_new > 0 else float("nan")

        G_m     = r.get("G_m", 0.0)
        f_m     = r.get("f_m_GHz", 0.0)

        # Updated score
        G_MIN = 0.20; LAMBDA_OPT = 20.0; LAMBDA_MECH = 20.0
        LAMBDA_OPT_F = 2.0; LAMBDA_MECH_F = 15.0
        MECH_F_MIN = 5e9; MECH_F_MAX = 10e9
        pen_g  = (LAMBDA_OPT  * max(0.0, G_MIN - G_o_new)**2
                  + LAMBDA_MECH * max(0.0, G_MIN - G_m)**2)
        pen_of = LAMBDA_OPT_F * ((f_o_new - OPT_F_TARGET) / OPT_F_TARGET)**2
        f_m_c  = f_m * 1e9
        pen_mf = 0.0
        if f_m_c < MECH_F_MIN:
            pen_mf = LAMBDA_MECH_F * ((MECH_F_MIN - f_m_c) / MECH_F_MIN)**2
        elif f_m_c > MECH_F_MAX:
            pen_mf = LAMBDA_MECH_F * ((f_m_c - MECH_F_MAX) / MECH_F_MAX)**2
        score_new = G_o_new + G_m - pen_g - pen_of - pen_mf

        p = r["params"]
        print(f"{p['a']*1e9:>6.0f} {p['w']*1e9:>6.0f} {p['hx']*1e9:>6.0f} {p['hy']*1e9:>6.0f} "
              f"{r.get('G_o',0)*100:>7.1f}% {G_o_new*100:>7.1f}% {f_o_nm:>8.0f} "
              f"{G_m*100:>7.1f}% {f_m:>8.2f}  {score_new:>+10.4f}")

        r_new = dict(r)
        r_new["G_o"] = G_o_new
        r_new["f_o_THz"] = f_o_new * 1e-12
        r_new["f_o_nm"] = f_o_nm
        r_new["optical_gap"] = G_o_new
        r_new["optical_center_frequency"] = float(f_o_new)
        r_new["score"] = score_new
        updated.append(r_new)

    # Print top-10 by new score (success only)
    successful = [r for r in updated if r.get("status") == "success" and "_freqs_o" in r]
    if successful:
        top = sorted(successful, key=lambda r: r["score"], reverse=True)[:10]
        print("\n=== TOP 10 by updated score ===")
        for r in top:
            p = r["params"]
            print(f"  a={p['a']*1e9:.0f}  w={p['w']*1e9:.0f}  "
                  f"hx={p['hx']*1e9:.0f}  hy={p['hy']*1e9:.0f} nm  "
                  f"G_o={r['G_o']*100:.1f}%  G_m={r.get('G_m',0)*100:.1f}%  "
                  f"f_o={r.get('f_o_nm',float('nan')):.0f}nm  "
                  f"f_m={r.get('f_m_GHz',0):.2f}GHz  score={r['score']:+.4f}")

    # Also print zone-edge band frequencies for top-5 for debugging
    print("\n=== Zone-edge band frequencies (top 5) ===")
    for r in top[:5]:
        a = r["params"]["a"]
        freqs_o_all = np.array(r["_freqs_o"])
        edge = freqs_o_all[-1, :]
        f_light = C0 / (2.0 * a)
        print(f"  a={a*1e9:.0f}nm  f_light={f_light*1e-12:.1f}THz "
              f"({C0/f_light*1e9:.0f}nm):")
        for j, f in enumerate(edge):
            if np.isfinite(f) and f > 0:
                guided = "guided" if f < f_light else "LEAKY"
                print(f"    band{j}: {f*1e-12:.1f}THz = {C0/f*1e9:.0f}nm  [{guided}]")

    # Save updated JSON
    out_path = args.json.replace(".json", "_reanalyzed.json")
    with open(out_path, "w") as fh:
        json.dump(updated, fh, indent=2)
    print(f"\nSaved updated results to {out_path}")


if __name__ == "__main__":
    main()
