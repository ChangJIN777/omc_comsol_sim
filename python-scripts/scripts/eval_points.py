#!/usr/bin/env python3
"""Evaluate specific targeted u-vectors and append to opt_results.json.

Usage:
    python scripts/eval_points.py --n-k 9
"""
import argparse
import json
import os
import sys
import warnings
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from geometry import u_to_geometry, check_feasibility, load_bounds
from comsol_client import get_model
from database import save_result

# Import helpers from run_opt_comsol
sys.path.insert(0, os.path.dirname(__file__))
from run_opt_comsol import evaluate, _make_id, _save_json

_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "omc_unitcell_iso.mph")

# Physics-motivated points: large hx (strong optical modulation) + a tuned for 1550nm
# Format: [u_a, u_w, u_hx, u_hy]  (all in [0,1], mapping via bounds.yaml)
# a in [500,900]nm: u_a=(a-500)/400
# w in [450,900]nm: u_w=(w-450)/450
# hx in [80,380]nm: u_hx=(hx-80)/300
# hy in [100,400]nm: u_hy=(hy-100)/300
#
# Target: a~840-870nm for 1550nm optical gap, hx close to feasibility limit,
#         w~700-800nm for breathing mode near 7-9 GHz
TARGETED_POINTS = [
    # a=860, w=720, hx=360, hy=220  bridge=140, sidewall=140, hx_frac=0.84
    [0.900, 0.600, 0.933, 0.400],
    # a=850, w=730, hx=350, hy=230  bridge=150, sidewall=135, hx_frac=0.82
    [0.875, 0.622, 0.900, 0.433],
    # a=840, w=700, hx=340, hy=200  bridge=160, sidewall=150, hx_frac=0.81
    [0.850, 0.556, 0.867, 0.333],
    # a=870, w=750, hx=365, hy=240  bridge=140, sidewall=135, hx_frac=0.84
    [0.925, 0.667, 0.950, 0.467],
    # a=830, w=720, hx=340, hy=240  bridge=150, sidewall=120, hx_frac=0.82
    [0.825, 0.600, 0.867, 0.467],
    # a=850, w=800, hx=350, hy=270  bridge=150, sidewall=130, hx_frac=0.82
    [0.875, 0.778, 0.900, 0.567],
    # a=860, w=760, hx=360, hy=250  bridge=140, sidewall=130, hx_frac=0.84
    [0.900, 0.689, 0.933, 0.500],
    # a=855, w=740, hx=355, hy=230  bridge=145, sidewall=140
    [0.888, 0.644, 0.917, 0.433],
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-k", type=int, default=9)
    ap.add_argument("--n-bands-mech", type=int, default=8)
    ap.add_argument("--n-bands-opt", type=int, default=4)
    ap.add_argument("--study-mech", default="mech sym")
    ap.add_argument("--study-opt", default="opt TE")
    ap.add_argument("--out-json", default="results/opt_results.json")
    args = ap.parse_args()

    bounds = load_bounds()
    model = get_model(_TEMPLATE)
    print(f"Evaluating {len(TARGETED_POINTS)} physics-seeded large-hx points")
    print(f"(These target a≈840-870nm + large holes for 1550nm optical gap)\n")

    # Load existing results to append to
    existing = []
    if os.path.exists(args.out_json):
        with open(args.out_json) as fh:
            existing = json.load(fh)
        print(f"Loaded {len(existing)} prior results from {args.out_json}")

    results = list(existing)
    for i, u in enumerate(TARGETED_POINTS):
        g = u_to_geometry(u, bounds)
        ok, reasons = check_feasibility(g, bounds)
        print(f"\n[{i}] a={g.a*1e9:.0f}  w={g.w*1e9:.0f}  "
              f"hx={g.hx*1e9:.0f}  hy={g.hy*1e9:.0f} nm")
        if not ok:
            print(f"  INFEASIBLE: {reasons}")
            continue

        bridge = (g.a - 2*g.hx)*1e9
        sidewall = (g.w/2 - g.hy)*1e9
        fill = np.pi*g.hx*g.hy/(g.a*g.w)
        print(f"  bridge={bridge:.0f}nm  sidewall={sidewall:.0f}nm  fill={fill:.1%}")

        try:
            rec = evaluate(model, u, g,
                           n_k=args.n_k,
                           n_bands_mech=args.n_bands_mech,
                           n_bands_opt=args.n_bands_opt,
                           study_mech=args.study_mech,
                           study_opt=args.study_opt)
        except Exception as e:
            print(f"  [FAIL] {e}")
            continue

        save_result(rec)
        results.append(rec)
        _save_json(results, args.out_json)

        print(f"  G_o={rec['G_o']*100:.1f}%  G_m={rec['G_m']*100:.1f}%  "
              f"f_o={rec.get('f_o_nm', float('nan')):.0f}nm  "
              f"f_m={rec.get('f_m_GHz', 0.0):.2f}GHz  "
              f"score={rec['score']:+.4f}")

    print("\n=== Done ===")
    successful = [r for r in results if r.get("status") == "success"]
    if successful:
        best = max(successful, key=lambda r: r["score"])
        print(f"Overall best: score={best['score']:+.4f}  "
              f"G_o={best['G_o']*100:.1f}%  G_m={best['G_m']*100:.1f}%  "
              f"f_o={best.get('f_o_nm', float('nan')):.0f}nm  "
              f"f_m={best.get('f_m_GHz', 0.0):.2f}GHz")
        p = best['params']
        print(f"  a={p['a']*1e9:.0f}  w={p['w']*1e9:.0f}  "
              f"hx={p['hx']*1e9:.0f}  hy={p['hy']*1e9:.0f} nm")


if __name__ == "__main__":
    main()
