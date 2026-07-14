#!/usr/bin/env python3
"""Evaluate targeted points around the 1550 nm zone-edge sweet spot.

Zone-edge optical gap condition for t=220nm diamond nanobeam:
  f_zone_edge = C0 / (2 * n_eff * a),  n_eff ≈ 1.09 at zone edge
  → a ≈ 700 nm for 1550 nm gap center
  → f_light at zone edge = C0/(2a) = 214 THz (1401 nm)
  → max guided gap ≈ 2*(214-193)/193 ≈ 22%, so 20% target is achievable

Mechanical breathing gap also at a≈700nm: previous COMSOL runs found
  G_m≈35% at f_m≈6.9 GHz — excellent!

This script evaluates a grid of a in [660,740]nm × hx in [140,300]nm
at fixed w≈800nm, hy≈200nm, plus some w/hy variants.

Usage:
    python scripts/eval_zone_edge.py --n-k 9
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

sys.path.insert(0, os.path.dirname(__file__))
from run_opt_comsol import evaluate, _make_id, _save_json

_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "omc_unitcell_iso.mph")

# Bounds: a [500,900]nm, w [450,900]nm, hx [80,380]nm, hy [100,400]nm
# u = (x - x_min) / (x_max - x_min)
def to_u(a_nm, w_nm, hx_nm, hy_nm):
    return [
        (a_nm - 500) / 400,
        (w_nm - 450) / 450,
        (hx_nm - 80) / 300,
        (hy_nm - 100) / 300,
    ]

# Grid: a in [660,740]nm for zone-edge near 1530-1580nm
# hx in [140,280]nm for moderate optical/mechanical modulation
# w=800nm (good transverse confinement + mechanical)
# hy=200nm (moderate y-modulation)
TARGETED_POINTS = []

# --- Core grid: a × hx sweep at w=800, hy=200 ---
for a_nm in [660, 680, 700, 710, 720, 740]:
    for hx_nm in [140, 170, 200, 230, 260]:
        hy_nm = 200
        w_nm = 800
        # feasibility check: bridge = a - 2*hx > 60nm, sidewall = w/2 - hy > 60nm
        bridge = a_nm - 2 * hx_nm
        sidewall = w_nm / 2 - hy_nm
        if bridge < 70 or sidewall < 70:
            continue
        TARGETED_POINTS.append(to_u(a_nm, w_nm, hx_nm, hy_nm))

# --- Additional: larger hy for stronger y-modulation ---
for a_nm in [700, 720]:
    for hx_nm in [170, 200, 230]:
        for hy_nm in [240, 280]:
            w_nm = 850
            bridge = a_nm - 2 * hx_nm
            sidewall = w_nm / 2 - hy_nm
            if bridge < 70 or sidewall < 70:
                continue
            TARGETED_POINTS.append(to_u(a_nm, w_nm, hx_nm, hy_nm))

# --- Previous best geometries for re-evaluation with fixed gap detection ---
# a=697, w=875, hx=253, hy=179 → G_m=36.2% at 6.89 GHz
TARGETED_POINTS.append(to_u(697, 875, 253, 179))
# a=705, w=849, hx=258, hy=184 → G_m=35.2% at 6.91 GHz
TARGETED_POINTS.append(to_u(705, 849, 258, 184))
# a=690, w=888, hx=289, hy=215 → G_m=38.0% at 5.94 GHz
TARGETED_POINTS.append(to_u(690, 888, 289, 215))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-k", type=int, default=9)
    ap.add_argument("--n-bands-mech", type=int, default=8)
    ap.add_argument("--n-bands-opt", type=int, default=4)
    ap.add_argument("--study-mech", default="mech sym")
    ap.add_argument("--study-opt", default="opt TE")
    ap.add_argument("--out-json", default="results/opt_results.json")
    ap.add_argument("--max-points", type=int, default=None,
                    help="Limit number of points evaluated (for quick testing)")
    args = ap.parse_args()

    bounds = load_bounds()
    model = get_model(_TEMPLATE)

    pts = TARGETED_POINTS
    if args.max_points:
        pts = pts[:args.max_points]

    print(f"Zone-edge sweep: {len(pts)} candidate points")
    print("Target: a≈680-720nm for 1550nm zone-edge gap + G_m≈30-40% at 6-8 GHz\n")

    existing = []
    if os.path.exists(args.out_json):
        with open(args.out_json) as fh:
            existing = json.load(fh)
        print(f"Loaded {len(existing)} prior results from {args.out_json}\n")
    # Build set of already-evaluated IDs to skip duplicates
    seen_ids = {r["id"] for r in existing if "id" in r}

    results = list(existing)
    n_new = 0
    for i, u in enumerate(pts):
        g = u_to_geometry(u, bounds)
        uid = _make_id(u)
        if uid in seen_ids:
            print(f"[{i:2d}] a={g.a*1e9:.0f} hx={g.hx*1e9:.0f} -- already evaluated, skip")
            continue

        ok, reasons = check_feasibility(g, bounds)
        bridge = (g.a - 2*g.hx)*1e9
        sidewall = (g.w/2 - g.hy)*1e9
        fill = np.pi*g.hx*g.hy/(g.a*g.w)
        f_light = 299792458.0 / (2.0 * g.a)

        print(f"[{i:2d}] a={g.a*1e9:.0f}  w={g.w*1e9:.0f}  hx={g.hx*1e9:.0f}  "
              f"hy={g.hy*1e9:.0f} nm  bridge={bridge:.0f}  sw={sidewall:.0f}  "
              f"fill={fill:.1%}  f_light={f_light*1e-12:.1f}THz")

        if not ok:
            print(f"  INFEASIBLE: {reasons}")
            continue

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
        seen_ids.add(uid)
        _save_json(results, args.out_json)
        n_new += 1

        print(f"  G_o={rec['G_o']*100:.1f}%  G_m={rec['G_m']*100:.1f}%  "
              f"f_o={rec.get('f_o_nm', float('nan')):.0f}nm  "
              f"f_m={rec.get('f_m_GHz', 0.0):.2f}GHz  score={rec['score']:+.4f}")

    print(f"\n=== Done: {n_new} new evaluations ===")
    successful = [r for r in results if r.get("status") == "success"]
    if successful:
        best = max(successful, key=lambda r: r.get("score", -999))
        print(f"Overall best: score={best['score']:+.4f}  "
              f"G_o={best['G_o']*100:.1f}%  G_m={best['G_m']*100:.1f}%  "
              f"f_o={best.get('f_o_nm', float('nan')):.0f}nm  "
              f"f_m={best.get('f_m_GHz', 0.0):.2f}GHz")
        p = best["params"]
        print(f"  a={p['a']*1e9:.0f}  w={p['w']*1e9:.0f}  "
              f"hx={p['hx']*1e9:.0f}  hy={p['hy']*1e9:.0f} nm")


if __name__ == "__main__":
    main()
