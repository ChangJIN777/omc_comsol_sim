#!/usr/bin/env python3
"""Evaluate large-hx geometries at a≈840-880nm for strong optical modulation.

Physics target:
- a≈840-880nm: zone-edge gap near 1500-1600nm (closest to 1550nm target)
- hx/a≈35-47%: needed for G_o≥20% in diamond (n=2.4)
  - Diamond scales ≈80% of Si: need hx/a=45-50% vs Si's 40-47%
  - Bridge constraint: hx_max = (a-60nm)/2 → hx/a_max=47% at a=860nm
- G_m target: breathing gap in 5-10 GHz range

Usage:
    python scripts/eval_large_hx.py --n-k 9
"""
import argparse, json, os, sys, warnings
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from geometry import u_to_geometry, check_feasibility, load_bounds
from comsol_client import get_model
from database import save_result

sys.path.insert(0, os.path.dirname(__file__))
from run_opt_comsol import evaluate, _make_id, _save_json

_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol", "omc_unitcell_iso.mph")

def to_u(a_nm, w_nm, hx_nm, hy_nm):
    return [(a_nm-500)/400, (w_nm-450)/450, (hx_nm-80)/300, (hy_nm-100)/300]

# Targeted: large hx for strong optical gap, a~840-880nm for 1550nm target
# To get G_o=20% in diamond: need hx/a≈45-50% (vs Si at 40-47%)
# Bridge constraint: bridge = a-2*hx ≥ 100nm (soft: 60nm absolute minimum)
# Sidewall constraint: w/2-hy ≥ 100nm (soft: 60nm absolute minimum)
TARGETED_POINTS = [
    # a=840nm: zone-edge gap ~1490-1530nm, hx/a=30-45%
    to_u(840, 850, 252, 250),  # hx/a=30%, bridge=336, sw=175
    to_u(840, 850, 294, 250),  # hx/a=35%, bridge=252, sw=175
    to_u(840, 850, 336, 250),  # hx/a=40%, bridge=168, sw=175
    to_u(840, 850, 360, 300),  # hx/a=43%, bridge=120, sw=125
    to_u(840, 850, 294, 320),  # hx/a=35%, large hy
    to_u(840, 850, 336, 320),  # hx/a=40%, large hy

    # a=860nm: zone-edge gap ~1500-1550nm (sweet spot), hx/a=30-46%
    to_u(860, 850, 258, 250),  # hx/a=30%, bridge=344, sw=175
    to_u(860, 850, 301, 250),  # hx/a=35%, bridge=258, sw=175
    to_u(860, 850, 344, 250),  # hx/a=40%, bridge=172, sw=175
    to_u(860, 850, 380, 300),  # hx/a=44%, bridge=100, sw=125
    to_u(860, 850, 301, 320),  # hx/a=35%, large hy
    to_u(860, 850, 344, 320),  # hx/a=40%, large hy
    to_u(860, 850, 380, 320),  # hx/a=44%, large hy
    to_u(860, 900, 344, 350),  # hx/a=40%, wide beam, large hy
    to_u(860, 900, 380, 350),  # hx/a=44%, wide beam, large hy

    # a=880nm: zone-edge gap ~1550-1580nm, hx/a=30-46%
    to_u(880, 850, 264, 250),  # hx/a=30%, bridge=352
    to_u(880, 850, 308, 250),  # hx/a=35%, bridge=264
    to_u(880, 850, 352, 250),  # hx/a=40%, bridge=176
    to_u(880, 850, 380, 300),  # hx/a=43%, bridge=120
    to_u(880, 850, 308, 320),  # hx/a=35%, large hy
    to_u(880, 850, 352, 320),  # hx/a=40%, large hy
    to_u(880, 900, 352, 350),  # hx/a=40%, wide beam
    to_u(880, 900, 380, 350),  # hx/a=43%, wide beam

    # a=900nm: zone-edge gap ~1580-1620nm, hx/a=30-46%
    to_u(900, 850, 270, 250),  # hx/a=30%
    to_u(900, 850, 315, 250),  # hx/a=35%
    to_u(900, 850, 360, 300),  # hx/a=40%
    to_u(900, 850, 380, 320),  # hx/a=42%
    to_u(900, 900, 360, 350),  # hx/a=40%, wide beam
    to_u(900, 900, 380, 350),  # hx/a=42%, wide beam
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-k", type=int, default=7)
    ap.add_argument("--n-bands-mech", type=int, default=8)
    ap.add_argument("--n-bands-opt", type=int, default=4)
    ap.add_argument("--study-mech", default="mech sym")
    ap.add_argument("--study-opt", default="opt TE")
    ap.add_argument("--no-fy", action="store_true")
    ap.add_argument("--out-json", default="results/opt_results.json")
    ap.add_argument("--max-points", type=int, default=None)
    args = ap.parse_args()
    compute_fy = not args.no_fy

    bounds = load_bounds()
    model = get_model(_TEMPLATE)

    pts = TARGETED_POINTS[:args.max_points] if args.max_points else TARGETED_POINTS
    print(f"Large-hx sweep: {len(pts)} points (a=840-900nm, hx/a=30-35%)")
    print("Target: optical gap near 1550nm + mechanical G_m>15% at 5-10 GHz\n")

    existing = []
    if os.path.exists(args.out_json):
        with open(args.out_json) as fh:
            existing = json.load(fh)
        print(f"Loaded {len(existing)} prior results\n")
    seen_ids = {r["id"] for r in existing if "id" in r}

    results = list(existing)
    n_new = 0
    for i, u in enumerate(pts):
        g = u_to_geometry(u, bounds)
        uid = _make_id(u)
        if uid in seen_ids:
            print(f"[{i:2d}] a={g.a*1e9:.0f} hx={g.hx*1e9:.0f} -- skip (already done)")
            continue
        ok, reasons = check_feasibility(g, bounds)
        bridge = (g.a-2*g.hx)*1e9
        sidewall = (g.w/2-g.hy)*1e9
        fill = np.pi*g.hx*g.hy/(g.a*g.w)
        f_light = 299792458.0/(2*g.a)*1e-12
        print(f"[{i:2d}] a={g.a*1e9:.0f} w={g.w*1e9:.0f} hx={g.hx*1e9:.0f} hy={g.hy*1e9:.0f} "
              f"bridge={bridge:.0f} fill={fill:.1%} f_light={f_light:.1f}THz")
        if not ok:
            print(f"  INFEASIBLE: {reasons}"); continue
        try:
            rec = evaluate(model, u, g, n_k=args.n_k, n_bands_mech=args.n_bands_mech,
                           n_bands_opt=args.n_bands_opt, study_mech=args.study_mech,
                           study_opt=args.study_opt, compute_fy=compute_fy)
        except Exception as e:
            print(f"  [FAIL] {e}"); continue
        save_result(rec); results.append(rec); seen_ids.add(uid)
        _save_json(results, args.out_json); n_new += 1
        print(f"  G_o={rec['G_o']*100:.1f}%  G_m={rec['G_m']*100:.1f}%  "
              f"f_o={rec.get('f_o_nm',float('nan')):.0f}nm  "
              f"f_m={rec.get('f_m_GHz',0):.2f}GHz  score={rec['score']:+.4f}")

    print(f"\n=== Done: {n_new} new evaluations ===")
    if results:
        top = sorted([r for r in results if r.get('status')=='success'],
                     key=lambda r: r.get('score',-999), reverse=True)[:5]
        print("Top-5:")
        for r in top:
            p=r['params']
            print(f"  a={p['a']*1e9:.0f} w={p['w']*1e9:.0f} hx={p['hx']*1e9:.0f} hy={p['hy']*1e9:.0f} "
                  f"G_o={r['G_o']*100:.1f}% G_m={r.get('G_m',0)*100:.1f}% "
                  f"f_o={r.get('f_o_nm',float('nan')):.0f}nm f_m={r.get('f_m_GHz',0):.2f}GHz "
                  f"score={r['score']:+.4f}")


if __name__ == "__main__":
    main()
