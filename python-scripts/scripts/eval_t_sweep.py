#!/usr/bin/env python3
"""Sweep beam thickness t at a fixed near-optimal geometry to map G_o vs t.

Fixes a=900nm, w=850nm, hx=315nm, hy=250nm (best optical point from campaign 1)
and varies t from 220nm to 430nm in ~20nm steps.  This shows how strongly G_o
depends on t and guides the bounds for the 5-D campaign.

Usage:
    python -u scripts/eval_t_sweep.py [--no-fy] [--out results/t_sweep.json]
"""
import argparse
import hashlib
import json
import os
import sys
import warnings
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from geometry import Geometry, load_bounds
from comsol_client import get_model

C0 = 299_792_458.0
_TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "omc_unitcell_iso.mph")

# ── Reuse sweep helpers from run_opt_comsol ────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
from run_opt_comsol import (
    _opt_sweep, _mech_sweep, optical_gap_tracked, breathing_gap,
    _study_dataset, OPT_F_TARGET,
)


def _make_id_phys(a_nm, w_nm, hx_nm, hy_nm, t_nm):
    key = json.dumps([a_nm, w_nm, hx_nm, hy_nm, t_nm, "t_sweep"])
    return "ts_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fy",  action="store_true")
    ap.add_argument("--n-k",    type=int, default=7)
    ap.add_argument("--n-bands-mech", type=int, default=8)
    ap.add_argument("--n-bands-opt",  type=int, default=4)
    ap.add_argument("--study-mech",   default="mech sym")
    ap.add_argument("--study-opt",    default="opt TE")
    ap.add_argument("--out",    default="results/t_sweep.json")
    args = ap.parse_args()
    compute_fy = not args.no_fy

    # Fixed geometry (best G_o from campaign 1 at t=220nm)
    A_NM  = 900.0
    W_NM  = 850.0
    HX_NM = 315.0
    HY_NM = 250.0

    # t sweep: 220 → 430nm in ~25nm steps, plus 240nm for fine detail near bottom
    t_sweep_nm = [220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420]

    # Load existing results to skip already-evaluated points
    seen = {}
    if os.path.exists(args.out):
        with open(args.out) as fh:
            existing = json.load(fh)
        for r in existing:
            seen[r["id"]] = r
        print(f"Loaded {len(seen)} prior t-sweep results")

    model = get_model(_TEMPLATE)
    avail = set(model.studies())
    print(f"Available studies: {sorted(avail)}")

    results = list(seen.values())

    for t_nm in t_sweep_nm:
        rid = _make_id_phys(A_NM, W_NM, HX_NM, HY_NM, t_nm)
        if rid in seen:
            r = seen[rid]
            print(f"[t={t_nm:3.0f}nm] CACHED  G_o={r['G_o']*100:.1f}%  "
                  f"G_m={r['G_m']*100:.1f}%  f_o={r.get('f_o_nm',0):.0f}nm  "
                  f"f_m={r.get('f_m_GHz',0):.2f}GHz")
            continue

        g = Geometry(
            a  = A_NM  * 1e-9,
            w  = W_NM  * 1e-9,
            t  = t_nm  * 1e-9,
            hx = HX_NM * 1e-9,
            hy = HY_NM * 1e-9,
        )

        print(f"\n[t={t_nm:3.0f}nm] a={A_NM:.0f} w={W_NM:.0f} "
              f"hx={HX_NM:.0f} hy={HY_NM:.0f} nm")

        for name, val in g.as_dict().items():
            model.parameter(name, f"{val}[m]")

        try:
            k_m, freqs_m, fy_arr = _mech_sweep(model, g, args.study_mech,
                                                args.n_k, args.n_bands_mech,
                                                compute_fy=compute_fy)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gp_m = breathing_gap(freqs_m, fy_arr)

            k_o, freqs_o = _opt_sweep(model, g, args.study_opt,
                                      args.n_k, args.n_bands_opt)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gp_o = optical_gap_tracked(freqs_o)

            G_o   = gp_o.normalized_gap if gp_o.found else 0.0
            G_m   = gp_m.normalized_gap if gp_m.found else 0.0
            f_o_c = gp_o.f_center if gp_o.found else OPT_F_TARGET
            f_m_c = gp_m.f_center if gp_m.found else 0.0

            rec = dict(
                id=rid,
                t_nm=float(t_nm),
                params=g.as_dict(),
                G_o=G_o, G_m=G_m,
                f_o_nm=C0 / f_o_c * 1e9 if f_o_c > 0 else float("nan"),
                f_m_GHz=f_m_c * 1e-9,
                status="success",
            )
            print(f"  G_o={G_o*100:.1f}%  G_m={G_m*100:.1f}%  "
                  f"f_o={rec['f_o_nm']:.0f}nm  f_m={f_m_c*1e-9:.2f}GHz")
        except Exception as e:
            rec = dict(id=rid, t_nm=float(t_nm), params=g.as_dict(),
                       G_o=0.0, G_m=0.0, status="failed", error=str(e))
            print(f"  [FAIL] {e}")

        results.append(rec)
        seen[rid] = rec
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(results, fh, indent=2,
                      default=lambda x: float(x) if hasattr(x, "__float__") else str(x))

    print("\n=== t-sweep complete ===")
    successful = [r for r in results if r.get("status") == "success"]
    if successful:
        best = max(successful, key=lambda r: r.get("G_o", 0))
        print(f"Best G_o: {best['G_o']*100:.1f}%  at t={best['t_nm']:.0f}nm")
        print()
        for r in sorted(successful, key=lambda r: r["t_nm"]):
            print(f"  t={r['t_nm']:3.0f}nm  G_o={r['G_o']*100:.2f}%  "
                  f"G_m={r['G_m']*100:.1f}%  f_o={r.get('f_o_nm',0):.0f}nm  "
                  f"f_m={r.get('f_m_GHz',0):.2f}GHz")


if __name__ == "__main__":
    main()
