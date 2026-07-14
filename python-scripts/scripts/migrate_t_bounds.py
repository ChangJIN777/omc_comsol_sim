#!/usr/bin/env python3
"""Re-map u[4] (thickness) in stored results after changing t_max.

When t_max changes (e.g. 450nm → 500nm), old u4 values map to wrong physical t.
This script:
  1. Recovers physical t from the OLD bounds
  2. Re-encodes u4 under the NEW bounds
  3. Recomputes the result id (SHA1 of new u-vector)
  4. Re-scores every entry with the current scoring constants

Run ONCE after changing t_max in bounds.yaml, before the next optimization.

Usage:
    python scripts/migrate_t_bounds.py --t-max-old 450 --in results/opt_results_t5d.json
"""
import argparse
import hashlib
import json
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

C0 = 299_792_458.0
OPT_F_TARGET = C0 / 1550e-9
G_MIN = 0.20
LAMBDA_OPT   = 20.0
LAMBDA_MECH  = 5.0
LAMBDA_OPT_F = 3.0
LAMBDA_MECH_F = 0.0


def score_result(G_o, G_m, f_o_c, f_m_c):
    import math
    if not (math.isfinite(G_o) and math.isfinite(G_m)):
        return -5.0
    base   = G_o + G_m
    pen_g  = (LAMBDA_OPT  * max(0.0, G_MIN - G_o) ** 2
              + LAMBDA_MECH * max(0.0, G_MIN - G_m) ** 2)
    pen_of = LAMBDA_OPT_F * ((f_o_c - OPT_F_TARGET) / OPT_F_TARGET) ** 2
    return float(base - pen_g - pen_of)


def make_id(u):
    key = json.dumps([round(float(x), 6) for x in u] + ["comsol", "comsol"])
    return "opt_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-max-old", type=float, default=450.0,
                    help="Old t_max in nm (before the bounds change)")
    ap.add_argument("--t-min",     type=float, default=220.0,
                    help="t_min in nm (unchanged)")
    ap.add_argument("--t-max-new", type=float, default=500.0,
                    help="New t_max in nm (after the bounds change)")
    ap.add_argument("--in",  dest="infile",  default="results/opt_results_t5d.json")
    ap.add_argument("--out", dest="outfile", default=None,
                    help="Output file (default: overwrite --in)")
    args = ap.parse_args()

    outfile = args.outfile or args.infile

    t_min     = args.t_min * 1e-9
    t_max_old = args.t_max_old * 1e-9
    t_max_new = args.t_max_new * 1e-9
    span_old  = t_max_old - t_min
    span_new  = t_max_new - t_min

    with open(args.infile) as f:
        data = json.load(f)

    n_remapped = 0
    n_rescored = 0
    for r in data:
        u = r.get("u")
        if u is None:
            continue

        u = list(u)
        changed = False

        # Re-map u[4] if present and non-zero
        if len(u) == 5:
            u4_old = u[4]
            t_phys = t_min + u4_old * span_old       # recover physical t
            u4_new = (t_phys - t_min) / span_new     # re-encode with new span
            u4_new = max(0.0, min(1.0, u4_new))
            if abs(u4_new - u4_old) > 1e-9:
                u[4] = round(u4_new, 6)
                changed = True
        elif len(u) == 4:
            u = u + [0.0]   # legacy 4-D: t=t_min, u4=0 (same in both old and new)
            changed = True

        if changed:
            r["u"] = u
            r["id"] = make_id(u)
            n_remapped += 1

        # Re-score with current constants
        if r.get("status") == "success":
            G_o  = float(r.get("G_o", 0) or 0)
            G_m  = float(r.get("G_m", 0) or 0)
            f_oc = float(r.get("optical_center_frequency", OPT_F_TARGET) or OPT_F_TARGET)
            f_mc = float(r.get("mechanical_center_frequency", 0) or 0)
            if f_oc <= 0:
                f_oc = OPT_F_TARGET
            r["score"] = score_result(G_o, G_m, f_oc, f_mc)
            n_rescored += 1

    with open(outfile, "w") as f:
        json.dump(data, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))

    print(f"Migration complete:")
    print(f"  t range: {args.t_max_old:.0f}nm -> {args.t_max_new:.0f}nm  "
          f"(t_min={args.t_min:.0f}nm)")
    print(f"  u4 re-mapped : {n_remapped} entries")
    print(f"  re-scored    : {n_rescored} entries")
    print(f"  written to   : {outfile}")

    # Show top 5 by score after migration
    succ = [r for r in data if r.get("status") == "success" and r.get("G_o", 0) > 0]
    succ.sort(key=lambda r: r.get("score", -99), reverse=True)
    print(f"\nTop 5 after migration:")
    for r in succ[:5]:
        p = r.get("params", {})
        print(f"  a={p.get('a',0)*1e9:.0f} t={p.get('t',0.22e-6)*1e9:.0f}nm "
              f"hx={p.get('hx',0)*1e9:.0f} hy={p.get('hy',0)*1e9:.0f}  "
              f"G_o={r['G_o']*100:.1f}% G_m={r['G_m']*100:.1f}%  "
              f"score={r['score']:.4f}")


if __name__ == "__main__":
    main()
