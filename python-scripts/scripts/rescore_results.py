#!/usr/bin/env python3
"""Re-score existing result records with the CURRENT score_result() function
and constants (e.g. after changing thresholds or require_opt/require_mech),
without re-running COMSOL — score is a pure function of already-stored
G_o/G_m/f_o_c/f_m_c.

By default each record is rescored using its OWN stored require_opt/
require_mech/g_min_opt/g_min_mech (self-describing, per evaluate()) if
present, so a file mixing records from different requirement settings
rescores each consistently with how it was actually produced. Pass
--require-opt/--no-require-opt/--require-mech/--no-require-mech/--g-min-opt/
--g-min-mech to instead FORCE the same settings onto every record (e.g. to
see how an optical-only file WOULD have scored under the old co-optimized
objective, or vice versa).

The original score is preserved once under "score_prev" (never overwritten
on repeated re-scoring runs, so the very first pre-change value survives).

Usage:
    python scripts/rescore_results.py FILE [FILE ...]
    python scripts/rescore_results.py --no-require-mech FILE
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from run_opt_comsol import score_result, OPT_F_TARGET, G_MIN_OPT, G_MIN_MECH


def rescore_file(path, force_require_opt=None, force_require_mech=None,
                 force_g_min_opt=None, force_g_min_mech=None):
    with open(path) as f:
        data = json.load(f)

    n_changed = 0
    for r in data:
        if r.get("status") != "success":
            continue
        G_o = float(r.get("G_o", 0.0) or 0.0)
        G_m = float(r.get("G_m", 0.0) or 0.0)
        f_o_c = float(r.get("optical_center_frequency", OPT_F_TARGET) or OPT_F_TARGET)
        f_m_c = float(r.get("mechanical_center_frequency", 0.0) or 0.0)
        if f_o_c <= 0:
            f_o_c = OPT_F_TARGET
        params = r.get("params") or {}
        t_val = float(params["t"]) if params.get("t") is not None else None
        w_val = float(params["w"]) if params.get("w") is not None else None
        a_val = float(params["a"]) if params.get("a") is not None else None

        require_opt = (force_require_opt if force_require_opt is not None
                       else bool(r.get("require_opt", True)))
        require_mech = (force_require_mech if force_require_mech is not None
                        else bool(r.get("require_mech", True)))
        g_min_opt = (force_g_min_opt if force_g_min_opt is not None
                    else float(r.get("g_min_opt", G_MIN_OPT)))
        g_min_mech = (force_g_min_mech if force_g_min_mech is not None
                     else float(r.get("g_min_mech", G_MIN_MECH)))

        new_score = score_result(G_o, G_m, f_o_c, f_m_c, t=t_val, w=w_val, a=a_val,
                                 require_opt=require_opt,
                                 require_mech=require_mech, g_min_opt=g_min_opt,
                                 g_min_mech=g_min_mech)
        old_score = r.get("score")
        if "score_prev" not in r:
            r["score_prev"] = old_score
        if old_score != new_score:
            n_changed += 1
        r["score"] = new_score
        # Keep the record self-describing for whatever mode it was JUST
        # rescored under (matters when --force-* overrides are used).
        r["require_opt"], r["require_mech"] = require_opt, require_mech
        r["g_min_opt"], r["g_min_mech"] = g_min_opt, g_min_mech

    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"{path}: {len(data)} records, {n_changed} scores changed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--require-opt", dest="require_opt", action="store_true", default=None)
    ap.add_argument("--no-require-opt", dest="require_opt", action="store_false")
    ap.add_argument("--require-mech", dest="require_mech", action="store_true", default=None)
    ap.add_argument("--no-require-mech", dest="require_mech", action="store_false")
    ap.add_argument("--g-min-opt", type=float, default=None)
    ap.add_argument("--g-min-mech", type=float, default=None)
    args = ap.parse_args()
    for path in args.files:
        rescore_file(path, force_require_opt=args.require_opt,
                    force_require_mech=args.require_mech,
                    force_g_min_opt=args.g_min_opt,
                    force_g_min_mech=args.g_min_mech)


if __name__ == "__main__":
    main()
