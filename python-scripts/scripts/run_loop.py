#!/usr/bin/env python3
"""Closed-loop optimization: ask -> evaluate -> save -> tell.

Examples (pre-screen, runs anywhere with numpy):
    python scripts/run_loop.py --n-init 40 --n-iter 80 --optical surrogate --mech surrogate_stub

Full physics (on the Mac with MPB + COMSOL):
    python scripts/run_loop.py --n-init 60 --n-iter 200 --optical mpb --mech comsol
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from objective import evaluate_candidate          # noqa: E402
from database import save_result, load_completed_results, best_result  # noqa: E402
from optimizer import initialize_optimizer, ask, tell  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-init", type=int, default=40)
    ap.add_argument("--n-iter", type=int, default=80)
    ap.add_argument("--optical", default="surrogate",
                    choices=["surrogate", "mpb", "comsol"])
    ap.add_argument("--mech", default="surrogate_stub",
                    choices=["comsol", "surrogate_stub"])
    ap.add_argument("--backend", default="auto",
                    choices=["auto", "random", "optuna", "botorch"])
    args = ap.parse_args()

    data = load_completed_results()
    opt = initialize_optimizer(data, backend=args.backend, n_init=args.n_init)
    print(f"resuming with {len(data)} prior results; backend={type(opt).__name__}")

    total = args.n_init + args.n_iter
    for it in range(total):
        u = ask(opt)
        rec = evaluate_candidate(u, optical_backend=args.optical,
                                 mech_backend=args.mech)
        save_result(rec)
        tell(opt, u, rec["score"],
             constraints={"optical_gap": rec.get("optical_gap"),
                          "mechanical_gap": rec.get("mechanical_gap")})
        print(f"[{it:3d}] score={rec['score']:+.4f} "
              f"Go={rec.get('optical_gap')} Gm={rec.get('mechanical_gap')} "
              f"status={rec['status']}")

    b = best_result()
    if b:
        print("\nBEST:", {k: b[k] for k in
                          ("id", "score", "optical_gap", "mechanical_gap")})
        print("params (nm):", {k: round(v * 1e9, 1)
                               for k, v in b["params"].items()})


if __name__ == "__main__":
    main()
