#!/usr/bin/env python3
"""Closed-loop optimization: ask -> evaluate -> save -> tell.

Usage:
    python scripts/run_loop.py                             # uses configs/run_loop.yaml
    python scripts/run_loop.py --config path/to/run_loop.yaml

All settings live in the YAML config (default: configs/run_loop.yaml): n_init,
n_iter, optical/mech backends, and the optimizer backend. For a numpy-only
pre-screen keep optical: surrogate + mech: surrogate_stub; for full physics on
the Mac use optical: mpb + mech: comsol.
"""
import argparse
import os
import sys

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from objective import evaluate_candidate          # noqa: E402
from database import save_result, load_completed_results, best_result  # noqa: E402
from optimizer import initialize_optimizer, ask, tell  # noqa: E402


_DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                               "run_loop.yaml")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/run_loop.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    n_init = cfg["n_init"]
    n_iter = cfg["n_iter"]
    optical = cfg["optical"]
    mech = cfg["mech"]
    backend = cfg["backend"]

    data = load_completed_results()
    opt = initialize_optimizer(data, backend=backend, n_init=n_init)
    print(f"resuming with {len(data)} prior results; backend={type(opt).__name__}")

    total = n_init + n_iter
    for it in range(total):
        u = ask(opt)
        rec = evaluate_candidate(u, optical_backend=optical,
                                 mech_backend=mech)
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
