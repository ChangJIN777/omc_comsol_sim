#!/usr/bin/env python3
"""Closed-loop optimization of the 2D boomerang phononic cell: ask -> evaluate
-> save -> tell.

Usage:
    python scripts/run_loop_2d.py                             # uses configs/run_loop_2d.yaml
    python scripts/run_loop_2d.py --config path/to/config.yaml

All settings live in the YAML config (default: configs/run_loop_2d.yaml):
n_init, n_iter, mech/optical backends and require flags, and the optimizer
backend. Only the mechanical backend is implemented so far -- keep
require_opt: false / optical_backend: none until src/optical_comsol_2d.py is
built (see docs/optical_2d_plan.md).

The scored quantity is the COMPLETE gap (`gap_mode: "complete"` in
configs/targets_2d.yaml), so `parity=complete` on every line is expected; it
becomes `evenz`/`oddz` only under gap_mode: "symmetry".

The per-iteration line prints the scored gap; each family's own gap is in the
saved record (`mechanical_gap_{evenz,oddz}` and their edges -- see
objective2d._mechanical_gap), which is what to re-rank on later, or to explain a
narrow complete gap, rather than re-reading every `bands_npz`.

Expect a plateau early on: every candidate with no complete gap in the target
window scores an identical -1.40 (G_m = 0 and f_center = 0 together saturate the
frequency penalty). That is a property of the scoring function, not of the
sampler -- see configs/targets_2d.yaml's min_gap commentary.
"""
import argparse
import os
import sys

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from objective2d import evaluate_candidate                       # noqa: E402
from database import save_result, load_completed_results, best_result  # noqa: E402
from optimizer import initialize_optimizer, ask, tell             # noqa: E402


_DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                               "run_loop_2d.yaml")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/run_loop_2d.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    n_init = cfg["n_init"]
    n_iter = cfg["n_iter"]
    mech_backend = cfg.get("mech_backend", "comsol")
    optical_backend = cfg.get("optical_backend", "none")
    require_mech = cfg.get("require_mech", True)
    require_opt = cfg.get("require_opt", False)
    backend = cfg.get("backend", "auto")

    data = load_completed_results()
    opt = initialize_optimizer(data, backend=backend, n_init=n_init)
    print(f"resuming with {len(data)} prior results; backend={type(opt).__name__}")

    total = n_init + n_iter
    for it in range(total):
        u = ask(opt)
        rec = evaluate_candidate(u, mech_backend=mech_backend,
                                 optical_backend=optical_backend,
                                 require_mech=require_mech,
                                 require_opt=require_opt,
                                 save_bands=cfg.get("save_bands", True),
                                 bands_dir=cfg.get("bands_dir"))
        save_result(rec)
        tell(opt, u, rec["score"],
             constraints={"mechanical_gap": rec.get("mechanical_gap"),
                          "optical_gap": rec.get("optical_gap")})
        print(f"[{it:3d}] score={rec['score']:+.4f} "
              f"Gm={rec.get('mechanical_gap')} "
              f"parity={rec.get('mech_parity')} "
              f"Go={rec.get('optical_gap')} "
              f"status={rec['status']}")

    b = best_result()
    if b:
        print("\nBEST:", {k: b[k] for k in
                          ("id", "score", "mechanical_gap", "optical_gap")})
        print("params (nm):", {k: round(v * 1e9, 1)
                               for k, v in b["params"].items()})


if __name__ == "__main__":
    main()
