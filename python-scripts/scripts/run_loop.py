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
from geometry import n_dim                        # noqa: E402
from database import (save_result, load_completed_results, best_result,  # noqa: E402
                      db_path)
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
    # Optional keys, absent from hole configs (see configs/run_loop_rib.yaml).
    celltype = cfg.get("celltype", "hole")
    rib_options = cfg.get("rib") or None
    # `tidy3d:` section, read by optical: tidy3d only. Its max_credits is the
    # budget of this whole loop (src/optical_tidy3d.py).
    tidy3d_options = cfg.get("tidy3d") if optical == "tidy3d" else None
    xsect = cfg.get("xsect")     # rect (default) | wedge, see run_one.py
    db = db_path(celltype)       # rib results have their own sqlite file

    data = load_completed_results(db)
    opt = initialize_optimizer(data, backend=backend, n_init=n_init,
                               n_dim=n_dim(celltype))
    print(f"resuming with {len(data)} prior results; backend={type(opt).__name__}")

    # tidy3d.keep_hdf5: best -> only the best cell of THIS run keeps its .hdf5
    # files (optical_tidy3d.Hdf5Keeper; the best of an earlier run is not
    # known here and keeps its files).
    keeper = None
    if optical == "tidy3d":
        from optical_tidy3d import Hdf5Keeper
        keeper = Hdf5Keeper(tidy3d_options)

    total = n_init + n_iter
    for it in range(total):
        u = ask(opt)
        rec = evaluate_candidate(u, optical_backend=optical,
                                 mech_backend=mech, celltype=celltype,
                                 rib_options=rib_options,
                                 tidy3d_options=tidy3d_options, xsect=xsect)
        save_result(rec, db)
        if rec.get("optical_tidy3d_stop"):
            # Tidy3D did not run (cloud not allowed, dry run, credit budget
            # reached, ...) and would not for the next candidate either. The
            # record is saved; it is not a result, so the optimizer is not told.
            print(f"[{it:3d}] STOP ({rec['optical_tidy3d_stop']}): {rec['error']}")
            break
        tell(opt, u, rec["score"],
             constraints={"optical_gap": rec.get("optical_gap"),
                          "mechanical_gap": rec.get("mechanical_gap")})
        print(f"[{it:3d}] score={rec['score']:+.4f} "
              f"Go={rec.get('optical_gap')} Gm={rec.get('mechanical_gap')} "
              f"status={rec['status']}")
        if keeper is not None:
            keeper.update(rec)           # after the record is saved

    b = best_result(db)
    if b:
        print("\nBEST:", {k: b[k] for k in
                          ("id", "score", "optical_gap", "mechanical_gap")})
        print("params (nm):", {k: round(v * 1e9, 1)
                               for k, v in b["params"].items()})


if __name__ == "__main__":
    main()
