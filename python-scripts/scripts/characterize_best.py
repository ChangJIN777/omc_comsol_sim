"""Standalone characterization of the best result from opt_results_t5d.json.

Usage:
    python scripts/characterize_best.py [--mph PATH] [--json PATH] [--n-k N]
"""
import argparse, json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import mph

# Import characterization internals from the optimizer script
from scripts.run_opt_comsol import _run_characterization

DEFAULT_MPH  = "comsol/omc_unitcell_iso.mph"
DEFAULT_JSON = "results/opt_results_t5d_iso.json"
DEFAULT_OUT  = "results/figures"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mph",    default=DEFAULT_MPH)
    ap.add_argument("--json",   default=DEFAULT_JSON)
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    ap.add_argument("--n-k",    type=int, default=7)
    ap.add_argument("--n-bands-opt", type=int, default=4)
    ap.add_argument("--study-opt",   default="opt TE")
    ap.add_argument("--study-mech",  default="mech sym")
    ap.add_argument("--no-fy",    action="store_true")
    ap.add_argument("--no-cache", action="store_true",
                    help="Force fresh COMSOL sweeps, ignoring any cached NPZ data")
    ap.add_argument("--opt-mode-solnum", type=int, default=None,
                    help="Force this 1-based eigenmode index for the plotted "
                         "TE mode profile, instead of the 'highest mode below "
                         "the gap' heuristic. Inspect "
                         "results/mode_exports/*_all_modes_manifest.txt "
                         "(dumped automatically every run) to find which "
                         "solnum is the real physical dielectric-band mode "
                         "vs. a spurious scattering-boundary artifact.")
    args = ap.parse_args()

    # Allow relative paths from repo root
    root = os.path.join(os.path.dirname(__file__), "..")
    mph_path  = args.mph  if os.path.isabs(args.mph)  else os.path.join(root, args.mph)
    json_path = args.json if os.path.isabs(args.json) else os.path.join(root, args.json)
    out_dir   = args.out_dir if os.path.isabs(args.out_dir) else os.path.join(root, args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Fake args namespace for _run_characterization
    args.out_fig = os.path.join(out_dir, "opt_progress.png")
    args.skip_characterization = False

    # Load best result
    with open(json_path) as f:
        data = json.load(f)
    succ = [r for r in data if r.get("status") == "success" and r.get("G_o", 0) > 0]
    if not succ:
        print("No successful results found in JSON.")
        sys.exit(1)
    best = max(succ, key=lambda r: r.get("score", -99))
    print(f"Best result: G_o={best['G_o']*100:.1f}%  G_m={best['G_m']*100:.1f}%  "
          f"score={best['score']:+.4f}")

    # Start COMSOL and load model
    print(f"Starting COMSOL... loading {mph_path}")
    client = mph.start()
    model  = client.load(mph_path)

    try:
        out_path = _run_characterization(model, best, args, out_dir=out_dir,
                                         no_cache=args.no_cache)
        print(f"\nFigure saved to: {out_path}")
    finally:
        client.clear()


if __name__ == "__main__":
    main()
