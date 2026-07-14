#!/usr/bin/env python3
"""Evaluate ONE candidate and save the result.

Examples:
    python scripts/run_one.py --u 0.5 0.6 0.5 0.6 --optical surrogate --mech surrogate_stub
    python scripts/run_one.py --u 0.5 0.6 0.5 0.6 --optical mpb --mech comsol
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from objective import evaluate_candidate          # noqa: E402
from database import save_result                  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--u", nargs=4, type=float, required=True,
                    metavar=("ua", "uw", "uhx", "uhy"))
    ap.add_argument("--optical", default="surrogate",
                    choices=["surrogate", "mpb", "comsol"])
    ap.add_argument("--mech", default="surrogate_stub",
                    choices=["comsol", "surrogate_stub"])
    args = ap.parse_args()

    rec = evaluate_candidate(args.u, optical_backend=args.optical,
                             mech_backend=args.mech)
    save_result(rec)
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
