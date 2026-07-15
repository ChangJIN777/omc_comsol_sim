#!/usr/bin/env python3
"""Evaluate ONE candidate and save the result.

Examples:
    python scripts/run_one.py --u 0.5 0.6 0.5 0.6 0.5 --optical surrogate --mech surrogate_stub
    python scripts/run_one.py --u 0.5 0.6 0.5 0.6 0.5 --optical mpb --mech comsol

--u takes 5 normalized values (ua, uw, uhx, uhy, ut), each in [0,1], mapping
to (a, w, hx, hy, t) via src/geometry.py:u_to_geometry -- see configs/bounds.yaml
for the physical ranges. A legacy 4-value --u is still accepted for
backward compatibility and pins thickness t to its minimum (t=t_min).

Stage progress (feasibility / optical / mechanical / score) is printed to
stderr as each stage starts and finishes, so stdout stays pure JSON even when
redirected to a file. Pass --quiet to suppress it.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from objective import evaluate_candidate          # noqa: E402
from database import save_result                  # noqa: E402

_STAGES = ["feasibility", "optical", "mechanical", "score"]


def _fmt_time(seconds):
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    if h > 0:
        return f"{h}h{m:02d}m"
    if m > 0:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def make_stderr_reporter():
    """Build an on_stage(name, event, info) callback that prints
    "[i/n] stage ... running/done/skipped (elapsed)" lines to stderr."""
    t_start = {}

    def report(name, event, info=None):
        idx = _STAGES.index(name) + 1
        n = len(_STAGES)
        label = name if info in (None, "start", "done") else f"{name} ({info})"
        if event == "start":
            t_start[name] = time.time()
            print(f"[{idx}/{n}] {label} ... running", file=sys.stderr)
        elif event == "skip":
            print(f"[{idx}/{n}] {label} ... skipped", file=sys.stderr)
        elif event == "done":
            dt = time.time() - t_start.get(name, time.time())
            print(f"[{idx}/{n}] {name} ... done ({_fmt_time(dt)})", file=sys.stderr)
        elif event == "fail":
            dt = time.time() - t_start.get(name, time.time())
            print(f"[{idx}/{n}] {name} ... FAILED ({_fmt_time(dt)}): {info}",
                  file=sys.stderr)

    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--u", nargs="+", type=float, required=True,
                    metavar="UA UW UHX UHY [UT]",
                    help="5 normalized values (ua uw uhx uhy ut), each in "
                         "[0,1]. A legacy 4-value form is also accepted and "
                         "pins thickness t to its minimum.")
    ap.add_argument("--optical", default="surrogate",
                    choices=["surrogate", "mpb", "comsol"])
    ap.add_argument("--mech", default="surrogate_stub",
                    choices=["comsol", "surrogate_stub"])
    ap.add_argument("--quiet", action="store_true",
                    help="suppress stage progress on stderr")
    args = ap.parse_args()

    if len(args.u) not in (4, 5):
        ap.error(f"--u takes 4 (legacy) or 5 values (ua uw uhx uhy [ut]), got {len(args.u)}")

    on_stage = None if args.quiet else make_stderr_reporter()

    rec = evaluate_candidate(args.u, optical_backend=args.optical,
                             mech_backend=args.mech, on_stage=on_stage)

    if not args.quiet:
        t0 = time.time()
        print("[saving] writing result record ... running", file=sys.stderr)
    save_result(rec)
    if not args.quiet:
        print(f"[saving] writing result record ... done ({_fmt_time(time.time()-t0)})",
              file=sys.stderr)

    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
