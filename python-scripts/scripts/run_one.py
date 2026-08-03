#!/usr/bin/env python3
"""Evaluate ONE candidate and save the result.

Usage:
    python scripts/run_one.py                             # uses configs/run_one.yaml
    python scripts/run_one.py --config path/to/run_one.yaml

All settings live in the YAML config (default: configs/run_one.yaml). The `u`
key takes 5 normalized values (ua, uw, uhx, uhy, ut), each in [0,1], mapping to
(a, w, hx, hy, t) via src/geometry.py:u_to_geometry -- see configs/bounds.yaml
for the physical ranges. A legacy 4-value `u` is still accepted for
backward compatibility and pins thickness t to its minimum (t=t_min).

Stage progress (feasibility / optical / mechanical / score) is printed to
stderr as each stage starts and finishes, so stdout stays pure JSON even when
redirected to a file. Set `quiet: true` in the config to suppress it.
"""
import argparse
import json
import os
import sys
import time

import yaml

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


_DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                               "run_one.yaml")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/run_one.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    u = [float(x) for x in cfg["u"]]
    optical = cfg["optical"]
    mech = cfg["mech"]
    quiet = cfg["quiet"]

    if len(u) not in (4, 5):
        ap.error(f"config `u` takes 4 (legacy) or 5 values (ua uw uhx uhy [ut]), got {len(u)}")

    on_stage = None if quiet else make_stderr_reporter()

    rec = evaluate_candidate(u, optical_backend=optical,
                             mech_backend=mech, on_stage=on_stage)

    if not quiet:
        t0 = time.time()
        print("[saving] writing result record ... running", file=sys.stderr)
    save_result(rec)
    if not quiet:
        print(f"[saving] writing result record ... done ({_fmt_time(time.time()-t0)})",
              file=sys.stderr)

    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
