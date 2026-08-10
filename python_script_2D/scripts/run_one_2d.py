#!/usr/bin/env python3
"""Evaluate ONE 2D boomerang-cell candidate and save the result.

Usage:
    python scripts/run_one_2d.py                             # uses configs/run_one_2d.yaml
    python scripts/run_one_2d.py --config path/to/config.yaml

All settings live in the YAML config (default: configs/run_one_2d.yaml). The
`u` key takes 6 normalized values (ua, uw, ur, ur1, ur2, uth), each in [0,1],
mapping to (a, w, r, r1, r2, th) via src/geometry2d.py:u_to_geometry -- see
configs/bounds_2d.yaml for the physical ranges.

Stage progress (feasibility / mechanical / optical / score) is printed to
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

from objective2d import evaluate_candidate         # noqa: E402
from database import save_result                   # noqa: E402

_STAGES = ["feasibility", "mechanical", "optical", "score"]


def _fmt_time(seconds):
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    if h > 0:
        return f"{h}h{m:02d}m"
    if m > 0:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def make_stderr_reporter():
    """on_stage(name, event, info) callback -> "[i/n] stage ... running/done" to stderr."""
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
                               "run_one_2d.yaml")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/run_one_2d.yaml)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    u = [float(x) for x in cfg["u"]]
    if len(u) != 6:
        ap.error(f"config `u` takes 6 values (ua uw ur ur1 ur2 uth), got {len(u)}")

    quiet = cfg.get("quiet", False)
    on_stage = None if quiet else make_stderr_reporter()

    rec = evaluate_candidate(u, mech_backend=cfg.get("mech_backend", "comsol"),
                             optical_backend=cfg.get("optical_backend", "none"),
                             require_mech=cfg.get("require_mech", True),
                             require_opt=cfg.get("require_opt", False),
                             on_stage=on_stage)

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
