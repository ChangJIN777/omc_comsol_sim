"""Shared stderr stage-progress reporter for scripts/run_one.py and run_loop.py."""
import sys
import time

_STAGES = ["feasibility", "optical", "mechanical", "score"]


def _fmt_time(seconds):
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    if h > 0:
        return f"{h}h{m:02d}m"
    if m > 0:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def make_stderr_reporter(prefix="", stages=None):
    """Build an on_stage(name, event, info) callback that prints
    "[i/n] stage ... running/done/skipped (elapsed)" lines to stderr.

    stages: optional ordered stage-name list (defaults to _STAGES). Pass a
    custom order when the caller's pipeline runs stages in a different sequence
    (e.g. run_opt_comsol.py runs mechanical before optical)."""
    stages = stages if stages is not None else _STAGES
    t_start = {}

    def report(name, event, info=None):
        idx = stages.index(name) + 1
        n = len(stages)
        label = name if info in (None, "start", "done") else f"{name} ({info})"
        if event == "start":
            t_start[name] = time.time()
            print(f"{prefix}[{idx}/{n}] {label} ... running", file=sys.stderr)
        elif event == "skip":
            print(f"{prefix}[{idx}/{n}] {label} ... skipped", file=sys.stderr)
        elif event == "done":
            dt = time.time() - t_start.get(name, time.time())
            print(f"{prefix}[{idx}/{n}] {name} ... done ({_fmt_time(dt)})", file=sys.stderr)
        elif event == "fail":
            dt = time.time() - t_start.get(name, time.time())
            print(f"{prefix}[{idx}/{n}] {name} ... FAILED ({_fmt_time(dt)}): {info}",
                  file=sys.stderr)

    return report
