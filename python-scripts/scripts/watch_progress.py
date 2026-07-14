#!/usr/bin/env python3
"""Live progress monitor for a running optimization.

Polls the results JSON every --interval seconds and prints a compact dashboard:
  - ASCII progress bar (done / expected total)
  - Elapsed time, rate, ETA
  - Best result found so far (this session or all-time)
  - Pareto front summary (top-5 by score)

Usage:
    python scripts/watch_progress.py                        # defaults
    python scripts/watch_progress.py --n-iter 40           # tell it the total
    python scripts/watch_progress.py --file results/opt_results_t5d.json --n-iter 40
    python scripts/watch_progress.py --once                 # print once and exit
"""
import argparse
import json
import os
import sys
import time


# ── helpers ───────────────────────────────────────────────────────────────────

def _fmt_time(seconds):
    if seconds == float("inf") or seconds > 864000:
        return "?"
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    if h > 0:
        return f"{h}h{m:02d}m"
    if m > 0:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def _bar(done, total, width=30):
    frac = min(done / max(total, 1), 1.0)
    filled = int(frac * width)
    return f"[{'█' * filled}{'░' * (width - filled)}]"


def _read_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return []


def _render(data, n_iter, baseline_count, t_session_start):
    """Render one dashboard frame."""
    lines = []

    total = len(data)
    new_count = max(0, total - baseline_count)   # entries added this session

    succ_all = [r for r in data if r.get("status") == "success" and r.get("G_o", 0) > 0]
    succ_new = [r for r in data[baseline_count:] if r.get("status") == "success"]

    elapsed = time.time() - t_session_start
    rate_per_min = new_count / elapsed * 60 if elapsed > 0 else 0   # evals/min

    sec_per_eval = elapsed / new_count if new_count > 0 else float("inf")
    if n_iter and new_count > 0:
        remaining_s = (n_iter - new_count) * sec_per_eval if new_count < n_iter else 0
    else:
        remaining_s = float("inf")

    # ── header ────────────────────────────────────────────────────────────────
    lines.append("")
    lines.append("=" * 68)
    lines.append("  OMC OPTIMIZATION MONITOR")
    lines.append("=" * 68)

    # ── progress bar ─────────────────────────────────────────────────────────
    if n_iter:
        pct = new_count / n_iter * 100
        bar = _bar(new_count, n_iter)
        lines.append(f"  Progress : {bar} {new_count}/{n_iter} ({pct:.0f}%)")
    else:
        lines.append(f"  Progress : {new_count} new evaluations (total: {total})")

    # ── timing ────────────────────────────────────────────────────────────────
    eta = f"ETA ~{_fmt_time(remaining_s)}" if remaining_s < 1e8 else "ETA unknown"
    lines.append(f"  Elapsed  : {_fmt_time(elapsed)}    Rate: {rate_per_min:.2f} evals/min    {eta}")
    lines.append(f"  File     : {_json_path}  (last-modified: {_mtime_str(_json_path)})")

    # ── best result (new session) ──────────────────────────────────────────────
    lines.append("")
    if succ_new:
        best = max(succ_new, key=lambda r: r.get("score", -99))
        p = best.get("params", {})
        lines.append("  ┌── Best this session ─────────────────────────────────────")
        lines.append(f"  │  G_o = {best['G_o']*100:.1f}%   G_m = {best['G_m']*100:.1f}%   "
                     f"score = {best.get('score', 0):+.4f}")
        lines.append(f"  │  a={p.get('a',0)*1e9:.0f}nm  w={p.get('w',0)*1e9:.0f}nm  "
                     f"hx={p.get('hx',0)*1e9:.0f}nm  hy={p.get('hy',0)*1e9:.0f}nm  "
                     f"t={p.get('t',0.22e-6)*1e9:.0f}nm")
        lines.append(f"  │  f_o={best.get('f_o_nm',0):.0f}nm   "
                     f"f_m={best.get('f_m_GHz',0):.2f}GHz")
        lines.append("  └─────────────────────────────────────────────────────────")
    else:
        lines.append("  (no new successful evaluations yet)")

    # ── top 5 across all time ─────────────────────────────────────────────────
    lines.append("")
    lines.append("  Top-5 all-time by score:")
    lines.append(f"  {'a':>5} {'t':>5} {'hx':>5} {'hy':>5} {'G_o':>6} {'G_m':>6} "
                 f"{'f_o':>7} {'f_m':>8} {'score':>8}")
    lines.append("  " + "-" * 60)
    succ_all.sort(key=lambda r: r.get("score", -99), reverse=True)
    for r in succ_all[:5]:
        p = r.get("params", {})
        tag = " ← NEW" if r in data[baseline_count:] else ""
        lines.append(
            f"  {p.get('a',0)*1e9:5.0f} {p.get('t',0.22e-6)*1e9:5.0f} "
            f"{p.get('hx',0)*1e9:5.0f} {p.get('hy',0)*1e9:5.0f} "
            f"{r.get('G_o',0)*100:5.1f}% {r.get('G_m',0)*100:5.1f}% "
            f"{r.get('f_o_nm',0):7.0f} {r.get('f_m_GHz',0):7.2f}GHz "
            f"{r.get('score',0):+8.4f}{tag}"
        )

    # ── targets ───────────────────────────────────────────────────────────────
    lines.append("")
    lines.append("  Targets: G_o ≥ 20%   G_m ≥ 20%   f_o ≈ 1550nm   f_m ∈ 5-10GHz")
    lines.append("=" * 68)
    return "\n".join(lines)


def _mtime_str(path):
    try:
        return time.strftime("%H:%M:%S", time.localtime(os.path.getmtime(path)))
    except Exception:
        return "?"


# ── main ──────────────────────────────────────────────────────────────────────

ap = argparse.ArgumentParser(description="Monitor optimization progress.")
ap.add_argument("--file",     default="results/opt_results_t5d.json",
                help="Path to the results JSON being written by the optimizer")
ap.add_argument("--n-iter",   type=int, default=None,
                help="Total iterations expected (for progress bar and ETA)")
ap.add_argument("--interval", type=float, default=30.0,
                help="Refresh interval in seconds (default: 30)")
ap.add_argument("--once",     action="store_true",
                help="Print once and exit (no loop)")
args = ap.parse_args()

_json_path = args.file

# Baseline: how many entries existed before we started watching
baseline_data = _read_json(_json_path)
baseline_count = len(baseline_data)
t_start = time.time()

print(f"Watching {_json_path}  (baseline: {baseline_count} entries)")
if args.n_iter:
    print(f"Expecting {args.n_iter} new iterations.")
print("Press Ctrl-C to stop.\n")

try:
    while True:
        data = _read_json(_json_path)
        frame = _render(data, args.n_iter, baseline_count, t_start)
        # Clear screen and redraw (ANSI escape)
        if not args.once:
            print("\033[2J\033[H", end="")
        print(frame)
        sys.stdout.flush()
        if args.once:
            break
        time.sleep(args.interval)
except KeyboardInterrupt:
    print("\nMonitor stopped.")
