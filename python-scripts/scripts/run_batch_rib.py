#!/usr/bin/env python3
"""Evaluate a QUEUE of designs in one Python process / one COMSOL session.

Usage:
    python scripts/run_batch_rib.py                          # configs/run_batch_rib.yaml
    python scripts/run_batch_rib.py --config path/to/queue.yaml
    python scripts/run_batch_rib.py --config ... --dry-run   # list the queue, solve nothing

The explicit "list of designs" counterpart of run_loop.py. run_loop.py already
keeps one MPh session for all its iterations (src/comsol_client.py caches the
client and the loaded model process-wide), but its designs come from the
optimizer. Here they come from the config: a list of u vectors, a CSV of u
vectors, and/or several run_one-style configs. Every item is evaluated with
objective.evaluate_candidate and saved with save_result(rec, db_path(celltype)),
exactly as run_one.py does -- one record per item, failed ones included.

Why: every `run_one.py` call starts COMSOL, checks the licences out, and gives
them back at exit, so between two calls someone else can take the seat. Here
the seat obtained for the first solve is kept while queued work remains, and
only then.

Licence rules of this script:
  * COMSOL starts at the first item that needs a solver, not before (items
    that fail the feasibility / fabrication gates never start it).
  * No idle keep-alive: there is no waiting between items, and nothing keeps
    the session open once the queue is empty.
  * The session is released (comsol_client.release_client) as soon as no
    remaining item needs COMSOL, and in a `finally` block: at the end of the
    queue, on an error, and on Ctrl-C.
  * If a solve fails with a licence error the queue stops (the record of that
    item is saved, the rest are listed as not run) instead of failing every
    remaining item. Set stop_on_license_error: false to carry on.

To wait for a free seat first (nothing is held while waiting):
    scripts/wait_for_comsol_seat.sh -r 3 -- python scripts/run_batch_rib.py --config ...

Progress goes to stderr; stdout gets one JSON summary of the queue.
Exit status: 0 queue finished, 1 unexpected error, 3 stopped on a licence
error, 4 stopped because optical: tidy3d was not allowed to run (no cloud
opt-in, dry run, credit budget), 130 interrupted (Ctrl-C). Records of the items already evaluated are
in the database in every case.
"""
import argparse
import csv
import json
import os
import re
import sys
import time

import yaml

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(_ROOT, "src"))

import comsol_client                              # noqa: E402
from objective import evaluate_candidate          # noqa: E402
from geometry import n_dim                        # noqa: E402
from database import save_result, db_path         # noqa: E402
from cli_progress import make_stderr_reporter, _fmt_time  # noqa: E402


_DEFAULT_CONFIG = os.path.join(_ROOT, "configs", "run_batch_rib.yaml")

# Same wording as scripts/wait_for_comsol_seat.sh, so its -r retry sees it too.
_LICENCE_ERROR = re.compile(
    r"licen[sc]e error|licensed number of users already reached|"
    r"flexnet licensing error|\(-4,\d+|error: ?-4,|"
    r"could not (obtain|check ?out) licen[sc]e|"
    r"licen[sc]e.*(not available|unavailable)", re.I)


def _log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", file=sys.stderr,
          flush=True)


def _find(path, base):
    """A path from the config: as given, else relative to the config, else to
    python-scripts/."""
    for cand in (path, os.path.join(base, path), os.path.join(_ROOT, path)):
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(f"{path} (also tried under {base} and {_ROOT})")


def _read_csv(path):
    """u vectors from a CSV: one row each; blank lines, '#' comments and a
    non-numeric header row are skipped."""
    rows = []
    with open(path, newline="") as fh:
        for n, row in enumerate(csv.reader(fh), 1):
            cells = [c.strip() for c in row if c.strip()]
            if not cells or cells[0].startswith("#"):
                continue
            try:
                rows.append([float(c) for c in cells])
            except ValueError:
                if rows:
                    raise ValueError(f"{path}, line {n}: not numeric: {row}")
                continue                           # header row
    return rows


def _check_u(u, celltype, where):
    if celltype == "hole":
        if len(u) not in (4, 5):
            raise ValueError(f"{where}: `u` takes 4 (legacy) or 5 values for "
                             f"the hole cell, got {len(u)}")
    elif len(u) != n_dim(celltype):
        raise ValueError(f"{where}: `u` takes {n_dim(celltype)} values for "
                         f"celltype {celltype!r}, got {len(u)}")
    if not all(0.0 <= x <= 1.0 for x in u):
        raise ValueError(f"{where}: `u` values must be in [0, 1], got {u}")


def build_queue(cfg, base):
    """Resolve the config into a list of fully specified items.

    Top-level celltype / optical / mech / rib / xsect / tidy3d are the
    defaults of every item.
    A queue entry is
      [u...]                           a u vector with the defaults
      {u: [...], optical: ..., mech: ..., celltype: ..., rib: {...},
       xsect: ..., tidy3d: {...}}      a u vector with its own overrides
                                       (rib: and tidy3d: are merged over the
                                       default sections)
      {config: path}                   a run_one-style config (its own u,
                                       celltype, optical, mech, rib, xsect,
                                       tidy3d)
    As in run_one.py: xsect is rect (default) | wedge, and the tidy3d section
    is used by optical: tidy3d only (cloud opt-in and credit budget; one
    budget for the whole queue).
    queue_csv adds one item per row of a CSV, with the defaults.
    Everything is validated here, before any solver or licence is touched.
    """
    defaults = dict(celltype=cfg.get("celltype", "hole"),
                    optical=cfg.get("optical"), mech=cfg.get("mech"),
                    rib=cfg.get("rib") or {}, xsect=cfg.get("xsect"),
                    tidy3d=cfg.get("tidy3d") or {})
    entries = [(f"queue[{i}]", e) for i, e in enumerate(cfg.get("queue") or [])]
    if cfg.get("queue_csv"):
        path = _find(cfg["queue_csv"], base)
        entries += [(f"{os.path.basename(path)} row {i + 1}", u)
                    for i, u in enumerate(_read_csv(path))]
    items = []
    for where, e in entries:
        if isinstance(e, dict) and "config" in e:
            path = _find(e["config"], base)
            with open(path) as fh:
                one = yaml.safe_load(fh)
            item = dict(celltype=one.get("celltype", "hole"),
                        optical=one["optical"], mech=one["mech"],
                        rib=one.get("rib") or {}, u=one["u"],
                        xsect=one.get("xsect"),
                        tidy3d=one.get("tidy3d") or {})
            where = f"{where} ({os.path.basename(path)})"
        elif isinstance(e, dict):
            if "u" not in e:
                raise ValueError(f"{where}: needs `u` or `config`")
            unknown = set(e) - {"u", "celltype", "optical", "mech", "rib",
                                "label", "xsect", "tidy3d"}
            if unknown:
                raise ValueError(f"{where}: unknown key(s) {sorted(unknown)}")
            item = dict(defaults, **{k: e[k] for k in ("celltype", "optical",
                                                       "mech", "xsect")
                                     if k in e})
            item["rib"] = dict(defaults["rib"], **(e.get("rib") or {}))
            item["tidy3d"] = dict(defaults["tidy3d"], **(e.get("tidy3d") or {}))
            item["u"] = e["u"]
            if e.get("label"):
                where = f"{where} ({e['label']})"
        else:
            item = dict(defaults, u=e)
        if item["optical"] is None or item["mech"] is None:
            raise ValueError(f"{where}: `optical` and `mech` must be set, in "
                             "the entry or at the top level of the config")
        item["u"] = [float(x) for x in item["u"]]
        _check_u(item["u"], item["celltype"], where)
        if item["celltype"] != "rib":
            item["rib"] = {}
        # validated here, before any solver or licence: cross-section (and
        # that the optical engine can model it) and the tidy3d section
        from objective import check_xsect_backend
        try:
            item["xsect"] = check_xsect_backend(item.get("xsect"), item["optical"])
            if "rxtal_deg" in item["rib"]:     # 0 or 45 (mod 90) only
                from geometry import resolve_rxtal
                resolve_rxtal(item["rib"]["rxtal_deg"])
            if item["optical"] == "tidy3d":
                from optical_tidy3d import Tidy3DOptions
                Tidy3DOptions.from_config(item["tidy3d"])
            else:
                item["tidy3d"] = {}
        except ValueError as err:
            raise ValueError(f"{where}: {err}")
        item["where"] = where
        items.append(item)
    if not items:
        raise ValueError("the queue is empty: give `queue` and/or `queue_csv`")
    return items


def _needs_comsol(item):
    return item["optical"] == "comsol" or item["mech"] == "comsol"


def _release(why):
    if comsol_client.release_client():
        _log(f"[session] COMSOL session released ({why})")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=_DEFAULT_CONFIG,
                    help="Path to the YAML config (default: configs/run_batch_rib.yaml)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the resolved queue and exit; no solver, no licence")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    items = build_queue(cfg, os.path.dirname(os.path.abspath(args.config)))
    quiet = bool(cfg.get("quiet", False))
    stop_on_licence = bool(cfg.get("stop_on_license_error", True))

    if args.dry_run:
        print(json.dumps([{k: it[k] for k in
                           ("where", "celltype", "u", "optical", "mech", "rib",
                            "xsect", "tidy3d")}
                          for it in items], indent=2))
        return 0

    on_stage = None if quiet else make_stderr_reporter()
    n = len(items)
    _log(f"queue: {n} item(s), {sum(map(_needs_comsol, items))} need COMSOL")
    done, stopped, rc = [], None, 0
    t_start = time.time()
    # tidy3d.keep_hdf5 (top-level tidy3d section): `best` keeps the .hdf5 files
    # of the best-scoring Tidy3D item of this queue only; `none` is applied per
    # item by the adapter. See optical_tidy3d.Hdf5Keeper.
    keeper = None
    if any(it["optical"] == "tidy3d" for it in items):
        from optical_tidy3d import Hdf5Keeper
        keeper = Hdf5Keeper(cfg.get("tidy3d"))
    try:
        for i, it in enumerate(items):
            _log(f"[{i + 1}/{n}] {it['where']}: celltype={it['celltype']} "
                 f"optical={it['optical']} mech={it['mech']} u={it['u']}")
            t0 = time.time()
            rec = evaluate_candidate(it["u"], optical_backend=it["optical"],
                                     mech_backend=it["mech"], on_stage=on_stage,
                                     celltype=it["celltype"],
                                     rib_options=it["rib"] or None,
                                     tidy3d_options=it["tidy3d"] or None,
                                     xsect=it["xsect"], verbose=not quiet)
            save_result(rec, db_path(it["celltype"]))      # every item, as run_one.py
            if keeper is not None:
                keeper.update(rec)                         # after it is saved
            done.append(dict(where=it["where"], id=rec["id"],
                             status=rec["status"], score=rec.get("score"),
                             optical_gap=rec.get("optical_gap"),
                             mechanical_gap=rec.get("mechanical_gap"),
                             error=rec.get("error")))
            _log(f"[{i + 1}/{n}] saved {rec['id']}: status={rec['status']} "
                 f"score={rec.get('score')} Go={rec.get('optical_gap')} "
                 f"Gm={rec.get('mechanical_gap')} ({_fmt_time(time.time() - t0)})")
            if rec.get("optical_tidy3d_stop") not in (None, "dry_run",
                                                      "estimate_only"):
                # Tidy3D did not run (cloud not allowed, no or exhausted
                # credit budget, not installed): the record is saved and the
                # rest of the queue is listed as not run. A dry run or a cost
                # estimate is per item, so the queue carries on with those.
                stopped = (f"Tidy3D stop ({rec['optical_tidy3d_stop']}) at "
                           f"item {i + 1}: {rec['error']}")
                rc = 4
                break
            if (stop_on_licence and rec.get("error")
                    and _LICENCE_ERROR.search(str(rec["error"]))):
                stopped = f"licence error at item {i + 1}: {rec['error']}"
                rc = 3
                break
            # Nothing left that needs COMSOL: give the seat back now rather
            # than at the end of the queue.
            if (i + 1 < n and comsol_client.session_active()
                    and not any(map(_needs_comsol, items[i + 1:]))):
                _release("no remaining item needs COMSOL")
    except KeyboardInterrupt:
        stopped, rc = "interrupted (Ctrl-C)", 130
    except Exception as e:
        stopped, rc = f"unexpected error: {type(e).__name__}: {e}", 1
        raise
    finally:
        _release(stopped or "queue empty")
        not_run = [it["where"] for it in items[len(done):]]
        if stopped:
            # on stderr, so wait_for_comsol_seat.sh -r can match a licence error
            _log(f"queue stopped: {stopped}")
        if not_run:
            _log(f"not run ({len(not_run)}): " + "; ".join(not_run))
        print(json.dumps(dict(n_items=n, n_done=len(done), stopped=stopped,
                              not_run=not_run, elapsed_s=round(time.time() - t_start, 1),
                              results=done), indent=2))
    return rc


if __name__ == "__main__":
    sys.exit(main())
