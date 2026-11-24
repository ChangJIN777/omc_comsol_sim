"""Optical band structure via Tidy3D FDTD: alternative engine to COMSOL optics.

Selected with `optical: tidy3d` in a run config (both cell types). This module
is an ADAPTER around python-scripts/diamond_bandStructure_sim.py, which stays
a self-contained command line tool: it is imported from its path, not moved
or rewritten. tidy3d itself is imported lazily (inside load_bs), so the rest
of the pipeline and the tests work on a machine without tidy3d.

COST: Tidy3D runs on the Flexcompute CLOUD and charges FlexCredits. Nothing is
uploaded unless the config has a `tidy3d:` section with `allow_cloud: true`
AND a `max_credits` budget. See Tidy3DOptions and run_optical_tidy3d.

Frames and units
----------------
repo (COMSOL templates): beam z, width y, thickness x, lengths in METERS,
                         k_frac = k*a/pi (X at 1).
module (Tidy3D):         beam x, width y, thickness z, lengths in UM,
                         k in units of 2*pi/a (X at 0.5), so k_frac = 2*k.
Only the axis labels differ; y is the width in both.

Geometry mapping (band_params_kwargs)
-------------------------------------
rib  -> "rib_block": a -> a, t -> t (rib length along the beam),
        w -> w (TOTAL tip-to-tip width, = 2L in both), s -> spine_width (full
        width in both), th -> thickness, r1 -> fillet_r1 (concave root),
        r2 -> fillet_r2 (convex tip). COMSOL fillets are true arcs; the module
        uses `fillet_points` polygon vertices per arc.
hole -> "holey": a -> a, w -> beam_width, t -> thickness (hole cell: t IS the
        thickness), hole_x = 2*hx, hole_y = 2*hy. The repo's hx, hy are
        SEMI-axes (geometry.py: bridge = a - 2*hx, sidewall = w/2 - hy; the
        template's cylinder has radius hx); the module's hole_x, hole_y are
        FULL axes.

Cross-section (xsect, geometry.XSECTS)
--------------------------------------
"rect"  -> BandParams.xsect "rect": the slab only. COMSOL counterparts:
           omc_unitcell.mph (hole), omc_unitcell_rib.mph (rib).
"wedge" -> BandParams.xsect "wedge": slab + triangular underside with base
           angle geometry.WEDGE_ANGLE_DEG (30 deg), the definition of
           comsol/omc_unitcell_iso.mph: under a strip of width b the apex is
           b/(2*sqrt(3)) below the slab. Hole cell: one prism along the beam,
           b = w, the hole cut through both (exact). Rib cell: the roof of
           the FILLETED outline, depth = tan(30 deg) * distance to the
           outline (geometry.rib_roof_distance), so the wedge follows the
           root and tip fillets and is deepest under the junction. The module
           builds it as `wedge_layers` (12) stacked PolySlabs of exactly
           offset outlines: exact in plan view, a staircase in depth with an
           error of at most (deepest point)/(2*wedge_layers), about 6 nm for
           the default cell. The COMSOL template builds the same roof exactly.
           COMSOL counterparts: omc_unitcell_iso.mph (hole),
           omc_unitcell_rib_wedge.mph (rib, build_omc_unitcell_rib('wedge')).

Polarization
------------
The scored family is COMSOL study "opt TE": PEC on the y = 0 mirror plane,
i.e. tangential E = 0 there, so Ey (the component along the width) is EVEN in
y. The module's "TE" mode uses Ey dipoles and symmetry (0, -1, 1): Tidy3D
symmetry -1 on the y = 0 plane is a PEC plane (Ey even, same sector), and +1
on the mid-thickness plane is a PMC plane (in-plane E even across the
thickness, i.e. TE-like). "opt TM" (PMC on y = 0) maps to "TM" = Ez dipoles,
symmetry (0, 1, -1) in the same way.

Mid-thickness (z) mirror: OPTICS ONLY, set by the cross-section
  rect : the mirror exists. The module uses symmetry (0, -1, +1) for TE
         ((0, +1, -1) for TM), i.e. one thickness-parity family. The COMSOL
         templates have the full thickness and no symmetry plane, so the
         COMSOL drivers keep the same family by the parity of each solved
         field (comsol_optics.study_z_family / _modes_at_k). Both engines
         then score the thickness-even (TE-like) modes of the same cell.
  wedge: no mirror. The module uses symmetry (0, -1, 0) and spreads its
         dipoles and monitors over the slab thickness; COMSOL keeps all modes.
         Both engines then score all thickness parities. An FDTD run only
         finds modes its Ey dipoles excite, so a weakly excited (Ez-dominant)
         mode that COMSOL returns can be missing: gap(tidy3d) >= gap(COMSOL)
         is still possible here, and finalists should be confirmed in COMSOL.
Mechanics is not affected: full thickness, all modes.

Scoring (score_resonances)
--------------------------
The resonances of every k are fed to the gap functions used with COMSOL:
rib, k sweep     comsol_optics.guided_interval_gap (coverage = FDTD window)
hole (and any    comsol_optics.optical_gap_tracked on the guided resonances
x_only run)      (fundamental gap at X, tracked over k in a sweep)
The module's own "two lowest guided resonances at X" result is recorded next
to it (optical_tidy3d_detail["x_point"]).

Console output (tidy3d_console)
-------------------------------
Everything this engine prints goes to STDERR, so the stdout of the callers
(scripts/run_one.py: pure JSON) is not touched:
  * Tidy3D's rich console (web progress: upload, task status and progress
    bars, cost, download; log messages) is switched to stderr with
    tidy3d.log.set_logging_console(stderr=True) for the duration of one
    evaluation, and the previous console handler (and its level) is put back
    afterwards;
  * print() of diamond_bandStructure_sim.py and anything else written to
    sys.stdout is caught by contextlib.redirect_stdout(sys.stderr).
verbose=True (default) shows Tidy3D's progress (log level INFO, Batch
verbose) and sends our own step lines to `on_progress`; verbose=False
(`quiet: true` in the scripts) leaves errors only.
"""
from __future__ import annotations

import contextlib
import importlib
import importlib.util
import json
import os
import sys
import time
from dataclasses import dataclass, field, asdict

import numpy as np

from bandgap import Gap
from comsol_optics import (guided_interval_gap, guided_interval_gaps,
                           optical_gap_tracked, resolvable_step_hz, rib_k_grid,
                           select_gap_by_penalized_size)
from geometry import (WEDGE_ANGLE_DEG, diamond_index, resolve_xsect,
                      wedge_depths, z_mirror_symmetric)

ENGINE = "tidy3d"
C0 = 299_792_458.0
M_TO_UM = 1e6

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BS_PATH = os.path.join(_ROOT, "diamond_bandStructure_sim.py")
DEFAULT_RESULTS_DIR = os.path.join("results", "tidy3d")

# repo optical study -> module mode (see "Polarization" above)
STUDY_TO_MODE = {"opt TE": "TE", "opt TM": "TM"}
MODES = ("run", "dry_run", "estimate")
K_MODES = ("x_only", "sweep")
# tidy3d.keep_hdf5: what happens to the downloaded SimulationData files
#   all   keep them for every cell (default)
#   best  loops keep them for the best-scoring cell only (Hdf5Keeper)
#   none  delete them as soon as a cell is analyzed
# The small files of a cell (params.json, plan.json, bands.csv,
# evaluation.json) always stay, and a repeated cell is re-scored from them
# instead of being paid again (_load_cached_analysis).
KEEP_HDF5 = ("all", "best", "none")
NOTE_NAME = "hdf5_note.json"
# default k sampling per cell type: (k_mode, n_k for a sweep)
#   hole: the scored gap is the fundamental gap at X -> 1 simulation.
#   rib : guided_interval_gap needs guided modes at several k (None = the rib
#         option n_k_opt, i.e. as many k as the COMSOL sweep).
K_DEFAULTS = {"hole": ("x_only", 7), "rib": ("sweep", None)}
# BandParams numerics that `tidy3d: params:` may override (never the geometry,
# the material, the polarization or the target, which come from the repo).
PARAM_OVERRIDES = ("freq_span_thz", "min_steps_per_wvl", "num_dipoles",
                   "num_monitors", "dPML", "runtime_fwidth", "t_start_fwidth",
                   "seed", "minQ", "minamp", "maxerr", "fillet_points",
                   "hole_points", "wedge_layers")
_FREQ_SPAN_THZ = 100.0          # BandParams default, used when not overridden


class Tidy3DNotRun(RuntimeError):
    """The optical solve was NOT run (and nothing was charged).

    kind: "not_installed" | "cloud_not_allowed" | "no_budget" | "dry_run" |
          "estimate_only" | "budget_exceeded".
    Every kind repeats for the next candidate, so optimization loops stop
    cleanly on it (the record is saved with status "optical_not_run").
    """

    def __init__(self, kind, message, detail=None):
        super().__init__(message)
        self.tidy3d_stop = kind
        self.detail = dict(detail or {})


# ---- options ----------------------------------------------------------------

@dataclass
class Tidy3DOptions:
    """The `tidy3d:` section of a run config. Unknown keys are an error."""
    allow_cloud: bool = False       # must be exactly `true` to upload anything
    mode: str = "run"               # run | dry_run | estimate
    max_credits: float | None = None        # budget of THIS process (whole loop)
    max_credits_total: float | None = None  # budget over all runs (ledger file)
    k_mode: str | None = None       # x_only | sweep; None = K_DEFAULTS
    n_k: int | None = None          # k points of a sweep; None = K_DEFAULTS
    results_dir: str = DEFAULT_RESULTS_DIR  # relative to python-scripts/
    folder_name: str = "omc_unitcell"       # Tidy3D web folder of the tasks
    reuse_results: bool = True      # re-analyze downloaded results, no cloud
    keep_hdf5: str = "all"          # all | best | none (see KEEP_HDF5)
    index: float | None = None      # None = configs/materials.yaml
    params: dict = field(default_factory=dict)   # PARAM_OVERRIDES only

    @classmethod
    def from_config(cls, cfg) -> "Tidy3DOptions":
        if isinstance(cfg, cls):
            return cfg
        cfg = dict(cfg or {})
        known = set(cls.__dataclass_fields__)
        unknown = sorted(set(cfg) - known)
        if unknown:                     # a typo must not silently change cost
            raise ValueError(f"unknown key(s) under `tidy3d`: {unknown}; "
                             f"expected a subset of {sorted(known)}")
        opts = cls(**cfg)
        for name in ("allow_cloud", "reuse_results"):
            if not isinstance(getattr(opts, name), bool):
                raise ValueError(f"tidy3d.{name} must be true or false, got "
                                 f"{getattr(opts, name)!r}")
        if opts.mode not in MODES:
            raise ValueError(f"tidy3d.mode {opts.mode!r}: use one of {MODES}")
        if opts.keep_hdf5 not in KEEP_HDF5:
            raise ValueError(f"tidy3d.keep_hdf5 {opts.keep_hdf5!r}: use one of "
                             f"{KEEP_HDF5}")
        if opts.k_mode is not None and opts.k_mode not in K_MODES:
            raise ValueError(f"tidy3d.k_mode {opts.k_mode!r}: use one of {K_MODES}")
        for name in ("max_credits", "max_credits_total", "index"):
            val = getattr(opts, name)
            if val is not None:
                val = float(val)
                if not val > 0:
                    raise ValueError(f"tidy3d.{name} must be > 0, got {val}")
                setattr(opts, name, val)
        if opts.n_k is not None:
            opts.n_k = int(opts.n_k)
            if opts.n_k < 2:
                raise ValueError("tidy3d.n_k must be >= 2 (use k_mode: x_only "
                                 "for the X point alone)")
        opts.params = dict(opts.params or {})
        bad = sorted(set(opts.params) - set(PARAM_OVERRIDES))
        if bad:
            raise ValueError(f"tidy3d.params: {bad} cannot be overridden; "
                             f"allowed: {PARAM_OVERRIDES}")
        return opts

    def k_sampling(self, celltype, n_k_default=None):
        """(k_mode, n_k) for a cell type; n_k is None for x_only."""
        if celltype not in K_DEFAULTS:
            raise ValueError(f"unknown celltype {celltype!r}")
        k_mode = self.k_mode or K_DEFAULTS[celltype][0]
        if k_mode == "x_only":
            return k_mode, None
        n_k = self.n_k or n_k_default or K_DEFAULTS[celltype][1] or 13
        return k_mode, int(n_k)

    def results_root(self):
        return os.path.normpath(
            self.results_dir if os.path.isabs(self.results_dir)
            else os.path.join(_ROOT, self.results_dir))

    def index_value(self):
        return float(self.index) if self.index is not None else material_index()


def material_index(path=None) -> float:
    """Diamond refractive index of the repo (configs/materials.yaml): 2.386,
    the module's `diamond_telecom`. Same value as every other optical engine."""
    return diamond_index(path)


def id_settings(options, celltype, n_k_default=None, xsect=None) -> str:
    """What makes two Tidy3D evaluations of the same u different results.

    Goes into the record id next to optical_backend="tidy3d", so a Tidy3D
    record never replaces a COMSOL one, nor an x_only one a sweep, nor a wedge
    one a rectangular one, nor one with another index (records made with
    n = 2.40 keep their ids; n = 2.386 gives new ones). The cost-control keys
    (allow_cloud, mode, budgets, folders) are not part of it.
    """
    opts = Tidy3DOptions.from_config(options)
    k_mode, n_k = opts.k_sampling(celltype, n_k_default)
    return json.dumps([ENGINE, k_mode, n_k, round(opts.index_value(), 6),
                       sorted((k, float(v)) for k, v in opts.params.items()),
                       resolve_xsect(xsect)])


# ---- geometry and k mapping (no tidy3d needed) --------------------------------

def _um(x_m):
    return round(float(x_m) * M_TO_UM, 9)       # fm resolution, drops 1e-16 noise


def band_params_kwargs(g, *, index, f_target, optical_study="opt TE",
                       overrides=None, xsect=None) -> dict:
    """BandParams keyword arguments (um, THz) of a repo geometry (meters, Hz).

    See the module docstring for the mapping and the checks behind it. xsect
    "rect" | "wedge" goes to BandParams.xsect (with the repo's wedge angle),
    which also sets the z component of the module's symmetry tuple.
    """
    if optical_study not in STUDY_TO_MODE:
        raise ValueError(f"unknown optical_study {optical_study!r}; "
                         f"use one of {tuple(STUDY_TO_MODE)}")
    celltype = getattr(g, "celltype", "hole")
    kw = dict(index=float(index), mode=STUDY_TO_MODE[optical_study],
              target_freq_thz=round(float(f_target) * 1e-12, 6), a=_um(g.a),
              xsect=resolve_xsect(xsect), wedge_angle=WEDGE_ANGLE_DEG)
    if celltype == "rib":
        kw.update(geometry="rib_block", t=_um(g.t), w=_um(g.w),
                  spine_width=_um(g.s), thickness=_um(g.th),
                  fillet_r1=_um(g.r1), fillet_r2=_um(g.r2))
    elif celltype == "hole":
        kw.update(geometry="holey", beam_width=_um(g.w), thickness=_um(g.t),
                  hole_x=_um(2.0 * g.hx), hole_y=_um(2.0 * g.hy))
    else:
        raise ValueError(f"unknown celltype {celltype!r}")
    for name, val in (overrides or {}).items():
        if name not in PARAM_OVERRIDES:
            raise ValueError(f"{name!r} cannot be overridden")
        kw[name] = val
    return kw


@dataclass
class KPlan:
    Nk: int                 # BandParams.Nk: simulation i has k = i/(2*Nk)
    indices: list           # simulations to run ("sim_<i>"); the last is X
    k_frac: list            # i/Nk, repo convention (X at 1)

    @property
    def x_only(self):
        return self.indices == [self.Nk]

    def k_module(self):
        return [i / (2.0 * self.Nk) for i in self.indices]


def k_plan(celltype, a, *, k_mode, n_k=None, f_target, rel_tol=0.30) -> KPlan:
    """Which Bloch simulations to run.

    The module samples k = i/(2*Nk), i = 0..Nk. A sweep keeps the last n_k of
    them, with Nk chosen so that they start at or just above
      rib : the first k of comsol_optics.rib_k_grid (light line of the lowest
            scored frequency),
      hole: k_frac = 0.5 (the zone-edge half of comsol_optics._opt_sweep),
    and end at X. The spacing 1/Nk is never coarser than COMSOL's grid.
    """
    if k_mode == "x_only":
        return KPlan(1, [1], [1.0])
    if k_mode != "sweep":
        raise ValueError(f"k_mode {k_mode!r}: use one of {K_MODES}")
    n_k = int(n_k)
    if n_k < 2:
        raise ValueError("a sweep needs n_k >= 2")
    if celltype == "rib":
        k_start = float(rib_k_grid(a, 2, f_target, rel_tol)[0])
    else:
        k_start = 0.5
    Nk = max(n_k - 1, int(np.ceil((n_k - 1) / (1.0 - k_start) - 1e-9)))
    indices = list(range(Nk - (n_k - 1), Nk + 1))
    return KPlan(Nk, indices, [i / Nk for i in indices])


# ---- scoring (numpy only) ----------------------------------------------------------

def bands_to_arrays(bands):
    """Module resonances -> (k_frac, freqs [n_k, n_max] Hz, NaN padded).

    bands: list of {"k": k in 2*pi/a, "freq": Hz array}, as BandResult.bands.
    Rows are sorted by k, frequencies ascending.
    """
    bands = sorted(bands, key=lambda d: float(d["k"]))
    rows = [np.sort(np.atleast_1d(np.asarray(d["freq"], dtype=float)))
            for d in bands]
    n_max = max([r.size for r in rows] + [1])
    freqs = np.full((len(rows), n_max), np.nan)
    for i, r in enumerate(rows):
        freqs[i, :r.size] = r
    return np.array([2.0 * float(d["k"]) for d in bands]), freqs


def guided_only(k_frac, freqs, a):
    """Copy of freqs with everything at or above the light line set to NaN."""
    k_frac = np.asarray(k_frac, dtype=float).ravel()
    freqs = np.array(freqs, dtype=float)
    with np.errstate(invalid="ignore"):
        freqs[~(freqs < (C0 * k_frac / (2.0 * a))[:, None])] = np.nan
    return freqs


def score_resonances(celltype, k_frac, freqs, a, *, f_target, rel_tol=0.30,
                     k_mode="sweep", window_hz=None, optical_tol_hz=None,
                     out=None):
    """Gap of Tidy3D resonances with the gap logic used for COMSOL.

    optical_tol_hz (rib scoring, score version 2): when given, the +/- rel_tol
    centre rule is NOT applied. A sweep takes every guided interval gap
    (guided_interval_gaps) and returns the one that maximizes
    gap * exp(-((f_c - f_target)/optical_tol_hz)^2); x_only returns the
    fundamental gap at X wherever it is centred. `out` (a dict) then receives
    "freq_penalty" of the returned gap and "candidates" (every gap, as
    comsol_optics.select_gap_by_penalized_size lists them). Without
    optical_tol_hz (hole cell) nothing changes.

    Returns (Gap, method):
      "guided_interval_gap"   rib sweep. Same function as the COMSOL rib path,
                              with coverage = the FDTD frequency window.
      "optical_gap_tracked"   hole sweep (rows with k_frac >= 0.5, as
                              _opt_sweep), on the guided resonances.
      "optical_gap_tracked_x" x_only, either cell type: the fundamental gap at
                              X between the two lowest guided resonances, with
                              optical_gap_tracked's +/- rel_tol centre rule. For
                              the rib this is NOT the interval gap: bands that
                              cross the gap away from X are not seen.
    """
    k_frac = np.asarray(k_frac, dtype=float).ravel()
    freqs = np.asarray(freqs, dtype=float)
    if freqs.ndim != 2 or freqs.shape[0] != k_frac.size or k_frac.size == 0:
        raise ValueError("freqs must be [n_k, n], one row per k_frac")
    if np.any(np.diff(k_frac) <= 0) or abs(k_frac[-1] - 1.0) > 1e-6:
        raise ValueError("k_frac must increase and end at X (k_frac = 1)")
    if celltype == "rib" and k_mode == "sweep":
        if window_hz is None:
            raise ValueError("the rib interval gap needs the FDTD window")
        if optical_tol_hz is None:
            gp = guided_interval_gap(k_frac, freqs, a, f_target=f_target,
                                     rel_tol=rel_tol, coverage=tuple(window_hz))
            return gp, "guided_interval_gap"
        gaps = guided_interval_gaps(k_frac, freqs, a, coverage=tuple(window_hz))
        gp, pen, cands = select_gap_by_penalized_size(gaps, f_target,
                                                      optical_tol_hz)
        if out is not None:
            out.update(freq_penalty=pen, candidates=cands)
        return gp, "guided_interval_gap"
    keep = k_frac >= 0.5 - 1e-9
    gp = optical_gap_tracked(guided_only(k_frac, freqs, a)[keep],
                             f_target=f_target,
                             rel_tol=(rel_tol if optical_tol_hz is None
                                      else float("inf")))
    if optical_tol_hz is not None and out is not None:
        gp, pen, cands = select_gap_by_penalized_size(
            [gp] if gp.found else [], f_target, optical_tol_hz)
        out.update(freq_penalty=pen, candidates=cands)
    return gp, ("optical_gap_tracked_x" if k_frac[keep].size == 1
                else "optical_gap_tracked")


def _num(x):
    x = float(x)
    return x if np.isfinite(x) else None


def _gap_to_dict(gp):
    return dict(found=bool(gp.found), normalized_gap=float(gp.normalized_gap),
                f_lower=float(gp.f_lower), f_upper=float(gp.f_upper),
                f_center=float(gp.f_center))


# ---- module loading and the credit budget ------------------------------------

_bs = None
_budget = None


def have_tidy3d() -> bool:
    try:
        return importlib.util.find_spec("tidy3d") is not None
    except (ImportError, ValueError):
        return False


def load_bs():
    """Import diamond_bandStructure_sim.py from its path (imports tidy3d)."""
    global _bs
    if _bs is None:
        if not have_tidy3d():
            raise Tidy3DNotRun(
                "not_installed",
                "tidy3d is not installed in this Python environment "
                "(pip install tidy3d, then `tidy3d configure` for the API key). "
                "Use optical: comsol, or install tidy3d.")
        spec = importlib.util.spec_from_file_location(
            "diamond_bandStructure_sim", BS_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _bs = module
    return _bs


def reset_session():
    """Forget the loaded module and the credits counted so far (tests)."""
    global _bs, _budget
    _bs, _budget = None, None


def session_credits() -> float:
    """Summed cost estimates of the batches started by this process."""
    return float(_budget.estimated) if _budget is not None else 0.0


def _ledger_path(opts):
    return os.path.join(opts.results_root(), "credit_ledger.jsonl")


def ledger_total(opts) -> float:
    """Summed cost estimates of every batch ever started into results_dir."""
    path, total = _ledger_path(opts), 0.0
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    total += float(json.loads(line).get("estimated_credits", 0.0))
    return total


def _ledger_append(opts, entry):
    os.makedirs(opts.results_root(), exist_ok=True)
    with open(_ledger_path(opts), "a") as fh:
        fh.write(json.dumps(entry) + "\n")


def _session_budget(bs, opts):
    """The one Budget of this process, capped by both credit limits."""
    global _budget
    if _budget is None:
        _budget = bs.Budget(None)
    cap = opts.max_credits
    if opts.max_credits_total is not None:
        left = max(0.0, opts.max_credits_total - ledger_total(opts))
        cap = min(cap, _budget.estimated + left)
    _budget.max_credits = cap
    return _budget


def _submit(bs, opts, sims, run_dir, label, budget, estimate_only,
            verbose=True, on_progress=None):
    """THE ONLY PLACE that reaches the Tidy3D cloud (bs.run_cloud)."""
    if opts.allow_cloud is not True:
        raise Tidy3DNotRun("cloud_not_allowed", "tidy3d.allow_cloud is not true")
    return bs.run_cloud(sims, run_dir, label, opts.folder_name, budget,
                        estimate_only=estimate_only, verbose=bool(verbose),
                        on_stage=on_progress)


# ---- console output: everything to stderr ----------------------------------------

def _td_log_module():
    """The module tidy3d.log, or None (not installed, or without this API).

    importlib is needed: the attribute `tidy3d.log` is the Logger object,
    which hides the submodule of the same name.
    """
    try:
        module = importlib.import_module("tidy3d.log")
        module.log.handlers.get("console")
        if callable(module.set_logging_console):
            return module
    except Exception:
        pass
    return None


@contextlib.contextmanager
def tidy3d_console(verbose=True):
    """Send all Tidy3D / module console output to stderr while the block runs.

    1. tidy3d.log.set_logging_console(stderr=True): Tidy3D's one rich Console
       (web progress bars and status, console.log lines, log messages; all of
       tidy3d.web takes it from tidy3d.log.get_logging_console()) writes to
       stderr. Its level is INFO when verbose, ERROR otherwise.
    2. contextlib.redirect_stdout(sys.stderr) for plain print() calls and for
       a Tidy3D whose console could not be switched. A rich console finds
       sys.stdout when it writes, so the redirection reaches it too; what it
       cannot reach is a stream object that was stored before the block or
       output written through the OS file descriptor. Live displays
       (progress bars, spinners) are animated only when stderr is a terminal;
       into a file or a pipe rich prints their lines without the animation.
    On exit the previous console handler object is put back, which restores
    the stream and the level. td.config is not touched.
    """
    module = _td_log_module()
    handlers, saved = None, None
    if module is not None:
        handlers = module.log.handlers
        saved = ("console" in handlers, handlers.get("console"))
        try:
            module.set_logging_console(stderr=True)
            module.set_logging_level("INFO" if verbose else "ERROR")
        except Exception:
            pass                        # redirect_stdout below still holds
    try:
        with contextlib.redirect_stdout(sys.stderr):
            yield
    finally:
        if saved is not None:
            try:
                if saved[0]:
                    handlers["console"] = saved[1]
                else:
                    handlers.pop("console", None)
            except Exception:
                pass


# ---- one evaluation -------------------------------------------------------------

def _write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


def _same_params(saved, current):
    for key, val in current.items():
        if key == "Nk":                 # the k points are checked separately
            continue
        old = saved.get(key)
        if isinstance(val, float) and isinstance(old, (int, float)):
            if abs(float(old) - val) > 1e-9 * max(1.0, abs(val)):
                return False
        elif old != val:
            return False
    return True


def _load_downloaded(bs, run_dir, params, plan):
    """Results already in run_dir for exactly this cell and these k, or None."""
    path = os.path.join(run_dir, "params.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        if not _same_params(json.load(fh), params):
            return None
    try:
        found = bs.load_run(run_dir)
    except FileNotFoundError:
        return None
    wanted = ["k=%.6f" % k for k in plan.k_module()]
    if any(key not in found for key in wanted):
        return None
    return {key: found[key] for key in wanted}


def _load_cached_analysis(run_dir, params, plan, window, eval_id):
    """Resonances of this exact cell saved by an earlier analysis, or None.

    For a cell whose .hdf5 files are gone (tidy3d.keep_hdf5 best / none).
    evaluation.json holds what the scoring needs: the k points and every
    resonance frequency per k; the frequency window is in its detail. It is
    used only if params.json matches the current parameters (geometry, index,
    numerics, resonance filters), the id, the k points and the window agree.
    Returns (k_frac array, freqs [n_k, n] Hz with NaN padding).
    """
    try:
        with open(os.path.join(run_dir, "params.json")) as fh:
            if not _same_params(json.load(fh), params):
                return None
        with open(os.path.join(run_dir, "evaluation.json")) as fh:
            ev = json.load(fh)
        if str(ev.get("id")) != str(eval_id):
            return None
        k_frac = np.asarray(ev["k_frac"], dtype=float)
        rows = ev["freqs_hz"]
        win = (ev.get("detail") or {}).get("window_THz")
        if (k_frac.shape != (len(plan.k_frac),) or len(rows) != k_frac.size
                or np.max(np.abs(k_frac - np.asarray(plan.k_frac))) > 1e-6
                or win is None
                or abs(win[0] * 1e12 - window[0]) > 1e6
                or abs(win[1] * 1e12 - window[1]) > 1e6):
            return None
        n = max([len(r) for r in rows] + [1])
        freqs = np.full((k_frac.size, n), np.nan)
        for i, row in enumerate(rows):
            vals = [np.nan if x is None else float(x) for x in row]
            freqs[i, :len(vals)] = vals
        return k_frac, freqs
    except (OSError, ValueError, KeyError, TypeError):
        return None


# ---- tidy3d.keep_hdf5: pruning of downloaded .hdf5 files -------------------------

def hdf5_files(run_dir):
    """Names of the .hdf5 files directly in run_dir (batch.hdf5 included)."""
    try:
        return sorted(n for n in os.listdir(run_dir)
                      if n.endswith(".hdf5")
                      and os.path.isfile(os.path.join(run_dir, n)))
    except OSError:
        return []


def _read_note(run_dir):
    try:
        with open(os.path.join(run_dir, NOTE_NAME)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def mark_managed(run_dir, policy):
    """Note in the record folder: its .hdf5 files were downloaded under a
    keep_hdf5 policy and may be pruned. Folders without it are never pruned."""
    note = _read_note(run_dir)
    note.update(managed=True, policy=policy)
    note.setdefault("events", [])
    _write_json(os.path.join(run_dir, NOTE_NAME), note)


def is_managed(run_dir) -> bool:
    return _read_note(run_dir).get("managed") is True


def prune_hdf5(run_dir, results_root, reason, log=None):
    """Delete the .hdf5 files of ONE record folder. Returns (count, bytes).

    Refuses, and deletes nothing, unless ALL of these hold:
      * run_dir is <results_root>/<hole|rib>/<record id> (exactly that depth);
      * the folder has the note of mark_managed (so records downloaded
        without a keep_hdf5 policy, or by anything else, are never touched);
      * evaluation.json and params.json are there (the analysis is cached, so
        the cell can be re-scored and plotted without the files).
    Only regular files ending in .hdf5 directly in the folder are removed.
    Every deletion is logged (stderr by default) and added to the note.
    """
    log = log or (lambda m: print(m, file=sys.stderr, flush=True))
    root = os.path.realpath(results_root)
    folder = os.path.realpath(run_dir)
    parent = os.path.dirname(folder)
    if (os.path.dirname(parent) != root
            or os.path.basename(parent) not in K_DEFAULTS
            or not os.path.isdir(folder)):
        log(f"[keep_hdf5] not pruned (outside {root}/<celltype>/<id>): {run_dir}")
        return 0, 0
    if not is_managed(folder):
        log(f"[keep_hdf5] not pruned (not downloaded under keep_hdf5): {folder}")
        return 0, 0
    if not all(os.path.isfile(os.path.join(folder, n))
               for n in ("evaluation.json", "params.json")):
        log(f"[keep_hdf5] not pruned (analysis not saved yet): {folder}")
        return 0, 0
    count, freed, names = 0, 0, []
    for name in hdf5_files(folder):
        path = os.path.join(folder, name)
        size = os.path.getsize(path)
        os.remove(path)
        count, freed = count + 1, freed + size
        names.append(name)
        log(f"[keep_hdf5] deleted {path} ({size / 1e6:.2f} MB)")
    if count:
        note = _read_note(folder)
        note.setdefault("events", []).append(dict(
            ts=time.time(), reason=reason, deleted=names, freed_bytes=freed))
        note.update(hdf5_deleted=True)
        _write_json(os.path.join(folder, NOTE_NAME), note)
        log(f"[keep_hdf5] {os.path.basename(folder)}: {count} file(s), "
            f"{freed / 1e6:.2f} MB freed ({reason}); params.json, plan.json, "
            "bands.csv and evaluation.json are kept")
    return count, freed


class Hdf5Keeper:
    """`keep_hdf5: best` for a loop: keep the .hdf5 of the best cell only.

    Give it every record after it has been SAVED (update). A successful
    Tidy3D record that beats the best so far becomes the best and the files
    of the previous best are deleted; any other Tidy3D record has its files
    deleted. Records that are not Tidy3D results, and folders that were not
    downloaded under a keep_hdf5 policy (prune_hdf5), are left alone. For
    `all` and `none` it does nothing (`none` is handled per cell by
    run_optical_tidy3d). Nothing is deleted before update is called, so an
    interrupted run only ever lacks files that were already superseded.
    seed() tells it the best of earlier sessions without deleting anything.
    """

    def __init__(self, options, log=None):
        opts = Tidy3DOptions.from_config(options)
        self.policy = opts.keep_hdf5
        self.root = opts.results_root()
        self.log = log
        self.best = None                # (score, run_dir, id)
        self.freed_bytes = 0

    @property
    def active(self):
        return self.policy == "best"

    @staticmethod
    def _entry(rec):
        """(score, run_dir, id) of a successful Tidy3D record, else None."""
        detail = rec.get("optical_tidy3d_detail") or {}
        if (rec.get("optical_engine") != ENGINE or rec.get("status") != "success"
                or not detail.get("run_dir") or rec.get("score") is None):
            return None
        return float(rec["score"]), detail["run_dir"], rec.get("id")

    def seed(self, records):
        for rec in records or []:
            entry = self._entry(rec)
            if entry and (self.best is None or entry[0] > self.best[0]):
                self.best = entry

    def _prune(self, run_dir, reason):
        self.freed_bytes += prune_hdf5(run_dir, self.root, reason, self.log)[1]

    def update(self, rec):
        if not self.active or rec.get("optical_engine") != ENGINE:
            return
        run_dir = (rec.get("optical_tidy3d_detail") or {}).get("run_dir")
        if not run_dir or rec.get("optical_tidy3d_stop"):
            return                      # nothing was run for this record
        entry = self._entry(rec)
        same = lambda a, b: os.path.realpath(a) == os.path.realpath(b)  # noqa: E731
        if entry and (self.best is None or entry[0] > self.best[0]):
            old, self.best = self.best, entry
            if old and not same(old[1], run_dir):
                if hdf5_files(run_dir):
                    self._prune(old[1], f"superseded as best by {entry[2]}")
                # else: the new best was re-scored from its cached analysis
                # and has no files; the previous best keeps its own
        elif self.best and same(self.best[1], run_dir):
            return                      # the best cell again: keep
        else:
            self._prune(run_dir, "not the best cell (keep_hdf5: best)")


def _job_size(sims):
    """Grid cells and time steps of one simulation (local, no cloud)."""
    try:
        sim = next(iter(sims.values()))
        return dict(grid_cells=float(sim.num_cells),
                    time_steps=int(sim.num_time_steps))
    except Exception:
        return {}


def expected_symmetry(optical_study="opt TE", xsect=None):
    """Tidy3D symmetry tuple (x, y, z) the module must use; +1 PMC, -1 PEC.

    y: -1 for 'opt TE' (PEC on y = 0), +1 for 'opt TM' (PMC on y = 0).
    z: the mid-thickness mirror, +1 for TE (in-plane E even), -1 for TM, and
       0 for the wedge, which has no such mirror.
    """
    sy, sz = {"opt TE": (-1, 1), "opt TM": (1, -1)}[optical_study]
    return (0, sy, sz if z_mirror_symmetric(xsect) else 0)


def run_optical_tidy3d(g, *, celltype, eval_id, f_target, options,
                       rel_tol=0.30, optical_study="opt TE", n_k_default=None,
                       xsect=None, verbose=True, on_progress=None,
                       optical_tol_hz=None):
    """Optical gap of one cell with Tidy3D. Enforce feasibility BEFORE calling.

    optical_tol_hz: rib scoring (see score_resonances). rel_tol then only
    places the first k of a rib sweep (k_plan); it is no centre rule. The
    record additions then include optical_freq_penalty and
    optical_gap_candidates.

    verbose: show Tidy3D's web progress and our step lines; False = errors
    only. on_progress: optional callable(message) for the step lines
    (building, uploading, estimated cost, running, loading, analyzing); it is
    not called when verbose is False. Nothing is ever written to stdout (see
    tidy3d_console). The gates and the return value are those of
    _run_optical_tidy3d.
    """
    def say(message):
        if verbose and on_progress is not None:
            try:
                on_progress(message)
            except Exception:           # reporting must never break a solve
                pass

    with tidy3d_console(verbose):
        return _run_optical_tidy3d(
            g, celltype=celltype, eval_id=eval_id, f_target=f_target,
            options=options, rel_tol=rel_tol, optical_study=optical_study,
            n_k_default=n_k_default, xsect=xsect, verbose=bool(verbose),
            say=say, optical_tol_hz=optical_tol_hz)


def _run_optical_tidy3d(g, *, celltype, eval_id, f_target, options, rel_tol,
                        optical_study, n_k_default, xsect, verbose, say,
                        optical_tol_hz=None):
    """The evaluation behind run_optical_tidy3d (same arguments).

    g: geometry.Geometry or RibGeometry (meters). eval_id: record id, used as
    the results folder <results_dir>/<celltype>/<eval_id>/ and as task label.
    options: the `tidy3d:` config section (dict) or a Tidy3DOptions.

    Order of the safety gates (nothing is uploaded before the last one):
      1. tidy3d installed                      else Tidy3DNotRun not_installed
      2. simulations built and validated locally; plan.json written
      3. mode dry_run                          -> Tidy3DNotRun dry_run
      4. downloaded results of this exact cell are re-analyzed (free)
      5. allow_cloud is true                   else cloud_not_allowed
      6. mode estimate: drafts uploaded, priced, deleted -> estimate_only
      7. max_credits set and not yet used up   else no_budget / budget_exceeded
      8. batch uploaded and priced; started only if the estimate fits the
         budget (Budget in diamond_bandStructure_sim.run_cloud), else the
         drafts are deleted -> budget_exceeded

    Returns dict(gap=Gap, k_frac, freqs_hz, record=dict of record additions).
    """
    opts = Tidy3DOptions.from_config(options)
    xsect = resolve_xsect(xsect)
    k_mode, n_k = opts.k_sampling(celltype, n_k_default)
    plan = k_plan(celltype, g.a, k_mode=k_mode, n_k=n_k, f_target=f_target,
                  rel_tol=rel_tol)
    kw = band_params_kwargs(g, index=opts.index_value(), f_target=f_target,
                            optical_study=optical_study, overrides=opts.params,
                            xsect=xsect)
    span = float(kw.get("freq_span_thz", _FREQ_SPAN_THZ))
    window = ((kw["target_freq_thz"] - span) * 1e12,
              (kw["target_freq_thz"] + span) * 1e12)
    run_dir = os.path.join(opts.results_root(), celltype, str(eval_id))
    detail = dict(
        engine=ENGINE, mode=opts.mode, k_mode=k_mode, n_sims=len(plan.indices),
        k_frac=[float(k) for k in plan.k_frac], Nk=plan.Nk,
        index=kw["index"], polarization=kw["mode"], optical_study=optical_study,
        window_THz=[window[0] * 1e-12, window[1] * 1e-12],
        band_params={k: v for k, v in kw.items()}, run_dir=run_dir,
        xsect=xsect, z_symmetry=z_mirror_symmetric(xsect),
        wedge_depths_m=(wedge_depths(g) if xsect == "wedge" else None),
        estimated_credits=0.0, session_credits=session_credits())

    try:
        bs = load_bs()                                          # gate 1
    except Tidy3DNotRun as err:
        err.detail = detail
        raise
    say("building %d simulation(s) (%s, %s, xsect %s)"
        % (len(plan.indices), celltype, k_mode, xsect))
    p = bs.BandParams(Nk=plan.Nk, **kw)
    cell = bs.build_unit_cell(p)                                # gate 2
    built = bs.build_simulations(p, x_only=plan.x_only, cell=cell)
    sims = {"sim_%d" % i: built["sim_%d" % i] for i in plan.indices}
    params = asdict(p)
    symmetry = tuple(int(x) for x in p.symmetry)
    if symmetry != expected_symmetry(optical_study, xsect):
        raise RuntimeError(
            f"diamond_bandStructure_sim uses symmetry {symmetry} for "
            f"{optical_study} / xsect {xsect}, expected "
            f"{expected_symmetry(optical_study, xsect)}")
    detail.update(symmetry=list(symmetry), **_job_size(sims))
    _write_json(os.path.join(run_dir, "plan.json"), detail)
    say("plan written: %s" % os.path.join(run_dir, "plan.json"))

    def stop(kind, message, **more):
        detail.update(more, session_credits=session_credits())
        return Tidy3DNotRun(kind, message, detail)

    if opts.mode == "dry_run":                                  # gate 3
        try:
            with open(os.path.join(run_dir, "grid_report.txt"), "w") as fh:
                fh.write(bs.grid_report(p, cell, sims) + "\n")
        except Exception:
            pass
        raise stop("dry_run", f"dry run: {len(sims)} simulation(s) built in "
                   f"{run_dir}, nothing uploaded")

    sim_datas = (_load_downloaded(bs, run_dir, params, plan)    # gate 4
                 if opts.reuse_results else None)
    reused = sim_datas is not None
    cached = None
    if reused:
        say("re-using %d downloaded result(s) in %s (no cloud)"
            % (len(sim_datas), run_dir))
    elif opts.reuse_results:
        # the .hdf5 files are gone (keep_hdf5) but the analysis was saved:
        # score again from the saved resonances instead of paying again
        cached = _load_cached_analysis(run_dir, params, plan, window, eval_id)
        if cached is not None:
            say("re-using the saved analysis in %s (no .hdf5, no cloud)" % run_dir)
    if not reused and cached is None:
        if opts.allow_cloud is not True:                        # gate 5
            raise stop("cloud_not_allowed",
                       f"Tidy3D cloud not allowed: {len(sims)} simulation(s) "
                       "would be uploaded and charged. Set `tidy3d: "
                       "{allow_cloud: true, max_credits: <budget>}` in the "
                       "config to opt in (mode: dry_run builds without upload).")
        budget = _session_budget(bs, opts) if opts.max_credits else None
        label = "%s_%s" % (celltype, eval_id)
        if opts.mode == "estimate":                             # gate 6
            budget = budget or bs.Budget(None)
            _submit(bs, opts, sims, run_dir, label, budget, True,
                    verbose=verbose, on_progress=say)
            raise stop("estimate_only", "estimate only: %d simulation(s) would "
                       "cost at most %.3f FlexCredits; nothing was started"
                       % (len(sims), budget.last_estimate),
                       estimated_credits=float(budget.last_estimate))
        if budget is None:                                      # gate 7
            raise stop("no_budget", "tidy3d.max_credits is not set: refusing "
                       "to start paid simulations without a credit budget")
        if budget.estimated >= budget.max_credits:
            raise stop("budget_exceeded", "Tidy3D credit budget used up: "
                       "%.3f of %.3f FlexCredits; nothing was uploaded"
                       % (budget.estimated, budget.max_credits))
        before = budget.estimated
        entry = dict(ts=time.time(), id=str(eval_id), celltype=celltype,
                     n_sims=len(sims), run_dir=run_dir)
        budget.last_estimate = 0.0
        try:
            sim_datas = _submit(bs, opts, sims, run_dir, label, budget, False,
                                verbose=verbose, on_progress=say)
        except Exception as err:                                # gate 8
            cost = budget.estimated - before
            if cost > 0:        # started, then failed: the credits may be gone
                _ledger_append(opts, dict(entry, estimated_credits=cost,
                                          status="failed_after_start"))
                raise
            if (budget.last_estimate > 0
                    and before + budget.last_estimate > budget.max_credits):
                raise stop("budget_exceeded", "Tidy3D credit budget would be "
                           "exceeded: %.3f used + %.3f for this cell > %.3f "
                           "FlexCredits; the drafts were deleted, nothing was "
                           "started" % (before, budget.last_estimate,
                                        budget.max_credits),
                           estimated_credits=0.0,
                           refused_estimate=float(budget.last_estimate)
                           ) from err
            raise
        cost = budget.estimated - before
        _ledger_append(opts, dict(entry, estimated_credits=cost, status="run"))
        detail.update(estimated_credits=float(cost))
        _write_json(os.path.join(run_dir, "params.json"), params)
        if opts.keep_hdf5 != "all":     # these files may be pruned later
            mark_managed(run_dir, opts.keep_hdf5)

    if cached is not None:
        k_frac, freqs = cached
        result = None
        at_x = freqs[-1][np.isfinite(freqs[-1])]
        edges = bs.band_edges(at_x, p)  # (f_lower, f_upper, gap, mid, fitness)
    else:
        say("analyzing resonances of %d simulation(s)" % len(sim_datas))
        result = bs.analyze_bands(sim_datas, p)
        k_frac, freqs = bands_to_arrays(result.bands)
        if len(k_frac) != len(plan.k_frac) or np.max(
                np.abs(k_frac - np.asarray(plan.k_frac))) > 1e-6:
            raise RuntimeError(f"Tidy3D returned k_frac {k_frac.tolist()}, "
                               f"expected {plan.k_frac}")
        edges = (result.f_lower, result.f_upper, result.bandGap, result.midGap,
                 result.fitness)
    picked = {}
    gp, method = score_resonances(celltype, k_frac, freqs, g.a,
                                  f_target=f_target, rel_tol=rel_tol,
                                  k_mode=k_mode, window_hz=window,
                                  optical_tol_hz=optical_tol_hz, out=picked)
    detail.update(
        reused_download=bool(reused), gap_method=method,
        # where the resonances came from: a cloud run of this call, .hdf5
        # files already on disk, or the saved analysis of a pruned cell
        result_source=("cached_analysis" if cached is not None
                       else "hdf5" if reused else "cloud"),
        keep_hdf5=opts.keep_hdf5,
        session_credits=session_credits(),
        min_resolvable_Hz=_num(resolvable_step_hz(k_frac, g.a)),
        light_line_edge_Hz=float(C0 / (2.0 * g.a)),
        # the module's own figure of merit: two lowest guided resonances at X
        x_point=dict(f_lower_THz=_num(edges[0]),
                     f_upper_THz=_num(edges[1]),
                     bandGap_THz=_num(edges[2]),
                     midGap_THz=_num(edges[3]),
                     normalized_gap=(_num(edges[2] / edges[3])
                                     if np.isfinite(edges[3]) else None),
                     fitness=_num(edges[4])))
    freqs_list = [[_num(x) for x in row] for row in freqs]
    record = dict(optical_tidy3d_detail=detail,
                  optical_gap_detail=_gap_to_dict(gp),
                  _k_o=[float(x) for x in k_frac], _freqs_o=freqs_list)
    if optical_tol_hz is not None:      # rib: Gaussian-weighted gap choice
        record.update(optical_freq_penalty=float(picked.get("freq_penalty", 0.0)),
                      optical_gap_candidates=picked.get("candidates", []))
    if cached is None:
        # (a cell re-scored from its saved analysis keeps the files it has)
        try:
            bs.save_bands_csv(os.path.join(run_dir, "bands.csv"), result, p)
        except Exception:
            pass
        _write_json(os.path.join(run_dir, "evaluation.json"),
                    dict(id=str(eval_id), celltype=celltype, params_m=g.as_dict(),
                         gap=_gap_to_dict(gp), k_frac=record["_k_o"],
                         freqs_hz=freqs_list, detail=detail))
        if opts.keep_hdf5 == "none":    # analysis saved: the files can go
            prune_hdf5(run_dir, opts.results_root(), "keep_hdf5: none")
    return dict(gap=gp, k_frac=k_frac, freqs_hz=freqs, record=record)


def not_run_record(rec, err):
    """Record update for a Tidy3DNotRun: status "optical_not_run", nothing run.

    The id gets the stop kind appended, so this record can never replace a
    paid result of the same cell in the database (INSERT OR REPLACE on id).
    """
    kind = err.tidy3d_stop
    rec.update(id="%s_%s" % (rec["id"], kind), status="optical_not_run",
               error=str(err), optical_tidy3d_stop=kind,
               optical_tidy3d_detail=err.detail,
               optical_gap=0.0, mechanical_gap=0.0, score=-5.0,
               optical_center_frequency=0.0, mechanical_center_frequency=0.0)
    return rec
