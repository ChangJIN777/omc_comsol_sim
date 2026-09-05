"""Objective: 2D boomerang geometry -> (mechanical gap, [future] optical gap, score).

Adapted from python-scripts/src/objective.py (1D nanobeam project), with two
changes specific to this 2D phononic-crystal project:

1. Mechanical gap extraction is keyed on Z-PARITY (evenz/oddz), not a
   breathing-mode fy fraction -- see targets_2d.yaml's `gap_mode` and
   `_mechanical_gap` below. Both parities are always solved (the mechanical
   backend is cheap relative to optical) so "symmetry" vs "complete" gap_mode
   can be selected/re-scored after the fact without re-solving.

   `gap_mode` defaults to **"complete"**: `mechanical_gap` and `score` refer to
   a gap clear of BOTH parity families, which is what a phononic shield or a
   clamping-loss-free cavity actually requires. "symmetry" remains supported for
   the case where only one parity couples to the transducer, and the per-family
   diagnostics (`mechanical_gap_{evenz,oddz}` and their edges) are recorded in
   BOTH modes -- in complete mode they are what explains a small complete gap,
   so do not treat them as symmetry-mode-only.
2. Optical is NOT YET IMPLEMENTED (src/optical_comsol_2d.py is a documented
   stub). `optical_backend="none"` (the default) simply excludes it from the
   score -- no solve attempted. Set `require_opt=True` and
   `optical_backend="comsol"` once that module is implemented; until then
   doing so will record a clean "optical_failed" status rather than crashing
   the loop (see evaluate_candidate's try/except below), so a run config that
   flips it on early won't silently corrupt results.

Score (mirrors python-scripts/src/objective.py's formula, extended with the
require_opt/require_mech toggle from python-scripts/scripts/run_opt_comsol.py
so a mechanical-only run -- the only backend implemented so far -- doesn't
need a fake optical term):
  S = (G_m if require_mech else 0)
      - lm  * max(0, gmin_m - G_m)^2                [if require_mech]
      - lmf * ((f_m,c - f_m,target)/f_m,target)^2    [if require_mech]
      + (G_o if require_opt else 0)
      - lo  * max(0, gmin_o - G_o)^2                 [if require_opt]
      - lof * ((f_o,c - f_o,target)/f_o,target)^2    [if require_opt]
Infeasible geometry -> large negative score, no solve.

"Is a complete gap the OVERLAP of the evenz and oddz gaps?" -- not quite, and
the difference matters. See the EQUIVALENT FORMULATIONS section of
`_combined_bands`'s docstring below for the precise statement (short version:
stack-then-search returns a strict superset of the pairwise overlaps, and the
extra windows are real ones the overlap recipe cannot express).
"""
from __future__ import annotations

import os
import hashlib
import json

import numpy as np
import yaml

from geometry2d import u_to_geometry, check_feasibility, load_bounds
from bandgap import Gap, gap_near_frequency, largest_gap

_CFG = os.path.join(os.path.dirname(__file__), "..", "configs")
C0 = 299_792_458.0

# Where per-candidate band-structure .npz files go. results/ is git-ignored, so
# these stay local; override with OMC2D_BANDS_DIR (mirrors database.py's
# OMC2D_DB_PATH convention).
_BANDS_DIR = os.environ.get(
    "OMC2D_BANDS_DIR",
    os.path.join(os.path.dirname(__file__), "..", "results", "bands"))

# z-parity families, in TIE-BREAK ORDER. Both are always solved (see
# acoustic_comsol_2d.run_mechanical_comsol_2d's `parities` default), so a
# symmetry-mode score always has a loser worth recording.
PARITY_ORDER = ("evenz", "oddz")


def _load(name):
    with open(os.path.join(_CFG, name)) as fh:
        return yaml.safe_load(fh)


def _hash_u(u, backends):
    key = json.dumps([round(float(x), 6) for x in u] + list(backends))
    return "cand2d_" + hashlib.sha1(key.encode()).hexdigest()[:12]


def _combined_bands(freqs_evenz, freqs_oddz):
    """Stack both z-parity families, truncated to where the union is COMPLETE.

    Each family is solved for its own n_bands lowest modes, so at a given
    k-point the union of the two families is only complete up to
    min(evenz[k, -1], oddz[k, -1]) -- above that, at least one family has run
    out of computed modes. A gap must hold at EVERY k, so the usable ceiling is

        ceiling = min over k of min(evenz[k, -1], oddz[k, -1])

    Above it, any apparent gap is a truncation artifact rather than a real
    absence of states -- e.g. if oddz is only resolved to 13 GHz at one k-point,
    an oddz gap at 16-19 GHz looks "complete" because nothing was ever computed
    there to contradict it.

    NOTE: reducing with max-over-k instead of min-over-k is wrong and was the
    original bug here -- max_k is the frequency below which a family SOMETIMES
    has all its modes; a gap needs ALWAYS. Fuzzing 40k random spectra against a
    known ground truth: max-over-k admitted 419 false-positive complete gaps,
    min-over-k admitted 0.

    Entries above the ceiling are NaN'd and the trailing columns dropped.
    np.sort puts NaN at the end of each row, so the `keep` mask is a contiguous
    True-prefix -- that placement is what makes slicing to n_keep valid.

    EQUIVALENT FORMULATIONS -- "isn't a complete gap just the OVERLAP of the
    two families' gaps?" This question comes up every time; the answer is that
    stack-then-search (what this function feeds) returns a strict SUPERSET of
    the pairwise-intersection recipe, and the superset is the correct one.

    Write e_1(k) <= ... <= e_nE(k) for the evenz states at k and o_1..o_nO for
    oddz. `bandgap.all_gaps` reports a gap of a single family between its bands
    i and i+1 as (max_k e_i, min_k e_{i+1}), nonempty iff that is a real
    interval. Because e_i(k) <= e_{i+1}(k) pointwise, both endpoints of the
    band RANGES B_i = [min_k e_i, max_k e_i] are nondecreasing in i, so the
    ranges are ordered and the holes in their union are exactly those gaps.
    Hence:

        {f : no evenz state at any k}  =  [0, min_k e_1)                    (i)
                                       U  the evenz gaps                   (ii)
                                       U  (max_k e_nE, inf)                (iii)

    and likewise for oddz. A complete gap is a frequency free of BOTH families,
    i.e. the intersection of two sets of that shape. Expanding gives four kinds
    of term, only ONE of which is a pairwise intersection of named gaps:

        (evenz gap_i) ∩ (oddz gap_j)        <- the "overlap" recipe
        (evenz gap_i) ∩ [0, min_k o_1)      <- free because it is BELOW oddz
        [0, min_k e_1) ∩ (oddz gap_j)       <- free because it is BELOW evenz
        [0, min_k e_1) ∩ [0, min_k o_1)     <- below both; not reportable as a
                                               gap (no lower band brackets it)

    Term (iii) is the untrustworthy region the `ceiling` above exists to
    exclude: one family simply ran out of computed modes there.

    So the pairwise-overlap recipe MISSES any window that is clear of a family
    because it sits below that family's lowest band rather than inside one of
    that family's gaps. That is not a corner case: if one family has no gap at
    all, the overlap recipe returns nothing while a real complete gap can still
    exist below that family's first band. Measured on 4 000 random band pairs,
    84% produced at least one complete gap that is not any pairwise
    intersection.

    The containment direction holds exactly, with one qualifier. Let c_m(k) be
    the m-th smallest state of the union at k and N = n_keep. Then:

        every pairwise intersection (evenz gap_i ∩ oddz gap_j), restricted to
        below min_k c_N, is contained in exactly ONE gap returned by
        all_gaps(_combined_bands(...)[0]).

    Proof sketch: for f free of both families, m(f) = #{states < f} is constant
    in k (a change would require a state to cross f, making f an eigenvalue),
    so f sits in the single combined gap (max_k c_m, min_k c_{m+1}); it is
    retained iff m <= N-1, i.e. iff f < min_k c_N. Verified by
    test_complete_gap_contains_every_pairwise_intersection over 30 000 random
    band pairs (48 283 intersections, 0 misses, 0 multi-coverage), spanning flat
    and dispersive bands and unequal nE/nO.

    The qualifier is NOT the ceiling but min_k c_N, which can be well below it:
    a band straddling the ceiling is dropped by the `keep` mask, taking with it
    any gap whose upper edge was that band's bottom. Concretely, with evenz
    bands ranging 4-5, 6-6.4, 12-22 and oddz 4.5-5, 6.8-7.2, 24-30, the
    intersection 7.2->12 GHz is real but is NOT reported: the ceiling is 12.0
    (evenz's top band bottoms there), min_k c_N is 6.8, so nothing above
    6.8 GHz survives. Give both families a band above 12 and the 7.2->12 gap
    appears as expected. This is one-sided conservatism -- it can hide a real
    complete gap, never invent one -- and the fix is more bands (`neigs`), not
    a different reduction.

    Requires both arrays to have the same number of k-points; nE and nO may
    differ freely (the ceiling uses each family's own top band).
    """
    evenz = np.asarray(freqs_evenz, dtype=float)
    oddz = np.asarray(freqs_oddz, dtype=float)
    ceiling = float(np.nanmin(np.minimum(evenz[:, -1], oddz[:, -1])))
    combined = np.sort(np.hstack([evenz, oddz]), axis=1)
    combined = np.where(combined <= ceiling, combined, np.nan)
    keep = np.all(np.isfinite(combined), axis=0)
    n_keep = int(np.count_nonzero(keep))
    return combined[:, :n_keep], ceiling


def best_parity(gaps):
    """Pick the winning family from a {parity: Gap} mapping -> (parity, Gap).

    Ties break toward the EARLIER entry in PARITY_ORDER, i.e. evenz. That is a
    deliberate choice, not a side effect of a `>=`: an exact tie is entirely
    plausible for a structure close to z-symmetric, and a parity label that
    flipped between runs of the same design would be a miserable thing to
    debug. Shared with scripts/plot_bands_2d.py so a figure's legend and a
    result record can never disagree about which family won.

    Accepts a subset of PARITY_ORDER (the plotter is often given one family).
    """
    ordered = [(p, gaps[p]) for p in PARITY_ORDER if p in gaps]
    if not ordered:
        raise ValueError(f"no parity to choose from; got {sorted(gaps)}, "
                         f"expected some of {list(PARITY_ORDER)}")
    best_p, best_g = ordered[0]
    for p, g in ordered[1:]:
        if g.normalized_gap > best_g.normalized_gap:   # strict: see docstring
            best_p, best_g = p, g
    return best_p, best_g


def _mechanical_gap(freqs_evenz, freqs_oddz, f_target, gap_mode):
    """Gap selection keyed on z-parity, per targets_2d.yaml's gap_mode.

    Returns (Gap, info). The Gap is the SCORED one; `info` is merged straight
    into the result record and carries:

      mech_parity       "evenz" / "oddz" (symmetry mode) or "complete".
      mechanical_gap_{evenz,oddz}, mechanical_center_frequency_{evenz,oddz},
      mechanical_gap_lower_frequency_{evenz,oddz},
      mechanical_gap_upper_frequency_{evenz,oddz}
                        each family's own gap, ALWAYS populated in both modes.
      mechanical_gap_lower_frequency, mechanical_gap_upper_frequency
                        edges of the SCORED gap (the one `mechanical_gap` and
                        `mechanical_center_frequency` describe).
      truncation_ceiling_hz, n_bands_usable
                        complete mode only (None in symmetry mode) -- what
                        tells "no gap exists" apart from "not enough bands
                        were solved to say".

    All frequencies are Hz. Gap EDGES are stored, not just (centre, normalized
    width), even though the edges are algebraically recoverable as
    centre*(1 -+ G/2) -- that identity only holds for THIS project's
    normalization (G = (hi-lo)/mean(hi,lo), see bandgap.Gap), and "is
    Delta f / f_lower" is a common enough alternative convention that a reader
    reconstructing edges from the record could be quietly wrong. Storing them
    also makes the natural question -- "do the two families' gaps overlap, and
    where?" -- answerable by reading two numbers per family instead of doing
    algebra. A family with no gap reports 0.0 for all four of its fields.

    'symmetry': best single-family gap (evenz or oddz, whichever is larger
                near f_target) -- a quasi/symmetry-restricted gap.
    'complete': a gap with NO band of EITHER parity -- both families stacked,
                truncated to where the union is complete (see _combined_bands).

    The per-family numbers are recorded in BOTH modes because they are pure
    post-processing on data already in hand, and because the record otherwise
    collapses to one scalar with no way to tell which family earned it. The
    premise of gap_mode 'symmetry' is that a family-restricted gap is
    acceptable *because only one parity couples to the transducer* -- so
    without the label a scored candidate cannot be judged usable at all. In
    'complete' mode they answer "what would each family have given alone?",
    which is exactly the question when a complete gap comes back empty.

    Note the per-family gaps are computed on the UNTRUNCATED family arrays.
    That is correct and not an inconsistency with _combined_bands: within one
    family the n_bands lowest modes are all known, so any gap between
    consecutive bands of that family is real. The truncation ceiling exists
    only because *stacking* two families leaves the interleaving unknown above
    the point where either one runs out.
    """
    per = {"evenz": gap_near_frequency(freqs_evenz, f_target),
           "oddz": gap_near_frequency(freqs_oddz, f_target)}
    info = {}
    for p in PARITY_ORDER:
        # Gap.empty() carries normalized_gap=0.0 and f_lower/f_upper/f_center=0,
        # so a family with no gap reports 0.0 rather than NaN or a missing key.
        info[f"mechanical_gap_{p}"] = float(per[p].normalized_gap)
        info[f"mechanical_center_frequency_{p}"] = float(per[p].f_center)
        info[f"mechanical_gap_lower_frequency_{p}"] = float(per[p].f_lower)
        info[f"mechanical_gap_upper_frequency_{p}"] = float(per[p].f_upper)

    def _with_scored_edges(gap):
        info["mechanical_gap_lower_frequency"] = float(gap.f_lower)
        info["mechanical_gap_upper_frequency"] = float(gap.f_upper)
        return gap, info

    if gap_mode == "complete":
        combined, ceiling = _combined_bands(freqs_evenz, freqs_oddz)
        info.update(mech_parity="complete",
                    truncation_ceiling_hz=float(ceiling),
                    n_bands_usable=int(combined.shape[1]))
        if combined.shape[1] < 2:
            return _with_scored_edges(Gap.empty())  # nothing below the ceiling
        return _with_scored_edges(gap_near_frequency(combined, f_target))

    best_p, best = best_parity(per)
    info.update(mech_parity=best_p, truncation_ceiling_hz=None,
                n_bands_usable=None)
    return _with_scored_edges(best)


def evaluate_candidate(u, *, mech_backend="comsol", optical_backend="none",
                       require_mech=True, require_opt=False,
                       n_bands_mech=None, cache=None, bounds=None,
                       targets=None, on_stage=None,
                       save_bands=True, bands_dir=None):
    """Evaluate one normalized design vector u. Returns a result record dict.

    on_stage, if given, is called as on_stage(name, event, info=None) at each
    stage boundary (name in {"feasibility","mechanical","optical","score"},
    event in {"start","done","skip","fail"}) -- see scripts/run_one_2d.py for
    a stderr-reporting implementation. Never raises because of this hook.

    save_bands (default True) dumps the full [n_k, n_bands] band structure of
    every successful mechanical solve to `bands_dir`/<id>.npz and records the
    path as `bands_npz` in the result record. Defaulting to ON is deliberate:
    a mechanical solve is ~54 eigensolves, and the record otherwise keeps only
    two scalars (G_m, f_center), so the bands were the most expensive thing
    being thrown away. Keeping them also makes this module's "re-score
    'symmetry' vs 'complete' without re-solving" claim true instead of
    aspirational, and gives scripts/plot_bands_2d.py something to plot.
    Set save_bands=False for throwaway runs. Failing to write the .npz never
    fails the candidate -- it is recorded in `bands_error` and the solve stands.

    Mechanical fields in the record. `mechanical_gap` /
    `mechanical_center_frequency` are the SCORED gap and keep that meaning;
    alongside them, `mech_parity` names which family earned it ("evenz",
    "oddz", or "complete"), `mechanical_gap_lower_frequency` /
    `..._upper_frequency` give its edges, and the same four quantities suffixed
    `_evenz` / `_oddz` give both families' own gaps in either mode -- so
    "do the two families' gaps overlap, and where?" is answerable straight from
    a record. See _mechanical_gap for why the losing family is kept, and its
    EQUIVALENT FORMULATIONS note for how an overlap relates to a complete gap.
    """
    def _stage(name, event, info=None):
        if on_stage is not None:
            on_stage(name, event, info)

    bounds = bounds or load_bounds()
    targets = targets or _load("targets_2d.yaml")
    cid = _hash_u(u, (mech_backend, optical_backend))

    if cache is not None and cid in cache:
        return cache[cid]

    _stage("feasibility", "start")
    g = u_to_geometry(u, bounds)
    ok, reasons = check_feasibility(g, bounds)
    rec = dict(id=cid, u=[float(x) for x in u], params=g.as_dict(),
               mech_backend=mech_backend, optical_backend=optical_backend,
               status="pending", reasons=reasons)

    if not ok:
        _stage("feasibility", "fail", reasons)
        rec.update(status="infeasible", mechanical_gap=0.0, optical_gap=0.0,
                   score=-10.0, mechanical_center_frequency=0.0,
                   optical_center_frequency=0.0)
        if cache is not None:
            cache[cid] = rec
        return rec
    _stage("feasibility", "done")

    f_m_target = targets["mechanical"]["target_frequency_GHz"] * 1e9
    f_o_target = C0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    gap_mode = targets["mechanical"].get("gap_mode", "complete")
    n_bands_mech = n_bands_mech or targets["mechanical"].get("n_bands", 10)

    # ---- mechanical (always solved when require_mech, both z-parities) ----
    _stage("mechanical", "start" if require_mech else "skip", mech_backend)
    G_m, f_m_c = 0.0, 0.0
    mech_summary = None
    if require_mech:
        try:
            if mech_backend == "comsol":
                from acoustic_comsol_2d import run_mechanical_comsol_2d
                from acoustic_comsol_2d import save_bands as _write_bands
                d = run_mechanical_comsol_2d(g, n_bands=n_bands_mech)
                if save_bands:
                    dest = os.path.join(bands_dir or _BANDS_DIR, f"{cid}.npz")
                    try:
                        rec["bands_npz"] = os.path.abspath(_write_bands(d, dest))
                    except Exception as exc:      # disk full, bad path, ...
                        rec["bands_npz"] = None   # never lose a good solve to a
                        rec["bands_error"] = str(exc)   # failed side-effect
                gp, gap_info = _mechanical_gap(d["freqs_evenz"], d["freqs_oddz"],
                                               f_m_target, gap_mode)
                G_m, f_m_c = gp.normalized_gap, gp.f_center
                rec.update(gap_mode=gap_mode, **gap_info)
                mech_summary = (
                    f"{gap_info['mech_parity']}: "
                    f"{G_m*100:.1f}% @ {f_m_c/1e9:.2f} GHz | "
                    f"evenz {gap_info['mechanical_gap_evenz']*100:.1f}% "
                    f"oddz {gap_info['mechanical_gap_oddz']*100:.1f}%")
            else:
                raise ValueError(f"unknown mech_backend {mech_backend!r}")
        except Exception as e:  # solver failed -> record, don't crash loop
            _stage("mechanical", "fail", str(e))
            rec.update(status="mech_failed", error=str(e),
                       mechanical_gap=0.0, optical_gap=0.0, score=-5.0,
                       mechanical_center_frequency=0.0,
                       optical_center_frequency=0.0)
            if cache is not None:
                cache[cid] = rec
            return rec
        _stage("mechanical", "done", mech_summary)

    # ---- optical (not yet implemented; "none" = excluded from score) ------
    _stage("optical", "start" if require_opt else "skip", optical_backend)
    G_o, f_o_c = 0.0, f_o_target
    if require_opt:
        try:
            if optical_backend == "comsol":
                from optical_comsol_2d import run_optical_comsol_2d
                d = run_optical_comsol_2d(g)
                gp = gap_near_frequency(d["freqs_hz"], f_o_target)
                G_o, f_o_c = gp.normalized_gap, gp.f_center
            elif optical_backend == "none":
                raise ValueError(
                    "require_opt=True but optical_backend='none' -- set "
                    "optical_backend='comsol' (not yet implemented, see "
                    "src/optical_comsol_2d.py) or require_opt=False")
            else:
                raise ValueError(f"unknown optical_backend {optical_backend!r}")
        except Exception as e:
            _stage("optical", "fail", str(e))
            rec.update(status="optical_failed", error=str(e),
                       mechanical_gap=G_m, optical_gap=0.0,
                       score=-2.0 if not require_mech else None,
                       mechanical_center_frequency=f_m_c,
                       optical_center_frequency=0.0)
            if rec["score"] is None:
                # mechanical already succeeded -- keep it, just zero the score
                # contribution from optical rather than discarding the run.
                rec["score"] = _score(G_m, 0.0, f_m_c, f_o_target, f_m_target,
                                      f_o_target, targets,
                                      require_mech=require_mech,
                                      require_opt=False)
            if cache is not None:
                cache[cid] = rec
            return rec
        _stage("optical", "done")

    _stage("score", "start")
    score = _score(G_m, G_o, f_m_c, f_o_c, f_m_target, f_o_target, targets,
                   require_mech=require_mech, require_opt=require_opt)
    rec.update(status="success",
               mechanical_gap=float(G_m), optical_gap=float(G_o),
               mechanical_center_frequency=float(f_m_c),
               optical_center_frequency=float(f_o_c),
               require_mech=bool(require_mech), require_opt=bool(require_opt),
               score=float(score))
    _stage("score", "done")
    if cache is not None:
        cache[cid] = rec
    return rec


def _score(G_m, G_o, f_m_c, f_o_c, f_m_t, f_o_t, targets, *,
          require_mech, require_opt):
    s = targets["scoring"]
    base = 0.0
    pen = 0.0
    if require_mech:
        gmin_m = targets["mechanical"]["min_gap"]
        base += G_m
        pen += s["lambda_mech"] * max(0.0, gmin_m - G_m) ** 2
        pen += s["lambda_mech_freq"] * ((f_m_c - f_m_t) / f_m_t) ** 2
    if require_opt:
        gmin_o = targets["optical"]["min_gap"]
        base += G_o
        pen += s["lambda_optical"] * max(0.0, gmin_o - G_o) ** 2
        pen += s["lambda_optical_freq"] * ((f_o_c - f_o_t) / f_o_t) ** 2
    return base - pen


if __name__ == "__main__":
    # Smoke test with the surrogate-free path (feasibility only, no COMSOL):
    # exercises geometry + feasibility + score wiring without a solver.
    from geometry2d import u_to_geometry
    g = u_to_geometry([0.5, 0.5, 0.5, 0.3, 0.3, 0.5])
    ok, reasons = check_feasibility(g)
    print("geometry (nm):", {k: round(v * 1e9, 1) for k, v in g.as_dict().items()})
    print("feasible:", ok, reasons)
    print("(run acoustic_comsol_2d.py directly, with a COMSOL server up, "
          "for an actual band-structure solve)")
