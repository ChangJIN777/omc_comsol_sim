"""Bandgap extraction from band-structure data.

Backend-agnostic: every solver (MPB optical, COMSOL optical/mechanical, the
numpy surrogate) returns frequencies as an array of shape [n_k, n_bands] (Hz for
mechanical, Hz or normalized for optical). These helpers find gaps.

Two notions of gap:
  - complete gap : a frequency interval with NO band of ANY symmetry. Use when
                   `bands` already contains all bands.
  - symmetry gap : pass only the bands of one symmetry family (filter upstream
                   using mode parity), then call the same function. This is the
                   practical route for 1D nanobeams, whose mechanical spectrum
                   rarely has a complete gap but often has a clean gap for the
                   even-even "breathing" family.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class Gap:
    f_lower: float        # top of band below the gap
    f_upper: float        # bottom of band above the gap
    f_center: float
    gap_size: float       # absolute (same units as input)
    normalized_gap: float # (f_upper - f_lower) / f_center
    lower_band: int
    found: bool

    @staticmethod
    def empty() -> "Gap":
        return Gap(0, 0, 0, 0.0, 0.0, -1, False)


def largest_gap(freqs: np.ndarray, f_window: Optional[tuple] = None,
                nan_safe: bool = False) -> Gap:
    """Find the largest normalized gap between consecutive bands.

    Parameters
    ----------
    freqs : array [n_k, n_bands], frequencies along the Bloch path.
    f_window : optional (f_lo, f_hi). If given, only consider gaps whose center
        falls inside the window (keeps the gap in the target spectral range).
    nan_safe : if True, use nanmax/nanmin so NaN entries are ignored.  Use this
        when `freqs` has NaN for modes outside a symmetry family (e.g. only
        breathing modes kept, others set to NaN).  Both band edges must be
        finite for a gap to be accepted.

    A gap between band i and band i+1 exists iff
        max_k freqs[:, i]  <  min_k freqs[:, i+1].
    """
    freqs = np.asarray(freqs, dtype=float)
    if freqs.ndim != 2:
        raise ValueError("freqs must be [n_k, n_bands]")
    freqs = np.sort(freqs, axis=1)  # ensure ascending bands at each k; NaN → end
    n_bands = freqs.shape[1]

    _max = np.nanmax if nan_safe else np.max
    _min = np.nanmin if nan_safe else np.min

    best = Gap.empty()
    for i in range(n_bands - 1):
        top_lower = float(_max(freqs[:, i]))
        bot_upper = float(_min(freqs[:, i + 1]))
        if not (np.isfinite(top_lower) and np.isfinite(bot_upper)):
            continue  # one or both bands absent (all NaN)
        if bot_upper <= top_lower:
            continue  # bands overlap -> no gap here
        center = 0.5 * (top_lower + bot_upper)
        if f_window is not None and not (f_window[0] <= center <= f_window[1]):
            continue
        ng = (bot_upper - top_lower) / center
        if ng > best.normalized_gap:
            best = Gap(top_lower, bot_upper, center, bot_upper - top_lower,
                       ng, i, True)
    return best


def gap_near_frequency(freqs: np.ndarray, f_target: float,
                       rel_tol: float = 0.5) -> Gap:
    """Largest gap whose center lies within +/- rel_tol of f_target."""
    window = (f_target * (1 - rel_tol), f_target * (1 + rel_tol))
    return largest_gap(freqs, f_window=window)


def all_gaps(freqs: np.ndarray, f_window: Optional[tuple] = None,
            nan_safe: bool = False) -> list:
    """Find every gap between consecutive bands (not just the largest).

    Same gap criterion as `largest_gap` (max_k freqs[:,i] < min_k freqs[:,i+1]),
    but returns a Gap for EVERY qualifying (i, i+1) pair, sorted by
    normalized_gap descending. Use this to show/inspect all candidate gaps
    (e.g. to check whether the one picked by `largest_gap` is the physically
    relevant one) rather than silently keeping only the biggest.
    """
    freqs = np.asarray(freqs, dtype=float)
    if freqs.ndim != 2:
        raise ValueError("freqs must be [n_k, n_bands]")
    freqs = np.sort(freqs, axis=1)
    n_bands = freqs.shape[1]

    _max = np.nanmax if nan_safe else np.max
    _min = np.nanmin if nan_safe else np.min

    gaps = []
    for i in range(n_bands - 1):
        top_lower = float(_max(freqs[:, i]))
        bot_upper = float(_min(freqs[:, i + 1]))
        if not (np.isfinite(top_lower) and np.isfinite(bot_upper)):
            continue
        if bot_upper <= top_lower:
            continue
        center = 0.5 * (top_lower + bot_upper)
        if f_window is not None and not (f_window[0] <= center <= f_window[1]):
            continue
        ng = (bot_upper - top_lower) / center
        gaps.append(Gap(top_lower, bot_upper, center, bot_upper - top_lower,
                        ng, i, True))
    gaps.sort(key=lambda gp: gp.normalized_gap, reverse=True)
    return gaps


def trusted_window(band_sets) -> tuple:
    """Frequency window (f_lo, f_hi) covered by EVERY set at EVERY k-point.

    band_sets: sequence of [n_k, n_bands] arrays (NaN = mode not returned), one
    per symmetry sector. A shift-and-invert eigensolver returns the modes
    nearest its shift, so below the highest "lowest returned mode" and above
    the lowest "highest returned mode" some sector at some k was not computed.
    Returns (nan, nan) if a set has a k-point without any finite frequency.
    """
    f_lo, f_hi = -np.inf, np.inf
    for b in band_sets:
        b = np.asarray(b, dtype=float)
        if b.ndim != 2:
            raise ValueError("each band set must be [n_k, n_bands]")
        if b.size == 0 or not np.all(np.any(np.isfinite(b), axis=1)):
            return (float("nan"), float("nan"))
        f_lo = max(f_lo, float(np.max(np.nanmin(b, axis=1))))
        f_hi = min(f_hi, float(np.min(np.nanmax(b, axis=1))))
    return (f_lo, f_hi)


def interval_gap(band_sets, f_target: float, rel_tol: float = 0.5) -> Gap:
    """Largest gap near f_target by FREQUENCY INTERVAL, over one or more sectors.

    The widest interval, inside trusted_window(band_sets), that holds no
    eigenfrequency of ANY set at ANY k-point, with its center within
    +/- rel_tol of f_target. Pass one set for a symmetry-restricted gap, all
    sectors for a complete gap. Band indices are not used (lower_band = -1),
    because column j of a nearest-to-shift solve need not be the same band at
    every k or in every sector. Python port of gapReportUnion in
    comsol/omc_mechanical_livelink_rib.m.
    """
    f_lo, f_hi = trusted_window(band_sets)
    if not (np.isfinite(f_lo) and np.isfinite(f_hi) and f_hi > f_lo):
        return Gap.empty()
    f = np.concatenate([np.asarray(b, dtype=float).ravel() for b in band_sets])
    f = np.unique(f[np.isfinite(f) & (f >= f_lo) & (f <= f_hi)])
    if f.size < 2:
        return Gap.empty()
    lo, hi = f[:-1], f[1:]
    center = 0.5 * (lo + hi)
    ng = (hi - lo) / center
    ok = np.abs(center - f_target) < rel_tol * f_target
    if not np.any(ok):
        return Gap.empty()
    i = int(np.argmax(np.where(ok, ng, -np.inf)))
    return Gap(float(lo[i]), float(hi[i]), float(center[i]),
               float(hi[i] - lo[i]), float(ng[i]), -1, True)


if __name__ == "__main__":
    # toy demo: two bands with a clean gap
    k = np.linspace(0, np.pi, 11)
    b1 = 0.20 + 0.02 * np.sin(k)
    b2 = 0.30 + 0.02 * np.cos(k)
    g = largest_gap(np.column_stack([b1, b2]))
    print(g)
