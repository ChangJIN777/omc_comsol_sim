"""Ask/tell optimizer with a clean interface and graceful fallbacks.

Ported unchanged in structure from python-scripts/src/optimizer.py -- only
N_DIM changes (3 FREE design variables here: a, w, r; r1/r2/th are fixed in
configs/bounds_2d.yaml and are not part of u). See src/geometry2d.py:VARS for
the exact ordering -- N_DIM is asserted against it at import.

Stage 1 (always available, numpy-only): Sobol/LHS-ish space-filling + random.
Stage 2 (if installed): Optuna TPE.

The loop only uses: initialize_optimizer(data) -> opt; ask(opt) -> u;
tell(opt, u, score). Swap the backend without touching scripts/run_loop_2d.py.
"""
from __future__ import annotations

import numpy as np

N_DIM = 3  # a, w, r (normalized) -- see geometry2d.VARS

try:  # keep N_DIM and VARS from drifting apart; geometry2d is dependency-light
    from geometry2d import VARS as _VARS
    if len(_VARS) != N_DIM:
        raise AssertionError(
            f"optimizer.N_DIM={N_DIM} but geometry2d.VARS has {len(_VARS)} "
            f"entries {_VARS}. Update both together, and start a fresh results "
            f"DB -- stored u vectors are not comparable across a dimension "
            f"change.")
except ImportError:  # optimizer is usable standalone
    pass


# ----------------------------- Stage 1: numpy --------------------------------
class RandomSampler:
    """Space-filling fallback: deterministic Sobol-like Halton then random."""

    def __init__(self, n_init=50, seed=0):
        self.n_init = n_init
        self.rng = np.random.default_rng(seed)
        self.count = 0

    def _halton(self, i, base):
        f, r = 1.0, 0.0
        while i > 0:
            f /= base
            r += f * (i % base)
            i //= base
        return r

    def ask(self):
        bases = [2, 3, 5, 7, 11, 13]
        if self.count < self.n_init:
            u = [self._halton(self.count + 1, b) for b in bases[:N_DIM]]
        else:
            u = self.rng.random(N_DIM).tolist()
        self.count += 1
        return u

    def tell(self, u, score, **kw):
        pass


# ----------------------------- Stage 2: Optuna -------------------------------
class OptunaSampler:
    def __init__(self, n_init=20, seed=0):
        import optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        self._optuna = optuna
        self.study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(n_startup_trials=n_init, seed=seed))
        self._pending = {}

    def ask(self):
        trial = self.study.ask(
            {f"u{i}": self._optuna.distributions.FloatDistribution(0, 1)
             for i in range(N_DIM)})
        u = [trial.params[f"u{i}"] for i in range(N_DIM)]
        self._pending[tuple(round(x, 8) for x in u)] = trial
        return u

    def tell(self, u, score, **kw):
        key = tuple(round(x, 8) for x in u)
        trial = self._pending.pop(key, None)
        if trial is not None:
            self.study.tell(trial, score)


def initialize_optimizer(data=None, *, backend="auto", n_init=50, seed=0):
    """Return an optimizer object exposing .ask() and .tell().

    backend: 'auto' | 'random' | 'optuna'. 'auto' picks the best installed.
    """
    if backend in ("auto", "optuna"):
        try:
            opt = OptunaSampler(n_init=min(n_init, 20), seed=seed)
            _replay(opt, data)
            return opt
        except Exception:
            if backend == "optuna":
                raise
    opt = RandomSampler(n_init=n_init, seed=seed)
    return opt


def _replay(opt, data):
    """Seed an optimizer with previously completed results (resume)."""
    if not data:
        return
    for r in data:
        if r.get("status") == "success" and r.get("u") is not None:
            try:
                opt.tell(r["u"], r["score"])
            except Exception:
                pass


def ask(opt):
    return opt.ask()


def tell(opt, u, score, constraints=None):
    opt.tell(u, score, constraints=constraints)
