"""Ask/tell optimizer with a clean interface and graceful fallbacks.

Stage 1 (always available, numpy-only): Sobol/LHS-ish space-filling + random.
Stage 2 (if installed): Optuna TPE.
Stage 3 (if installed): BoTorch Bayesian optimization.

The loop only uses: initialize_optimizer(data) -> opt; ask(opt) -> u;
tell(opt, u, score). Swap the backend without touching scripts/run_loop.py.
"""
from __future__ import annotations

import numpy as np

N_DIM = 5  # a, w, hx, hy, t (normalized)


# ----------------------------- Stage 1: numpy --------------------------------
class RandomSampler:
    """Space-filling fallback: deterministic Sobol-like Halton then random."""

    _BASES = [2, 3, 5, 7, 11, 13, 17, 19]

    def __init__(self, n_init=50, seed=0, n_dim=N_DIM):
        if not 1 <= n_dim <= len(self._BASES):
            raise ValueError(f"n_dim must be 1..{len(self._BASES)}, got {n_dim}")
        self.n_init = n_init
        self.n_dim = n_dim
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
        if self.count < self.n_init:
            u = [self._halton(self.count + 1, b)
                 for b in self._BASES[:self.n_dim]]
        else:
            u = self.rng.random(self.n_dim).tolist()
        self.count += 1
        return u

    def tell(self, u, score, **kw):
        pass


# ----------------------------- Stage 2: Optuna -------------------------------
class OptunaSampler:
    def __init__(self, n_init=20, seed=0, n_dim=N_DIM):
        import optuna
        self.n_dim = n_dim
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        self._optuna = optuna
        self.study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(n_startup_trials=n_init, seed=seed))
        self._pending = {}

    def ask(self):
        trial = self.study.ask(
            {f"u{i}": self._optuna.distributions.FloatDistribution(0, 1)
             for i in range(self.n_dim)})
        u = [trial.params[f"u{i}"] for i in range(self.n_dim)]
        self._pending[tuple(round(x, 8) for x in u)] = trial
        return u

    def tell(self, u, score, **kw):
        key = tuple(round(x, 8) for x in u)
        trial = self._pending.pop(key, None)
        if trial is not None:
            self.study.tell(trial, score)


def initialize_optimizer(data=None, *, backend="auto", n_init=50, seed=0,
                         n_dim=N_DIM):
    """Return an optimizer object exposing .ask() and .tell().

    backend: 'auto' | 'random' | 'optuna' | 'botorch'.
    'auto' picks the best installed backend.
    n_dim: length of u (default 5 = hole cell; geometry.n_dim(celltype)).
    """
    if backend in ("auto", "optuna"):
        try:
            opt = OptunaSampler(n_init=min(n_init, 20), seed=seed, n_dim=n_dim)
            _replay(opt, data, n_dim)
            return opt
        except Exception:
            if backend == "optuna":
                raise
    # fallback
    opt = RandomSampler(n_init=n_init, seed=seed, n_dim=n_dim)
    return opt


def _replay(opt, data, n_dim=N_DIM):
    """Seed an optimizer with previously completed results (resume)."""
    if not data:
        return
    for r in data:
        if (r.get("status") == "success" and r.get("u") is not None
                and len(r["u"]) == n_dim):
            try:
                opt.tell(r["u"], r["score"])
            except Exception:
                pass


def ask(opt):
    return opt.ask()


def tell(opt, u, score, constraints=None):
    opt.tell(u, score, constraints=constraints)
