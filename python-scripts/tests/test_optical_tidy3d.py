"""Tests of the Tidy3D optical engine (src/optical_tidy3d.py).

NOTHING here can reach the Tidy3D cloud: the tests that exercise the module
diamond_bandStructure_sim.py replace `tidy3d` by a MagicMock in sys.modules
(so they also run, and stay offline, on a machine where tidy3d is installed).
Only test_real_tidy3d_builds_locally uses the real package, in dry_run mode
with the cloud function replaced by one that fails; it is skipped when tidy3d
is not installed.

Run:  python tests/test_optical_tidy3d.py
"""
import contextlib
import json
import os
import sys
import tempfile
from unittest import mock

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
sys.path.insert(0, os.path.join(_HERE, "..", "scripts"))

import numpy as np

import objective
import objective_rib
import optical_tidy3d as ot
from comsol_optics import C0, guided_interval_gap, optical_gap_tracked
from geometry import Geometry, RibGeometry, phys_to_u, u_to_geometry
from objective import evaluate_candidate

F_T = C0 / 1550e-9
RIB_CELL = dict(a=570e-9, w=1747e-9, t=391e-9, s=254e-9, th=250e-9,
                r1=40e-9, r2=40e-9)
HOLE_U = [0.5, 0.6, 0.5, 0.6, 0.5]
_FAKE_MODULES = ("tidy3d", "tidy3d.web", "tidy3d.plugins",
                 "tidy3d.plugins.resonance", "diamond_bandStructure_sim")


def _rib_u0():
    return phys_to_u(RIB_CELL, celltype="rib")


@contextlib.contextmanager
def fake_tidy3d(cost=0.6):
    """Load the REAL diamond_bandStructure_sim.py on top of a mocked tidy3d.

    Yields (td mock, module). The layout is the real one: `tidy3d.web` is a
    SUBMODULE that `import tidy3d` does not attach, so `td.web` raises
    AttributeError here exactly as with tidy3d 2.12 (a plain MagicMock would
    invent the attribute and hide a missing `import tidy3d.web`). The web
    mock is reachable for the tests as td.WEB; its Batch has an estimate_cost
    that returns `cost`, so run_cloud's own budget logic runs for real, offline.
    """
    td = mock.MagicMock(name="tidy3d")
    td.__name__ = "tidy3d"               # `import tidy3d.web as web` needs it
    td.constants.C_0 = C0 * 1e6          # um/s
    td.inf = float("inf")
    web = mock.MagicMock(name="tidy3d.web")
    web.Batch.return_value.estimate_cost.return_value = cost
    del td.web                           # not an attribute of the package
    td.WEB = web                         # test handle only
    saved = {name: sys.modules.get(name) for name in _FAKE_MODULES}
    sys.modules.update({"tidy3d": td, "tidy3d.web": web,
                        "tidy3d.plugins": td.plugins,
                        "tidy3d.plugins.resonance": td.plugins.resonance})
    sys.modules.pop("diamond_bandStructure_sim", None)
    real_have = ot.have_tidy3d
    ot.reset_session()
    ot.have_tidy3d = lambda: True
    try:
        yield td, ot.load_bs()
    finally:
        ot.have_tidy3d = real_have
        ot.reset_session()
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def _resonances(k_frac, a, kind="rib"):
    """Tidy3D-style resonances: one ragged list per k (Hz).

    Two guided bands below a 175-210 THz gap and an air band above it (the
    bands of tests/test_pipeline._fake_optical_bands), plus leaky resonances
    above the light line at some k only.
    """
    k = np.asarray(k_frac, dtype=float)
    x = (k - k[0]) / (1.0 - k[0]) if k.size > 1 else np.ones(1)
    rows = []
    for i, (kk, xx) in enumerate(zip(k, x)):
        ll = C0 * kk / (2 * a)
        row = [(120 + 20 * xx) * 1e12, (150 + 25 * xx ** 2) * 1e12,
               (230 - 20 * xx ** 2) * 1e12]
        if i % 2 == 0:
            row += [1.05 * ll, 1.2 * ll]
        rows.append(np.array(row))
    return rows


def _pad(rows):
    out = np.full((len(rows), max(len(r) for r in rows)), np.nan)
    for i, r in enumerate(rows):
        out[i, :len(r)] = np.sort(r)
    return out


def _install_fake_analysis(bs, a_m):
    """bs.analyze_bands -> synthetic resonances at the k of the simulations
    that were 'run' (read from the keys of sim_datas), real band_edges."""
    def analyze(sim_datas, p):
        idx = sorted(int(name.split("_")[1]) if name.startswith("sim_")
                     else int(round(float(name[2:]) * 2 * p.Nk))
                     for name in sim_datas)
        k_frac = [i / p.Nk for i in idx]
        rows = _resonances(k_frac, a_m)
        bands = [dict(k=kf / 2.0, freq=np.sort(r), Q=np.full(len(r), 1e4))
                 for kf, r in zip(k_frac, rows)]
        res = bs.BandResult(bands=bands, f_light_X=p.f_light_X, has_X=True)
        (res.f_lower, res.f_upper, res.bandGap, res.midGap,
         res.fitness) = bs.band_edges(bands[-1]["freq"], p)
        return res
    bs.analyze_bands = analyze


# ---- geometry mapping ---------------------------------------------------------

def test_geometry_mapping_units_and_fields():
    rib = RibGeometry(**RIB_CELL)
    kw = ot.band_params_kwargs(rib, index=2.386, f_target=F_T)
    assert kw["geometry"] == "rib_block" and kw["mode"] == "TE"
    assert kw["xsect"] == "rect" and kw["wedge_angle"] == 30.0
    # meters -> um, each repo field on the right BandParams field
    assert kw["a"] == 0.57 and kw["t"] == 0.391 and kw["w"] == 1.747
    assert kw["spine_width"] == 0.254 and kw["thickness"] == 0.25
    assert kw["fillet_r1"] == 0.04 and kw["fillet_r2"] == 0.04
    assert kw["index"] == 2.386 and abs(kw["target_freq_thz"] - 193.414) < 1e-3
    assert not {"beam_width", "hole_x", "hole_y"} & set(kw)

    hole = Geometry(a=700e-9, w=660e-9, t=325e-9, hx=230e-9, hy=280e-9)
    kw = ot.band_params_kwargs(hole, index=2.386, f_target=F_T)
    assert kw["geometry"] == "holey" and kw["a"] == 0.7
    assert kw["beam_width"] == 0.66
    assert kw["thickness"] == 0.325          # hole cell: t IS the thickness
    assert kw["hole_x"] == 0.46 and kw["hole_y"] == 0.56   # FULL axes = 2*semi
    assert not {"t", "w", "spine_width"} & set(kw)

    assert ot.band_params_kwargs(rib, index=2.4, f_target=F_T,
                                 optical_study="opt TM")["mode"] == "TM"
    for bad in (dict(optical_study="opt XX"), dict(overrides={"a": 1.0})):
        try:
            ot.band_params_kwargs(rib, index=2.4, f_target=F_T, **bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted {bad}")
    # the index is the repo's material config: 2.386 for every engine
    assert ot.material_index() == 2.386


def test_mapping_against_the_module_geometry():
    """The mapped numbers mean the same thing in diamond_bandStructure_sim."""
    with fake_tidy3d() as (td, bs):
        rib = RibGeometry(**RIB_CELL)
        p = bs.BandParams(**ot.band_params_kwargs(rib, index=2.386, f_target=F_T))
        assert p.n == 2.386 and abs(p.f_light_X - C0 / (2 * rib.a)) < 1.0
        assert p.n == bs.MATERIAL_INDEX["diamond_telecom"]
        # "opt TE" (PEC on y = 0, Ey even) -> Ey dipoles, PEC y plane (-1),
        # PMC mid-thickness plane (+1)
        assert p.polarization == "Ey" and tuple(p.symmetry) == (0, -1, 1)
        v = bs.rib_block_vertices(0.0, p.a, p.t, p.w, p.spine_width,
                                  p.fillet_r1, p.fillet_r2, p.fillet_points)
        um = 1e6
        assert abs(v[:, 1].max() - rib.w / 2 * um) < 1e-9       # w = tip to tip
        assert abs(v[:, 0].max() - rib.a / 2 * um) < 1e-9
        at_edge = v[np.isclose(v[:, 0], rib.a / 2 * um)]
        assert abs(at_edge[:, 1].max() - rib.s / 2 * um) < 1e-9  # s = full spine
        tip = v[v[:, 1] > rib.w / 2 * um - 1e-9]
        assert abs(np.ptp(tip[:, 0]) - (rib.t - 2 * rib.r2) * um) < 1e-9  # t along the beam
        # area of the polygon vs the analytic filleted cell (repo meaning)
        corner = (1 - np.pi / 4)
        area = (rib.a * rib.s + rib.t * (rib.w - rib.s)
                - 4 * corner * rib.r2 ** 2 + 4 * corner * rib.r1 ** 2) * um ** 2
        assert abs(bs._polygon_area(v) - area) / area < 2e-3    # 8-point arcs

        tm = bs.BandParams(**ot.band_params_kwargs(
            rib, index=2.386, f_target=F_T, optical_study="opt TM"))
        assert tm.polarization == "Ez" and tuple(tm.symmetry) == (0, 1, -1)

        hole = u_to_geometry(HOLE_U)
        p = bs.BandParams(**ot.band_params_kwargs(hole, index=2.386, f_target=F_T))
        cell = bs.build_unit_cell(p)
        beam, hole_poly = cell.outline
        assert abs(np.ptp(beam[:, 1]) - hole.w * um) < 1e-9
        assert abs(hole_poly[:, 0].max() - hole.hx * um) < 1e-9  # semi-axis hx
        assert abs(hole_poly[:, 1].max() - hole.hy * um) < 1e-9  # semi-axis hy
        fill = 1 - np.pi * hole.hx * hole.hy / (hole.a * hole.w)
        assert abs(cell.fill_factor - fill) < 1e-3


def test_k_plan():
    plan = ot.k_plan("hole", 700e-9, k_mode="x_only", f_target=F_T)
    assert (plan.Nk, plan.indices, plan.k_frac, plan.x_only) == (1, [1], [1.0], True)
    assert plan.k_module() == [0.5]
    # hole sweep = the zone-edge half of comsol_optics._opt_sweep
    plan = ot.k_plan("hole", 700e-9, k_mode="sweep", n_k=7, f_target=F_T)
    assert np.allclose(plan.k_frac, np.linspace(0.5, 1.0, 7)) and not plan.x_only
    # rib sweep: n_k points from the light line of the lowest scored frequency
    from comsol_optics import rib_k_grid
    for a in (300e-9, 570e-9, 900e-9):
        for n_k in (2, 5, 13):
            plan = ot.k_plan("rib", a, k_mode="sweep", n_k=n_k, f_target=F_T)
            ref = rib_k_grid(a, n_k, F_T)
            assert len(plan.k_frac) == n_k and plan.k_frac[-1] == 1.0
            assert plan.k_frac[0] >= ref[0] - 1e-12
            assert np.diff(plan.k_frac).max() <= np.diff(ref).max() + 1e-12
            assert plan.indices[-1] == plan.Nk
            assert np.allclose(np.array(plan.k_module()) * 2, plan.k_frac)


# ---- scoring ------------------------------------------------------------------

def test_rib_scoring_is_the_comsol_interval_gap():
    a = RIB_CELL["a"]
    plan = ot.k_plan("rib", a, k_mode="sweep", n_k=13, f_target=F_T)
    rows = _resonances(plan.k_frac, a)
    k_frac, freqs = ot.bands_to_arrays(
        [dict(k=kf / 2, freq=r) for kf, r in zip(plan.k_frac, rows)][::-1])
    assert np.allclose(k_frac, plan.k_frac)             # sorted by k again
    window = (93.4e12, 293.4e12)
    gp, method = ot.score_resonances("rib", k_frac, freqs, a, f_target=F_T,
                                     k_mode="sweep", window_hz=window)
    assert method == "guided_interval_gap" and gp.found
    assert abs(gp.f_lower - 175e12) < 1e6 and abs(gp.f_upper - 210e12) < 1e6
    # same numbers as the COMSOL-style call on a full [n_k, n_bands] array
    full = np.column_stack([_pad(rows)[:, :3],
                            C0 * k_frac / (2 * a) * 1.05])
    ref = guided_interval_gap(k_frac, full, a, f_target=F_T)
    assert (gp.f_lower, gp.f_upper, gp.normalized_gap) == (
        ref.f_lower, ref.f_upper, ref.normalized_gap)
    # FDTD: a k without any resonance is allowed (COMSOL: untrusted -> none)
    hollow = freqs.copy()
    hollow[0, :] = np.nan
    assert ot.score_resonances("rib", k_frac, hollow, a, f_target=F_T,
                               k_mode="sweep", window_hz=window)[0].found
    assert not guided_interval_gap(k_frac, hollow, a, f_target=F_T).found
    # the gap must lie inside the simulated frequency window: with the window
    # cut at 200 THz the 175-210 THz gap has no trusted upper edge
    cut = ot.score_resonances("rib", k_frac, freqs, a, f_target=F_T,
                              k_mode="sweep", window_hz=(93.4e12, 200e12))[0]
    assert cut.f_upper <= 200e12 and cut.normalized_gap < 0.1
    assert not ot.score_resonances("rib", k_frac, freqs, a, f_target=F_T,
                                   k_mode="sweep",
                                   window_hz=(155e12, 200e12))[0].found
    # a band crossing the gap closes it, as with COMSOL
    x = (k_frac - k_frac[0]) / (1 - k_frac[0])
    steep = np.column_stack([freqs, (130 + 110 * x) * 1e12])
    assert not ot.score_resonances("rib", k_frac, steep, a, f_target=F_T,
                                   k_mode="sweep", window_hz=window)[0].found
    # x_only on the rib: X-point gap only, labelled as such
    gx, method = ot.score_resonances("rib", k_frac[-1:], freqs[-1:], a,
                                     f_target=F_T, k_mode="x_only")
    assert method == "optical_gap_tracked_x" and gx.found
    assert abs(gx.f_lower - 140e12) < 1e6 and abs(gx.f_upper - 175e12) < 1e6
    for bad_k in (k_frac[:-1], k_frac[::-1]):           # must end at X
        try:
            ot.score_resonances("rib", bad_k, freqs[:len(bad_k)], a,
                                f_target=F_T, window_hz=window)
        except ValueError:
            pass
        else:
            raise AssertionError("k grid not ending at X accepted")


def test_hole_scoring_is_the_comsol_tracked_gap():
    a = 500e-9                                   # light line at X: 299.8 THz
    ll = C0 / (2 * a)
    # X only: two guided resonances + a leaky one
    row = np.array([[170e12, 215e12, 1.1 * ll]])
    gp, method = ot.score_resonances("hole", [1.0], row, a, f_target=F_T,
                                     k_mode="x_only")
    assert method == "optical_gap_tracked_x" and gp.found
    assert abs(gp.normalized_gap - 45.0 / 192.5) < 1e-9
    ref = optical_gap_tracked(np.array([[170e12, 215e12]]), f_target=F_T)
    assert gp.normalized_gap == ref.normalized_gap
    # one guided resonance only: the leaky one is NOT used as the air band
    one = np.array([[170e12, 1.02 * ll, 1.1 * ll]])
    assert not ot.score_resonances("hole", [1.0], one, a, f_target=F_T,
                                   k_mode="x_only")[0].found
    # fundamental gap far from the target is rejected (+/- 30 %)
    far = np.array([[100e12, 120e12]])
    assert not ot.score_resonances("hole", [1.0], far, a, f_target=F_T,
                                   k_mode="x_only")[0].found
    # sweep: the gap is tracked over k in [0.5, 1]; rows below 0.5 are ignored
    k = np.array([0.25, 0.5, 0.75, 1.0])
    freqs = np.array([[60e12, 191e12, np.nan],        # ignored (k < 0.5)
                      [110e12, np.nan, np.nan],       # air band not guided here
                      [172e12, 200e12, np.nan],       # narrowest: 28 THz
                      [170e12, 215e12, 1.1 * ll]])
    gp, method = ot.score_resonances("hole", k, freqs, a, f_target=F_T,
                                     k_mode="sweep")
    assert method == "optical_gap_tracked" and gp.found
    assert abs(gp.normalized_gap - 28.0 / 192.5) < 1e-9


def test_module_x_point_gap_is_recorded_next_to_the_score():
    with fake_tidy3d() as (td, bs):
        a = 600e-9
        p = bs.BandParams(geometry="holey", a=a * 1e6, target_freq_thz=193.414)
        ll = C0 / (2 * a)
        f_lo, f_up, gap, mid, fit = bs.band_edges([215e12, 170e12, 1.1 * ll], p)
        assert abs(f_lo - 170) < 1e-9 and abs(f_up - 215) < 1e-9
        gp, _ = ot.score_resonances("hole", [1.0],
                                    np.array([[170e12, 215e12, 1.1 * ll]]), a,
                                    f_target=F_T, k_mode="x_only")
        assert abs(gp.normalized_gap - gap / mid) < 1e-12   # same two edges


# ---- engine switch, ids --------------------------------------------------------

def _fake_run(calls):
    def run(g, *, celltype, eval_id, f_target, options, rel_tol=0.30,
            optical_study="opt TE", n_k_default=None, xsect=None,
            verbose=True, on_progress=None, optical_tol_hz=None):
        opts = ot.Tidy3DOptions.from_config(options)
        k_mode, n_k = opts.k_sampling(celltype, n_k_default)
        calls.append((celltype, type(g).__name__, eval_id, k_mode, n_k))
        from bandgap import Gap
        gp = Gap(175e12, 210e12, 192.5e12, 35e12, 35.0 / 192.5, -1, True)
        return dict(gap=gp, k_frac=np.array([1.0]),
                    freqs_hz=np.array([[175e12, 210e12]]),
                    record=dict(
                        optical_tidy3d_detail=dict(engine="tidy3d", k_mode=k_mode,
                                                   min_resolvable_Hz=None),
                        optical_gap_detail=ot._gap_to_dict(gp),
                        _k_o=[1.0], _freqs_o=[[175e12, 210e12]]))
    return run


def test_engine_switch_and_distinct_ids():
    import optical_comsol
    calls = []
    real_run, real_comsol = ot.run_optical_tidy3d, optical_comsol.run_optical_comsol
    real_rib = objective_rib._solve_optical_bands
    k15 = np.linspace(1e-3, 1.0, 15)
    ot.run_optical_tidy3d = _fake_run(calls)
    optical_comsol.run_optical_comsol = lambda g, **kw: dict(
        freqs_hz=np.column_stack([150e12 + 25e12 * k15, 230e12 - 20e12 * k15]),
        k_frac=k15, a=g.a)
    objective_rib._solve_optical_bands = (
        lambda g, study, k_frac, n_bands, f_target, **kw:
        _pad(_resonances(k_frac, g.a)))
    stub = dict(mech_backend="surrogate_stub")
    try:
        # hole: one key switches the engine
        h_com = evaluate_candidate(HOLE_U, optical_backend="comsol", **stub)
        h_td = evaluate_candidate(HOLE_U, optical_backend="tidy3d", **stub)
        h_sw = evaluate_candidate(HOLE_U, optical_backend="tidy3d",
                                  tidy3d_options=dict(k_mode="sweep"), **stub)
        h_sur = evaluate_candidate(HOLE_U, optical_backend="surrogate", **stub)
        # rib
        r_com = evaluate_candidate(_rib_u0(), celltype="rib",
                                   optical_backend="comsol", **stub)
        r_td = evaluate_candidate(_rib_u0(), celltype="rib",
                                  optical_backend="tidy3d", **stub)
        r_x = evaluate_candidate(_rib_u0(), celltype="rib",
                                 optical_backend="tidy3d",
                                 tidy3d_options=dict(k_mode="x_only"), **stub)
        r_n7 = evaluate_candidate(_rib_u0(), celltype="rib",
                                  optical_backend="tidy3d",
                                  tidy3d_options=dict(n_k=7), **stub)
    finally:
        ot.run_optical_tidy3d = real_run
        optical_comsol.run_optical_comsol = real_comsol
        objective_rib._solve_optical_bands = real_rib

    recs = [h_com, h_td, h_sw, h_sur, r_com, r_td, r_x, r_n7]
    assert all(r["status"] == "success" for r in recs), [r["status"] for r in recs]
    assert len({r["id"] for r in recs}) == len(recs)        # no id collides
    assert [r["optical_engine"] for r in recs] == [
        "comsol", "tidy3d", "tidy3d", "surrogate",
        "comsol", "tidy3d", "tidy3d", "tidy3d"]
    assert all(r["optical_engine"] == r["optical_backend"] for r in recs)
    # every record with optics carries index, cross-section (and, for COMSOL,
    # the parity family) in its id: none of them equals an id made before the
    # index / xsect change (n = 2.40, no xsect), so old records are kept
    assert h_com["id"] == objective._hash_u(
        HOLE_U, ("comsol", "surrogate_stub", "n=2.3860|xsect=rect|z=even"))
    assert h_sur["id"] == objective._hash_u(
        HOLE_U, ("surrogate", "surrogate_stub", "n=2.3860|xsect=rect"))
    old = ("comsol", "surrogate_stub", "complete", 15, 20, "opt TE", 13, 10, 0.30)
    assert r_com["id"] == objective_rib._rib_id(
        _rib_u0(), old + ("n=2.3860|xsect=rect|z=even", "score=v2|otol=2THz"))
    assert h_com["id"] != objective._hash_u(HOLE_U, ("comsol", "surrogate_stub"))
    assert r_com["id"] != objective_rib._rib_id(_rib_u0(), old)
    # defaults per cell type: hole x_only, rib sweep with n_k_opt points
    assert [c[3:] for c in calls] == [("x_only", None), ("sweep", 7),
                                      ("sweep", 13), ("x_only", None),
                                      ("sweep", 7)]
    assert calls[0][:2] == ("hole", "Geometry") and calls[0][2] == h_td["id"]
    assert calls[2][:2] == ("rib", "RibGeometry") and calls[2][2] == r_td["id"]
    # same gap -> same score as the same gap from COMSOL (scoring is shared)
    for r in (h_td, r_td):
        assert abs(r["optical_gap"] - 35.0 / 192.5) < 1e-12
        assert r["optical_tidy3d_detail"]["engine"] == "tidy3d"
        assert r["optical_gap_detail"]["found"] and r["mechanical_gap"] is None
        json.dumps(r, allow_nan=False)
    assert abs(r_td["score"] - r_com["score"]) < 1e-9       # both 175-210 THz
    assert "optical_tidy3d_detail" not in h_com and "celltype" not in h_td
    # mechanics is never tidy3d; surrogate / mpb stay hole-only
    for kw in (dict(mech_backend="tidy3d", optical_backend="surrogate"),
               dict(celltype="rib", optical_backend="mpb",
                    mech_backend="surrogate_stub")):
        try:
            rec = evaluate_candidate(_rib_u0() if "celltype" in kw else HOLE_U, **kw)
        except ValueError:
            continue
        assert rec["status"] == "mech_failed", rec["status"]


def test_infeasible_cells_never_reach_tidy3d():
    calls = []
    real_run = ot.run_optical_tidy3d
    ot.run_optical_tidy3d = _fake_run(calls)
    try:
        rec = evaluate_candidate([0.0, 0.0, 1.0, 1.0, 0.5],
                                 optical_backend="tidy3d",
                                 mech_backend="surrogate_stub",
                                 tidy3d_options=dict(allow_cloud=True,
                                                     max_credits=5.0))
        bad = phys_to_u({**RIB_CELL, "t": RIB_CELL["a"]}, celltype="rib")
        rib = evaluate_candidate(bad, celltype="rib", optical_backend="tidy3d",
                                 mech_backend="surrogate_stub",
                                 tidy3d_options=dict(allow_cloud=True,
                                                     max_credits=5.0))
    finally:
        ot.run_optical_tidy3d = real_run
    assert rec["status"] == "infeasible" and rib["status"] == "infeasible"
    assert calls == []                      # min feature size before the solver


# ---- cost safety ---------------------------------------------------------------

def test_options_are_strict():
    assert ot.Tidy3DOptions.from_config(None).allow_cloud is False
    assert ot.Tidy3DOptions.from_config({}).max_credits is None
    for bad in (dict(alow_cloud=True), dict(allow_cloud="true"),
                dict(allow_cloud=1), dict(mode="go"), dict(k_mode="full"),
                dict(max_credits=0), dict(max_credits=-1), dict(n_k=1),
                dict(params=dict(a=0.5)), dict(params=dict(Nk=3))):
        try:
            ot.Tidy3DOptions.from_config(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted {bad}")
    opts = ot.Tidy3DOptions.from_config(dict(allow_cloud=True, max_credits="2"))
    assert opts.allow_cloud is True and opts.max_credits == 2.0
    assert opts.k_sampling("hole") == ("x_only", None)
    assert opts.k_sampling("rib", 13) == ("sweep", 13)
    assert ot.Tidy3DOptions.from_config(dict(n_k=5)).k_sampling("rib", 13) == ("sweep", 5)
    # cost-control keys are not part of the record id, physics keys are
    base = ot.id_settings(None, "rib", 13)
    assert ot.id_settings(dict(allow_cloud=True, max_credits=3, mode="dry_run",
                               results_dir="x"), "rib", 13) == base
    assert ot.id_settings(dict(index=2.40), "rib", 13) != base
    assert ot.id_settings(None, "rib", 13, xsect="wedge") != base
    assert ot.id_settings(None, "rib", 13, xsect="rect") == base
    assert ot.id_settings(dict(params=dict(min_steps_per_wvl=30)), "rib", 13) != base


def test_without_tidy3d_the_pipeline_still_works():
    if ot.have_tidy3d():
        print("  (skipped: tidy3d is installed)")
        return
    ot.reset_session()
    assert "tidy3d" not in sys.modules      # importing the adapter is lazy
    cache = {}
    with tempfile.TemporaryDirectory() as tmp:
        rec = evaluate_candidate(HOLE_U, optical_backend="tidy3d",
                                 mech_backend="surrogate_stub", cache=cache,
                                 tidy3d_options=dict(allow_cloud=True,
                                                     max_credits=1.0,
                                                     results_dir=tmp))
    assert rec["status"] == "optical_not_run"
    assert rec["optical_tidy3d_stop"] == "not_installed"
    assert rec["id"].endswith("_not_installed") and rec["score"] == -5.0
    assert rec["optical_engine"] == "tidy3d" and cache == {}
    assert rec["optical_tidy3d_detail"]["n_sims"] == 1
    json.dumps(rec, allow_nan=False)


def _eval(u, tmp, celltype="hole", cache=None, **tidy):
    return evaluate_candidate(u, celltype=celltype, optical_backend="tidy3d",
                              mech_backend="surrogate_stub", cache=cache,
                              tidy3d_options=dict(results_dir=tmp, **tidy))


def test_allow_cloud_guard_blocks_every_upload():
    with fake_tidy3d() as (td, bs), tempfile.TemporaryDirectory() as tmp:
        _install_fake_analysis(bs, u_to_geometry(HOLE_U).a)
        real_cloud = bs.run_cloud
        bs.run_cloud = mock.MagicMock(side_effect=AssertionError("cloud reached"))
        # 1) default: no tidy3d section at all / allow_cloud false
        for tidy in ({}, dict(allow_cloud=False, max_credits=100.0)):
            rec = _eval(HOLE_U, tmp, **tidy)
            assert rec["status"] == "optical_not_run"
            assert rec["optical_tidy3d_stop"] == "cloud_not_allowed"
            assert rec["id"].endswith("_cloud_not_allowed")
            assert "allow_cloud" in rec["error"]
        # 2) opted in, but no credit budget
        rec = _eval(HOLE_U, tmp, allow_cloud=True)
        assert rec["optical_tidy3d_stop"] == "no_budget"
        # 3) dry run: built locally, plan written, nothing uploaded
        rec = _eval(HOLE_U, tmp, allow_cloud=True, max_credits=5.0, mode="dry_run")
        assert rec["optical_tidy3d_stop"] == "dry_run"
        d = rec["optical_tidy3d_detail"]
        assert d["n_sims"] == 1 and d["k_frac"] == [1.0]
        assert d["symmetry"] == [0, -1, 1] and d["index"] == 2.386
        assert d["xsect"] == "rect" and d["z_symmetry"] is True
        assert os.path.exists(os.path.join(d["run_dir"], "plan.json"))
        assert d["run_dir"].startswith(tmp) and os.sep + "hole" + os.sep in d["run_dir"]
        # rib dry run: the default is a k sweep of n_k_opt simulations
        rec = _eval(_rib_u0(), tmp, celltype="rib", mode="dry_run")
        d = rec["optical_tidy3d_detail"]
        assert rec["optical_tidy3d_stop"] == "dry_run" and d["n_sims"] == 13
        assert d["band_params"]["geometry"] == "rib_block"
        # none of the above touched the cloud
        assert bs.run_cloud.call_count == 0
        assert td.WEB.Batch.call_count == 0 and td.WEB.run.call_count == 0
        assert ot.session_credits() == 0.0
        assert not os.path.exists(os.path.join(tmp, "credit_ledger.jsonl"))
        # 4) _submit itself refuses without the opt-in (defence in depth)
        bs.run_cloud = real_cloud
        try:
            ot._submit(bs, ot.Tidy3DOptions(), {}, tmp, "x", bs.Budget(1.0), False)
        except ot.Tidy3DNotRun as err:
            assert err.tidy3d_stop == "cloud_not_allowed"
        else:
            raise AssertionError("_submit ran without allow_cloud")
        assert td.WEB.Batch.call_count == 0


def test_estimate_mode_starts_nothing():
    with fake_tidy3d(cost=0.37) as (td, bs), tempfile.TemporaryDirectory() as tmp:
        rec = _eval(HOLE_U, tmp, mode="estimate")           # not opted in
        assert rec["optical_tidy3d_stop"] == "cloud_not_allowed"
        assert td.WEB.Batch.call_count == 0
        rec = _eval(HOLE_U, tmp, mode="estimate", allow_cloud=True)
        batch = td.WEB.Batch.return_value
        assert rec["optical_tidy3d_stop"] == "estimate_only"
        assert rec["optical_tidy3d_detail"]["estimated_credits"] == 0.37
        assert batch.upload.call_count == 1 and batch.delete.call_count == 1
        assert batch.run.call_count == 0                    # never started
        assert ot.session_credits() == 0.0
        assert not os.path.exists(os.path.join(tmp, "credit_ledger.jsonl"))


def test_credit_budget_is_cumulative_and_stops_cleanly():
    u2 = [0.55, 0.6, 0.5, 0.6, 0.5]
    with fake_tidy3d(cost=0.6) as (td, bs), tempfile.TemporaryDirectory() as tmp:
        _install_fake_analysis(bs, u_to_geometry(HOLE_U).a)
        batch = td.WEB.Batch.return_value
        tidy = dict(allow_cloud=True, max_credits=1.0)
        cache = {}
        rec = _eval(HOLE_U, tmp, cache=cache, **tidy)        # 0.6 <= 1.0: runs
        assert rec["status"] == "success", rec.get("error")
        d = rec["optical_tidy3d_detail"]
        assert d["estimated_credits"] == 0.6 and d["session_credits"] == 0.6
        assert d["gap_method"] == "optical_gap_tracked_x" and not d["reused_download"]
        # two lowest guided resonances at X (140, 175 THz): the module's own
        # result and the scored optical_gap_tracked gap agree for x_only
        assert abs(d["x_point"]["f_lower_THz"] - 140.0) < 1e-6
        assert abs(d["x_point"]["f_upper_THz"] - 175.0) < 1e-6
        assert abs(rec["optical_gap"] - d["x_point"]["normalized_gap"]) < 1e-9
        assert abs(rec["optical_gap"] - 35.0 / 157.5) < 1e-9
        assert batch.run.call_count == 1
        assert batch.run.call_args.kwargs["path_dir"] == d["run_dir"]
        tasks = td.WEB.Batch.call_args.kwargs["simulations"]
        assert list(tasks) == ["hole_%s_sim_1" % rec["id"]]
        for name in ("plan.json", "params.json", "evaluation.json", "bands.csv"):
            assert os.path.exists(os.path.join(d["run_dir"], name)), name
        assert cache[rec["id"]] is rec and rec["_k_o"] == [1.0]
        json.dumps(rec, allow_nan=False)

        # second cell: 0.6 + 0.6 > 1.0 -> priced, drafts deleted, NOT started
        rec2 = _eval(u2, tmp, cache=cache, **tidy)
        assert rec2["status"] == "optical_not_run"
        assert rec2["optical_tidy3d_stop"] == "budget_exceeded"
        assert rec2["optical_tidy3d_detail"]["refused_estimate"] == 0.6
        assert batch.run.call_count == 1 and batch.delete.call_count == 1
        assert ot.session_credits() == 0.6 and len(cache) == 1
        assert rec2["id"] != rec["id"] and rec2["id"].endswith("_budget_exceeded")
        with open(os.path.join(tmp, "credit_ledger.jsonl")) as fh:
            ledger = [json.loads(line) for line in fh]
        assert [e["estimated_credits"] for e in ledger] == [0.6]
        assert ledger[0]["id"] == rec["id"] and ledger[0]["status"] == "run"

        # budget already used up: stop before any upload
        n_batches = td.WEB.Batch.call_count
        rec3 = _eval(u2, tmp, allow_cloud=True, max_credits=0.5)
        assert rec3["optical_tidy3d_stop"] == "budget_exceeded"
        assert td.WEB.Batch.call_count == n_batches

        # a new process (session reset) still honours max_credits_total
        ot._budget = None
        rec4 = _eval(u2, tmp, allow_cloud=True, max_credits=5.0,
                     max_credits_total=1.0)
        assert rec4["optical_tidy3d_stop"] == "budget_exceeded"
        assert batch.run.call_count == 1

        # downloaded results of the same cell are re-analyzed, never paid twice,
        # even with the cloud switched off
        bs.load_run = lambda run_dir: {"k=0.500000": object()}
        n_batches = td.WEB.Batch.call_count
        again = _eval(HOLE_U, tmp, allow_cloud=False)
        assert again["status"] == "success" and again["id"] == rec["id"]
        assert again["optical_tidy3d_detail"]["reused_download"] is True
        assert again["optical_tidy3d_detail"]["estimated_credits"] == 0.0
        assert again["optical_gap"] == rec["optical_gap"]
        assert td.WEB.Batch.call_count == n_batches and batch.run.call_count == 1
        # ... but not when the stored parameters differ (e.g. bounds changed)
        path = os.path.join(d["run_dir"], "params.json")
        with open(path) as fh:
            saved = json.load(fh)
        with open(path, "w") as fh:
            json.dump(dict(saved, beam_width=saved["beam_width"] + 0.01), fh)
        other = _eval(HOLE_U, tmp, allow_cloud=False)
        assert other["optical_tidy3d_stop"] == "cloud_not_allowed"


def test_failure_after_start_is_charged_to_the_ledger():
    with fake_tidy3d(cost=0.4) as (td, bs), tempfile.TemporaryDirectory() as tmp:
        td.WEB.Batch.return_value.run.side_effect = OSError("connection lost")
        rec = _eval(HOLE_U, tmp, allow_cloud=True, max_credits=1.0)
        assert rec["status"] == "optical_failed" and "connection" in rec["error"]
        assert "optical_tidy3d_stop" not in rec       # a real failure: loop goes on
        with open(os.path.join(tmp, "credit_ledger.jsonl")) as fh:
            ledger = [json.loads(line) for line in fh]
        assert ledger[0]["status"] == "failed_after_start"
        assert ledger[0]["estimated_credits"] == 0.4 and ot.session_credits() == 0.4
        # an upload failure is not mistaken for a budget stop
        td.WEB.Batch.return_value.upload.side_effect = OSError("no network")
        rec = _eval([0.55, 0.6, 0.5, 0.6, 0.5], tmp, allow_cloud=True,
                    max_credits=0.5)
        assert rec["status"] == "optical_failed"


def test_rib_sweep_end_to_end_offline():
    with fake_tidy3d(cost=0.1) as (td, bs), tempfile.TemporaryDirectory() as tmp:
        _install_fake_analysis(bs, RIB_CELL["a"])
        rec = _eval(_rib_u0(), tmp, celltype="rib", allow_cloud=True,
                    max_credits=2.0)
        assert rec["status"] == "success", rec.get("error")
        d = rec["optical_tidy3d_detail"]
        assert d["gap_method"] == "guided_interval_gap" and d["n_sims"] == 13
        assert len(td.WEB.Batch.call_args.kwargs["simulations"]) == 13
        assert abs(rec["optical_gap"] - 35.0 / 192.5) < 1e-6   # 175-210 THz
        assert abs(d["x_point"]["f_lower_THz"] - 140.0) < 1e-6  # two lowest at X
        assert rec["optical_gap"] != d["x_point"]["normalized_gap"]
        assert rec["optical_min_resolvable_Hz"] > 0 and len(rec["_k_o"]) == 13
        assert rec["celltype"] == "rib" and rec["optical_engine"] == "tidy3d"
        json.dumps(rec, allow_nan=False)


def test_rib_tidy3d_gap_choice_has_no_centre_window():
    """Rib, score version 2: Gaussian-weighted choice instead of the window,
    for a sweep and for x_only; the hole path keeps its window."""
    a = RIB_CELL["a"]
    plan = ot.k_plan("rib", a, k_mode="sweep", n_k=13, f_target=F_T)
    k = np.array(plan.k_frac)
    flat = lambda cols: np.column_stack(                       # noqa: E731
        [np.full_like(k, f * 1e12) for f in cols])
    window = (93.4e12, 293.4e12)
    # two gaps: 150-185 THz (20.9 %, far) and 188-200 THz (6.2 %, on target)
    freqs = flat([130.0, 150.0, 185.0, 188.0, 200.0])
    out = {}
    gp, method = ot.score_resonances("rib", k, freqs, a, f_target=F_T,
                                     k_mode="sweep", window_hz=window,
                                     optical_tol_hz=2e12, out=out)
    assert method == "guided_interval_gap" and gp.found
    assert (round(gp.f_lower / 1e12), round(gp.f_upper / 1e12)) == (188, 200)
    assert abs(out["freq_penalty"] - np.exp(-((194e12 - F_T) / 2e12) ** 2)) < 1e-12
    assert len(out["candidates"]) == 3 and out["candidates"][0]["chosen"]
    # the old windowed call picks the widest gap inside +/- 30 %
    old, _ = ot.score_resonances("rib", k, freqs, a, f_target=F_T,
                                 k_mode="sweep", window_hz=window)
    assert round(old.f_center / 1e12, 1) == 167.5
    # the FDTD window is still the trusted range
    cut, _ = ot.score_resonances("rib", k, freqs, a, f_target=F_T, k_mode="sweep",
                                 window_hz=(93.4e12, 190e12), optical_tol_hz=2e12)
    assert cut.f_upper <= 190e12
    # x_only: the fundamental gap at X wherever it is centred (100-120 THz)
    row = np.array([[100e12, 120e12]])
    out = {}
    gx, method = ot.score_resonances("rib", [1.0], row, a, f_target=F_T,
                                     k_mode="x_only", optical_tol_hz=2e12, out=out)
    assert method == "optical_gap_tracked_x" and gx.found
    assert abs(gx.f_center - 110e12) < 1.0 and out["freq_penalty"] < 1e-30
    assert len(out["candidates"]) == 1
    assert not ot.score_resonances("rib", [1.0], row, a, f_target=F_T,
                                   k_mode="x_only")[0].found       # windowed
    assert not ot.score_resonances("hole", [1.0], row, 600e-9, f_target=F_T,
                                   k_mode="x_only")[0].found       # hole: unchanged
    # end to end (mocked tidy3d): the record carries choice, weight, version
    with fake_tidy3d(cost=0.1) as (td, bs), tempfile.TemporaryDirectory() as tmp:
        _install_fake_analysis(bs, a)
        rec = _eval(_rib_u0(), tmp, celltype="rib", allow_cloud=True, max_credits=2.0)
        assert rec["status"] == "success", rec.get("error")
        pen = np.exp(-((192.5e12 - F_T) / 2e12) ** 2)              # 175-210 THz
        assert abs(rec["optical_freq_penalty"] - pen) < 1e-9
        assert abs(rec["optical_gap_effective"] - 35.0 / 192.5 * pen) < 1e-9
        assert rec["optical_gap"] == rec["optical_gap_unpenalized"]
        assert rec["score_version"] == 2 and rec["optical_tol_thz"] == 2.0
        assert any(c["chosen"] for c in rec["optical_gap_candidates"])
        g_eff = rec["optical_gap_effective"]
        assert abs(rec["score"] - (g_eff - 10 * (0.2 - g_eff) ** 2)) < 1e-9
        json.dumps(rec, allow_nan=False)
        rx = _eval(_rib_u0(), tmp, celltype="rib", allow_cloud=True,
                   max_credits=2.0, k_mode="x_only")
        assert rx["status"] == "success" and rx["id"] != rec["id"]
        # two lowest guided at X: 140 / 175 THz, centre 157.5 THz: found, weight ~ 0
        assert abs(rx["optical_center_frequency"] - 157.5e12) < 1e6
        assert rx["optical_freq_penalty"] < 1e-30 and len(rx["optical_gap_candidates"]) == 1


def test_rib_loop_saves_everything_and_stops():
    """scripts/run_opt_comsol_rib.main_rib with optical: tidy3d and no opt-in."""
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol_rib as loop
    import database
    with fake_tidy3d() as (td, bs), tempfile.TemporaryDirectory() as tmp:
        cfg = dict(celltype="rib", n_init=5, n_iter=60, backend="random",
                   optical="tidy3d", mech="surrogate_stub",
                   tidy3d=dict(allow_cloud=False, max_credits=1.0,
                               results_dir=os.path.join(tmp, "t3d")),
                   out_json=os.path.join(tmp, "res.json"),
                   out_fig=os.path.join(tmp, "fig.png"),
                   db_path=os.path.join(tmp, "runs.sqlite3"),
                   skip_characterization=True)
        loop.main_rib(cfg, "unused")
        with open(cfg["out_json"]) as fh:
            results = json.load(fh)
        saved = database.load_completed_results(cfg["db_path"])
        assert td.WEB.Batch.call_count == 0
    statuses = [r["status"] for r in results]
    assert statuses[-1] == "optical_not_run" and len(results) < 60
    assert results[-1]["optical_tidy3d_stop"] == "cloud_not_allowed"
    assert set(statuses[:-1]) <= {"infeasible", "unfabricable"}
    assert len(saved) == len(results)                 # every candidate saved


def test_hole_comsol_loop_refuses_another_optical_engine():
    import subprocess
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "cfg.yaml")
        with open(path, "w") as fh:
            fh.write("optical: tidy3d\n")
        proc = subprocess.run(
            [sys.executable, os.path.join(_HERE, "..", "scripts",
                                          "run_opt_comsol.py"), "--config", path],
            capture_output=True, text=True)
    assert proc.returncode == 2 and "not available" in proc.stdout, (
        proc.returncode, proc.stdout[-300:], proc.stderr[-300:])


# ---- the real tidy3d package (only where it is installed), strictly offline ------

def test_web_api_is_imported_explicitly():
    """`import tidy3d` does not load the web API: td.web does not exist until
    tidy3d.web is imported (tidy3d 2.12). Nothing may rely on `td.web`."""
    import re
    root = os.path.join(_HERE, "..")
    with open(os.path.join(root, "diamond_bandStructure_sim.py")) as fh:
        src = fh.read()
    code = "\n".join(line.split("#")[0] for line in src.splitlines())
    assert re.search(r"^import tidy3d\.web as web\s*$", code, re.M)
    assert "web.Batch(" in code
    for path in ("diamond_bandStructure_sim.py",
                 os.path.join("src", "optical_tidy3d.py")):
        with open(os.path.join(root, path)) as fh:
            body = "\n".join(line.split("#")[0] for line in fh.read().splitlines())
        assert not re.search(r"\b(td|tidy3d)\.web\.\w+\(", body), path
    # the in-process fake has the real layout: no td.web attribute ...
    with fake_tidy3d() as (td, bs):
        try:
            td.web
        except AttributeError:
            pass
        else:
            raise AssertionError("the fake tidy3d exposes td.web")
        # ... and the module got its web API from the submodule
        assert bs.web is sys.modules["tidy3d.web"] is td.WEB
        with tempfile.TemporaryDirectory() as tmp:
            bs.run_cloud({"sim_1": object()}, tmp, "x", estimate_only=True,
                         verbose=False)
        assert td.WEB.Batch.call_count == 1
    # ... and so has the on-disk fake of the subprocess tests
    import subprocess
    with tempfile.TemporaryDirectory() as tmp:
        _write_fake_tidy3d_package(tmp)
        probe = ("import tidy3d as td\n"
                 "assert not hasattr(td, 'web'), 'td.web exists without import'\n"
                 "import tidy3d.web as web\n"
                 "assert web.Batch and td.web is web\n"
                 "print('ok')\n")
        proc = subprocess.run([sys.executable, "-c", probe], capture_output=True,
                              text=True, env=dict(os.environ, PYTHONPATH=tmp))
        assert proc.stdout.strip() == "ok", proc.stderr[-500:]


@contextlib.contextmanager
def _network_blocked():
    """Any attempt to open a connection fails loudly while the block runs."""
    import socket

    def refuse(*args, **kw):
        raise AssertionError("network access attempted in an offline test")

    saved = (socket.socket.connect, socket.socket.connect_ex,
             socket.create_connection)
    socket.socket.connect = refuse
    socket.socket.connect_ex = refuse
    socket.create_connection = refuse
    try:
        yield
    finally:
        (socket.socket.connect, socket.socket.connect_ex,
         socket.create_connection) = saved


def _real_bs():
    """The module on top of the REAL tidy3d, or None when it is not installed."""
    if not ot.have_tidy3d():
        print("  (skipped: tidy3d is not installed)")
        return None
    ot.reset_session()
    return ot.load_bs()


def test_real_tidy3d_api_matches_what_the_module_calls():
    """Names and signatures only: nothing is constructed that could reach the
    cloud (no web.Batch instance, no API key, no account call)."""
    import inspect
    with _network_blocked():
        bs = _real_bs()
        if bs is None:
            return
        try:
            import tidy3d as td
            import tidy3d.web as web
            assert bs.web is web and inspect.isclass(web.Batch)
            fields = getattr(web.Batch, "model_fields", None) or web.Batch.__fields__
            assert {"simulations", "folder_name", "verbose"} <= set(fields)
            assert list(inspect.signature(web.Batch.upload).parameters) == ["self"]
            assert "verbose" in inspect.signature(web.Batch.estimate_cost).parameters
            assert "path_dir" in inspect.signature(web.Batch.run).parameters
            assert list(inspect.signature(web.Batch.delete).parameters) == ["self"]
            assert bs.ResonanceFinder.__module__.startswith("tidy3d.plugins.resonance")
            rf = getattr(bs.ResonanceFinder, "model_fields", None) or bs.ResonanceFinder.__fields__
            assert "freq_window" in rf
            assert "signals" in inspect.signature(bs.ResonanceFinder.run).parameters
            assert td.Boundary.bloch(0.25).plus.bloch_vec == 0.25
            auto = inspect.signature(td.GridSpec.auto).parameters
            assert {"min_steps_per_wvl", "snapping_points"} <= set(auto)
            ps = getattr(td.PolySlab, "model_fields", None) or td.PolySlab.__fields__
            assert {"vertices", "slab_bounds", "axis", "sidewall_angle",
                    "reference_plane"} <= set(ps)
            assert "fname" in inspect.signature(td.SimulationData.from_file).parameters
            assert abs(bs.C_0 - C0 * 1e6) < 1.0              # um/s
            # the console switch of optical_tidy3d.tidy3d_console
            logmod = ot._td_log_module()
            assert logmod is not None
            assert "stderr" in inspect.signature(logmod.set_logging_console).parameters
            assert "level" in inspect.signature(logmod.set_logging_level).parameters
            before = logmod.log.handlers.get("console")
            with ot.tidy3d_console(verbose=False):
                pass
            assert logmod.log.handlers.get("console") is before
        finally:
            ot.reset_session()


def _material_depth(sim, th, x, y):
    """Top and depth below the slab of the dielectric under (x, y), read from
    the structures of a real td.Simulation (later structures override)."""
    zz = np.linspace(-th / 2 - 0.45, th / 2 + 0.05, 1001)
    solid = np.zeros(zz.size, bool)
    for s in sim.structures:
        hit = np.asarray(s.geometry.inside(np.full(zz.size, x), np.full(zz.size, y), zz))
        solid = np.where(hit, s.medium.permittivity > 1.5, solid)
    if not solid.any():
        return None, -1.0
    return zz[solid].max(), -th / 2 - zz[solid].min()


def test_real_tidy3d_builds_locally():
    """Hole and rib, rect and wedge, built and validated by the real tidy3d
    (pydantic models, grid, symmetry), and the dielectric they contain compared
    point by point with the repo's geometry. No upload: mode dry_run, the cloud
    function replaced by one that fails, the network blocked."""
    from geometry import rib_roof_distance, wedge_depth
    with _network_blocked():
        bs = _real_bs()
        if bs is None:
            return
        real_cloud = bs.run_cloud
        bs.run_cloud = mock.MagicMock(side_effect=AssertionError("cloud reached"))
        hole_u = phys_to_u(dict(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9))
        rng = np.random.default_rng(5)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                for celltype, u, n in (("hole", hole_u, 1), ("rib", _rib_u0(), 13)):
                    for xs in ("rect", "wedge"):
                        rec = evaluate_candidate(
                            u, celltype=celltype, optical_backend="tidy3d",
                            mech_backend="surrogate_stub", xsect=xs,
                            tidy3d_options=dict(results_dir=tmp, mode="dry_run"))
                        d = rec["optical_tidy3d_detail"]
                        assert rec["optical_tidy3d_stop"] == "dry_run", rec.get("error")
                        assert d["n_sims"] == n and d["grid_cells"] > 1e4
                        assert d["time_steps"] > 1000 and d["index"] == 2.386
                        assert d["symmetry"] == [0, -1, 1 if xs == "rect" else 0]
                        with open(os.path.join(d["run_dir"], "grid_report.txt")) as fh:
                            report = fh.read()
                        assert "n = 2.3860" in report
                        assert ("wedge" in report) == (xs == "wedge")

                        # the structures themselves
                        g = u_to_geometry(u, celltype=celltype)
                        p = bs.BandParams(Nk=d["Nk"], **d["band_params"])
                        sims = bs.build_simulations(p, x_only=(n == 1))
                        assert len(sims) == (1 if n == 1 else d["Nk"] + 1)
                        sim = sims["sim_%d" % d["Nk"]]
                        assert bs.bloch_k(sim) == 0.5                # X point
                        assert tuple(sim.symmetry) == tuple(d["symmetry"])
                        assert abs(sim.size[0] - g.a * 1e6) < 1e-9   # one period
                        z_src = [abs(s.center[2]) for s in sim.sources]
                        assert (max(z_src) > 0) == (xs == "wedge")
                        th = (g.th if celltype == "rib" else g.t) * 1e6
                        tan30 = np.tan(np.radians(30.0))
                        worst = 0.0
                        for _ in range(60):
                            x = rng.uniform(-g.a / 2, g.a / 2) * 1e6
                            y = rng.uniform(-g.w / 2, g.w / 2) * 1.05e6
                            if celltype == "rib":
                                dist = rib_roof_distance(g, x * 1e-6, y * 1e-6) * 1e6
                                want = dist * tan30
                            else:
                                dist = min(g.w / 2 * 1e6 - abs(y),
                                           (np.hypot(x / (g.hx * 1e6), y / (g.hy * 1e6)) - 1)
                                           * min(g.hx, g.hy) * 1e6)
                                want = (g.w / 2 * 1e6 - abs(y)) * tan30
                            if abs(dist) < 0.012:
                                continue                 # too close to an edge
                            top, depth = _material_depth(sim, th, x, y)
                            assert (top is not None) == (dist > 0), (celltype, xs, x, y)
                            if top is None:
                                continue
                            assert abs(top - th / 2) < 1e-3
                            worst = max(worst, abs(depth - (want if xs == "wedge" else 0.0)))
                        # sampling step 0.5 nm; the rib wedge is a 12-layer staircase
                        layer = (d["wedge_depths_m"]["max"] * 1e6 / 24
                                 if (celltype, xs) == ("rib", "wedge") else 0.0)
                        assert worst <= layer + 1.5e-3, (celltype, xs, worst, layer)
                        if (celltype, xs) == ("hole", "wedge"):
                            assert abs(d["wedge_depths_m"]["beam"] - wedge_depth(g.w)) < 1e-15
            assert bs.run_cloud.call_count == 0
        finally:
            bs.run_cloud = real_cloud
            ot.reset_session()


def test_real_tidy3d_run_and_analysis_with_a_fake_batch():
    """The whole `mode: run` path on real tidy3d objects: run_cloud, the real
    ResonanceFinder, SimulationData, scoring, ledger. The ONLY replaced piece
    is web.Batch (a local class that fabricates the monitor data); the real
    one is never constructed and the network is blocked."""
    with _network_blocked():
        bs = _real_bs()
        if bs is None:
            return
        import tidy3d as td
        f_true, q_true = (150e12, 195e12, 260e12), (5e3, 4e3, 60.0)
        made = []

        class FakeBatch:
            def __init__(self, simulations, folder_name="default", verbose=True):
                self.simulations = simulations
                made.append(self)

            def upload(self):
                self.uploaded = True

            def estimate_cost(self, verbose=True):
                return 0.05

            def delete(self):
                self.deleted = True

            def run(self, path_dir="."):
                out, rng = {}, np.random.default_rng(0)
                for name, sim in self.simulations.items():
                    assert isinstance(sim, td.Simulation)
                    datas = []
                    for m in sim.monitors:
                        t = np.linspace(m.start, sim.run_time, 5000)
                        sig = sum(a * np.exp(-np.pi * f * t / q) * np.cos(2 * np.pi * f * t + ph)
                                  for f, q, a, ph in zip(f_true, q_true, rng.uniform(5, 10, 3),
                                                         rng.uniform(0, 6, 3)))
                        arr = td.ScalarFieldTimeDataArray(
                            sig[None, None, None, :],
                            coords=dict(x=[m.center[0]], y=[m.center[1]],
                                        z=[m.center[2]], t=t))
                        datas.append(td.FieldTimeData(
                            monitor=m, Ey=arr, symmetry=(0, 0, 0),
                            symmetry_center=(0, 0, 0),
                            grid_expanded=sim.discretize_monitor(m)))
                    out[name] = td.SimulationData(simulation=sim, data=tuple(datas))
                    out[name].to_file(os.path.join(path_dir, name + ".hdf5"))
                return out

        real_batch = bs.web.Batch
        bs.web.Batch = FakeBatch
        try:
            with tempfile.TemporaryDirectory() as tmp:
                rec = _eval(HOLE_U, tmp, allow_cloud=True, max_credits=1.0)
                assert rec["status"] == "success", rec.get("error")
                d = rec["optical_tidy3d_detail"]
                # a = 700 nm: light line at X 214 THz, so 150 and 195 THz are the
                # two guided resonances and 260 THz is leaky
                assert abs(d["x_point"]["f_lower_THz"] - 150.0) < 0.05
                assert abs(d["x_point"]["f_upper_THz"] - 195.0) < 0.05
                assert abs(rec["optical_gap"] - 45.0 / 172.5) < 1e-3
                assert d["gap_method"] == "optical_gap_tracked_x"
                assert len(made) == 1 and d["estimated_credits"] == 0.05
                assert len(rec["_freqs_o"][0]) == 3
                json.dumps(rec, allow_nan=False)
                with open(os.path.join(tmp, "credit_ledger.jsonl")) as fh:
                    assert json.loads(fh.readline())["estimated_credits"] == 0.05
                # the downloaded file is found and re-analyzed, no second batch
                again = _eval(HOLE_U, tmp, allow_cloud=False)
                assert again["status"] == "success", again.get("error")
                assert again["optical_tidy3d_detail"]["reused_download"] is True
                assert abs(again["optical_gap"] - rec["optical_gap"]) < 1e-12
                assert len(made) == 1
        finally:
            bs.web.Batch = real_batch
            ot.reset_session()


# ---- console output: Tidy3D progress on stderr, stdout stays JSON --------------
# The run_one.py tests below start a subprocess whose `tidy3d` is the FAKE
# package written by _write_fake_tidy3d_package (first on PYTHONPATH): its
# web.Batch only prints, it has no network code at all.

_FAKE_TD_INIT = r'''
"""FAKE tidy3d for tests/test_optical_tidy3d.py: offline, prints only."""
import atexit, json, os, sys, types
from unittest import mock as _mock

from .log import log, set_logging_console, get_logging_console  # noqa: F401
# as in the real package, the attribute `log` (the Logger) hides the submodule

_m = _mock.MagicMock(name="tidy3d")
inf = float("inf")
constants = types.SimpleNamespace(C_0=299792458.0e6)


def __getattr__(name):                  # Medium, PolySlab, GridSpec, ...
    # `web` is a submodule and, as in the real package, NOT loaded by
    # `import tidy3d`: td.web fails until someone imports tidy3d.web.
    if name.startswith("__") or name == "web":
        raise AttributeError("module 'tidy3d' has no attribute %r" % name)
    return getattr(_m, name)


class Boundary:
    def __init__(self, k=0.0):
        self.plus = types.SimpleNamespace(bloch_vec=k)

    @classmethod
    def bloch(cls, k):
        return cls(k)

    @classmethod
    def pml(cls):
        return cls()


class BoundarySpec:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


class Simulation:
    num_cells, num_time_steps = 1.0e6, 1000

    def __init__(self, **kw):
        self.__dict__.update(kw)


@atexit.register
def _write_state():                     # console state left behind by the run
    path = os.environ.get("FAKE_TD_STATE")
    if path:
        h = log.handlers.get("console")
        with open(path, "w") as fh:
            json.dump(dict(console_stderr=h.console.stderr if h else None,
                           level=h.level if h else None,
                           switched_to_stderr=log.switched_to_stderr), fh)
'''

_FAKE_TD_WEB = r'''
"""FAKE tidy3d.web: a submodule, as in the real package. Prints only."""
from .log import log, get_logging_console


class _SimData:
    def __init__(self, sim):
        self.simulation, self.data = sim, ()


class Batch:
    def __init__(self, simulations, folder_name="default", verbose=True):
        self.simulations, self.verbose = simulations, verbose

    def _say(self, text):
        if self.verbose:
            get_logging_console().log(text)

    def upload(self):
        self._say("FAKE-TD console: Created task, uploading 100%")
        log.info("FAKE-TD log.info: upload finished")
        log.warning("FAKE-TD log.warning: something to know")
        log.error("FAKE-TD log.error: always shown")
        if self.verbose:
            print("FAKE-TD raw print() during upload")

    def estimate_cost(self, verbose=True):
        return 0.25

    def delete(self):
        self._say("FAKE-TD console: drafts deleted")

    def run(self, path_dir):
        for status in ("queued", "preprocess", "running 100%", "postprocess",
                       "success"):
            self._say("FAKE-TD console: status = " + status)
        self._say("FAKE-TD console: downloading data 100%")
        return {name: _SimData(sim) for name, sim in self.simulations.items()}


def run(*args, **kw):
    raise AssertionError("fake tidy3d: web.run must not be called")
'''

_FAKE_TD_LOG = r'''
"""FAKE tidy3d.log: the API optical_tidy3d.tidy3d_console relies on."""
import sys

_LEVELS = dict(DEBUG=10, INFO=20, WARNING=30, ERROR=40)


class Console:                          # stands in for rich.console.Console
    def __init__(self, stderr=False):
        self.stderr = stderr

    @property
    def file(self):                     # looked up at write time, as rich does
        return sys.stderr if self.stderr else sys.stdout

    def log(self, text):
        print(text, file=self.file)


class LogHandler:
    def __init__(self, console, level):
        self.console, self.level = console, level


class Logger:
    def __init__(self):
        self.handlers, self.switched_to_stderr = {}, 0

    def _log(self, level, text):
        h = self.handlers.get("console")
        if h is not None and level >= h.level:
            h.console.log(text)

    def info(self, text):
        self._log(20, text)

    def warning(self, text):
        self._log(30, text)

    def error(self, text):
        self._log(40, text)


log = Logger()


def set_logging_level(level="WARNING"):
    if "console" in log.handlers:
        log.handlers["console"].level = _LEVELS[level]


def set_logging_console(stderr=False):
    level = log.handlers["console"].level if "console" in log.handlers else 30
    log.switched_to_stderr += bool(stderr)
    log.handlers["console"] = LogHandler(Console(stderr=stderr), level)


def get_logging_console():
    if "console" not in log.handlers:
        set_logging_console()
    return log.handlers["console"].console


set_logging_console()                   # default: stdout, WARNING
'''

_FAKE_TD_RESONANCE = r'''
import numpy as np


class _A(np.ndarray):
    def to_numpy(self):
        return np.asarray(self)


class _Data:                            # the bit of xarray the module uses
    def __init__(self):
        self.freq = np.array([175e12, 210e12, 380e12]).view(_A)
        self.Q = np.full(3, 1e4).view(_A)
        self.amplitude = np.full(3, 10.0).view(_A)
        self.error = np.full(3, 1e-3).view(_A)

    def where(self, cond, drop=True):
        return self


class ResonanceFinder:
    def __init__(self, freq_window=None):
        pass

    def run(self, signals=None):
        return _Data()
'''


def _write_fake_tidy3d_package(root):
    for rel, text in (("tidy3d/__init__.py", _FAKE_TD_INIT),
                      ("tidy3d/log.py", _FAKE_TD_LOG),
                      ("tidy3d/web.py", _FAKE_TD_WEB),
                      ("tidy3d/plugins/__init__.py", ""),
                      ("tidy3d/plugins/resonance.py", _FAKE_TD_RESONANCE)):
        path = os.path.join(root, *rel.split("/"))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)


def _run_one_subprocess(tmp, quiet=False, **tidy):
    """scripts/run_one.py (hole, optical: tidy3d, mech: surrogate_stub) in a
    subprocess with the fake tidy3d. Returns (record, stdout, stderr, state)."""
    import subprocess
    import yaml
    fake_root = os.path.join(tmp, "fake_site")
    if not os.path.isdir(fake_root):
        _write_fake_tidy3d_package(fake_root)
    state = os.path.join(tmp, "fake_td_state.json")
    if os.path.exists(state):
        os.remove(state)
    cfg = dict(u=HOLE_U, optical="tidy3d", mech="surrogate_stub", quiet=quiet,
               tidy3d=dict(results_dir=os.path.join(tmp, "td"), **tidy))
    path = os.path.join(tmp, "run_one.yaml")
    with open(path, "w") as fh:
        yaml.safe_dump(cfg, fh)
    env = dict(os.environ, PYTHONPATH=fake_root, FAKE_TD_STATE=state,
               OMC_DB_PATH=os.path.join(tmp, "results.sqlite"),
               OMC_DB_PATH_RIB=os.path.join(tmp, "results_rib.sqlite"))
    proc = subprocess.run(
        [sys.executable, os.path.join(_HERE, "..", "scripts", "run_one.py"),
         "--config", path], capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr[-2000:]
    rec = json.loads(proc.stdout)       # stdout is PURE JSON
    assert os.path.exists(state), "the fake tidy3d was not the one imported"
    with open(state) as fh:
        return rec, proc.stdout, proc.stderr, json.load(fh)


def _assert_console_restored(state):
    assert state["switched_to_stderr"] >= 1     # the console was re-routed
    assert state["console_stderr"] is False     # ... and put back afterwards
    assert state["level"] == 30                 # WARNING, as before


def test_run_one_tidy3d_progress_goes_to_stderr():
    with tempfile.TemporaryDirectory() as tmp:
        rec, out, err, state = _run_one_subprocess(
            tmp, allow_cloud=True, max_credits=1.0)
        assert rec["status"] == "success", rec.get("error")
        assert rec["optical_engine"] == "tidy3d"
        assert abs(rec["optical_gap"] - 35.0 / 192.5) < 1e-6
        assert rec["optical_tidy3d_detail"]["estimated_credits"] == 0.25
        # Tidy3D's own output: console lines, log.info, raw print()
        for text in ("FAKE-TD console: Created task, uploading 100%",
                     "FAKE-TD console: status = queued",
                     "FAKE-TD console: status = running 100%",
                     "FAKE-TD console: downloading data 100%",
                     "FAKE-TD log.info: upload finished",
                     "FAKE-TD log.error: always shown",
                     "FAKE-TD raw print() during upload",
                     # our own step lines, through the on_stage reporter
                     "[2/4] optical (tidy3d) ... running",
                     "[2/4] optical: building 1 simulation(s)",
                     "[2/4] optical: plan written: ",
                     "[2/4] optical: uploading 1 simulation(s)",
                     "[2/4] optical: estimated maximum cost of 1 "
                     "simulation(s): 0.250 FlexCredits",
                     "[2/4] optical: running 1 task(s)",
                     "[2/4] optical: loading task 1/1",
                     "[2/4] optical: analyzing resonances",
                     "[2/4] optical ... done"):
            assert text in err, (text, err)
        assert err.count("estimated maximum cost") == 1     # not printed twice
        order = [err.index(t) for t in ("building 1", "uploading 1 sim",
                                        "FAKE-TD console: Created",
                                        "estimated maximum cost", "running 1 task",
                                        "status = queued", "loading task 1/1",
                                        "analyzing resonances")]
        assert order == sorted(order), err
        assert "FAKE-TD" not in out and "FlexCredits" not in out.split("{", 1)[0]
        assert out.lstrip().startswith("{")
        _assert_console_restored(state)


def test_run_one_quiet_silences_tidy3d_progress():
    with tempfile.TemporaryDirectory() as tmp:
        rec, out, err, state = _run_one_subprocess(
            tmp, quiet=True, allow_cloud=True, max_credits=1.0)
        assert rec["status"] == "success", rec.get("error")
        assert "FAKE-TD" not in out
        for text in ("FAKE-TD console", "FAKE-TD log.info", "FAKE-TD raw",
                     "FAKE-TD log.warning", "optical", "FlexCredits"):
            assert text not in err, (text, err)
        assert "FAKE-TD log.error: always shown" in err     # errors remain
        _assert_console_restored(state)


def test_run_one_dry_run_and_estimate_show_their_steps():
    with tempfile.TemporaryDirectory() as tmp:
        rec, out, err, state = _run_one_subprocess(tmp, mode="dry_run")
        assert rec["optical_tidy3d_stop"] == "dry_run"
        for text in ("[2/4] optical: building 1 simulation(s)",
                     "[2/4] optical: plan written: ", "FAILED", "dry run"):
            assert text in err, (text, err)
        assert "uploading" not in err and "FAKE-TD" not in err + out
        _assert_console_restored(state)

        rec, out, err, state = _run_one_subprocess(
            tmp, mode="estimate", allow_cloud=True)
        assert rec["optical_tidy3d_stop"] == "estimate_only"
        assert rec["optical_tidy3d_detail"]["estimated_credits"] == 0.25
        for text in ("[2/4] optical: uploading 1 draft(s)",
                     "FAKE-TD console: Created task, uploading 100%",
                     "[2/4] optical: estimated maximum cost of 1 "
                     "simulation(s): 0.250 FlexCredits",
                     "FAKE-TD console: drafts deleted",
                     "[2/4] optical: drafts deleted, nothing was started"):
            assert text in err, (text, err)
        assert "running 1 task" not in err and "status = queued" not in err
        assert "FAKE-TD" not in out
        _assert_console_restored(state)

        rec, out, err, state = _run_one_subprocess(tmp, quiet=True,
                                                   mode="dry_run")
        assert rec["optical_tidy3d_stop"] == "dry_run" and err.strip() == "", err


def test_stage_callback_gets_the_tidy3d_steps_in_process():
    """Both cell types send "info" events; verbose=False sends none; the
    Batch gets verbose from the caller; nothing lands on stdout."""
    import io
    with fake_tidy3d(cost=0.1) as (td, bs), tempfile.TemporaryDirectory() as tmp:
        _install_fake_analysis(bs, RIB_CELL["a"])
        for verbose in (True, False):
            events, out = [], io.StringIO()
            with contextlib.redirect_stdout(out), \
                    contextlib.redirect_stderr(io.StringIO()):
                rec = evaluate_candidate(
                    _rib_u0(), celltype="rib", optical_backend="tidy3d",
                    mech_backend="surrogate_stub", verbose=verbose,
                    on_stage=lambda *a: events.append(a),
                    tidy3d_options=dict(results_dir=os.path.join(tmp, str(verbose)),
                                        allow_cloud=True, max_credits=5.0))
            assert rec["status"] == "success", rec.get("error")
            assert out.getvalue() == ""
            assert td.WEB.Batch.call_args.kwargs["verbose"] is verbose
            info = [a[2] for a in events if a[:2] == ("optical", "info")]
            if not verbose:
                assert info == []
                continue
            assert info[0].startswith("building 13 simulation(s)")
            assert any(m.startswith("uploading 13 simulation(s)") for m in info)
            assert any("0.100 FlexCredits" in m for m in info)
            assert sum(m.startswith("loading task") for m in info) == 13
            assert any(m.startswith("loading task 13/13") for m in info)
            assert info[-1].startswith("analyzing resonances")
        # the command line tool's defaults are unchanged: verbose Batch, print
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            bs.run_cloud({"sim_1": object()}, tmp, "cli", "default",
                         bs.Budget(None), estimate_only=True)
        assert td.WEB.Batch.call_args.kwargs["verbose"] is True
        assert "estimated maximum cost of 1 simulation(s): 0.100" in out.getvalue()


def test_tidy3d_console_restores_the_global_state():
    import io
    import types
    with tempfile.TemporaryDirectory() as tmp:
        _write_fake_tidy3d_package(tmp)
        spec = __import__("importlib.util").util.spec_from_file_location(
            "tidy3d.log", os.path.join(tmp, "tidy3d", "log.py"))
        logmod = __import__("importlib.util").util.module_from_spec(spec)
        spec.loader.exec_module(logmod)
    saved = {n: sys.modules.get(n) for n in ("tidy3d", "tidy3d.log")}
    sys.modules["tidy3d"] = types.ModuleType("tidy3d")
    sys.modules["tidy3d.log"] = logmod
    try:
        before = logmod.log.handlers["console"]
        for verbose, level in ((True, 20), (False, 40)):
            err = io.StringIO()
            try:
                with contextlib.redirect_stderr(err), ot.tidy3d_console(verbose):
                    h = logmod.log.handlers["console"]
                    assert h is not before and h.console.stderr and h.level == level
                    logmod.get_logging_console().log("to the console")
                    print("plain print")
                    raise KeyError("boom")          # restored on errors too
            except KeyError:
                pass
            assert err.getvalue() == "to the console\nplain print\n"
            assert logmod.log.handlers["console"] is before
            assert before.level == 30 and before.console.stderr is False
        # no console handler before -> none afterwards
        del logmod.log.handlers["console"]
        with ot.tidy3d_console(True):
            assert logmod.log.handlers["console"].console.stderr
        assert "console" not in logmod.log.handlers
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
    # without tidy3d (or with an unknown API): plain stdout redirection
    sys.modules["tidy3d"] = None        # import tidy3d.log -> ImportError
    try:
        err, out = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), \
                ot.tidy3d_console(False):
            print("still redirected")
        assert out.getvalue() == "" and err.getvalue() == "still redirected\n"
    finally:
        if saved["tidy3d"] is None:
            sys.modules.pop("tidy3d", None)
        else:
            sys.modules["tidy3d"] = saved["tidy3d"]


TESTS = [v for k, v in sorted(globals().items(), key=lambda kv: kv[1].__code__.co_firstlineno
                              if callable(kv[1]) and hasattr(kv[1], "__code__") else 0)
         if k.startswith("test_") and callable(v)]

if __name__ == "__main__":
    for test in TESTS:
        print(test.__name__)
        test()
    print("all %d tests passed" % len(TESTS))
