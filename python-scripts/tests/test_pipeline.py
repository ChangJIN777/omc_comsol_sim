"""Minimal tests that run with numpy only (no MPB/COMSOL)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import json
import sqlite3
import tempfile

from geometry import (u_to_geometry, check_feasibility, check_fabricability,
                      phys_to_u, n_dim, rib_search_vars, RibGeometry)
from bandgap import largest_gap, gap_near_frequency, interval_gap, trusted_window
from optical_surrogate import optical_gap_surrogate
from objective import evaluate_candidate


def test_geometry_mapping_and_bounds():
    g = u_to_geometry([0.0, 0.0, 0.0, 0.0])
    assert g.a > 0 and g.w > 0 and g.hx > 0 and g.hy > 0
    g2 = u_to_geometry([1, 1, 1, 1])
    assert g2.a > g.a and g2.w > g.w


def test_feasibility_rejects_oversized_hole():
    # force a hole bigger than the beam -> infeasible
    g = u_to_geometry([0.0, 0.0, 1.0, 1.0])
    ok, reasons = check_feasibility(g)
    assert ok is False and len(reasons) > 0


def test_largest_gap_detects_clean_gap():
    k = np.linspace(0, np.pi, 11)
    b1 = 0.20 + 0.01 * np.sin(k)
    b2 = 0.30 + 0.01 * np.cos(k)
    g = largest_gap(np.column_stack([b1, b2]))
    assert g.found and g.normalized_gap > 0.3


def test_surrogate_returns_telecom_gap():
    g = u_to_geometry([0.5, 0.6, 0.5, 0.6])
    r = optical_gap_surrogate(g)
    assert r["found"]
    assert 0.0 < r["normalized_gap"] < 1.0
    assert 1000 < r["center_wavelength_nm"] < 2400


def test_objective_end_to_end_surrogate():
    rec = evaluate_candidate([0.5, 0.6, 0.5, 0.6],
                             optical_backend="surrogate",
                             mech_backend="surrogate_stub")
    assert rec["status"] in ("success",)
    assert rec["optical_gap"] > 0
    assert "score" in rec


# ---- rib cell (celltype="rib"); no COMSOL needed ----------------------------
RIB_CELL = dict(a=570e-9, w=1747e-9, t=391e-9, s=254e-9, th=250e-9,
                r1=40e-9, r2=40e-9)
FAB_LIMITS = dict(min_feature=50e-9, min_gap=50e-9, min_thickness=50e-9,
                  max_etch_aspect=5.0, max_solid_aspect=5.0,
                  max_rib_slenderness=10.0, max_rib_overhang=10.0,
                  max_release_width=None)


def _rib(**change):
    return RibGeometry(**{**RIB_CELL, **change})


def _fab_codes(g, **limits):
    bounds = {"rib": {"fabrication": {**FAB_LIMITS, **limits}}}
    ok, failures = check_fabricability(g, bounds)
    assert ok == (len(failures) == 0)
    return [f["code"] for f in failures]


def _rib_u0():
    return phys_to_u(RIB_CELL, celltype="rib")


def test_rib_u_mapping_roundtrip():
    names = rib_search_vars()
    assert n_dim("hole") == 5 and n_dim("rib") == len(names)
    g = u_to_geometry(_rib_u0(), celltype="rib")
    assert set(g.as_dict()) == set(RIB_CELL)
    assert all(abs(g.as_dict()[k] - RIB_CELL[k]) < 1e-12 for k in names)
    lo = u_to_geometry([0.0] * len(names), celltype="rib").as_dict()
    hi = u_to_geometry([1.0] * len(names), celltype="rib").as_dict()
    assert all(hi[k] > lo[k] for k in names)
    for bad in ([0.5] * (len(names) + 1), [0.5] * (len(names) - 1)):
        try:
            u_to_geometry(bad, celltype="rib")
        except ValueError:
            pass
        else:
            raise AssertionError("wrong-length rib u accepted")


def test_rib_template_validity_rules():
    assert check_feasibility(_rib()) == (True, [])
    cases = {"GEOM_RIBS_MERGE": dict(t=570e-9),
             "GEOM_NO_RIB": dict(s=1747e-9),
             "GEOM_TIP_FILLET": dict(r2=200e-9),
             "GEOM_ROOT_FILLET": dict(r1=95e-9),
             "GEOM_FILLETS_MEET": dict(w=400e-9),
             "GEOM_NONPOSITIVE": dict(r1=0.0)}
    for code, change in cases.items():
        ok, reasons = check_feasibility(_rib(**change))
        assert ok is False and any(r.startswith(code) for r in reasons), code


def test_rib_fabricability_rules():
    assert _fab_codes(_rib()) == []
    assert "FAB_MIN_RIB_WIDTH" in _fab_codes(_rib(t=45e-9))
    assert "FAB_MIN_SPINE_WIDTH" in _fab_codes(_rib(s=40e-9))
    assert "FAB_MIN_RIB_ARM" in _fab_codes(_rib(w=340e-9))
    assert "FAB_MIN_GAP" in _fab_codes(_rib(t=530e-9))
    assert "FAB_MIN_THICKNESS" in _fab_codes(_rib(th=40e-9))
    assert _fab_codes(_rib(th=1000e-9)) == ["FAB_ETCH_ASPECT"]
    assert _fab_codes(_rib(t=70e-9)) == ["FAB_RIB_SLENDERNESS"]
    assert _fab_codes(_rib(th=60e-9)) == ["FAB_RIB_OVERHANG"]
    assert _fab_codes(_rib(s=45e-9, t=600e-9, a=900e-9)) == [
        "FAB_MIN_SPINE_WIDTH", "FAB_SOLID_ASPECT"]
    # thresholds are configurable; null switches a rule off
    assert _fab_codes(_rib(t=70e-9), max_rib_slenderness=None) == []
    assert _fab_codes(_rib(), max_release_width=300e-9) == ["FAB_RELEASE_WIDTH"]
    assert _fab_codes(_rib(), min_gap=200e-9) == ["FAB_MIN_GAP"]
    # the hole cell has no rules here
    assert check_fabricability(u_to_geometry([0.5, 0.6, 0.5, 0.6])) == (True, [])


def test_hole_cell_unchanged_by_rib_support():
    g = u_to_geometry([0.5, 0.6, 0.5, 0.6])
    assert set(g.as_dict()) == {"a", "w", "t", "hx", "hy"}
    rec = evaluate_candidate([0.5, 0.6, 0.5, 0.6],
                             optical_backend="surrogate",
                             mech_backend="surrogate_stub")
    assert "celltype" not in rec and set(rec["params"]) == set(g.as_dict())


def _fake_mech_bands():
    k = np.linspace(0, 1, 9)
    sym = np.column_stack([5e9 + 0.5e9 * k, 9e9 - 0.5e9 * k])
    anti = np.column_stack([4e9 + 0 * k, 7e9 + 0.2e9 * k, 12e9 + 0 * k])
    return k, sym, anti


def test_interval_gap_single_and_union():
    _, sym, anti = _fake_mech_bands()
    g = interval_gap([sym], 8e9)
    assert g.found and abs(g.f_lower - 5.5e9) < 1 and abs(g.f_upper - 8.5e9) < 1
    g = interval_gap([sym, anti], 8e9)       # 5.5-7.0 (24 %) beats 7.2-8.5
    assert g.found and abs(g.f_lower - 5.5e9) < 1 and abs(g.f_upper - 7.0e9) < 1
    assert trusted_window([sym, anti]) == (5.5e9, 8.5e9)
    missing = sym.copy()
    missing[3, :] = np.nan                   # a k-point with no modes
    assert not interval_gap([missing], 8e9).found
    assert not interval_gap([], 8e9).found


def _fake_optical_bands(k_frac):
    """6 'modes' per k: pseudo-mode, two guided bands below a 175-210 THz gap,
    an air band above it, two radiation-like modes above the light line."""
    from comsol_optics import C0
    k = np.asarray(k_frac, dtype=float)
    x = (k - k[0]) / (1.0 - k[0])
    thz = 1e12
    rad = C0 * k / (2 * 570e-9) * 1.05
    return np.column_stack([np.full_like(k, 2 * thz), (120 + 20 * x) * thz,
                            (150 + 25 * x ** 2) * thz, (230 - 20 * x ** 2) * thz,
                            rad, 1.2 * rad])


def test_guided_interval_gap_on_synthetic_bands():
    from comsol_optics import guided_interval_gap, rib_k_grid, optical_gap_tracked
    assert callable(optical_gap_tracked)     # hole function, moved here
    a = 570e-9
    k = rib_k_grid(a, 13)
    assert 0.4 < k[0] < 0.5 and k[-1] == 1.0
    assert abs(rib_k_grid(300e-9, 13)[0] - 0.244) < 0.01    # small a: start lower
    bands = _fake_optical_bands(k)
    g = guided_interval_gap(k, bands, a)
    assert g.found and abs(g.f_lower - 175e12) < 1e6 and abs(g.f_upper - 210e12) < 1e6
    assert abs(g.normalized_gap - 35.0 / 192.5) < 1e-9
    # a steep band crossing the gap leaves only unresolvable intervals
    x = (k - k[0]) / (1.0 - k[0])
    steep = np.column_stack([bands, (130 + 110 * x) * 1e12])
    assert not guided_interval_gap(k, steep, a).found
    # solver stopped below the light line at X: upper side untrusted
    cut = bands.copy()
    cut[-1, 3:] = np.nan
    assert not guided_interval_gap(k, cut, a).found
    # target above the zone-edge light line (a = 900 nm): no guided gap
    assert not guided_interval_gap(rib_k_grid(900e-9, 13), bands, 900e-9).found


def test_rib_objective_gates_and_no_solver():
    kw = dict(celltype="rib", optical_backend="none", mech_backend="surrogate_stub")
    rec = evaluate_candidate(_rib_u0(), **kw)
    assert rec["status"] == "no_solver" and rec["celltype"] == "rib"
    assert rec["optical_gap"] is None and rec["score"] == 0.0
    assert rec["mech_gap_mode"] == "complete"
    names = rib_search_vars()
    # t = a: the template cannot build it
    u_bad = phys_to_u({**RIB_CELL, "t": RIB_CELL["a"]}, celltype="rib")
    rec = evaluate_candidate(u_bad, **kw)
    assert rec["status"] == "infeasible" and rec["score"] == -10.0
    assert any(r.startswith("GEOM_") for r in rec["reasons"])
    # buildable but too slender for the shipped fabrication limits (if any)
    u_fab = phys_to_u({**RIB_CELL, "t": 100e-9, "w": 2900e-9}, celltype="rib")
    rec = evaluate_candidate(u_fab, **kw)
    if rec["status"] == "unfabricable":      # depends on bounds.yaml thresholds
        assert rec["score"] == -10.0 and rec["fab_failures"]
        assert all(f["code"].startswith("FAB_") for f in rec["fab_failures"])
        assert len(rec["reasons"]) == len(rec["fab_failures"])
    assert len(rec["u"]) == len(names)
    for backend in ("surrogate", "mpb"):     # hole-only optical backends
        try:
            evaluate_candidate(_rib_u0(), celltype="rib",
                               optical_backend=backend,
                               mech_backend="surrogate_stub")
        except ValueError:
            pass
        else:
            raise AssertionError(f"{backend} accepted for the rib cell")


def test_rib_objective_with_fake_comsol():
    import acoustic_comsol
    import objective_rib
    k, sym, anti = _fake_mech_bands()
    bands = {"mech sym": sym, "mech antisym": anti}
    calls = []

    def fake_mech(g, *, n_k=15, n_bands=12, template=None, study="mech sym", **kw):
        calls.append((study, os.path.basename(template), n_bands,
                      type(g).__name__))
        return dict(freqs_hz=bands[study], k_frac=k, a=g.a)

    def fake_opt(g, study, k_frac, n_bands, f_target, **kw):
        # kw: xsect / info, added with the cross-section switch
        calls.append((study, len(k_frac), n_bands))
        return _fake_optical_bands(k_frac)

    real_mech = acoustic_comsol.run_mechanical_comsol
    real_opt = objective_rib._solve_optical_bands
    acoustic_comsol.run_mechanical_comsol = fake_mech
    objective_rib._solve_optical_bands = fake_opt
    try:
        rec = evaluate_candidate(_rib_u0(), celltype="rib",
                                 optical_backend="comsol", mech_backend="comsol")
        rec_sym = evaluate_candidate(_rib_u0(), celltype="rib",
                                     optical_backend="none", mech_backend="comsol",
                                     rib_options=dict(mech_gap_mode="sym"))
    finally:
        acoustic_comsol.run_mechanical_comsol = real_mech
        objective_rib._solve_optical_bands = real_opt

    rib = ("omc_unitcell_rib.mph", 20, "RibGeometry")
    assert calls == [("opt TE", 13, 10), ("mech sym",) + rib,
                     ("mech antisym",) + rib, ("mech sym",) + rib]
    assert rec["status"] == "success"
    assert abs(rec["mechanical_gap"] - 0.24) < 1e-9           # complete gap
    assert sorted(rec["mech_gaps"]) == ["antisym", "complete", "sym"]
    assert abs(rec["optical_gap"] - 35.0 / 192.5) < 1e-6
    # score version 2: the optical gap (175-210 THz, centre 192.5 THz, 0.914 THz
    # below the 193.414 THz target) is weighted by exp(-(0.914/2)^2) = 0.811
    pen = np.exp(-((192.5e12 - 299_792_458.0 / 1550e-9) / 2e12) ** 2)
    g_eff = 35.0 / 192.5 * pen
    assert abs(rec["optical_freq_penalty"] - pen) < 1e-12
    assert abs(rec["optical_gap_effective"] - g_eff) < 1e-12
    assert abs(rec["score"] - (g_eff - 10.0 * (0.20 - g_eff) ** 2
                               - ((6.25 - 8.0) / 8.0) ** 2)) < 1e-9
    assert rec["score_version"] == 2 and rec_sym["score_version"] == 1
    assert rec_sym["optical_gap"] is None
    assert abs(rec_sym["mechanical_gap"] - 3.0 / 7.0) < 1e-9
    assert rec["id"] != rec_sym["id"]
    json.dumps(rec, allow_nan=False)                           # strict JSON


def test_optimizer_n_dim():
    from optimizer import RandomSampler
    assert len(RandomSampler().ask()) == 5
    assert len(RandomSampler(n_dim=4).ask()) == 4


def test_rib_results_use_their_own_database():
    import database
    assert database.db_path("hole") == database._DB
    assert database.db_path("rib") == database._DB_RIB != database._DB
    if "OMC_DB_PATH_RIB" not in os.environ:
        assert os.path.basename(database.db_path("rib")) == "runs_rib.sqlite3"
    rib = evaluate_candidate(_rib_u0(), celltype="rib", optical_backend="none",
                             mech_backend="surrogate_stub")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "runs_rib.sqlite3")
        database.save_result(rib, path)
        loaded = database.load_completed_results(path)
        assert [r["id"] for r in loaded] == [rib["id"]]
        assert loaded[0]["celltype"] == "rib"
        c = sqlite3.connect(path)
        cols = [row[1] for row in c.execute("PRAGMA table_info(candidates)")]
        c.close()
        assert "celltype" not in cols and len(cols) == 13      # schema unchanged


# ---- rib optical scoring, score version 2: no centre window, Gaussian weight ----
F_T = 299_792_458.0 / 1550e-9                 # 193.414 THz


def _two_gap_bands(k_frac, a=570e-9, upper=(215.0, 240.0)):
    """Guided bands with resolvable gaps 130-150, 150-185 THz (centre 167.5)
    and upper[0]-upper[1] THz; flat bands, so each band is one frequency."""
    from comsol_optics import C0
    k = np.asarray(k_frac, dtype=float)
    cols = [130.0, 150.0, 185.0, upper[0], upper[1]]
    rad = C0 * k / (2 * a) * 1.05              # one mode above the light line
    return np.column_stack([np.full_like(k, f * 1e12) for f in cols] + [rad])


def test_freq_penalty_numbers():
    from comsol_optics import freq_penalty
    assert freq_penalty(F_T, F_T, 2e12) == 1.0                       # on target
    assert abs(freq_penalty(F_T + 2e12, F_T, 2e12) - np.exp(-1)) < 1e-15
    assert abs(freq_penalty(F_T - 2e12, F_T, 2e12) - np.exp(-1)) < 1e-15
    assert abs(freq_penalty(F_T + 4e12, F_T, 2e12) - np.exp(-4)) < 1e-15
    assert abs(freq_penalty(F_T + 1e12, F_T, 2e12) - np.exp(-0.25)) < 1e-15
    assert freq_penalty(0.0, F_T, 2e12) == 0.0                       # no gap
    for bad in (0.0, -1.0, float("nan")):
        try:
            freq_penalty(F_T, F_T, bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"tolerance {bad} accepted")


def test_gap_centre_window_is_gone_and_best_weighted_gap_is_chosen():
    from comsol_optics import (guided_interval_gap, guided_interval_gaps,
                               rib_k_grid, select_gap_by_penalized_size)
    a = 570e-9
    k = rib_k_grid(a, 13)
    # 1) gaps far from the target only: found now, weight ~ 0
    far = _two_gap_bands(k, upper=(186.0, 187.0))      # 185-186-187: unresolvable
    gaps = guided_interval_gaps(k, far, a)
    assert [(round(g.f_lower / 1e12), round(g.f_upper / 1e12)) for g in gaps] == [
        (130, 150), (150, 185)]
    gp, pen, cands = select_gap_by_penalized_size(gaps, F_T, 2e12)
    assert gp.found and pen < 1e-30 and len(cands) == 2
    # ... the windowed function is unchanged: found only inside its window
    assert guided_interval_gap(k, far, a, rel_tol=0.30).found
    assert not guided_interval_gap(k, far, a, rel_tol=0.05).found
    # 2) the other rules stay: wider than the resolvable step (11.4 THz here)
    assert all(g.gap_size > 11e12 for g in gaps)
    # 3) two gaps: the WIDER one (150-185, 20.9 %) loses to the narrower one
    #    near the target (188-200, 6.2 %, centre 194 THz) because of the weight
    two = _two_gap_bands(k, upper=(188.0, 200.0))
    gaps = guided_interval_gaps(k, two, a)
    widest = max(gaps, key=lambda g: g.normalized_gap)
    assert round(widest.f_center / 1e12, 1) == 167.5
    gp, pen, cands = select_gap_by_penalized_size(gaps, F_T, 2e12)
    assert (round(gp.f_lower / 1e12), round(gp.f_upper / 1e12)) == (188, 200)
    assert abs(pen - np.exp(-((194e12 - F_T) / 2e12) ** 2)) < 1e-12
    assert [c["chosen"] for c in cands] == [True] + [False] * (len(cands) - 1)
    assert cands[0]["penalized_gap"] == max(c["penalized_gap"] for c in cands)
    assert abs(cands[0]["penalized_gap"] - gp.normalized_gap * pen) < 1e-15
    # a very wide tolerance gives back "the widest gap"
    wide, pen_w, _ = select_gap_by_penalized_size(gaps, F_T, 1e15)
    assert round(wide.f_center / 1e12, 1) == 167.5 and pen_w > 0.99
    # no gap at all
    gp, pen, cands = select_gap_by_penalized_size([], F_T, 2e12)
    assert not gp.found and pen == 0.0 and cands == []


def _rib_with_fakes(optical_bands, **kw):
    """evaluate_candidate (rib, optical comsol) on fake solvers."""
    import acoustic_comsol
    import comsol_client
    import objective_rib
    k, sym, anti = _fake_mech_bands()
    bands = {"mech sym": sym, "mech antisym": anti}
    real = (acoustic_comsol.run_mechanical_comsol, objective_rib._solve_optical_bands,
            comsol_client.template_path)
    acoustic_comsol.run_mechanical_comsol = (
        lambda g, *, n_k=15, n_bands=12, template=None, study="mech sym", **o:
        dict(freqs_hz=bands[study], k_frac=k, a=g.a))
    objective_rib._solve_optical_bands = (
        lambda g, study, k_frac, n_bands, f_target, **o: optical_bands(k_frac))
    comsol_client.template_path = lambda *a, **o: "omc_unitcell_rib.mph"
    try:
        return evaluate_candidate(_rib_u0(), celltype="rib", optical_backend="comsol",
                                  **{"mech_backend": "comsol", **kw})
    finally:
        (acoustic_comsol.run_mechanical_comsol, objective_rib._solve_optical_bands,
         comsol_client.template_path) = real


def test_rib_score_uses_the_weighted_optical_gap():
    import inspect
    import objective_rib
    targets = dict(optical=dict(min_gap=0.20), mechanical=dict(min_gap=0.20),
                   scoring=dict(lambda_optical=10.0, lambda_mech=10.0,
                                lambda_optical_freq=1.0, lambda_mech_freq=1.0))
    score = objective_rib._score_rib
    # both gaps: min(G_o_eff, G_m) - size penalties - mechanical centre penalty
    assert abs(score(0.30, 0.25, 8e9, 8e9, targets) - 0.25) < 1e-15
    assert abs(score(0.10, 0.25, 8e9, 8e9, targets) - (0.10 - 10 * 0.10 ** 2)) < 1e-15
    assert abs(score(0.30, 0.25, 7.2e9, 8e9, targets) - (0.25 - 0.1 ** 2)) < 1e-15
    # the optical centre penalty of version 1 (lambda_optical_freq) is gone:
    # the function does not even take the optical centre any more
    assert list(inspect.signature(score).parameters) == [
        "G_o_eff", "G_m", "f_m_c", "f_m_t", "targets"]
    nan = float("nan")                                   # one gap only / none
    assert abs(score(0.30, nan, nan, 8e9, targets) - 0.30) < 1e-15
    assert abs(score(nan, 0.25, 8e9, 8e9, targets) - 0.25) < 1e-15
    assert score(nan, nan, nan, 8e9, targets) == 0.0
    assert abs(score(0.0, nan, nan, 8e9, targets) + 10 * 0.2 ** 2) < 1e-15

    # end to end: an off-target gap is found, recorded, and weighted to ~ 0
    rec = _rib_with_fakes(lambda k: _two_gap_bands(k, upper=(186.0, 187.0)))
    assert rec["status"] == "success", rec.get("error")
    assert abs(rec["optical_gap"] - 35.0 / 167.5) < 1e-9         # physical gap
    assert rec["optical_gap"] == rec["optical_gap_unpenalized"]
    assert abs(rec["optical_center_frequency"] - 167.5e12) < 1.0
    assert rec["optical_freq_penalty"] < 1e-30 and rec["optical_gap_effective"] < 1e-30
    assert abs(rec["score"] - (0.0 - 10 * 0.2 ** 2 - ((6.25 - 8) / 8) ** 2)) < 1e-9
    assert rec["optical_gap_detail"]["found"] and len(rec["optical_gap_candidates"]) == 2
    assert rec["optical_tol_thz"] == 2.0 and rec["score_version"] == 2
    assert rec["optical_rel_tol"] is None                        # deprecated key kept
    json.dumps(rec, allow_nan=False)
    # two gaps near and far: the one near the target is scored, all are listed
    near = lambda k: _two_gap_bands(k, upper=(188.0, 200.0))     # noqa: E731
    rec2 = _rib_with_fakes(near)
    pen = np.exp(-((194e12 - F_T) / 2e12) ** 2)
    assert abs(rec2["optical_center_frequency"] - 194e12) < 1.0
    assert abs(rec2["optical_freq_penalty"] - pen) < 1e-12
    assert abs(rec2["optical_gap_effective"] - 12.0 / 194.0 * pen) < 1e-12
    assert abs(rec2["optical_gap"] - 12.0 / 194.0) < 1e-12
    g_eff = 12.0 / 194.0 * pen
    assert abs(rec2["score"] - (g_eff - 10 * (0.2 - g_eff) ** 2
                                - ((6.25 - 8) / 8) ** 2)) < 1e-9
    assert sum(c["chosen"] for c in rec2["optical_gap_candidates"]) == 1
    assert len(rec2["optical_gap_candidates"]) == 3
    # a wider tolerance: another objective, another id, and the wide gap wins
    rec3 = _rib_with_fakes(near, rib_options=dict(optical_tol_thz=1000.0))
    assert rec3["id"] != rec2["id"] and rec3["optical_tol_thz"] == 1000.0
    assert abs(rec3["optical_center_frequency"] - 167.5e12) < 1.0
    for bad in (0.0, -2.0):
        try:
            _rib_with_fakes(near, rib_options=dict(optical_tol_thz=bad))
        except ValueError:
            pass
        else:
            raise AssertionError(f"optical_tol_thz {bad} accepted")


def test_deprecated_optical_rel_tol_warns_and_is_ignored():
    import warnings
    near = lambda k: _two_gap_bands(k, upper=(188.0, 200.0))     # noqa: E731
    ref = _rib_with_fakes(near)
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        old = _rib_with_fakes(near, rib_options=dict(optical_rel_tol=0.001))
    assert any(issubclass(w.category, FutureWarning)
               and "optical_rel_tol" in str(w.message) for w in seen)
    # ignored: a +/- 0.1 % window would have rejected every gap here
    assert old["id"] == ref["id"] and old["score"] == ref["score"]
    assert old["optical_gap"] == ref["optical_gap"] > 0
    with warnings.catch_warnings(record=True) as seen:           # no key: no warning
        warnings.simplefilter("always")
        _rib_with_fakes(near)
    assert not [w for w in seen if "optical_rel_tol" in str(w.message)]
    # the loop script still accepts the old key in a config
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol_rib as loop
    assert loop._rib_options(dict(rib=dict(optical_rel_tol=0.3, n_k=9))) == dict(
        optical_rel_tol=0.3, n_k=9)
    assert "optical_tol_thz" in loop._RIB_OPTION_KEYS
    # every shipped rib config has the new key and not the old one
    import yaml
    cfg_dir = os.path.join(os.path.dirname(__file__), "..", "configs")
    for name in ("run_one_rib", "run_loop_rib", "run_opt_comsol_rib", "run_batch_rib",
                 "run_one_rib_tidy3d", "run_opt_comsol_rib_tidy3d"):
        with open(os.path.join(cfg_dir, name + ".yaml")) as fh:
            rib = yaml.safe_load(fh)["rib"]
        assert rib["optical_tol_thz"] == 2.0 and "optical_rel_tol" not in rib, name


def test_score_version_keeps_ids_and_warm_starts_apart():
    import objective_rib
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol_rib as loop
    assert objective_rib.SCORE_VERSION == 2
    assert objective_rib.score_version("none") == 1
    assert objective_rib.score_version("comsol") == objective_rib.score_version("tidy3d") == 2
    rec = _rib_with_fakes(lambda k: _two_gap_bands(k, upper=(188.0, 200.0)))
    # the id of a version-1 record of the same cell and settings (the old
    # tuple ended with optical_rel_tol = 0.30 and had no score tag)
    base = ("comsol", "comsol", "complete", 15, 20, "opt TE", 13, 10, 0.30,
            "n=2.3860|xsect=rect|z=even")
    assert rec["id"] != objective_rib._rib_id(_rib_u0(), base)
    assert rec["id"] == objective_rib._rib_id(_rib_u0(), base + ("score=v2|otol=2THz",))
    # without optics the score did not change: version 1, the old id
    none = evaluate_candidate(_rib_u0(), celltype="rib", optical_backend="none",
                              mech_backend="surrogate_stub")
    assert none["score_version"] == 1 and none["id"] == objective_rib._rib_id(
        _rib_u0(), ("none", "surrogate_stub", "complete", 15, 20, "opt TE",
                    13, 10, 0.30))
    # warm start filter
    same = loop._same_scoring
    old_rec = dict(status="success", score=0.5, params=rec["params"])   # no marker: v1
    assert same(old_rec, 1, 2.0) and not same(old_rec, 2, 2.0)
    assert same(rec, 2, 2.0) and not same(rec, 1, 2.0)
    assert not same(rec, 2, 5.0)                                  # other tolerance
    assert not same(dict(score_version=2), 2, 2.0)                # tolerance unknown
    try:
        import optuna  # noqa: F401
    except ImportError:
        print("  (warm-start seeding skipped: optuna not installed)")
        return
    from geometry import load_bounds

    class FakeStudy:
        def __init__(self):
            self.trials = []

        def add_trial(self, trial):
            self.trials.append(trial)

    prior = [dict(rec, score=0.11), dict(rec, score=0.12, optical_tol_thz=5.0),
             old_rec, dict(rec, score=0.13, status="infeasible")]
    opt = type("O", (), {})()
    opt.study = FakeStudy()
    n = loop._warm_start(opt, prior, load_bounds(), len(rec["u"]), version=2, tol_thz=2.0)
    assert n == 1 and [t.value for t in opt.study.trials] == [0.11]
    opt.study = FakeStudy()
    n = loop._warm_start(opt, prior, load_bounds(), len(rec["u"]), version=1, tol_thz=2.0)
    assert n == 1 and [t.value for t in opt.study.trials] == [0.5]


def test_hole_scoring_is_unchanged_by_the_rib_score_change():
    import objective
    targets = objective._load("targets.yaml")
    f_t = 299_792_458.0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    # objective._score still has the quadratic optical centre penalty
    s = objective._score(0.25, 0.25, 0.9 * f_t, 8e9, f_t, 8e9, targets)
    assert abs(s - (0.25 - targets["scoring"]["lambda_optical_freq"] * 0.1 ** 2)) < 1e-12
    rec = evaluate_candidate([0.5, 0.6, 0.5, 0.6], optical_backend="surrogate",
                             mech_backend="surrogate_stub")
    for key in ("score_version", "optical_tol_thz", "optical_freq_penalty",
                "optical_gap_effective", "optical_gap_candidates"):
        assert key not in rec, key
    # the hole gap function keeps its centre window
    from comsol_optics import optical_gap_tracked
    assert not optical_gap_tracked(np.array([[100e12, 120e12]])).found
    assert optical_gap_tracked(np.array([[170e12, 215e12]])).found


if __name__ == "__main__":
    test_geometry_mapping_and_bounds()
    test_feasibility_rejects_oversized_hole()
    test_largest_gap_detects_clean_gap()
    test_surrogate_returns_telecom_gap()
    test_objective_end_to_end_surrogate()
    test_rib_u_mapping_roundtrip()
    test_rib_template_validity_rules()
    test_rib_fabricability_rules()
    test_hole_cell_unchanged_by_rib_support()
    test_interval_gap_single_and_union()
    test_guided_interval_gap_on_synthetic_bands()
    test_rib_objective_gates_and_no_solver()
    test_rib_objective_with_fake_comsol()
    test_optimizer_n_dim()
    test_rib_results_use_their_own_database()
    test_freq_penalty_numbers()
    test_gap_centre_window_is_gone_and_best_weighted_gap_is_chosen()
    test_rib_score_uses_the_weighted_optical_gap()
    test_deprecated_optical_rel_tol_warns_and_is_ignored()
    test_score_version_keeps_ids_and_warm_starts_apart()
    test_hole_scoring_is_unchanged_by_the_rib_score_change()
    print("all tests passed")
