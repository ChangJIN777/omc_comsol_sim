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
    # score version 4: G_o + G_m - penalties. Optical gap 175-210 THz (18.2 %,
    # centre 192.5 THz): size penalty 10*(0.2 - 0.1818)^2, flat-top detuning
    # 9*r/(1 + r), r = (0.914/5)^8 (~1e-5). Mechanical gap 24 % centred at
    # 6.25 GHz: inside the 5-10 GHz window, and no 8 GHz centre term.
    g_o, f_t = 35.0 / 192.5, 299_792_458.0 / 1550e-9
    r = (abs(192.5e12 - f_t) / 5e12) ** 8
    lof = __import__("yaml").safe_load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "targets.yaml")))["scoring_rib"]["lambda_optical_freq"]   # shipped weight
    pens = 10.0 * (0.20 - g_o) ** 2 + lof * r / (1 + r)
    assert abs(rec["score"] - (g_o + 0.24 - pens)) < 1e-9
    assert rec["pen_window_m"] == 0.0 and rec["pen_lightcone_o"] == 0.0
    assert rec["pen_center_m"] == 0.0
    assert rec["score_version"] == 4 and rec_sym["score_version"] == 4
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
        return evaluate_candidate(_rib_u0(), celltype="rib",
                                  **{"optical_backend": "comsol",
                                     "mech_backend": "comsol", **kw})
    finally:
        (acoustic_comsol.run_mechanical_comsol, objective_rib._solve_optical_bands,
         comsol_client.template_path) = real


def test_rib_score_v4_terms_and_sum():
    """Score version 4: G_o + G_m minus five separate penalties (flat-top
    detuning, no 8 GHz centre term); version-3 weights still score as v3."""
    import objective_rib
    targets = dict(optical=dict(min_gap=0.20), mechanical=dict(min_gap=0.20),
                   scoring=dict(lambda_optical=10.0, lambda_mech=10.0,
                                lambda_optical_freq=1.0, lambda_mech_freq=1.0))
    w = objective_rib.rib_score_weights(targets)
    assert w == dict(gmin_o=0.2, gmin_m=0.2, lambda_optical=10.0, lambda_mech=10.0,
                     lambda_optical_freq=9.0, optical_freq_tol_THz=5.0,
                     optical_freq_power=8.0,
                     lambda_light_cone=2000.0, lambda_mech_window=2000.0,
                     mech_window_GHz=[5.0, 10.0])          # rib defaults; hole's lof = 1 unused
    assert "lambda_mech_freq" not in w                     # no 8 GHz centre term
    assert objective_rib.SCORE_VERSION == 4
    w3 = objective_rib.rib_score_weights(targets, version=3)   # an older record's set
    assert w3 == dict(gmin_o=0.2, gmin_m=0.2, lambda_optical=10.0, lambda_mech=10.0,
                      lambda_mech_freq=1.0, lambda_optical_freq=9.0,
                      lambda_light_cone=2000.0, lambda_mech_window=2000.0,
                      mech_window_GHz=[5.0, 10.0])
    a = 570e-9                                             # light line at X: 263 THz
    kw = dict(a=a, f_o_target=F_T, f_m_target=8e9, weights=w)
    score = lambda *args: objective_rib.rib_score(*args, **kw)   # noqa: E731

    # the two gaps ADD (versions 1 and 2 took the smaller one)
    s, terms = score(0.30, F_T, True, 0.25, 8e9, True)
    assert abs(s - 0.55) < 1e-15 and all(v == 0.0 for v in terms.values())
    assert s > min(0.30, 0.25)
    # size penalties, each on its own gap
    s, terms = score(0.10, F_T, True, 0.15, 8e9, True)
    assert abs(terms["pen_size_o"] - 10 * 0.10 ** 2) < 1e-15
    assert abs(terms["pen_size_m"] - 10 * 0.05 ** 2) < 1e-15
    assert abs(s - (0.25 - 0.1 - 0.025)) < 1e-15
    # detuning: flat-top, subtracted outside the gap term
    #   lof * [1 - 1/(1 + (|df|/f_tol)^8)], lof = 9, f_tol = 5 THz
    det = lambda df: score(0.30, F_T + df, True, 0.25, 8e9, True)[1]["pen_detune_o"]  # noqa: E731
    assert det(0.0) == 0.0
    assert det(1e12) < 3e-5 and det(2.5e12) < 0.04         # flat inside the tolerance
    assert abs(det(5e12) - 4.5) < 1e-9                     # lof/2 at f_tol
    assert abs(det(10e12) - 9.0 * (1 - 1 / 257.0)) < 1e-9  # saturated: 8.965
    assert abs(det(4e12) - 9.0 * (1 - 1 / (1 + 0.8 ** 8))) < 1e-9
    assert abs(det(100e12) - 9.0) < 1e-8 and det(100e12) <= 9.0
    for df in (0.5e12, 3e12, 5e12, 7e12, 20e12):           # symmetric in +/- df
        assert abs(det(df) - det(-df)) < 1e-9
    xs = [det(i * 0.25e12) for i in range(120)]
    assert all(b >= a_ for a_, b in zip(xs, xs[1:]))       # monotonic in |df|
    s, terms = score(0.30, F_T + 10e12, True, 0.25, 8e9, True)
    assert abs(s - (0.55 - terms["pen_detune_o"])) < 1e-12
    w2thz = dict(w, optical_freq_tol_THz=2.0)              # the tolerance is a weight
    assert abs(objective_rib.optical_terms(0.3, F_T - 2e12, F_T, a, w2thz)["pen_detune_o"]
               - 4.5) < 1e-9
    # version-3 weights: the quadratic of that version
    assert abs(objective_rib.optical_terms(0.3, F_T + 10e12, F_T, a, w3)["pen_detune_o"]
               - 9.0 * (10e12 / F_T) ** 2) < 1e-15
    # light cone: centre above c/2a
    ll = 299_792_458.0 / (2 * a)
    _, terms = score(0.30, 1.05 * ll, True, 0.25, 8e9, True)
    assert abs(terms["pen_lightcone_o"] - 2000 * 0.05 ** 2) < 1e-9          # = 5
    assert score(0.30, ll, True, 0.25, 8e9, True)[1]["pen_lightcone_o"] == 0.0
    assert score(0.30, 0.9 * ll, True, 0.25, 8e9, True)[1]["pen_lightcone_o"] == 0.0
    # mechanical window [5, 10] GHz: quadratic in the fractional excess
    win = lambda f: score(0.30, F_T, True, 0.25, f, True)[1]["pen_window_m"]   # noqa: E731
    assert win(5e9) == 0.0 and win(7.3e9) == 0.0 and win(10e9) == 0.0
    assert abs(win(4e9) - 2000 * 0.2 ** 2) < 1e-9          # (5-4)/5  -> 80
    assert abs(win(12e9) - 2000 * 0.2 ** 2) < 1e-9         # (12-10)/10 -> 80
    assert abs(win(4.5e9) - 2000 * 0.1 ** 2) < 1e-9 and abs(win(11e9) - 2000 * 0.01) < 1e-9
    # no 8 GHz centre term: anywhere in the window scores the same,
    # whatever mechanical.target_frequency_GHz says
    for f_m in (5e9, 6e9, 8e9, 10e9):
        for f_m_t in (8e9, 6e9):
            s, terms = objective_rib.rib_score(0.30, F_T, True, 0.25, f_m, True, a=a,
                                               f_o_target=F_T, f_m_target=f_m_t, weights=w)
            assert terms["pen_center_m"] == 0.0 and abs(s - 0.55) < 1e-15
    # ... while a version-3 record's weights still give that term
    s3, t3 = objective_rib.rib_score(0.30, F_T, True, 0.25, 6e9, True, a=a,
                                     f_o_target=F_T, f_m_target=8e9, weights=w3)
    assert abs(t3["pen_center_m"] - 0.25 ** 2) < 1e-15 and abs(s3 - (0.55 - 0.0625)) < 1e-15
    # a gap that was searched for but not found
    s, terms = score(0.0, 0.0, False, 0.25, 8e9, True)     # optical
    assert abs(terms["pen_size_o"] - 0.4) < 1e-15 and terms["pen_detune_o"] == 0.0
    assert terms["pen_lightcone_o"] == 0.0 and abs(s - (0.25 - 0.4)) < 1e-15
    big = dict(kw, a=1000e-9)                              # target above c/2a = 150 THz
    t2 = objective_rib.rib_score(0.0, 0.0, False, 0.25, 8e9, True, **big)[1]
    assert abs(t2["pen_lightcone_o"] - 2000 * ((F_T - 149.896229e12) / 149.896229e12) ** 2) < 1e-3
    s, terms = score(0.30, F_T, True, 0.0, 0.0, False)     # mechanical
    assert abs(terms["pen_size_m"] - 0.4) < 1e-15 and terms["pen_window_m"] == 0.0
    assert terms["pen_center_m"] == 0.0 and abs(s - (0.30 - 0.4)) < 1e-15
    # physics that was not solved: only the other gap's terms
    nan = float("nan")
    s, terms = score(0.30, F_T, True, nan, nan, False)
    assert abs(s - 0.30) < 1e-15 and terms["pen_size_m"] is None and terms["pen_window_m"] is None
    s, terms = score(nan, nan, False, 0.25, 4e9, True)
    assert abs(s - (0.25 - 80)) < 1e-9 and terms["pen_detune_o"] is None
    assert score(nan, nan, False, nan, nan, False) == (0.0, dict.fromkeys(terms))

    # the same light-cone and window numbers as scripts/run_opt_comsol.py
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol as hole
    assert hole.LAMBDA_LIGHTCONE == 2000.0
    assert (hole.MECH_F_MIN, hole.MECH_F_MAX) == (5e9, 10e9)
    for f_c, aa in ((F_T + 10e12, 570e-9), (F_T - 25e12, 700e-9), (240e12, 700e-9),
                    (200e12, 900e-9)):
        o = objective_rib.optical_terms(0.30, f_c, F_T, aa, w)
        ref_lc = hole.LAMBDA_LIGHTCONE * max(
            0.0, (f_c - 299_792_458.0 / (2 * aa)) / (299_792_458.0 / (2 * aa))) ** 2
        assert abs(o["pen_lightcone_o"] - ref_lc) < 1e-12
    # on target the whole score is the hole loop's (its detuning term is the
    # quadratic of version 3; its size penalties are not triggered here)
    for f_m in (5.5e9, 8e9, 9.5e9):           # in the window (its weight there is 0)
        ref = hole.score_result(0.30, 0.25, hole.OPT_F_TARGET, f_m, a=570e-9)
        mine = objective_rib.rib_score(0.30, F_T, True, 0.25, f_m, True, a=570e-9,
                                       f_o_target=F_T, f_m_target=8e9, weights=w)[0]
        assert abs(mine - ref) < 1e-9, f_m
    # weights: the rib block of targets.yaml, validated
    t2 = dict(targets, scoring_rib=dict(lambda_optical_freq=3.0, mech_window_GHz=[6, 9]))
    w2 = objective_rib.rib_score_weights(t2)
    assert w2["lambda_optical_freq"] == 3.0 and w2["mech_window_GHz"] == [6.0, 9.0]
    assert w2["lambda_light_cone"] == 2000.0 and w2["optical_freq_tol_THz"] == 5.0
    t3 = dict(targets, scoring_rib=dict(optical_freq_tol_THz=3.0))
    assert objective_rib.rib_score_weights(t3)["optical_freq_tol_THz"] == 3.0
    for bad in (dict(lambda_optcal_freq=1.0), dict(mech_window_GHz=[10, 5]),
                dict(lambda_light_cone=-1.0), dict(optical_freq_tol_THz=0.0)):
        try:
            objective_rib.rib_score_weights(dict(targets, scoring_rib=bad))
        except ValueError:
            pass
        else:
            raise AssertionError(f"scoring_rib {bad} accepted")
    import yaml
    with open(os.path.join(os.path.dirname(__file__), "..", "configs", "targets.yaml")) as fh:
        shipped = yaml.safe_load(fh)
    # every rib key is spelled out there (the values are the user's to tune)
    assert set(shipped["scoring_rib"]) == set(objective_rib.RIB_SCORING_DEFAULTS)
    objective_rib.rib_score_weights(shipped)
    assert shipped["scoring"]["lambda_optical_freq"] == 1.0          # hole: untouched


def test_rib_gap_choice_and_record_under_v4():
    """The candidate with the largest optical score contribution is scored;
    the record carries the terms."""
    import objective_rib
    from comsol_optics import guided_interval_gaps, rib_k_grid, select_gap_by_value
    a = 570e-9
    k = rib_k_grid(a, 13)
    import objective as obj
    w = objective_rib.rib_score_weights(obj._load("targets.yaml"))
    value = lambda gp: objective_rib.optical_terms(   # noqa: E731
        gp.normalized_gap, gp.f_center, F_T, a, w)["value"]
    # 150-185 THz (20.9 %, 26 THz below target) against 188-200 THz (6.2 %, on target)
    gaps = guided_interval_gaps(k, _two_gap_bands(k, upper=(188.0, 200.0)), a)
    gp, cands = select_gap_by_value(gaps, value, F_T)
    lof = w["lambda_optical_freq"]                                   # shipped weight
    flat = lambda df, lof=lof: lof * (1 - 1 / (1 + (abs(df) / 5e12) ** 8))   # noqa: E731
    v_wide = 35 / 167.5 - flat(167.5e12 - F_T)                         # 0.209 - lof
    v_near = 12 / 194.0 - 10 * (0.2 - 12 / 194.0) ** 2 - flat(194e12 - F_T)
    assert abs(cands[0]["value"] - max(v_wide, v_near)) < 1e-12
    assert v_near > v_wide                # the on-target gap wins: flat-top penalty
    assert abs(flat(167.5e12 - F_T) - lof) < 1e-4 and flat(194e12 - F_T) < 1e-6
    assert round(gp.f_center / 1e12, 1) == 194.0 and cands[0]["chosen"]
    assert [c["chosen"] for c in cands] == [True, False, False]
    assert all("value" in c and "freq_penalty" not in c for c in cands)
    # with a weak detuning weight the wide off-target gap is chosen instead
    w_weak = dict(w, lambda_optical_freq=0.05)
    gp9, _ = select_gap_by_value(gaps, lambda g_: objective_rib.optical_terms(
        g_.normalized_gap, g_.f_center, F_T, a, w_weak)["value"], F_T)
    assert round(gp9.f_center / 1e12, 1) == 167.5
    # ... as with the weights of a version-3 record (quadratic detuning)
    w3 = objective_rib.rib_score_weights(obj._load("targets.yaml"), version=3)
    gp3, _ = select_gap_by_value(gaps, lambda g_: objective_rib.optical_terms(
        g_.normalized_gap, g_.f_center, F_T, a, w3)["value"], F_T)
    assert round(gp3.f_center / 1e12, 1) == 167.5
    assert select_gap_by_value([], value, F_T)[1] == []

    rec = _rib_with_fakes(lambda kk: _two_gap_bands(kk, upper=(188.0, 200.0)))
    assert rec["status"] == "success", rec.get("error")
    assert rec["score_version"] == 4 and rec["score_weights"] == w
    assert abs(rec["optical_center_frequency"] - 194e12) < 1.0
    assert abs(rec["optical_gap"] - 12 / 194.0) < 1e-12
    assert rec["optical_gap_effective"] == rec["optical_gap"] == rec["optical_gap_unpenalized"]
    det = flat(194e12 - F_T)
    assert abs(rec["pen_detune_o"] - det) < 1e-12 and rec["optical_freq_penalty"] == rec["pen_detune_o"]
    size = 10 * (0.2 - 12 / 194.0) ** 2
    assert abs(rec["pen_size_o"] - size) < 1e-12 and rec["pen_lightcone_o"] == 0.0
    # mechanics of the fake: complete gap 5.5-7.0 GHz (24 %), centre 6.25 GHz:
    # in the window, and no penalty for being 1.75 GHz from 8 GHz
    assert rec["pen_size_m"] == 0.0 and rec["pen_window_m"] == 0.0
    assert rec["pen_center_m"] == 0.0
    assert abs(rec["score"] - (12 / 194.0 + 0.24 - size - det)) < 1e-9
    assert rec["optical_tol_thz"] is None and rec["optical_rel_tol"] is None
    assert len(rec["optical_gap_candidates"]) == 3
    json.dumps(rec, allow_nan=False)
    # no optical gap at all: size penalty only (centre taken as the target)
    none = _rib_with_fakes(lambda kk: _two_gap_bands(kk)[:, :1])   # one band only
    assert none["optical_gap"] == 0.0 and abs(none["pen_size_o"] - 0.4) < 1e-12
    assert none["pen_detune_o"] == 0.0 and none["pen_lightcone_o"] == 0.0
    # mechanics only: the window term applies there too
    mech = _rib_with_fakes(lambda kk: None, optical_backend="none")
    assert mech["optical_gap"] is None and mech["pen_size_o"] is None
    assert mech["score_version"] == 4
    assert abs(mech["score"] - 0.24) < 1e-9


def test_deprecated_rib_options_warn_and_are_ignored():
    import warnings
    near = lambda k: _two_gap_bands(k, upper=(188.0, 200.0))     # noqa: E731
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        ref = _rib_with_fakes(near)
    assert not [w for w in seen if "deprecated" in str(w.message)]
    for name, val in (("optical_tol_thz", 2.0), ("optical_tol_thz", 500.0),
                      ("optical_rel_tol", 0.001)):
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter("always")
            old = _rib_with_fakes(near, rib_options={name: val})
        assert any(issubclass(w.category, FutureWarning) and name in str(w.message)
                   and "ignored" in str(w.message) for w in seen), name
        # ignored: same id, same gap, same score
        assert old["id"] == ref["id"] and old["score"] == ref["score"], name
        assert old["optical_center_frequency"] == ref["optical_center_frequency"]
        assert old["optical_tol_thz"] is None
    # the loop script still accepts both keys in a config
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol_rib as loop
    assert loop._rib_options(dict(rib=dict(optical_rel_tol=0.3, optical_tol_thz=5.0,
                                           n_k=9)))["n_k"] == 9
    # the shipped rib configs carry neither key any more
    import yaml
    cfg_dir = os.path.join(os.path.dirname(__file__), "..", "configs")
    for name in ("run_one_rib", "run_loop_rib", "run_opt_comsol_rib", "run_batch_rib",
                 "run_one_rib_tidy3d", "run_opt_comsol_rib_tidy3d"):
        with open(os.path.join(cfg_dir, name + ".yaml")) as fh:
            rib = yaml.safe_load(fh)["rib"]
        assert "optical_tol_thz" not in rib and "optical_rel_tol" not in rib, name


def test_score_version_keeps_ids_and_warm_starts_apart():
    import objective_rib
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol_rib as loop
    assert objective_rib.score_version("none") == objective_rib.score_version("comsol") == 4
    rec = _rib_with_fakes(lambda k: _two_gap_bands(k, upper=(188.0, 200.0)))
    base = ("comsol", "comsol", "complete", 15, 20, "opt TE", 13, 10, 0.30,
            "n=2.3860|xsect=rect|z=even")
    v1, v2 = objective_rib._rib_id(_rib_u0(), base), objective_rib._rib_id(
        _rib_u0(), base + ("score=v2|otol=2THz",))
    v3 = objective_rib._rib_id(_rib_u0(), base + ("score=v3",))
    assert rec["id"] not in (v1, v2, v3)                   # never replaces v1 / v2 / v3
    assert rec["id"] == objective_rib._rib_id(_rib_u0(), base + ("score=v4",))
    # mechanics-only records changed score too (window term): new id as well
    none = evaluate_candidate(_rib_u0(), celltype="rib", optical_backend="none",
                              mech_backend="surrogate_stub")
    legacy = ("none", "surrogate_stub", "complete", 15, 20, "opt TE", 13, 10, 0.30)
    assert none["score_version"] == 4
    assert none["id"] == objective_rib._rib_id(_rib_u0(), legacy + ("score=v4",))
    assert none["id"] != objective_rib._rib_id(_rib_u0(), legacy)
    assert none["id"] != objective_rib._rib_id(_rib_u0(), legacy + ("score=v3",))
    # warm start filter: version and weights
    w = rec["score_weights"]
    same = loop._same_scoring
    v1_rec = dict(status="success", score=0.5, params=rec["params"])        # no marker
    v2_rec = dict(v1_rec, score_version=2, optical_tol_thz=5.0, score=0.4)
    import objective as obj
    w3 = objective_rib.rib_score_weights(obj._load("targets.yaml"), version=3)
    v3_rec = dict(rec, score_version=3, score_weights=w3, score=0.3)
    assert same(rec, 4, w) and not same(v1_rec, 4, w) and not same(v2_rec, 4, w)
    assert not same(v3_rec, 4, w) and not same(rec, 3, w3) and same(v3_rec, 3, w3)
    assert not same(dict(v3_rec, score_version=4), 4, w)   # v3 weights are not v4's
    assert same(v1_rec, 1, w) and same(v2_rec, 2, w)
    other = dict(rec, score_weights=dict(w, lambda_optical_freq=1.0))
    assert not same(other, 4, w)                           # other weights: other objective
    assert not same(dict(rec, score_weights=dict(w, optical_freq_tol_THz=2.0)), 4, w)
    assert not same(dict(rec, score_weights=dict(w, optical_freq_power=4.0)), 4, w)
    assert not same(dict(rec, score_weights=None), 4, w)
    assert same(json.loads(json.dumps(rec)), 4, w)         # survives the JSON round trip
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

    prior = [dict(rec, score=0.11), dict(other, score=0.12), v1_rec, v2_rec, v3_rec,
             dict(rec, score=0.13, status="infeasible")]
    opt = type("O", (), {})()
    opt.study = FakeStudy()
    n = loop._warm_start(opt, prior, load_bounds(), len(rec["u"]), version=4, weights=w)
    assert n == 1 and [t.value for t in opt.study.trials] == [0.11]


def test_hole_scoring_is_unchanged_by_the_rib_score_change():
    """objective._score and run_opt_comsol.score_result against independent
    re-statements of their formulas, and no rib key in a hole record."""
    import objective
    targets = objective._load("targets.yaml")
    f_t = 299_792_458.0 / (targets["optical"]["target_wavelength_nm"] * 1e-9)
    s = targets["scoring"]
    assert (s["lambda_optical"], s["lambda_mech"], s["lambda_optical_freq"],
            s["lambda_mech_freq"]) == (10.0, 10.0, 1.0, 1.0)
    rng = np.random.default_rng(7)
    for _ in range(200):
        G_o, G_m = rng.uniform(0, 0.4, 2)
        f_o, f_m = f_t * rng.uniform(0.6, 1.4), 8e9 * rng.uniform(0.3, 2.0)
        ref = (min(G_o, G_m) - 10 * max(0, 0.2 - G_o) ** 2 - 10 * max(0, 0.2 - G_m) ** 2
               - ((f_o - f_t) / f_t) ** 2 - ((f_m - 8e9) / 8e9) ** 2)
        assert abs(objective._score(G_o, G_m, f_o, f_m, f_t, 8e9, targets) - ref) < 1e-12
        ref1 = G_o - 10 * max(0, 0.2 - G_o) ** 2 - ((f_o - f_t) / f_t) ** 2
        assert abs(objective._score(G_o, float("nan"), f_o, f_m, f_t, 8e9, targets) - ref1) < 1e-12
    rec = evaluate_candidate([0.5, 0.6, 0.5, 0.6], optical_backend="surrogate",
                             mech_backend="surrogate_stub")
    for key in ("score_version", "score_weights", "optical_tol_thz", "optical_freq_penalty",
                "optical_gap_effective", "optical_gap_candidates", "pen_size_o",
                "pen_detune_o", "pen_lightcone_o", "pen_size_m", "pen_window_m",
                "pen_center_m"):
        assert key not in rec, key
    # the hole loop: constants and formula as they were
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol as hole
    assert (hole.G_MIN_OPT, hole.G_MIN_MECH, hole.LAMBDA_OPT, hole.LAMBDA_MECH,
            hole.LAMBDA_OPT_F, hole.LAMBDA_MECH_F, hole.LAMBDA_LIGHTCONE) == (
        0.15, 0.20, 20.0, 5.0, 9.0, 0.0, 2000.0)
    for _ in range(200):
        G_o, G_m = rng.uniform(0, 0.4, 2)
        f_o, a = f_t * rng.uniform(0.6, 1.4), rng.uniform(500e-9, 900e-9)
        ll = 299_792_458.0 / (2 * a)
        ref = (G_o + G_m - 20 * max(0, 0.15 - G_o) ** 2 - 5 * max(0, 0.20 - G_m) ** 2
               - 9 * ((f_o - f_t) / f_t) ** 2 - 2000 * max(0, (f_o - ll) / ll) ** 2)
        assert abs(hole.score_result(G_o, G_m, f_o, 3e9, a=a) - ref) < 1e-10
    assert hole.score_result(float("nan"), 0.2, f_t, 8e9) == -5.0
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
    test_rib_score_v4_terms_and_sum()
    test_rib_gap_choice_and_record_under_v4()
    test_deprecated_rib_options_warn_and_are_ignored()
    test_score_version_keeps_ids_and_warm_starts_apart()
    test_hole_scoring_is_unchanged_by_the_rib_score_change()
    print("all tests passed")
