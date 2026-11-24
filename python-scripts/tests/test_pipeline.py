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
    assert 0.10 < rec["score"] < 0.19                          # min(G_o, G_m) - penalties
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
    print("all tests passed")
