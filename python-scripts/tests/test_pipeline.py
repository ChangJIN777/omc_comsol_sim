"""Minimal tests that run with numpy only (no MPB/COMSOL)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from geometry import u_to_geometry, check_feasibility
from bandgap import largest_gap, gap_near_frequency
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


if __name__ == "__main__":
    test_geometry_mapping_and_bounds()
    test_feasibility_rejects_oversized_hole()
    test_largest_gap_detects_clean_gap()
    test_surrogate_returns_telecom_gap()
    test_objective_end_to_end_surrogate()
    print("all tests passed")
