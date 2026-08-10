"""Minimal tests that run with numpy only (no MPh/COMSOL server needed).

Mirrors python-scripts/tests/test_pipeline.py's style. `acoustic_comsol_2d.py`
is safe to import without MPh installed (see comsol_client.py's try/except),
so bz_path() -- pure math, no COMSOL calls -- is testable here too.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from geometry2d import u_to_geometry, check_feasibility
from bandgap import largest_gap
from objective2d import evaluate_candidate, _mechanical_gap
from acoustic_comsol_2d import bz_path


def test_geometry_mapping_and_bounds():
    g = u_to_geometry([0.0] * 6)
    assert g.a > 0 and g.w > 0 and g.r > 0 and g.th > 0
    g2 = u_to_geometry([1.0] * 6)
    assert g2.a > g.a and g2.w > g.w and g2.th > g.th


def test_feasibility_rejects_oversized_fillet():
    # r1 fillet forced to the top of its range while w is forced to the
    # bottom -> r1 > w/2, which check_feasibility must reject.
    g = u_to_geometry([0.5, 0.0, 0.5, 1.0, 0.0, 0.5])
    ok, reasons = check_feasibility(g)
    assert ok is False and len(reasons) > 0


def test_feasibility_rejects_thin_web():
    # small a with large r/w pushes the hexagon "web" below min_feature.
    g = u_to_geometry([0.0, 1.0, 1.0, 0.0, 0.0, 0.5])
    ok, reasons = check_feasibility(g)
    assert ok is False and any("web" in r for r in reasons)


def test_largest_gap_detects_clean_gap():
    k = np.linspace(0, np.pi, 11)
    b1 = 0.20 + 0.01 * np.sin(k)
    b2 = 0.30 + 0.01 * np.cos(k)
    gp = largest_gap(np.column_stack([b1, b2]))
    assert gp.found and gp.normalized_gap > 0.3


def test_bz_path_hexagonal_shape_and_closure():
    a = 480e-9
    k_norm, kx, ky = bz_path(a, n_per_segment=5, unitcell="hexagonal")
    n_expected = 3 * 5 + 1   # 3 segments * n_per_segment, plus closing Gamma point
    assert k_norm.shape == kx.shape == ky.shape == (n_expected,)
    assert k_norm[0] == 0.0 and k_norm[-1] == 3.0
    # path starts and ends at Gamma (kx=ky=0)
    assert abs(kx[0]) < 1e-9 and abs(ky[0]) < 1e-9
    assert abs(kx[-1]) < 1e-9 and abs(ky[-1]) < 1e-9


def test_mechanical_gap_symmetry_vs_complete():
    # evenz has a clean gap between band 0 and 1; oddz has bands filling
    # that same window -> "symmetry" mode should still find the evenz gap,
    # "complete" mode (which stacks both families) must not.
    k = np.linspace(0, np.pi, 9)
    evenz = np.column_stack([7.0e9 + 0.05e9 * np.sin(k), 9.0e9 + 0.05e9 * np.cos(k)])
    oddz = np.column_stack([7.9e9 + 0.02e9 * np.sin(k), 12.0e9 + 0.02e9 * np.cos(k)])

    gp_sym = _mechanical_gap(evenz, oddz, 8.0e9, "symmetry")
    assert gp_sym.found and gp_sym.normalized_gap > 0.15

    gp_complete = _mechanical_gap(evenz, oddz, 8.0e9, "complete")
    assert (not gp_complete.found) or gp_complete.normalized_gap < gp_sym.normalized_gap


def test_objective_end_to_end_no_backends():
    # require_mech=False, require_opt=False -> no solver called at all;
    # exercises feasibility gating + score wiring in isolation.
    rec = evaluate_candidate([0.5, 0.5, 0.5, 0.3, 0.3, 0.5],
                             require_mech=False, require_opt=False)
    assert rec["status"] == "success"
    assert rec["score"] == 0.0


if __name__ == "__main__":
    test_geometry_mapping_and_bounds()
    test_feasibility_rejects_oversized_fillet()
    test_feasibility_rejects_thin_web()
    test_largest_gap_detects_clean_gap()
    test_bz_path_hexagonal_shape_and_closure()
    test_mechanical_gap_symmetry_vs_complete()
    test_objective_end_to_end_no_backends()
    print("all tests passed")
