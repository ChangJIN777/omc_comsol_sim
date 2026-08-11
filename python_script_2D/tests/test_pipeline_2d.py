"""Minimal tests that run with numpy only (no MPh/COMSOL server needed).

Mirrors python-scripts/tests/test_pipeline.py's style. `acoustic_comsol_2d.py`
is safe to import without MPh installed (see comsol_client.py's try/except),
so bz_path() -- pure math, no COMSOL calls -- is testable here too.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from geometry2d import (FIXED_VARS, REFERENCE, VARS, Geometry2D, u_to_geometry,
                        check_feasibility, clearances, load_bounds, phys_to_u)
from bandgap import largest_gap, gap_near_frequency
from objective2d import evaluate_candidate, _mechanical_gap, _combined_bands
from acoustic_comsol_2d import bz_path
import optimizer

N = len(VARS)


def test_geometry_mapping_and_bounds():
    g = u_to_geometry([0.0] * N)
    assert g.a > 0 and g.w > 0 and g.r > 0 and g.th > 0
    g2 = u_to_geometry([1.0] * N)
    assert g2.a > g.a and g2.w > g.w and g2.r > g.r


def test_free_and_fixed_variable_split():
    """a, w, r are free; r1, r2, th are fixed and absent from u."""
    assert VARS == ["a", "w", "r"]
    assert FIXED_VARS == ["r1", "r2", "th"]
    assert optimizer.N_DIM == N
    b = load_bounds()
    assert set(b["variables"]) == set(VARS)
    assert set(b["fixed"]) == set(FIXED_VARS)


def test_fixed_params_are_constant_across_u():
    """Sweeping u must never move r1, r2 or th."""
    b = load_bounds()
    for u in ([0.0] * N, [0.5] * N, [1.0] * N, [0.2, 0.9, 0.4]):
        g = u_to_geometry(u)
        for name in FIXED_VARS:
            assert getattr(g, name) == float(b["fixed"][name]), name


def test_phys_to_u_rejects_conflicting_fixed_value():
    """Silently dropping a fixed param would make the round-trip lossy."""
    bad = dict(REFERENCE, th=300e-9)   # th is fixed at 220 nm
    try:
        phys_to_u(bad)
    except ValueError as e:
        assert "th" in str(e) and "FIXED" in str(e)
    else:
        raise AssertionError("phys_to_u accepted a conflicting fixed value")


def test_u_to_geometry_rejects_wrong_length():
    for bad in ([0.5] * 6, [0.5] * 2):
        try:
            u_to_geometry(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted u of length {len(bad)}")


def test_selw_matches_the_comsol_template():
    """bounds_2d.yaml:fillet_select_halfwidth and the .m's `selw` are the same
    builder constant stored twice. check_feasibility gates on the YAML copy, so
    drift means Python accepts candidates COMSOL will mis-fillet."""
    import re
    yaml_val = load_bounds()["feasibility"]["fillet_select_halfwidth"]
    mfile = os.path.join(os.path.dirname(__file__), "..", "comsol",
                         "trusty_boomerang_script.m")
    with open(mfile) as fh:
        m = re.search(r"param\.set\(\s*'selw'\s*,\s*'([\d.]+)\[(\w+)\]'", fh.read())
    assert m, "could not find the selw parameter in trusty_boomerang_script.m"
    scale = {"nm": 1e-9, "um": 1e-6, "m": 1.0}[m.group(2)]
    assert np.isclose(float(m.group(1)) * scale, yaml_val), (
        f".m selw={m.group(0)} vs bounds_2d.yaml {yaml_val*1e9:.1f} nm")


def test_reference_design_is_feasible():
    """The MATLAB reference design (test_Boomerang.m) must pass the gate.

    This is the anchor for check_feasibility. It is the only design known to
    have been built and solved in the MATLAB pipeline, so if the Python gate
    rejects it, the gate is wrong -- not the design.
    """
    ok, reasons = check_feasibility(Geometry2D(**REFERENCE))
    assert ok, f"reference design rejected: {reasons}"


def test_reference_design_clearances_are_pinned():
    """Pin the FORMULAS (independently of bounds_2d.yaml's thresholds).

    Values hand-derived from buildBoomerangUnitCell.m for a=480, w=140,
    r=177 nm. If any of these move, the geometry model changed -- check it
    against COMSOL's mphgeom before updating the numbers here.
    """
    c = clearances(Geometry2D(**REFERENCE))
    assert np.isclose(c["in_radius"], 207.85e-9, atol=0.05e-9)   # a*sqrt(3)/4
    assert np.isclose(c["d_tip_face"], 177.0e-9, atol=0.05e-9)   # r
    assert np.isclose(c["d_tip_corner"], 190.33e-9, atol=0.05e-9)  # hypot(r,w/2)
    assert np.isclose(c["d_junction"], 70.0e-9, atol=0.05e-9)    # w/2
    assert np.isclose(c["web"], 17.52e-9, atol=0.05e-9)


def test_in_radius_is_quarter_not_half():
    """Guard against reintroducing the a*sqrt(3)/2 in-radius (2x too large).

    The primitive cell is a RHOMBUS of side a, whose center-to-edge distance
    is a*sqrt(3)/4 -- not a hexagon of side a, whose in-radius is a*sqrt(3)/2.
    """
    g = Geometry2D(**REFERENCE)
    assert np.isclose(clearances(g)["in_radius"], np.sqrt(3) / 4 * g.a)


def test_feasibility_rejects_oversized_fillet():
    # r1/r2 are fixed at 10 nm, so drive w below 2*r1 directly.
    g = Geometry2D(a=480e-9, w=15e-9, r=177e-9, r1=10e-9, r2=10e-9, th=220e-9)
    ok, reasons = check_feasibility(g)
    assert ok is False and any("fillet r1" in r for r in reasons)


def test_feasibility_rejects_thin_web():
    # small a with large r/w drives the hole's tip corner through the cell wall.
    g = u_to_geometry([0.0, 1.0, 1.0])
    ok, reasons = check_feasibility(g)
    assert ok is False and any("web" in r for r in reasons)


def test_fillet_check_accounts_for_junction_displacement():
    """h_fil1 runs before h_disksel2 is evaluated, pushing the junction
    vertices out by ~sqrt(3)*r1. Comparing the PRE-fillet distance w/2 would
    pass this case; the post-fillet distance must reject it.

    w=60, r=100, r1=35 nm: pre-fillet junction 30 nm < r-sw = 75 nm (passes),
    post-fillet 30 + sqrt(3)*35 = 90.6 nm >= 75 nm (must fail).
    """
    g = Geometry2D(a=600e-9, w=60e-9, r=100e-9, r1=35e-9, r2=10e-9, th=220e-9)
    assert g.w / 2 < g.r - 25e-9, "premise: the naive pre-fillet check passes"
    ok, reasons = check_feasibility(g)
    assert ok is False and any("post-fillet junction" in r for r in reasons)


def test_feasibility_rejects_broken_fillet_selection():
    """w=220, r=100 nm is inside bounds_2d.yaml but breaks the MATLAB fillet
    annuli: the tip corner (hypot(100,110)=148.7 nm) is beyond r+25=125 nm, so
    h_disksel2 misses the tips entirely, and r-25=75 nm < w/2=110 nm, so it
    grabs the junctions instead."""
    g = Geometry2D(a=600e-9, w=220e-9, r=100e-9, r1=10e-9, r2=10e-9, th=220e-9)
    ok, reasons = check_feasibility(g)
    assert ok is False and any("fillet selection" in r for r in reasons)


def test_reference_design_roundtrips_through_u():
    u = phys_to_u(REFERENCE)
    assert all(0.0 <= x <= 1.0 for x in u)
    g = u_to_geometry(u)
    for k, v in REFERENCE.items():
        assert np.isclose(getattr(g, k), v, rtol=1e-9), f"{k} did not roundtrip"


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

    gp_sym, info = _mechanical_gap(evenz, oddz, 8.0e9, "symmetry")
    assert gp_sym.found and gp_sym.normalized_gap > 0.15
    assert info["truncation_ceiling_hz"] is None   # not meaningful in this mode

    gp_complete, _ = _mechanical_gap(evenz, oddz, 8.0e9, "complete")
    assert (not gp_complete.found) or gp_complete.normalized_gap < gp_sym.normalized_gap


def test_complete_mode_ignores_bands_above_truncation_ceiling():
    """A gap that exists only because one family ran out of computed bands is
    NOT a complete gap.

    evenz tops out at 9 GHz, oddz at 40 GHz. The 12->40 GHz window contains no
    evenz band purely because evenz was never solved that high -- reporting it
    as "complete" is a truncation artifact. The ceiling is min(9, 40) = 9 GHz,
    so nothing above it may be used as a gap edge.
    """
    k = np.linspace(0, np.pi, 9)
    evenz = np.column_stack([7.0e9 + 0.05e9 * np.sin(k), 9.0e9 + 0.05e9 * np.cos(k)])
    oddz = np.column_stack([8.0e9 + 0.02e9 * np.sin(k), 12.0e9 + 0.02e9 * np.cos(k),
                            40.0e9 + 0.02e9 * np.sin(k)])

    combined, ceiling = _combined_bands(evenz, oddz)
    assert combined.size == 0 or np.nanmax(combined) <= 9.1e9, (
        "bands above the truncation ceiling leaked into the complete-gap search")

    # The spurious 12->40 GHz gap must not be reported at a 20 GHz target.
    gp, _ = _mechanical_gap(evenz, oddz, 20.0e9, "complete")
    assert not gp.found, f"truncation artifact reported as a complete gap: {gp}"


def test_ceiling_uses_min_over_k_not_max():
    """Regression for the original bug: the ceiling must be the MINIMUM over k
    of each family's top resolved band, not the maximum.

    max-over-k is the frequency below which a family SOMETIMES has all its
    modes; a gap has to hold at EVERY k. Here oddz's top band is only 13.4 GHz
    at k=0 but 33.7 GHz at k=2, so the two reductions differ by 2.5x and the
    max-over-k version admits bands in a window where oddz was never solved.
    The earlier truncation test cannot catch this -- its top bands barely vary.
    """
    evenz = np.array([[5.0, 8.0, 20.0, 40.0],
                      [5.1, 8.1, 20.1, 40.1],
                      [5.2, 8.2, 20.2, 40.2]]) * 1e9
    oddz = np.array([[3.12, 4.02, 6.79, 13.36],
                     [3.13, 4.02, 14.91, 16.44],
                     [1.30, 9.51, 22.53, 33.74]]) * 1e9

    combined, ceiling = _combined_bands(evenz, oddz)
    wrong = min(np.nanmax(evenz[:, -1]), np.nanmax(oddz[:, -1]))   # the old bug
    assert np.isclose(ceiling, 13.36e9), f"ceiling {ceiling/1e9:.2f} GHz"
    assert wrong > 2 * ceiling, "premise: the two reductions must differ here"
    assert combined.size == 0 or np.nanmax(combined) <= ceiling


def test_objective_end_to_end_no_backends():
    # require_mech=False, require_opt=False -> no solver called at all;
    # exercises feasibility gating + score wiring in isolation. Anchored on the
    # MATLAB reference design: the midpoint of the bounds box is NOT feasible
    # (14 nm web), so a hard-coded [0.5]*6 would test the infeasible path here.
    rec = evaluate_candidate(phys_to_u(REFERENCE),
                             require_mech=False, require_opt=False)
    assert rec["status"] == "success", rec.get("reasons")
    assert rec["score"] == 0.0


def test_objective_records_infeasible_without_solving():
    # Midpoint of the bounds box: web = 14.0 nm < min_web. Must be rejected
    # before any backend is touched (mech_backend is deliberately bogus --
    # reaching a solver at all would raise).
    rec = evaluate_candidate([0.5] * N, mech_backend="no-such-backend",
                             require_mech=True, require_opt=False)
    assert rec["status"] == "infeasible"
    assert any("web" in r for r in rec["reasons"])


_BASELINE = os.path.join(os.path.dirname(__file__), "data",
                         "baseline_oddz_reference.csv")


def load_baseline_oddz():
    """[27, 10] odd-z eigenfrequencies (Hz) of the pre-parameterization model.

    See the header of tests/data/baseline_oddz_reference.csv for provenance.
    Compare a rebuilt template's odd-z solve against this at the reference
    design to confirm the geometry survived being parameterized.
    """
    rows = []
    with open(_BASELINE) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("k_norm"):
                continue
            if line.strip():
                rows.append([float(x) for x in line.split(",")])
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1:]


def test_baseline_fixture_is_intact():
    """The fixture is the only surviving record of the pre-parameterization
    frequencies -- the source .mph was deleted. Guard its shape and contents so
    a bad edit or a truncated file is caught rather than silently weakening any
    comparison made against it."""
    k, F = load_baseline_oddz()
    assert F.shape == (27, 10)
    assert np.isclose(k[0], 0.0) and np.isclose(k[-1], 3 - 1 / 9)
    assert np.all(np.diff(F, axis=1) >= 0), "bands must be ascending at each k"
    # Band 0 is the acoustic branch, which goes to zero at Gamma -- but only to
    # within the eigensolver's numerical floor. COMSOL returns 593 Hz here
    # (shift=0 with the filter real(freq)+1e-6>0), i.e. ~5e-8 of band 1. Assert
    # "negligible", not "exactly zero".
    assert F[0, 0] < 1e-5 * F[0, 1], (
        f"band 0 at Gamma is {F[0,0]:.3e} Hz vs band 1 {F[0,1]:.3e} Hz -- "
        f"expected an acoustic branch near zero")
    # Spot-check Gamma against the values quoted in the fixture header.
    assert np.allclose(F[0, 1:4] / 1e9, [11.994, 12.291, 18.558], atol=1e-3)


def test_baseline_reproduces_the_known_gaps():
    """Pins what the baseline says about the reference design, so the
    target-frequency discussion in docs/TODO.md stays anchored to data.

    The largest odd-z gap is 12.29 -> 16.64 GHz (30%), NOT at the 8 GHz
    placeholder in targets_2d.yaml. A second, smaller gap at 5.29 -> 6.42 GHz
    (19%) DOES fall inside the +/-50% window around 8 GHz, so the scorer picks
    that one -- reporting a plausible 19% near-miss instead of the real 30%
    gap. This test exists to keep that trap visible.
    """
    _, F = load_baseline_oddz()

    best = largest_gap(F)
    assert np.isclose(best.f_lower / 1e9, 12.291, atol=1e-3)
    assert np.isclose(best.f_upper / 1e9, 16.636, atol=1e-3)
    assert 0.29 < best.normalized_gap < 0.31

    at8 = gap_near_frequency(F, 8.0e9)
    assert at8.found, "the 8 GHz window is not empty -- it finds the WRONG gap"
    assert np.isclose(at8.f_center / 1e9, 5.855, atol=1e-3)
    assert at8.normalized_gap < best.normalized_gap


def test_bz_path_hits_the_high_symmetry_points():
    """Pin the port of runBands_2D.m:63-64, not just the array shape.

    Hexagonal BZ: M at k=1 is (-pi/(sqrt(3)a), pi/a); K at k=2 is
    (0, 4pi/(3a)); and |M|/|K| = sqrt(3)/2 for a triangular lattice.
    """
    a = 480e-9
    n = 6                                   # so k=1 and k=2 land on samples
    k_norm, kx, ky = bz_path(a, n_per_segment=n, unitcell="hexagonal")
    i_m, i_k = list(k_norm).index(1.0), list(k_norm).index(2.0)

    assert np.isclose(kx[i_m], -np.pi / (np.sqrt(3) * a))
    assert np.isclose(ky[i_m], np.pi / a)
    assert np.isclose(kx[i_k], 0.0)
    assert np.isclose(ky[i_k], 4 * np.pi / (3 * a))

    mag = lambda i: np.hypot(kx[i], ky[i])  # noqa: E731
    assert np.isclose(mag(i_m) / mag(i_k), np.sqrt(3) / 2)


def test_bz_path_is_continuous_across_segment_joins():
    """The piecewise kx/ky branches must agree at k=1 and k=2."""
    a = 480e-9
    k_norm, kx, ky = bz_path(a, n_per_segment=200, unitcell="hexagonal")
    step = np.hypot(np.diff(kx), np.diff(ky))
    assert step.max() < 3 * np.median(step), "discontinuity at a segment join"


if __name__ == "__main__":
    for name, fn in sorted(list(globals().items())):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok  {name}")
    print("all tests passed")
