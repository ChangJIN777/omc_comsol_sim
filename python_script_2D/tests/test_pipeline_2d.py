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


def _two_parity_fixture():
    """(evenz, oddz) where ODDZ is the clear winner near 8 GHz.

    evenz: a narrow 7.90->8.10 GHz gap (2.5%).
    oddz:  a wide   7.00->9.00 GHz gap (25%).
    A third evenz band sits high so both families have the same n_bands.
    """
    k = np.linspace(0, np.pi, 9)
    evenz = np.column_stack([
        7.90e9 - 0.02e9 * np.sin(k),        # band 0 tops out at 7.90
        8.10e9 + 0.02e9 * np.cos(k),        # band 1 bottoms at 8.10  -> 2.5%
        30.0e9 + 0.02e9 * np.sin(k)])
    oddz = np.column_stack([
        7.00e9 - 0.02e9 * np.sin(k),        # tops out at 7.00
        9.00e9 + 0.02e9 * np.cos(k),        # bottoms at 9.00        -> 25%
        31.0e9 + 0.02e9 * np.sin(k)])
    return evenz, oddz


def test_mech_parity_names_the_winning_family():
    """symmetry mode must report WHICH family was scored, plus both families'
    own gaps. Collapsing to one scalar makes a scored candidate impossible to
    judge: gap_mode 'symmetry' is only acceptable because a single parity
    couples to the transducer, so the label is load-bearing."""
    evenz, oddz = _two_parity_fixture()
    gp, info = _mechanical_gap(evenz, oddz, 8.0e9, "symmetry")

    assert info["mech_parity"] == "oddz"
    assert np.isclose(gp.normalized_gap, info["mechanical_gap_oddz"])
    assert 0.24 < info["mechanical_gap_oddz"] < 0.26
    assert 0.02 < info["mechanical_gap_evenz"] < 0.03
    assert np.isclose(info["mechanical_center_frequency_oddz"] / 1e9, 8.0, atol=0.05)
    assert np.isclose(info["mechanical_center_frequency_evenz"] / 1e9, 8.0, atol=0.05)
    # symmetry mode has no stacking, so no truncation ceiling
    assert info["truncation_ceiling_hz"] is None
    assert info["n_bands_usable"] is None


def test_losing_family_with_no_gap_reports_zero_not_nan():
    """A family with no gap in the window must come back 0.0 -- present and
    numeric. NaN would poison any downstream comparison silently, and a missing
    key would make the record's shape depend on the data."""
    evenz, oddz = _two_parity_fixture()
    # flatten evenz into a single dispersive band pair with no gap at all
    k = np.linspace(0, np.pi, evenz.shape[0])
    evenz = np.column_stack([7.0e9 + 2.0e9 * np.sin(k),
                             7.5e9 + 2.0e9 * np.cos(k),
                             30.0e9 + 0.0 * k])
    gp, info = _mechanical_gap(evenz, oddz, 8.0e9, "symmetry")

    assert info["mech_parity"] == "oddz"
    assert info["mechanical_gap_evenz"] == 0.0
    assert info["mechanical_center_frequency_evenz"] == 0.0
    assert not np.isnan(info["mechanical_gap_evenz"])
    assert info["mechanical_gap_oddz"] > 0.2


def test_complete_mode_still_reports_both_families():
    """'complete' sets mech_parity='complete' but keeps the per-family numbers.

    They answer "what would each family have given alone?", which is exactly
    the question when a complete gap comes back empty. They are computed on the
    UNTRUNCATED family arrays -- within one family every mode below its top
    band is known, so a single-family gap needs no ceiling.
    """
    evenz, oddz = _two_parity_fixture()
    gp, info = _mechanical_gap(evenz, oddz, 8.0e9, "complete")

    assert info["mech_parity"] == "complete"
    assert 0.24 < info["mechanical_gap_oddz"] < 0.26     # unchanged by mode
    assert 0.02 < info["mechanical_gap_evenz"] < 0.03
    assert info["truncation_ceiling_hz"] is not None
    assert info["n_bands_usable"] >= 2
    # the union of the two families has no 8 GHz gap: evenz band 1 (8.10) sits
    # inside oddz's 7.00->9.00 window, so the complete gap must be smaller than
    # either family's own -- or absent.
    assert gp.normalized_gap < info["mechanical_gap_oddz"]


def test_targets_yaml_scores_the_complete_gap():
    """The shipped config must score the COMPLETE gap, and the code paths that
    read gap_mode must agree on that default.

    Pinned because the flip from "symmetry" touches four readers (objective2d's
    fallback, the plotter's fallback, and both run scripts' documented
    behaviour); a partial revert would silently score a family-restricted gap
    while the docs claimed otherwise.
    """
    import yaml
    import objective2d
    with open(os.path.join(os.path.dirname(__file__), "..", "configs",
                           "targets_2d.yaml")) as fh:
        targets = yaml.safe_load(fh)
    assert targets["mechanical"]["gap_mode"] == "complete"

    # the in-code fallback must match the shipped config, so a targets file
    # missing the key cannot silently switch modes
    src = open(objective2d.__file__).read()
    assert 'get("gap_mode", "complete")' in src, (
        "objective2d's gap_mode fallback disagrees with targets_2d.yaml")

    # and _mechanical_gap must actually take the complete branch for that value
    evenz, oddz = _two_parity_fixture()
    _, info = _mechanical_gap(evenz, oddz, 8.0e9,
                              targets["mechanical"]["gap_mode"])
    assert info["mech_parity"] == "complete"
    assert info["truncation_ceiling_hz"] is not None, (
        "complete mode must report a truncation ceiling")
    # the per-family diagnostics survive the flip -- they are what explains a
    # narrow complete gap, so losing them would be a regression, not tidying
    for p in ("evenz", "oddz"):
        for key in (f"mechanical_gap_{p}", f"mechanical_center_frequency_{p}",
                    f"mechanical_gap_lower_frequency_{p}",
                    f"mechanical_gap_upper_frequency_{p}"):
            assert key in info, f"{key} missing in complete mode"


def test_default_plot_config_still_renders_after_the_flip():
    """configs/plot_bands_2d.yaml must keep working on the odd-z-only baseline.

    Its default input is single-parity, and complete mode deliberately raises
    on that -- so the config pins gap_mode: symmetry. If that pin is removed
    while targets_2d.yaml says complete, the only no-COMSOL demo in the project
    dies via a guard that is behaving correctly. Assert the pin is present AND
    that the config actually renders.
    """
    import tempfile
    import yaml
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import plot_bands_2d as pb

    cfg_path = os.path.join(os.path.dirname(__file__), "..", "configs",
                            "plot_bands_2d.yaml")
    cfg = yaml.safe_load(open(cfg_path))
    assert cfg.get("kind", "baseline-csv") == "baseline-csv"
    assert cfg.get("gap_mode") == "symmetry", (
        "the default plot config must pin gap_mode: symmetry while its input is "
        "the odd-z-only baseline CSV")

    k_norm, families, a = pb.load_baseline_csv(cfg.get("csv"))
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "default_cfg.png")
        pb.render(k_norm, families, out, a=a, gap_mode=cfg["gap_mode"],
                  f_target_hz=8.0e9)
        assert os.path.getsize(out) > 5000


def _bands_from_ranges(ranges, n_k=13):
    """Band matrix whose column [min, max] ranges are exactly `ranges` (in Hz).

    Each band is a raised cosine across the path, so max_k / min_k land on the
    requested endpoints and largest_gap/all_gaps see the intended gaps.
    """
    k = np.linspace(0, np.pi, n_k)
    cols = [lo + (hi - lo) * 0.5 * (1 - np.cos(k)) for lo, hi in ranges]
    return np.sort(np.column_stack(cols), axis=1)


def _pairwise_intersections(evenz, oddz, cap):
    """All (evenz gap_i) ∩ (oddz gap_j), clipped above at `cap`."""
    from bandgap import all_gaps
    out = []
    for a in all_gaps(evenz):
        for b in all_gaps(oddz):
            lo = max(a.f_lower, b.f_lower)
            hi = min(a.f_upper, b.f_upper, cap)
            if hi > lo:
                out.append((lo, hi))
    return out


def test_complete_gap_contains_every_pairwise_intersection():
    """The complete-gap search returns a SUPERSET of the pairwise overlaps.

    Encodes the claim in _combined_bands' EQUIVALENT FORMULATIONS note rather
    than one example: over random band pairs (flat and dispersive, unequal band
    counts), every intersection of an evenz gap with an oddz gap that lies below
    min_k c_N must be contained in exactly ONE gap of the combined spectrum.
    The min_k c_N qualifier is load-bearing and is NOT the ceiling -- a band
    straddling the ceiling is dropped, taking any gap that ended at its bottom
    with it (see test_ceiling_can_hide_a_real_overlap).
    """
    from bandgap import all_gaps
    from objective2d import _combined_bands

    rng = np.random.default_rng(20260812)
    checked = 0
    for _ in range(400):
        n_k = int(rng.integers(3, 9))
        flat = bool(rng.integers(0, 2))
        def fam(n_b):
            base = np.cumsum(rng.uniform(0.3, 2.0, n_b)) * 1e9
            if flat:
                F = np.tile(base, (n_k, 1))
            else:
                F = base * (1 + 0.25 * rng.standard_normal((n_k, n_b)))
            return np.sort(F, axis=1)
        evenz, oddz = fam(int(rng.integers(2, 7))), fam(int(rng.integers(2, 7)))

        combined, _ceiling = _combined_bands(evenz, oddz)
        if combined.shape[1] < 2:
            continue
        cap = float(combined[:, -1].min())          # min_k c_N
        found = [(g.f_lower, g.f_upper) for g in all_gaps(combined)]

        for lo, hi in _pairwise_intersections(evenz, oddz, cap):
            checked += 1
            covering = [(x, y) for x, y in found
                        if x <= lo * (1 + 1e-9) and hi <= y * (1 + 1e-9)]
            assert covering, (
                f"pairwise overlap ({lo/1e9:.4f}, {hi/1e9:.4f}) GHz is not in "
                f"any combined gap {[(round(x/1e9,4), round(y/1e9,4)) for x, y in found]}")
            assert len(covering) == 1, (
                f"overlap ({lo/1e9:.4f}, {hi/1e9:.4f}) GHz spans {len(covering)} "
                f"combined gaps; the bracketing band index must be unique")
    assert checked > 200, f"fixture generated too few overlaps to be meaningful ({checked})"


def test_complete_gap_is_a_STRICT_superset_of_the_overlaps():
    """The extra windows are real, and the overlap recipe cannot express them.

    A window can be free of a family because it lies BELOW that family's lowest
    band, not only because it is inside one of that family's gaps. Here evenz
    has no gap at all and starts at 10 GHz, so the pairwise recipe yields
    nothing -- yet 3 -> 5 GHz is genuinely free of both families.
    """
    from bandgap import all_gaps
    from objective2d import _combined_bands

    evenz = _bands_from_ranges([(10e9, 14e9), (14e9, 18e9), (18e9, 40e9)])
    oddz = _bands_from_ranges([(1e9, 3e9), (5e9, 8e9), (8e9, 41e9)])

    assert all_gaps(evenz) == [], "fixture broken: evenz must have no gap"
    combined, _ = _combined_bands(evenz, oddz)
    assert _pairwise_intersections(evenz, oddz, np.inf) == []

    found = all_gaps(combined)
    assert len(found) == 1
    assert np.isclose(found[0].f_lower / 1e9, 3.0, atol=1e-6)
    assert np.isclose(found[0].f_upper / 1e9, 5.0, atol=1e-6)


def test_an_extra_window_can_BEAT_the_overlap():
    """The superset is not a technicality: an "extra" can be the LARGEST
    complete gap, so the overlap recipe can badly understate the answer.

    evenz gaps 7.90->8.08 (2.3%), oddz gaps 7.00->8.98 (25%). Their overlap is
    7.90->8.08 -- only 2.3%. But 7.00->7.88 is also free of both: it is inside
    oddz's gap and below evenz's LOWEST band (which bottoms at 7.88). At 11.8%
    that is 5x the overlap, and it is what the code scores.
    """
    from bandgap import all_gaps
    from objective2d import _combined_bands

    evenz, oddz = _two_parity_fixture()
    combined, _ = _combined_bands(evenz, oddz)
    edges = {(round(g.f_lower / 1e9, 2), round(g.f_upper / 1e9, 2)): g
             for g in all_gaps(combined)}

    overlap = (7.90, 8.08)
    extra = (7.00, 7.88)
    assert overlap in edges, f"the overlap must still be found; got {sorted(edges)}"
    assert extra in edges, f"the below-evenz window must be found; got {sorted(edges)}"
    assert edges[extra].normalized_gap > 4 * edges[overlap].normalized_gap

    scored, info = _mechanical_gap(evenz, oddz, 8.0e9, "complete")
    assert np.isclose(scored.f_lower / 1e9, extra[0], atol=1e-2)
    assert np.isclose(scored.f_upper / 1e9, extra[1], atol=1e-2)
    # ... and the overlap is exactly what the recorded per-family edges give
    ov_lo = max(info["mechanical_gap_lower_frequency_evenz"],
                info["mechanical_gap_lower_frequency_oddz"])
    ov_hi = min(info["mechanical_gap_upper_frequency_evenz"],
                info["mechanical_gap_upper_frequency_oddz"])
    assert np.isclose(ov_lo / 1e9, overlap[0], atol=1e-2)
    assert np.isclose(ov_hi / 1e9, overlap[1], atol=1e-2)


def test_ceiling_can_hide_a_real_overlap():
    """One-sided conservatism: truncation can drop an overlap, never invent one.

    evenz bands 4-5, 6-6.4, 12-22 and oddz 4.5-5, 6.8-7.2, 24-30 give a genuine
    7.2 -> 12 GHz overlap. But evenz's top band bottoms at 12, so the ceiling is
    12 and the `keep` mask drops every column above 6.8 GHz -- the overlap is
    NOT reported. Adding a band above it to each family brings it back. The fix
    for a missing high gap is more `neigs`, not a different reduction.
    """
    from bandgap import all_gaps
    from objective2d import _combined_bands

    lo_e = [(4e9, 5e9), (6e9, 6.4e9)]
    lo_o = [(4.5e9, 5e9), (6.8e9, 7.2e9)]

    hidden, ceil_h = _combined_bands(
        _bands_from_ranges(lo_e + [(12e9, 22e9)]),
        _bands_from_ranges(lo_o + [(24e9, 30e9)]))
    assert np.isclose(ceil_h / 1e9, 12.0, atol=1e-6)
    edges = [(round(g.f_lower / 1e9, 2), round(g.f_upper / 1e9, 2))
             for g in all_gaps(hidden)]
    assert (7.2, 12.0) not in edges, f"expected the overlap to be hidden, got {edges}"
    assert (5.0, 6.0) in edges and (6.4, 6.8) in edges   # the low ones survive

    shown, ceil_s = _combined_bands(
        _bands_from_ranges(lo_e + [(12e9, 18e9), (30e9, 32e9)]),
        _bands_from_ranges(lo_o + [(24e9, 26e9), (31e9, 33e9)]))
    assert ceil_s > ceil_h
    edges = [(round(g.f_lower / 1e9, 2), round(g.f_upper / 1e9, 2))
             for g in all_gaps(shown)]
    assert (7.2, 12.0) in edges, f"overlap should reappear, got {edges}"


def test_gap_edges_are_recorded_for_the_scored_gap_and_both_families():
    """Edges in Hz, so an overlap is readable straight from a record.

    Also pins the algebraic identity that made them merely 'reconstructible':
    lower/upper == centre*(1 -+ G/2) for THIS project's normalization.
    """
    evenz, oddz = _two_parity_fixture()
    gp, info = _mechanical_gap(evenz, oddz, 8.0e9, "symmetry")

    assert np.isclose(info["mechanical_gap_lower_frequency"], gp.f_lower)
    assert np.isclose(info["mechanical_gap_upper_frequency"], gp.f_upper)
    for p in ("evenz", "oddz"):
        lo = info[f"mechanical_gap_lower_frequency_{p}"]
        hi = info[f"mechanical_gap_upper_frequency_{p}"]
        c = info[f"mechanical_center_frequency_{p}"]
        G = info[f"mechanical_gap_{p}"]
        assert hi > lo > 0
        assert np.isclose(lo, c * (1 - G / 2)) and np.isclose(hi, c * (1 + G / 2))

    # the point of storing them: the overlap is now direct arithmetic
    ov_lo = max(info["mechanical_gap_lower_frequency_evenz"],
                info["mechanical_gap_lower_frequency_oddz"])
    ov_hi = min(info["mechanical_gap_upper_frequency_evenz"],
                info["mechanical_gap_upper_frequency_oddz"])
    assert ov_hi > ov_lo, "fixture's two families should overlap near 8 GHz"

    # a family with no gap reports 0.0 across all four of its fields
    k = np.linspace(0, np.pi, evenz.shape[0])
    gapless = np.column_stack([7.0e9 + 2.0e9 * np.sin(k),
                               7.5e9 + 2.0e9 * np.cos(k), 30.0e9 + 0.0 * k])
    _, info = _mechanical_gap(gapless, oddz, 8.0e9, "symmetry")
    for key in ("mechanical_gap_evenz", "mechanical_center_frequency_evenz",
                "mechanical_gap_lower_frequency_evenz",
                "mechanical_gap_upper_frequency_evenz"):
        assert info[key] == 0.0, f"{key} should be 0.0, got {info[key]!r}"


def test_parity_tie_breaks_to_evenz_deterministically():
    """An exact tie is plausible for a near-z-symmetric structure. Pin the
    winner so the label cannot flip between runs of the same design."""
    from objective2d import best_parity, PARITY_ORDER
    evenz, _ = _two_parity_fixture()
    gp, info = _mechanical_gap(evenz, evenz.copy(), 8.0e9, "symmetry")
    assert info["mech_parity"] == "evenz"
    assert info["mechanical_gap_evenz"] == info["mechanical_gap_oddz"]
    assert PARITY_ORDER[0] == "evenz", "tie-break order is part of the contract"
    # and the helper itself, on a subset
    only_odd = {"oddz": gap_near_frequency(evenz, 8.0e9)}
    assert best_parity(only_odd)[0] == "oddz"


def test_plotter_and_record_agree_on_the_winning_parity():
    """The figure legend and the result record must never name different
    families. Both go through objective2d.best_parity; assert that holds for
    data where evenz wins, which is the reverse of the other fixture."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import plot_bands_2d as pb

    oddz, evenz = _two_parity_fixture()          # swapped: evenz now wins
    _, info = _mechanical_gap(evenz, oddz, 8.0e9, "symmetry")
    sel = pb.select_gaps({"evenz": evenz, "oddz": oddz}, 8.0e9, 0.5, "symmetry")

    assert info["mech_parity"] == "evenz"
    assert sel["scored_label"].startswith(info["mech_parity"])
    assert np.isclose(sel["scored"].normalized_gap,
                      info[f"mechanical_gap_{info['mech_parity']}"])


def test_save_bands_keeps_scalar_lattice_constant():
    """save_bands must not drop `a`.

    run_mechanical_comsol_2d returns `a` as a plain Python float, and an
    earlier version filtered the dict on isinstance(v, np.ndarray), so the
    lattice constant vanished from every .npz. scripts/plot_bands_2d.py reads
    it for the BZ annotation, and a missing key is invisible until something
    downstream raises.
    """
    import tempfile
    from acoustic_comsol_2d import save_bands

    k_norm, kx, ky = bz_path(480e-9, n_per_segment=3)
    data = dict(k_norm=k_norm, kx=kx, ky=ky, a=480e-9,
                freqs_oddz=np.zeros((len(k_norm), 4)),
                skipped_because_not_numeric={"nope": 1})
    with tempfile.TemporaryDirectory() as td:
        # nested path exercises the makedirs branch
        path = save_bands(data, os.path.join(td, "sub", "bands.npz"))
        assert os.path.isfile(path)
        with np.load(path) as d:          # no allow_pickle: nothing may pickle
            assert "a" in d, f"lattice constant dropped; keys={list(d.keys())}"
            assert np.isclose(float(d["a"]), 480e-9)
            assert d["freqs_oddz"].shape == (len(k_norm), 4)
            assert "skipped_because_not_numeric" not in d


def test_plot_bands_2d_renders_the_baseline_csv():
    """Smoke-test the plotter end to end on real data, with no COMSOL.

    The baseline CSV is odd-z only, so this also covers the single-family path.
    Asserts a non-empty PNG and that the reported gaps are the two the baseline
    is known to contain -- so the figure cannot silently start shading nothing.
    """
    import tempfile
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import plot_bands_2d as pb

    k_norm, families, a = pb.load_baseline_csv()
    assert sorted(families) == ["oddz"] and np.isclose(a, 480e-9)

    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "figs", "bands_2d.png")
        info = pb.render(k_norm, families, out, a=a, f_target_hz=8.0e9,
                         rel_tol=0.5, gap_mode="symmetry")
        assert os.path.isfile(out)
        assert os.path.getsize(out) > 5000, "PNG is suspiciously small"

    # scored gap is the one inside the [4, 12] GHz window (5.86 GHz, ~19%);
    # the largest gap anywhere is the 30% one at 14.46 GHz, outside it. Both
    # must be reported, and they must be different -- that difference is the
    # trap the figure exists to show.
    assert np.isclose(info["scored"].f_center / 1e9, 5.855, atol=1e-3)
    assert np.isclose(info["best"].f_center / 1e9, 14.464, atol=1e-3)
    assert info["best"].normalized_gap > info["scored"].normalized_gap
    assert info["ceiling_hz"] is None       # symmetry mode has no ceiling


def test_plot_bands_2d_renders_complete_mode_with_ceiling():
    """Cover the two-parity 'complete' path, including the ceiling line.

    Synthesizes an even-z family from the baseline by shifting it, so both
    parities exist. The point is that render() reports a finite truncation
    ceiling and draws it -- a gap above that line is an artifact of finite
    n_bands, and the line is what makes it recognizable in the figure.
    """
    import tempfile
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import plot_bands_2d as pb

    k_norm, families, a = pb.load_baseline_csv()
    oddz = families["oddz"]
    families = {"oddz": oddz, "evenz": oddz * 0.83 + 0.4e9}

    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "complete.png")
        info = pb.render(k_norm, families, out, a=a, f_target_hz=8.0e9,
                         gap_mode="complete")
        assert os.path.isfile(out) and os.path.getsize(out) > 5000

    assert info["ceiling_hz"] is not None
    assert info["n_bands_usable"] >= 2
    ceiling = info["ceiling_hz"]
    # the ceiling must be min over k of min(top_evenz, top_oddz) -- never above
    # either family's lowest top band, or it would license artifact gaps
    assert ceiling <= min(oddz[:, -1].min(), families["evenz"][:, -1].min()) + 1.0
    for key in ("scored", "best"):
        if info[key] is not None:
            assert info[key].f_upper <= ceiling + 1.0, (
                f"{key} gap edge {info[key].f_upper:.3e} is above the "
                f"truncation ceiling {ceiling:.3e}")


def test_plot_bands_2d_refuses_complete_mode_on_one_parity():
    """'complete' with a single family must fail loudly, not fall back.

    A complete gap is one no band of EITHER parity crosses; evaluating it from
    one family would silently report a symmetry gap under the wrong name.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import plot_bands_2d as pb

    _, families, _ = pb.load_baseline_csv()
    try:
        pb.select_gaps(families, 8.0e9, 0.5, "complete")
    except SystemExit as exc:
        assert "both parities" in str(exc)
    else:
        raise AssertionError("complete mode accepted single-parity data")


def test_plot_bands_2d_closes_the_bz_loop():
    """The COMSOL sweep stops at 3 - 1/kpts; the plot must reach Gamma.

    runBands_2D.m:556 closes the path by copying row 0. Duplicating an existing
    row cannot change any max_k/min_k, so gaps are unaffected -- assert that.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import plot_bands_2d as pb

    k_norm, families, _ = pb.load_baseline_csv()
    assert k_norm[-1] < 3.0
    k2, f2 = pb.close_bz_loop(k_norm, families)
    assert np.isclose(k2[-1], 3.0) and len(k2) == len(k_norm) + 1
    assert np.allclose(f2["oddz"][-1], families["oddz"][0])
    assert np.isclose(largest_gap(f2["oddz"]).normalized_gap,
                      largest_gap(families["oddz"]).normalized_gap)
    # idempotent once already closed
    k3, f3 = pb.close_bz_loop(k2, f2)
    assert len(k3) == len(k2)


if __name__ == "__main__":
    for name, fn in sorted(list(globals().items())):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok  {name}")
    print("all tests passed")
