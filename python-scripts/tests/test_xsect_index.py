"""Tests of the cross-section switch (xsect: rect | wedge), the optical
thickness-parity rule that follows from it, and the diamond index 2.386.

No COMSOL, no MPh session and no Tidy3D cloud: COMSOL models are fakes, and
`tidy3d` is a MagicMock (tests/test_optical_tidy3d.fake_tidy3d).

Run:  python tests/test_xsect_index.py
"""
import json
import os
import re
import sys
import tempfile
from unittest import mock

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

import numpy as np

import comsol_client
import comsol_optics
import geometry
import objective
import objective_rib
import optical_tidy3d as ot
from geometry import Geometry, RibGeometry, phys_to_u, u_to_geometry
from objective import evaluate_candidate
from test_optical_tidy3d import (fake_tidy3d, _fake_run, _pad, _resonances,
                                 RIB_CELL, HOLE_U, F_T, C0)

N = 2.386
SQ3 = np.sqrt(3.0)


def _rib_u0():
    return phys_to_u(RIB_CELL, celltype="rib")


def _raises(exc, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except exc as err:
        return err
    raise AssertionError(f"{fn.__name__} did not raise {exc.__name__}")


# ---- cross-section: definition, templates, feasibility ---------------------------

def test_wedge_definition():
    assert geometry.XSECTS == ("rect", "wedge") and geometry.WEDGE_ANGLE_DEG == 30.0
    assert geometry.resolve_xsect(None) == "rect"
    _raises(ValueError, geometry.resolve_xsect, "triangle")
    # apex b/(2*sqrt(3)) below the slab: Polygon 1 of omc_unitcell_iso.mph
    assert abs(geometry.wedge_depth(810e-9) - 810e-9 / (2 * SQ3)) < 1e-18
    hole = Geometry(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9)
    assert geometry.wedge_depths(hole) == {"beam": geometry.wedge_depth(hole.w)}
    rib = RibGeometry(**RIB_CELL)
    d = geometry.wedge_depths(rib)
    assert abs(d["spine"] - rib.s / (2 * SQ3)) < 1e-18       # plain triangles
    assert abs(d["rib"] - rib.t / (2 * SQ3)) < 1e-18
    # the roof is deepest under the junction centre (cone about the root arc)
    junction = (np.hypot(rib.t / 2 + rib.r1, rib.s / 2 + rib.r1) - rib.r1) / SQ3
    assert abs(d["junction"] - junction) < 1e-18 and d["max"] == d["junction"]
    assert d["max"] > d["rib"] > d["spine"]
    # the mirror plane, and with it the optical parity rule
    assert geometry.z_mirror_symmetric("rect") and geometry.z_mirror_symmetric(None)
    assert not geometry.z_mirror_symmetric("wedge")


def test_template_choice_by_celltype_and_xsect():
    name = lambda *a: os.path.basename(comsol_client.template_path(*a))  # noqa: E731
    assert name("hole", "rect") == "omc_unitcell.mph"
    assert name("hole", "wedge") == "omc_unitcell_iso.mph"
    assert name("rib", "rect") == name("rib") == "omc_unitcell_rib.mph"
    assert name("rib", "wedge") == "omc_unitcell_rib_wedge.mph"
    assert name("hole") == "omc_unitcell.mph"            # default xsect = rect
    _raises(ValueError, comsol_client.template_path, "hole", "round")
    _raises(ValueError, comsol_client.template_path, "disk", "rect")
    assert comsol_client.THICKNESS_PARAM == {"hole": "t", "rib": "th"}
    # a template that is not built: clear error with the build step
    real = comsol_client._COMSOL_DIR
    comsol_client._COMSOL_DIR = os.path.join(real, "no_such_dir")
    try:
        err = _raises(FileNotFoundError, comsol_client.template_path, "rib",
                      "wedge", True)
        assert "build_omc_unitcell_rib('wedge')" in str(err)
        assert "omc_unitcell_rib_wedge.mph" in str(err)
        comsol_client.template_path("rib", "wedge")       # no check: no error
    finally:
        comsol_client._COMSOL_DIR = real
    # a rib wedge template of the first (straight-prism) builder is refused
    import zipfile
    with tempfile.TemporaryDirectory() as tmp:
        comsol_client._COMSOL_DIR = tmp
        try:
            path = os.path.join(tmp, "omc_unitcell_rib_wedge.mph")
            with zipfile.ZipFile(path, "w") as z:
                z.writestr("dmodel.xml", 's("wd") s("max(s,t)/2*tan(wedge_angle)")')
            assert comsol_client.rib_wedge_is_roof(path) is False
            err = _raises(RuntimeError, comsol_client.template_path, "rib",
                          "wedge", True)
            assert "rebuild" in str(err) and "roof" in str(err)
            os.utime(path, (1, 1))                        # a new file: re-read
            with zipfile.ZipFile(path, "w") as z:
                z.writestr("dmodel.xml", 's("wdj") s("tan(wedge_angle)*(...)")')
            assert comsol_client.rib_wedge_is_roof(path) is True
            assert comsol_client.template_path("rib", "wedge", True) == path
            with open(path, "w") as fh:
                fh.write("not a zip")
            os.utime(path, (5, 5))
            assert comsol_client.rib_wedge_is_roof(path) is True   # unreadable: allowed
        finally:
            comsol_client._COMSOL_DIR = real


def test_hole_templates_are_what_the_table_says():
    """Read the shipped .mph files (zip): which one has the wedge polygon."""
    import zipfile
    seen = {}
    for xsect in ("rect", "wedge"):
        path = comsol_client.template_path("hole", xsect)
        if not os.path.exists(path):
            print(f"  (skipped: {os.path.basename(path)} not present)")
            return
        with zipfile.ZipFile(path) as z:
            seen[xsect] = z.read("dmodel.xml").decode("utf-8", "replace")
    assert "w/2/sqrt(3)" in seen["wedge"] and "pol1" in seen["wedge"]
    assert "sqrt(3)" not in seen["rect"]
    assert "/feature/pol1" not in seen["rect"] and "/feature/ext1" not in seen["rect"]


def test_wedge_feasibility_rule():
    ok_cell = Geometry(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9)
    assert geometry.check_feasibility(ok_cell, xsect="wedge") == (True, [])
    blind = Geometry(a=891e-9, w=900e-9, t=150e-9, hx=218e-9, hy=334e-9)
    assert geometry.check_feasibility(blind)[0] is True          # rect: fine
    ok, reasons = geometry.check_feasibility(blind, xsect="wedge")
    assert ok is False and reasons[0].startswith("WEDGE_HOLE_BLIND")
    # exactly the limit t = w/(2*sqrt(3))
    edge = Geometry(a=891e-9, w=900e-9, t=900e-9 / (2 * SQ3) + 1e-12,
                    hx=218e-9, hy=334e-9)
    assert geometry.check_feasibility(edge, xsect="wedge")[0] is True
    # the rib wedge adds no rule
    rib = RibGeometry(**RIB_CELL)
    assert geometry.check_feasibility(rib, xsect="wedge") == (True, [])
    # ... and the objective applies it before any solver
    u = phys_to_u(blind.as_dict())
    rec = evaluate_candidate(u, optical_backend="comsol", mech_backend="comsol",
                             xsect="wedge")
    assert rec["status"] == "infeasible" and rec["xsect"] == "wedge"
    assert any(r.startswith("WEDGE_HOLE_BLIND") for r in rec["reasons"])


def test_mpb_and_surrogate_refuse_wedge_before_any_solve():
    import optical_surrogate
    calls = []
    real = optical_surrogate.optical_gap_surrogate
    optical_surrogate.optical_gap_surrogate = lambda g, **kw: calls.append(g)
    try:
        for backend in ("surrogate", "mpb"):
            err = _raises(ValueError, evaluate_candidate, HOLE_U,
                          optical_backend=backend, mech_backend="surrogate_stub",
                          xsect="wedge")
            assert "rectangular" in str(err) and backend in str(err)
        _raises(ValueError, evaluate_candidate, HOLE_U,
                optical_backend="surrogate", mech_backend="surrogate_stub",
                xsect="trapezoid")
    finally:
        optical_surrogate.optical_gap_surrogate = real
    assert calls == []
    assert objective.check_xsect_backend(None, "mpb") == "rect"
    assert objective.check_xsect_backend("wedge", "comsol") == "wedge"
    assert objective.check_xsect_backend("wedge", "tidy3d") == "wedge"
    rec = evaluate_candidate(HOLE_U, optical_backend="surrogate",
                             mech_backend="surrogate_stub", xsect="rect")
    assert rec["status"] == "success" and rec["xsect"] == "rect"


def test_xsect_in_record_and_id():
    import acoustic_comsol
    import optical_comsol
    seen = []
    k15 = np.linspace(1e-3, 1.0, 15)

    def fake_optics(g, **kw):
        seen.append(("opt", kw.get("xsect")))
        return dict(freqs_hz=np.column_stack([150e12 + 25e12 * k15,
                                              230e12 - 20e12 * k15]),
                    k_frac=k15, a=g.a, z_family=None, freqs_all_hz=[[1.0]],
                    z_even_fraction=[[0.9]])

    def fake_mech(g, *, n_k=15, n_bands=12, template=None, study="mech sym", **kw):
        seen.append(("mech", os.path.basename(template)))
        k = np.linspace(0, 1, 5)
        return dict(freqs_hz=np.column_stack([5e9 + 0 * k, 9e9 + 0 * k]),
                    k_frac=k, a=g.a)

    real = (optical_comsol.run_optical_comsol, acoustic_comsol.run_mechanical_comsol,
            ot.run_optical_tidy3d, objective_rib._solve_optical_bands)
    optical_comsol.run_optical_comsol = fake_optics
    acoustic_comsol.run_mechanical_comsol = fake_mech
    ot.run_optical_tidy3d = _fake_run([])
    objective_rib._solve_optical_bands = (
        lambda g, study, k_frac, n_bands, f_target, **kw:
        seen.append(("rib-opt", kw.get("xsect"))) or _pad(_resonances(k_frac, g.a)))
    # a hole cell whose hole goes through the wedge
    u = phys_to_u(dict(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9))
    exists = os.path.exists
    try:
        # templates: pretend every one is there and current (no file is read)
        with mock.patch("os.path.exists", lambda p: True if str(p).endswith(".mph")
                        else exists(p)), \
                mock.patch.object(comsol_client, "rib_wedge_is_roof", lambda p: True):
            recs = {}
            for xs in ("rect", "wedge"):
                recs["hole", "comsol", xs] = evaluate_candidate(
                    u, optical_backend="comsol", mech_backend="comsol", xsect=xs)
                recs["hole", "tidy3d", xs] = evaluate_candidate(
                    u, optical_backend="tidy3d", mech_backend="surrogate_stub",
                    xsect=xs)
                recs["rib", "comsol", xs] = evaluate_candidate(
                    _rib_u0(), celltype="rib", optical_backend="comsol",
                    mech_backend="comsol", xsect=xs,
                    rib_options=dict(mech_gap_mode="sym"))
                recs["rib", "none", xs] = evaluate_candidate(
                    _rib_u0(), celltype="rib", optical_backend="none",
                    mech_backend="surrogate_stub", xsect=xs)
            default = evaluate_candidate(u, optical_backend="comsol",
                                         mech_backend="comsol")
    finally:
        (optical_comsol.run_optical_comsol, acoustic_comsol.run_mechanical_comsol,
         ot.run_optical_tidy3d, objective_rib._solve_optical_bands) = real

    assert all(r["status"] in ("success", "no_solver") for r in recs.values()), {
        k: (r["status"], r.get("error")) for k, r in recs.items()}
    # the key is in every record, and rect / wedge never share an id
    for (cell, engine, xs), r in recs.items():
        assert r["xsect"] == xs and r["optical_z_symmetry"] == (xs == "rect")
        json.dumps(r, allow_nan=False)
    ids = [r["id"] for r in recs.values()]
    assert len(set(ids)) == len(ids)
    assert default["id"] == recs["hole", "comsol", "rect"]["id"]   # default = rect
    # one device: optics AND mechanics use the template of the cross-section
    assert seen[:2] == [("opt", "rect"), ("mech", "omc_unitcell.mph")]
    assert ("opt", "wedge") in seen and ("mech", "omc_unitcell_iso.mph") in seen
    assert ("rib-opt", "wedge") in seen
    assert ("mech", "omc_unitcell_rib.mph") in seen
    assert ("mech", "omc_unitcell_rib_wedge.mph") in seen
    # per-mode parity goes into the optical detail
    d = recs["hole", "comsol", "rect"]["optical_comsol_detail"]
    assert d["z_even_fraction"] == [[0.9]] and d["freqs_all_hz"] == [[1.0]]
    assert recs["rib", "comsol", "rect"]["optical_comsol_detail"]["z_family"] == "even"
    assert recs["rib", "comsol", "wedge"]["optical_comsol_detail"]["z_family"] is None
    assert recs["rib", "comsol", "wedge"]["wedge_depths"]["rib"] > 0
    # a rectangular rib cell without optics: the old tuple plus the score
    # version (version 3 changed the mechanical scoring too)
    legacy = objective_rib._rib_id(_rib_u0(), (
        "none", "surrogate_stub", "complete", 15, 20, "opt TE", 13, 10, 0.30,
        "score=v4"))
    assert recs["rib", "none", "rect"]["id"] == legacy
    assert recs["rib", "none", "wedge"]["id"] != legacy
    assert recs["rib", "none", "rect"]["optical_index"] is None
    # the id tag itself
    tag = objective.physics_id_tag
    assert tag("none", "rect") is None and tag("none", "wedge") == "xsect=wedge"
    assert tag("comsol", "rect", "even") == "n=2.3860|xsect=rect|z=even"
    assert tag("comsol", "wedge", None) == "n=2.3860|xsect=wedge|z=all"
    assert tag("mpb", "rect") == "n=2.3860|xsect=rect"
    assert tag("tidy3d", "wedge", None, "T") == "T"


def test_hole_comsol_loop_ids_and_default():
    import matplotlib
    matplotlib.use("Agg")
    import run_opt_comsol as loop
    p = dict(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9)
    legacy = loop._make_id(p)                      # what the loop always wrote
    assert loop._make_id(p, "wedge", 2.40) == legacy
    ids = {loop._make_id(p, xs, n) for xs in ("rect", "wedge") for n in (2.40, N)}
    assert len(ids) == 4 and legacy in ids
    assert comsol_client.LEGACY_XSECT == {"hole": "wedge", "rib": "rect"}

    # evaluate(): the parity family follows xsect; mechanics is untouched
    calls = []

    def fake_opt_sweep(model, g, study, n_k, n_bands, min_THz=10.0,
                       z_family=None, thickness=None, info=None):
        calls.append((study, z_family, thickness))
        info.setdefault("z_even_fraction", []).append([0.97, 0.02])
        info.setdefault("freqs_all_hz", []).append([170e12, 215e12])
        return (np.array([1.0]), np.array([[170e12, 215e12]]))

    def fake_mech_sweep(model, g, study, n_k, n_bands, compute_fy=True):
        k = np.linspace(0, 1, 3)
        f = np.column_stack([5e9 + 0 * k, 9e9 + 0 * k])
        return k, f, np.ones_like(f)

    model = mock.MagicMock()
    real = (loop._opt_sweep, loop._mech_sweep)
    loop._opt_sweep, loop._mech_sweep = fake_opt_sweep, fake_mech_sweep
    try:
        g = Geometry(**p)
        kw = dict(n_k=3, n_bands_mech=2, n_bands_opt=2, study_mech="mech sym",
                  study_opt="opt TE")
        r_legacy = loop.evaluate(model, [0.5] * 5, g, **kw)
        r_rect = loop.evaluate(model, [0.5] * 5, g, xsect="rect", **kw)
        r_wedge = loop.evaluate(model, [0.5] * 5, g, xsect="wedge", **kw)
    finally:
        loop._opt_sweep, loop._mech_sweep = real
    assert calls == [("opt TE", None, "t"), ("opt TE", "even", "t"),
                     ("opt TE", None, "t")]
    assert r_legacy["xsect"] == "wedge" and r_legacy["id"] == r_wedge["id"]
    assert r_rect["id"] != r_wedge["id"] and r_wedge["id"] != legacy  # n changed
    assert r_wedge["id"] == loop._make_id(p, "wedge", N)
    assert r_rect["optical_z_family"] == "even" and r_wedge["optical_z_family"] is None
    assert r_rect["optical_index"] == N and r_rect["_z_even_fraction"] == [[0.97, 0.02]]
    assert r_rect["mechanical_gap"] == r_wedge["mechanical_gap"] > 0


# ---- Tidy3D: symmetry tuple and wedge geometry --------------------------------

def test_tidy3d_symmetry_tuple_per_xsect_and_mode():
    # +1 = even = PMC, -1 = odd = PEC, 0 = none (Tidy3D convention)
    exp = ot.expected_symmetry
    assert exp("opt TE", "rect") == (0, -1, 1)     # TE-like: even across thickness
    assert exp("opt TM", "rect") == (0, 1, -1)
    assert exp("opt TE", "wedge") == (0, -1, 0)    # no mid-thickness mirror
    assert exp("opt TM", "wedge") == (0, 1, 0)
    assert exp("opt TE") == (0, -1, 1)             # default xsect = rect
    with fake_tidy3d() as (td, bs):
        rib = RibGeometry(**RIB_CELL)
        hole = u_to_geometry(HOLE_U)
        for g in (rib, hole):
            for study in ("opt TE", "opt TM"):
                for xs in ("rect", "wedge"):
                    kw = ot.band_params_kwargs(g, index=N, f_target=F_T,
                                               optical_study=study, xsect=xs)
                    p = bs.BandParams(**kw)
                    assert tuple(p.symmetry) == exp(study, xs), (study, xs)
                    assert p.xsect == xs and p.wedge_angle == 30.0
        # the module default and its command line stay rectangular
        assert bs.BandParams().xsect == "rect"
        assert tuple(bs.BandParams().symmetry) == (0, -1, 1)
        args = bs.build_parser().parse_args(["--xsect", "wedge", "--mode", "TM"])
        p = bs.params_from_args(args)
        assert p.xsect == "wedge" and tuple(p.symmetry) == (0, 1, 0)
        assert bs.params_from_args(bs.build_parser().parse_args([])).xsect == "rect"
        _raises(ValueError, bs.build_unit_cell,
                bs.BandParams(geometry="rib_cos6", xsect="wedge"))
        _raises(ValueError, bs.build_unit_cell, bs.BandParams(xsect="round"))


def _polyslabs(td):
    return [c.kwargs for c in td.PolySlab.call_args_list]


def test_tidy3d_wedge_geometry_uses_the_repo_definition():
    um = 1e6
    with fake_tidy3d() as (td, bs):
        # ---- hole: one prism along the beam, base w; hole through both ----
        hole = Geometry(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9)
        depth = geometry.wedge_depth(hole.w) * um
        p = bs.BandParams(**ot.band_params_kwargs(hole, index=N, f_target=F_T,
                                                  xsect="wedge"))
        td.PolySlab.reset_mock()
        cell = bs.build_unit_cell(p)
        slabs = _polyslabs(td)
        prism = [s for s in slabs if s["axis"] == 0]
        assert len(prism) == 1 and abs(cell.depth_below - depth) < 1e-9
        v = np.array(prism[0]["vertices"])                    # (y, z)
        assert np.allclose(sorted(v[:, 0]), [-0.405, 0.0, 0.405])
        assert abs(v[:, 1].max() + hole.t / 2 * um) < 1e-9    # base = slab bottom
        assert abs(v[:, 1].min() + hole.t / 2 * um + depth) < 1e-9
        assert abs(depth - 0.81 / (2 * SQ3)) < 1e-9           # the .mph polygon
        cut = [s for s in slabs if s["axis"] == 2][0]
        assert cut["slab_bounds"][0] < -hole.t / 2 * um - depth   # through the wedge
        assert abs(cut["slab_bounds"][1] - hole.t / 2 * um) < 1e-9
        assert len(cell.structures) == 3                      # slab, wedge, hole

        td.Simulation.reset_mock(); td.PointDipole.reset_mock()
        bs.build_simulations(p, x_only=True, cell=cell)
        sim = td.Simulation.call_args.kwargs
        assert tuple(sim["symmetry"]) == (0, -1, 0)
        spacing = (sim["size"][1] - cell.width_y) / 2
        assert abs(sim["size"][2] - (p.thickness + 2 * depth + 2 * spacing)) < 1e-9
        z = [c.kwargs["center"][2] for c in td.PointDipole.call_args_list]
        assert max(abs(x) for x in z) > 0 and max(abs(x) for x in z) <= 0.3 * p.thickness

        # rect: unchanged (mid-plane sources, z symmetry, no prism)
        p0 = bs.BandParams(**ot.band_params_kwargs(hole, index=N, f_target=F_T))
        td.PolySlab.reset_mock(); td.Simulation.reset_mock(); td.PointDipole.reset_mock()
        cell0 = bs.build_unit_cell(p0)
        bs.build_simulations(p0, x_only=True, cell=cell0)
        assert all(s["axis"] == 2 for s in _polyslabs(td)) and cell0.depth_below == 0
        assert tuple(td.Simulation.call_args.kwargs["symmetry"]) == (0, -1, 1)
        assert all(c.kwargs["center"][2] == 0 for c in td.PointDipole.call_args_list)
        assert abs(td.Simulation.call_args.kwargs["size"][2]
                   - (p0.thickness + 2 * spacing)) < 1e-9

        # ---- rib: the roof of the filleted outline, as stacked offset outlines ----
        rib = RibGeometry(**RIB_CELL)
        pr = bs.BandParams(**ot.band_params_kwargs(rib, index=N, f_target=F_T,
                                                   xsect="wedge"))
        assert pr.wedge_layers == 12
        td.PolySlab.reset_mock()
        cell = bs.build_unit_cell(pr)
        slabs = _polyslabs(td)
        assert all(s["axis"] == 2 and s["sidewall_angle"] == 0 for s in slabs)
        z0 = -rib.th / 2 * um
        layers = [s for s in slabs if s["slab_bounds"][1] <= z0 + 1e-12]
        assert len(slabs) - len(layers) == 3                 # the slab, 3 periods
        d_max = geometry.wedge_depths(rib)["max"] * um
        assert abs(cell.depth_below - d_max) < 1e-9          # deepest: the junction
        tops = sorted({round(s["slab_bounds"][1], 9) for s in layers}, reverse=True)
        bots = sorted({round(s["slab_bounds"][0], 9) for s in layers}, reverse=True)
        assert len(tops) == 12 and abs(tops[0] - z0) < 1e-9  # starts at the slab
        assert np.allclose(tops[1:], bots[:-1])              # no gap between layers
        assert abs(bots[-1] - (z0 - d_max)) < 1e-9
        assert np.allclose(np.diff(tops), -d_max / 12)
        # each layer is the outline offset by the roof distance of its middle
        mid = [s for s in layers if abs(s["slab_bounds"][1] - z0) < 1e-12]
        assert len(mid) == 3                                 # one per period
        v = np.array([s for s in mid if abs(np.mean(np.array(s["vertices"])[:, 0])) < 1e-6
                      ][0]["vertices"])
        d = bs.rib_roof_distance(v[:, 0], v[:, 1], pr.a, pr.t, pr.w, pr.spine_width,
                                 pr.fillet_r1, pr.fillet_r2)
        assert np.allclose(d, d_max * SQ3 / 24, atol=2e-4)   # level delta = d_max_d/(2*12)

        # the adapter refuses a module whose symmetry does not follow the rule
        with tempfile.TemporaryDirectory() as tmp:
            rec = evaluate_candidate(_rib_u0(), celltype="rib",
                                     optical_backend="tidy3d",
                                     mech_backend="surrogate_stub", xsect="wedge",
                                     tidy3d_options=dict(mode="dry_run",
                                                         results_dir=tmp))
            d = rec["optical_tidy3d_detail"]
            assert rec["optical_tidy3d_stop"] == "dry_run", rec.get("error")
            assert d["symmetry"] == [0, -1, 0] and d["xsect"] == "wedge"
            assert d["z_symmetry"] is False and d["band_params"]["xsect"] == "wedge"
            assert abs(d["wedge_depths_m"]["rib"] - geometry.wedge_depth(
                rec["params"]["t"])) < 1e-15
            assert d["wedge_depths_m"]["max"] > d["wedge_depths_m"]["rib"]
            assert d["band_params"].get("wedge_layers") is None    # module default
            real = bs.BandParams.symmetry
            bs.BandParams.symmetry = property(lambda self: (0, -1, 1))
            try:
                bad = evaluate_candidate(_rib_u0(), celltype="rib",
                                         optical_backend="tidy3d",
                                         mech_backend="surrogate_stub",
                                         xsect="wedge",
                                         tidy3d_options=dict(mode="dry_run",
                                                             results_dir=tmp))
            finally:
                bs.BandParams.symmetry = real
            assert bad["status"] == "optical_failed" and "symmetry" in bad["error"]


# ---- rib wedge: the roof of the filleted outline -------------------------------

def _seg_dist(P, A, B):
    """Distance from each point of P to the nearest of the segments A[i]-B[i]."""
    AB = B - A
    s = np.clip(((P[:, None, :] - A[None]) * AB[None]).sum(-1)
                / np.maximum((AB * AB).sum(-1), 1e-30)[None], 0, 1)
    C = A[None] + s[..., None] * AB[None]
    return np.sqrt(((P[:, None, :] - C) ** 2).sum(-1)).min(1)


# rib cells in um: the default one, a wide spine, a short stubby rib (the tip
# plane cuts the junction cone), sharp corners, unequal fillets
ROOF_CELLS = [dict(a=.57, t=.391, w=1.747, s=.254, r1=.04, r2=.04),
              dict(a=.9, t=.2, w=1.0, s=.6, r1=.03, r2=.05),
              dict(a=.8, t=.5, w=.62, s=.2, r1=.05, r2=.05),
              dict(a=.6, t=.3, w=1.2, s=.25, r1=0., r2=0.),
              dict(a=.6, t=.3, w=1.2, s=.3, r1=.06, r2=.02)]


def _rib_m(c):
    return RibGeometry(th=250e-9, **{k: v * 1e-6 for k, v in c.items()})


def test_rib_roof_analytic_values():
    """Triangle on straight strips, cones under the fillets, hip at the tip."""
    g = RibGeometry(**RIB_CELL)
    d, depth = geometry.rib_roof_distance, geometry.rib_wedge_depth
    a, t, w, s, r1, r2 = g.a, g.t, g.w, g.s, g.r1, g.r2
    tol = 1e-15
    # spine, between two ribs: the triangle with base s, for any z there
    for z in (a / 2, a / 2 - 20e-9, -a / 2):
        for y in (0.0, 30e-9, -90e-9):
            assert abs(d(g, z, y) - (s / 2 - abs(y))) < tol
    assert abs(depth(g, a / 2, 0.0) - s / (2 * SQ3)) < tol        # apex b/(2*sqrt(3))
    # rib arm, away from junction and tip: the triangle with base t
    y_arm = (s / 2 + r1 + w / 2 - r2) / 2
    for z in (0.0, 50e-9, -150e-9):
        assert abs(d(g, z, y_arm) - (t / 2 - abs(z))) < tol
        assert abs(d(g, z, -y_arm) - (t / 2 - abs(z))) < tol
    assert abs(depth(g, 0.0, y_arm) - t / (2 * SQ3)) < tol
    # zero on the outline: edges, root arc, tip arc
    cz, cy = t / 2 + r1, s / 2 + r1                               # root arc centre
    for th in np.linspace(np.pi, 1.5 * np.pi, 7):
        assert abs(d(g, cz + r1 * np.cos(th), cy + r1 * np.sin(th))) < 1e-14
    tz, ty = t / 2 - r2, w / 2 - r2                               # tip arc centre
    for th in np.linspace(0, np.pi / 2, 7):
        assert abs(d(g, tz + r2 * np.cos(th), ty + r2 * np.sin(th))) < 1e-14
    assert abs(d(g, a / 2, s / 2)) < tol and abs(d(g, t / 2, y_arm)) < tol
    assert abs(d(g, 0.0, w / 2)) < tol
    # concave cone about the root arc centre: d = rho - r1, also over the junction
    for rho, th in ((r1 + 15e-9, 1.25 * np.pi), (r1 + 60e-9, 1.1 * np.pi),
                    (r1 + 60e-9, 1.4 * np.pi)):
        assert abs(d(g, cz + rho * np.cos(th), cy + rho * np.sin(th)) - (rho - r1)) < 1e-14
    assert abs(d(g, 0.0, 0.0) - (np.hypot(cz, cy) - r1)) < tol    # junction centre
    assert d(g, 0.0, 0.0) > max(t, s) / 2                         # deeper than both
    assert abs(depth(g, 0.0, 0.0) - geometry.wedge_depths(g)["max"]) < tol
    # convex cone inside the tip arc: d = r2 - rho, apex r2*tan(30) below the centre
    for rho, th in ((0.0, 0.0), (0.4 * r2, 0.3), (0.9 * r2, 1.2)):
        assert abs(d(g, tz + rho * np.cos(th), ty + rho * np.sin(th)) - (r2 - rho)) < 1e-14
    assert abs(depth(g, tz, ty) - r2 / SQ3) < tol
    # hip: near the tip the tip plane w/2 - |y| takes over from the rib triangle
    assert abs(d(g, 0.0, w / 2 - 50e-9) - 50e-9) < tol
    assert abs(d(g, 0.0, w / 2 - t / 2 - 1e-9) - t / 2) < tol
    # outside the cell: not positive; periodic along the beam
    assert d(g, a / 2, s / 2 + 10e-9) < 0 and d(g, 0.0, w / 2 + 1e-9) < 0
    assert d(g, cz - 0.5 * r1, cy - 0.5 * r1) < 0                 # air inside the root arc
    assert depth(g, a / 2, s / 2 + 10e-9) == 0.0
    assert abs(d(g, 0.123e-6 + 3 * a, 0.2e-6) - d(g, 0.123e-6, 0.2e-6)) < 1e-14


def test_rib_roof_is_the_distance_to_the_filleted_outline():
    """Closed form vs brute force: distance to the filleted polygon the slab
    is made of (3 periods), for five cells; and geometry.py == the module."""
    from matplotlib.path import Path
    rng = np.random.default_rng(1)
    with fake_tidy3d() as (td, bs):
        for c in ROOF_CELLS:
            outline = bs.rib_block_vertices(0.0, c["a"], c["t"], c["w"], c["s"],
                                            c["r1"], c["r2"], 200)
            P = np.vstack([outline, outline[:1]])
            A, B = P[:-1], P[1:]
            cut = ((np.abs(np.abs(A[:, 0]) - c["a"] / 2) < 1e-9)
                   & (np.abs(np.abs(B[:, 0]) - c["a"] / 2) < 1e-9))   # not outline
            A = np.vstack([A[~cut] + [dx, 0] for dx in (-c["a"], 0, c["a"])])
            B = np.vstack([B[~cut] + [dx, 0] for dx in (-c["a"], 0, c["a"])])
            pts = np.column_stack([rng.uniform(-c["a"] / 2, c["a"] / 2, 1500),
                                   rng.uniform(-c["w"] / 2, c["w"] / 2, 1500)])
            inside = Path(outline).contains_points(pts)
            d = bs.rib_roof_distance(pts[:, 0], pts[:, 1], **c)
            assert inside.sum() > 200
            assert np.abs(d[inside] - _seg_dist(pts[inside], A, B)).max() < 2e-6, c
            assert np.all(d[~inside] < 1e-6) and np.all(d[inside] > -1e-6)
            assert abs(d.max() - bs.rib_roof_max_distance(c["t"], c["w"], c["s"],
                                                          c["r1"])) < 0.02
            assert abs(bs.rib_roof_distance(0.0, 0.0, **c)
                       - bs.rib_roof_max_distance(c["t"], c["w"], c["s"], c["r1"])) < 1e-12
            # the same function in geometry.py (meters)
            g = _rib_m(c)
            ref = np.array([geometry.rib_roof_distance(g, x * 1e-6, y * 1e-6)
                            for x, y in pts[:300]]) * 1e6
            assert np.abs(ref - d[:300]).max() < 1e-9
            assert abs(geometry.wedge_depths(g)["max"] * 1e6
                       - bs.rib_roof_max_distance(c["t"], c["w"], c["s"], c["r1"])
                       / SQ3) < 1e-9


def test_tidy3d_layer_stack_reproduces_the_roof():
    """Offset outlines are level sets of the roof distance in every topology
    (spanning, island, short rib), and the stack is the roof to half a layer."""
    from matplotlib.path import Path
    rng = np.random.default_rng(2)
    with fake_tidy3d() as (td, bs):
        for c in ROOF_CELLS:
            d_max = bs.rib_roof_max_distance(c["t"], c["w"], c["s"], c["r1"])
            pts = np.column_stack([rng.uniform(-c["a"] / 2, c["a"] / 2, 4000),
                                   rng.uniform(-c["w"] / 2, c["w"] / 2, 4000)])
            d = bs.rib_roof_distance(pts[:, 0], pts[:, 1], **c)
            kinds = set()
            for frac in (0.0, 0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 0.999):
                delta = frac * d_max
                polys = bs.rib_roof_level_polygons(0.0, c["a"], c["t"], c["w"], c["s"],
                                                   c["r1"], c["r2"], delta, 48)
                assert len(polys) == 1, (c, frac)
                v = polys[0]
                kinds.add("span" if np.abs(v[:, 0]).max() > c["a"] / 2 - 1e-9 else "island")
                # vertices lie on the level set (cusps on the axis included) ...
                on_cut = np.abs(np.abs(v[:, 0]) - c["a"] / 2) < 1e-9
                dv = bs.rib_roof_distance(v[:, 0], v[:, 1], **c)
                # (1e-4 um: a tip cut through the junction arc lands on a chord)
                assert np.abs(dv[~on_cut] - delta).max() < 1e-4, (c, frac)
                # ... and the polygon is the set d >= delta (chord error only)
                wrong = Path(v).contains_points(pts) != (d >= delta)
                assert not wrong.any() or np.abs(d[wrong] - delta).max() < 3e-4, (c, frac)
            assert kinds == {"span", "island"}, (c, kinds)
            assert bs.rib_roof_level_polygons(0.0, c["a"], c["t"], c["w"], c["s"],
                                              c["r1"], c["r2"], 1.001 * d_max) == []
            # delta = 0 is the filleted outline of the slab itself
            v0 = bs.rib_roof_level_polygons(0.0, c["a"], c["t"], c["w"], c["s"],
                                            c["r1"], c["r2"], 0.0, 48)[0]
            ref = bs.rib_block_vertices(0.0, c["a"], c["t"], c["w"], c["s"],
                                        c["r1"], c["r2"], 49)
            assert abs(bs._polygon_area(v0) - bs._polygon_area(ref)) < 1e-9

            # the stack built by rib_block_cell: depth at a point = bottom of the
            # deepest layer whose polygon holds it
            p = bs.BandParams(geometry="rib_block", xsect="wedge", a=c["a"], t=c["t"],
                              w=c["w"], spine_width=c["s"], fillet_r1=c["r1"],
                              fillet_r2=c["r2"], thickness=0.25, fillet_points=24)
            td.PolySlab.reset_mock()
            cell = bs.build_unit_cell(p)
            z0 = -p.thickness / 2
            stack = np.zeros(len(pts))
            for call in td.PolySlab.call_args_list:
                lo, hi = call.kwargs["slab_bounds"]
                if hi > z0 + 1e-12:
                    continue                                  # the slab itself
                hit = Path(np.asarray(call.kwargs["vertices"])).contains_points(pts)
                stack[hit] = np.maximum(stack[hit], z0 - lo)
            roof = np.maximum(d, 0) / SQ3
            half_layer = cell.depth_below / (2 * p.wedge_layers)
            assert abs(cell.depth_below - d_max / SQ3) < 1e-12
            assert np.abs(stack - roof).max() <= half_layer + 5e-4, (
                c, np.abs(stack - roof).max(), half_layer)
            assert np.abs(stack - roof).mean() < 0.6 * half_layer
            # more layers: smaller error (the documented bound)
            p2 = bs.BandParams(**{**bs.asdict(p), "wedge_layers": 48})
            td.PolySlab.reset_mock()
            bs.build_unit_cell(p2)
            n_layers = {round(k.kwargs["slab_bounds"][1], 9)
                        for k in td.PolySlab.call_args_list} - {round(-z0, 9)}
            assert len(n_layers) == 48
        _raises(ValueError, bs.build_unit_cell,
                bs.BandParams(xsect="wedge", wedge_layers=0))
        # no block: the plain strip keeps the exact prism
        td.PolySlab.reset_mock()
        cell = bs.build_unit_cell(bs.BandParams(xsect="wedge", t=0.0))
        assert sum(k.kwargs["axis"] == 0 for k in td.PolySlab.call_args_list) == 1
        assert abs(cell.depth_below - 0.2 / (2 * SQ3)) < 1e-12


def _builder_wedge_depth(c, z, y, tan_a=1 / SQ3):
    """Depth of the wedge that comsol/build_omc_unitcell_rib.m ('wedge') builds,
    evaluated piece by piece exactly as its geometry sequence combines them:
    (S | A | J) & tip block, minus the corner boxes, plus the quarter cones.
    y >= 0 half; lengths in um. A solid given by `depth <= h(z, y)` inside a
    footprint has depth h there and -1 (absent) elsewhere."""
    a, t, w, s, r1, r2 = (c[k] for k in ("a", "t", "w", "s", "r1", "r2"))
    az = np.abs(z)
    none = -np.ones_like(z)
    # S: spine prism, extruded over z in [-a/2, -(t/2+r1)] and [t/2+r1, a/2]
    S = np.where((az >= t / 2 + r1) & (az <= a / 2) & (y <= s / 2), tan_a * (s / 2 - y), none)
    # A: rib prism from y = s/2 + r1 over the distance w/2 - s/2 - r1
    A = np.where((az <= t / 2) & (y >= s / 2 + r1) & (y <= w / 2), tan_a * (t / 2 - az), none)
    # J: box minus cone frustum (radius r1 at the slab, + depth/tan below)
    wdj = tan_a * (np.hypot(t / 2 + r1, s / 2 + r1) - r1)
    rho = np.hypot(az - (t / 2 + r1), y - (s / 2 + r1))
    inbox = (az <= t / 2 + r1) & (y <= s / 2 + r1)
    J = np.where(inbox & (rho >= r1), np.minimum(tan_a * (rho - r1), 1.25 * wdj), none)
    U = np.maximum(np.maximum(S, A), J)                       # union
    tip = tan_a * (w / 2 - y)                                 # tip block
    H = np.where(U >= 0, np.minimum(U, tip), none)            # intersection
    corner = (az >= t / 2 - r2) & (az <= t / 2) & (y >= w / 2 - r2) & (y <= w / 2)
    H = np.where(corner, none, H)                             # minus corner boxes
    rc = np.hypot(az - (t / 2 - r2), y - (w / 2 - r2))        # cone, apex down
    C = np.where(corner & (rc <= r2), tan_a * (r2 - rc), none)
    return np.maximum(np.maximum(H, C), 0.0), wdj


def test_comsol_builder_construction_is_the_same_roof():
    """The exact plane / cone construction of the rib wedge template equals
    depth = tan(30 deg) * (distance to the filleted outline) everywhere."""
    rng = np.random.default_rng(3)
    with fake_tidy3d() as (td, bs):
        for c in ROOF_CELLS:
            if min(c["r1"], c["r2"]) <= 0:
                continue                 # the template needs r1, r2 > 0 (Fillet)
            z = rng.uniform(-c["a"] / 2, c["a"] / 2, 20000)
            y = rng.uniform(0, c["w"] / 2 * 1.02, 20000)
            built, wdj = _builder_wedge_depth(c, z, y)
            roof = np.maximum(bs.rib_roof_distance(z, y, **c), 0) / SQ3
            assert np.abs(built - roof).max() < 1e-12, (c, np.abs(built - roof).max())
            # the junction box (1.25*wdj deep) never truncates the roof, and the
            # parameter wd of the template is the deepest point
            assert roof.max() < 1.25 * wdj
            wd = min(wdj, c["w"] / 2 / SQ3)
            assert abs(wd - bs.rib_roof_max_distance(c["t"], c["w"], c["s"], c["r1"])
                       / SQ3) < 1e-12
    # the builder source holds these very pieces and parameters
    with open(os.path.join(_ROOT, "comsol", "build_omc_unitcell_rib.m")) as fh:
        src = fh.read()
    for needle in ("'wdj', 'tan(wedge_angle)*(sqrt((t/2+r1)^2+(s/2+r1)^2)-r1)'",
                   "'wd', 'min(wdj,tan(wedge_angle)*w/2)'",
                   "ws.setIndex('distance', 'a/2-t/2-r1', 0);",
                   "wpr.set('quicky', 's/2+r1');",
                   "wr.setIndex('distance', 'w/2-s/2-r1', 0);",
                   "blk.set('size', {hJ, 's/2+r1', 't/2+r1'});",
                   "con.set('r', ['r1+' hJ '/' tanA]);", "con.set('rtop', 'r1');",
                   "[xBot '-w/2*' tanA], '0'});",
                   "con.set('h', ['r2*' tanA]);", "con.set('rtop', '0');",
                   "{'ribslab', 'wedge_hc', 'wedge_c1', 'wedge_c2'}",
                   "xBeamLo = '-th/2-wd'"):
        assert needle in src, needle
    # order of the sequence: filleted outline and slab first, wedge after
    order = [src.index(k) for k in ("'fil2_u', 'Fillet'", "'fil1_u', 'Fillet'",
                                    "'ribslab', 'Extrude'", "'wedge_s%d'",
                                    "'wedge_u', 'Union'", "'wedge_h', 'Intersection'",
                                    "'ribcell', 'Union'")]
    assert order == sorted(order)


# ---- COMSOL: thickness parity read from the field ----------------------------

def _slab_fields(parity, th=250e-9, n=4001, kx=None):
    """Field of a slab mode sampled across the thickness (x), slab at |x| < th/2.

    even: in-plane E_y, E_z ~ cos(kx x), E_x ~ sin(kx x)   (TE-like)
    odd : in-plane E_y, E_z ~ sin(kx x), E_x ~ cos(kx x)   (TM-like)
    """
    kx = kx or 0.8 * np.pi / th
    x = np.linspace(-1.5 * th, 1.5 * th, n)
    env = np.exp(-np.clip(np.abs(x) - th / 2, 0, None) / (0.4 * th))
    c, s = np.cos(kx * np.clip(x, -th / 2, th / 2)), np.sin(kx * np.clip(x, -th / 2, th / 2))
    if parity == "even":
        return x, 0.6 * s * env, 1.0 * c * env, 0.3j * c * env, th
    return x, 1.0 * c * env, 0.5 * s * env, 0.2j * s * env, th


def test_parity_fraction_on_synthetic_fields():
    F = comsol_optics.z_even_fraction_from_fields
    x, Ex, Ey, Ez, th = _slab_fields("even")
    f_even = F(x, Ex, Ey, Ez, th)
    x, Ex, Ey, Ez, th = _slab_fields("odd")
    f_odd = F(x, Ex, Ey, Ez, th)
    assert f_even > 0.95 and f_odd < 0.05, (f_even, f_odd)
    # robust to an offset mesh and to unequal weights
    x, Ex, Ey, Ez, th = _slab_fields("even", n=37)
    assert F(x, Ex, Ey, Ez, th, weights=np.gradient(x)) > 0.95
    # a 50/50 mixture (broken mirror) sits near the threshold, not at an end
    xe, Exe, Eye, Eze, th = _slab_fields("even")
    xo, Exo, Eyo, Ezo, _ = _slab_fields("odd")
    mixed = F(xe, Exe + Exo, Eye + Eyo, Eze + Ezo, th)
    assert 0.3 < mixed < 0.8
    assert np.isnan(F(np.array([1.0]), [1], [1], [1], 250e-9))    # no sample in layer
    # family selection, threshold 0.5; NaN belongs to neither family
    frac = np.array([0.97, 0.03, 0.5, np.nan])
    sel = comsol_optics.select_z_family
    assert sel(frac, "even").tolist() == [True, False, True, False]
    assert sel(frac, "odd").tolist() == [False, True, False, False]
    _raises(ValueError, sel, frac, "both")
    assert comsol_optics.Z_PARITY_THRESHOLD == 0.5
    # which family a study scores
    fam = comsol_optics.study_z_family
    assert fam("opt TE", "rect") == "even" and fam("opt TM", "rect") == "odd"
    assert fam("opt TE", "wedge") is None and fam("opt TM", "wedge") is None
    _raises(ValueError, fam, "opt XX", "rect")
    # the COMSOL expressions are the same ratio, on the thickness parameter
    num, den = comsol_optics.z_parity_exprs("th")
    assert num == "intop1((abs(x)<0.125*th)*(abs(ewfd.Ey)^2+abs(ewfd.Ez)^2))"
    assert "abs(ewfd.Ex)^2" in den and "(abs(x)<0.125*th)" in den
    assert comsol_optics.z_parity_exprs("t")[0].count("0.125*t)") == 1


class FakeModel:
    """Minimal MPh model: one optical study, modes with a known parity."""

    def __init__(self, params=("a", "w", "t", "s", "th", "r1", "r2", "kF"),
                 freqs=(2e12, 150e12, 160e12, 200e12, 230e12),
                 frac=(0.5, 0.98, 0.03, 0.96, 0.01)):
        self._params = {p: "1" for p in params}
        self.freqs, self.frac = np.array(freqs), np.array(frac)
        self.solved, self.exprs = [], []
        self.java = mock.MagicMock()
        mats = self.java.component.return_value.material
        mats.return_value.tags.return_value = ["mat1", "mat2"]
        mats.return_value.label.return_value = "Diamond (100)"

    def parameters(self):
        return dict(self._params)

    def parameter(self, name, value):
        self._params[name] = value

    def datasets(self):
        return ["opt TE//Solution 3"]

    def solve(self, study):
        self.solved.append(study)

    def evaluate(self, expr, dataset=None):
        self.exprs.append(expr)
        if expr == "freq":
            return self.freqs + 1e-3j
        if "abs(ewfd.Ex)^2" in expr:                          # denominator
            return np.ones_like(self.freqs)
        return self.frac                                      # numerator

    def eps_set(self):
        pg = self.java.component.return_value.material.return_value.propertyGroup
        return [(c.args[0], [str(x) for x in c.args[1]])
                for c in pg.return_value.set.call_args_list]


def test_comsol_sweep_keeps_one_parity_family_and_sets_the_index():
    g = RibGeometry(**RIB_CELL)
    k = [0.5, 1.0]
    # all modes (wedge / legacy): only the pseudo-mode is dropped
    m = FakeModel()
    _, f_all = comsol_optics.opt_sweep_at(m, g, "opt TE", k, 6)
    assert np.allclose(f_all[0, :4], [150e12, 160e12, 200e12, 230e12])
    assert not any("intop1" in e for e in m.exprs)            # no parity query
    # rect, opt TE: thickness-even modes only, parity recorded per mode
    m, info = FakeModel(), {}
    _, f_even = comsol_optics.opt_sweep_at(m, g, "opt TE", k, 6, z_family="even",
                                           thickness="th", info=info)
    assert np.allclose(f_even[0, :2], [150e12, 200e12]) and np.isnan(f_even[0, 2])
    assert len(info["z_even_fraction"]) == 2 == len(info["freqs_all_hz"])
    assert info["freqs_all_hz"][0] == [2e12, 150e12, 160e12, 200e12, 230e12]
    assert info["z_even_fraction"][0] == [0.5, 0.98, 0.03, 0.96, 0.01]
    assert any("(abs(x)<0.125*th)" in e for e in m.exprs)
    # opt TM on a rectangular cell: the odd family
    m = FakeModel()
    _, f_odd = comsol_optics.opt_sweep_at(m, g, "opt TM", k, 6, z_family="odd",
                                          thickness="th")
    assert np.allclose(f_odd[0, :2], [160e12, 230e12])
    # the hole sweeps take the same arguments
    hole = Geometry(a=891e-9, w=810e-9, t=253e-9, hx=218e-9, hy=334e-9)
    m = FakeModel()
    _, f = comsol_optics._opt_sweep(m, hole, "opt TE", 2, 6, z_family="even",
                                    thickness="t")
    assert np.allclose(f[-1, :2], [150e12, 200e12])
    assert any("(abs(x)<0.125*t)" in e for e in m.exprs)
    m = FakeModel()
    _, f = comsol_optics._full_opt_sweep(m, hole, "opt TE", n_k=2, n_bands=6)
    assert np.allclose(f[-1, :4], [150e12, 160e12, 200e12, 230e12])
    # a filter without the thickness parameter, or a size mismatch, is an error
    _raises(ValueError, comsol_optics.opt_sweep_at, FakeModel(), g, "opt TE", k, 6,
            z_family="even")
    bad = FakeModel(frac=(0.9, 0.1))
    _raises(RuntimeError, comsol_optics.opt_sweep_at, bad, g, "opt TE", k, 6,
            z_family="even", thickness="th")
    # twice as many eigenvalues are requested with the filter
    assert comsol_optics._prepare_optics(FakeModel(), 10, 6, None) == 10
    assert comsol_optics._prepare_optics(FakeModel(), 10, 6, "even") == 20
    assert comsol_optics._prepare_optics(FakeModel(), 4, 6, "even") == 12

    # every sweep wrote n = 2.386 into the template's diamond material
    for model in (m, FakeModel()):
        comsol_optics.opt_sweep_at(model, g, "opt TE", k, 6)
        name, tensor = model.eps_set()[-1]
        assert name == "relpermittivity" and len(tensor) == 9
        assert abs(float(tensor[0]) - N ** 2) < 1e-12 and float(tensor[1]) == 0
        assert tensor[0] == tensor[4] == tensor[8]
    # a template with the parameter n_diamond gets the parameter instead
    m = FakeModel(params=("a", "kF", "n_diamond"))
    assert comsol_client.set_diamond_index(m) == N
    assert float(m.parameters()["n_diamond"]) == N and m.eps_set() == []
    assert comsol_client.set_diamond_index(m, 2.0) == 2.0
    # and a failure is loud: never solve with the template's own index
    m = FakeModel()
    m.java.component.side_effect = RuntimeError("no such component")
    err = _raises(RuntimeError, comsol_client.set_diamond_index, m)
    assert "refractive index" in str(err)
    _raises(RuntimeError, comsol_optics.opt_sweep_at, m, g, "opt TE", k, 6)


def test_rib_objective_drives_comsol_with_xsect():
    """objective_rib._solve_optical_bands: template, parity family, index."""
    import acoustic_comsol
    loaded = []
    m = FakeModel()
    real = (comsol_client.get_model, comsol_client.template_path,
            comsol_client.rib_air_is_disk)

    def fake_template(celltype="hole", xsect="rect", must_exist=False):
        return os.path.join("comsol", comsol_client.TEMPLATES[celltype, xsect])

    comsol_client.get_model = lambda t: loaded.append(os.path.basename(t)) or m
    comsol_client.template_path = fake_template
    comsol_client.rib_air_is_disk = lambda p: True     # no file is read
    try:
        g = RibGeometry(**RIB_CELL)
        info = {}
        f = objective_rib._solve_optical_bands(g, "opt TE", [0.5, 1.0], 6, F_T,
                                               xsect="rect", info=info)
        assert np.allclose(f[0, :2], [150e12, 200e12])          # even only
        f = objective_rib._solve_optical_bands(g, "opt TE", [0.5, 1.0], 6, F_T,
                                               xsect="wedge")
        assert np.allclose(f[0, :4], [150e12, 160e12, 200e12, 230e12])
    finally:
        (comsol_client.get_model, comsol_client.template_path,
         comsol_client.rib_air_is_disk) = real
    assert loaded == ["omc_unitcell_rib.mph", "omc_unitcell_rib_wedge.mph"]
    assert abs(float(m.eps_set()[-1][1][0]) - N ** 2) < 1e-12
    assert info["z_even_fraction"][0][1] == 0.98
    # mechanics has no parity filter: run_mechanical_comsol has no such option
    import inspect
    assert "z_family" not in inspect.signature(
        acoustic_comsol.run_mechanical_comsol).parameters


# ---- diamond index 2.386 everywhere ----------------------------------------------

def test_index_2386_propagates_to_every_engine():
    import yaml
    with open(os.path.join(_ROOT, "configs", "materials.yaml")) as fh:
        assert yaml.safe_load(fh)["diamond"]["refractive_index_1550nm"] == N
    assert geometry.diamond_index() == N == ot.material_index()
    # surrogate
    import optical_surrogate
    g = u_to_geometry(HOLE_U)
    assert optical_surrogate._effective_indices(g) == \
        optical_surrogate._effective_indices(g, N)
    assert optical_surrogate._effective_indices(g) != \
        optical_surrogate._effective_indices(g, 2.40)
    assert optical_surrogate.optical_gap_surrogate(g)["n_hi"] == \
        optical_surrogate._effective_indices(g, N)[0]
    # MPB (importable without meep)
    import optical_mpb
    assert optical_mpb.N_DIAMOND == N and optical_mpb.EPS_DIAMOND == N ** 2
    # Tidy3D adapter and module
    assert ot.Tidy3DOptions().index_value() == N
    assert json.loads(ot.id_settings(None, "hole"))[3] == N
    with fake_tidy3d() as (td, bs):
        assert bs.MATERIAL_INDEX["diamond_telecom"] == N
        p = bs.BandParams(**ot.band_params_kwargs(g, index=ot.material_index(),
                                                  f_target=F_T))
        assert p.n == N
    # records
    rec = evaluate_candidate(HOLE_U, optical_backend="surrogate",
                             mech_backend="surrogate_stub")
    assert rec["optical_index"] == N
    # comsol batch (standard library reader of materials.yaml)
    import comsol_batch
    assert comsol_batch.repo_diamond_index() == N
    over = {}
    note = comsol_batch.diamond_index_override("t.mph", {"n_diamond": ("2.386", N)},
                                               ["std2"], over, N, False)
    assert over == {"n_diamond": "2.386"} and "pname" in note
    assert "not used" in comsol_batch.diamond_index_override(
        "t.mph", {}, ["std1", "std3"], {}, N, False)
    err = _raises(comsol_batch.BatchError, comsol_batch.diamond_index_override,
                  "no_such.mph", {"a": ("1", 1.0)}, ["std2"], {}, N, False)
    assert err.code == 2 and "n_diamond" in str(err)
    # a hole template: unpatched (a number in the material) vs patched in
    # place by comsol/patch_hole_templates_index.m (parameter n_diamond)
    import zipfile

    def blob(name, value, real, descr):
        return ('{"description":"%s","name":"%s","scalarImag":"0.0",'
                '"scalarReal":"%s","value":"%s"}' % (descr, name, real, value))

    with tempfile.TemporaryDirectory() as tmp:
        old, new = os.path.join(tmp, "old.mph"), os.path.join(tmp, "new.mph")
        base = [blob("a", "891[nm]", "8.91E-7", "lattice"), blob("kF", "pi/a", "3.5E6", "k"),
                blob("relpermittivity", "1", "1.0", "Relative permittivity")]
        with zipfile.ZipFile(old, "w") as z:
            z.writestr("smodel.json", "[" + ",".join(base + [blob(
                "relpermittivity", "5.76", "5.76", "Relative permittivity")]) + "]")
        with zipfile.ZipFile(new, "w") as z:
            z.writestr("smodel.json", "[" + ",".join(base + [
                blob("n_diamond", "2.386", "2.386", "diamond refractive index"),
                blob("relpermittivity", "n_diamond^2", "5.692996",
                     "Relative permittivity")]) + "]")
        tp_old = comsol_batch.template_parameters(old)
        tp_new = comsol_batch.template_parameters(new)
        assert "n_diamond" not in tp_old and tp_new["n_diamond"] == ("2.386", N)
        assert comsol_batch.template_diamond_eps(old) == 5.76
        assert abs(comsol_batch.template_diamond_eps(new) - N ** 2) < 1e-9
        err = _raises(comsol_batch.BatchError, comsol_batch.diamond_index_override,
                      old, tp_old, ["std2"], {}, N, False)
        assert "patch_hole_templates_index" in str(err) and err.code == 2
        over = {}
        note = comsol_batch.diamond_index_override(new, tp_new, ["std2"], over, N, False)
        assert over == {"n_diamond": "2.386"} and "pname" in note   # not refused
        over = {"n_diamond": "2.40"}                                # user's -p wins
        comsol_batch.diamond_index_override(new, tp_new, ["std2"], over, N, False)
        assert over == {"n_diamond": "2.40"}
    # the MPh drivers: parameter when the template has it (patched), material
    # property otherwise (unpatched) -- checked in
    # test_comsol_sweep_keeps_one_parity_family_and_sets_the_index
    with open(os.path.join(_ROOT, "comsol", "patch_hole_templates_index.m")) as fh:
        psrc = fh.read()
    for needle in ("PARAM = 'n_diamond';",
                   "DESCR = 'diamond refractive index, telecom (configs/materials.yaml)';",
                   "nDiamond = 2.386;", ".set('relpermittivity', {[PARAM '^2']});",
                   "copyfile(file, r.backup)", "mphsave(model, file);",
                   "'omc_unitcell.mph'", "'omc_unitcell_iso.mph'"):
        assert needle in psrc, needle
    assert "model.study(" not in psrc and ".run" not in psrc      # nothing is solved
    assert psrc.index("copyfile(file, r.backup)") < psrc.index("mphsave(model, file);")
    # the rib template builder
    with open(os.path.join(_ROOT, "comsol", "build_omc_unitcell_rib.m")) as fh:
        src = fh.read()
    assert "'n_diamond', '2.386'" in src and "{'n_diamond^2'}" in src
    assert "'diamond refractive index, telecom (configs/materials.yaml)'" in src
    assert "{'5.76'}" not in src and "omc_unitcell_rib_wedge.mph" in src
    assert "tan(wedge_angle)" in src and "'30[deg]'" in src


def test_matlab_diamond_index_is_2386():
    """omc-comsol-chang: no diamond index other than 2.386 left (static)."""
    root = os.path.join(_ROOT, "..", "omc-comsol-chang")
    if not os.path.isdir(root):
        print("  (skipped: omc-comsol-chang not present)")
        return
    def read(*parts):
        with open(os.path.join(root, *parts), errors="replace") as fh:
            return fh.read()
    stale = re.compile(r"2\.40[0-9]{1,2}\b")
    files = ["LoadMaterialParams.m", "runOpticalBand_1D.m", "runOpticalBand_2D.m",
             "runOpticalBand_3D.m"]
    nb = os.path.join(root, "nanobeam")
    files += [os.path.join("nanobeam", f) for f in os.listdir(nb) if f.endswith(".m")]
    for f in files:
        hits = [ln for ln in read(f).splitlines()
                if stale.search(ln.split("%")[0])]            # code, not comments
        assert not hits, (f, hits[:2])
    both = [re.search(r"P\.nbeam = ([0-9.]+);", read(*p)).group(1)
            for p in (("LoadMaterialParams.m",), ("nanobeam", "LoadMaterialParams.m"))]
    assert both == ["2.386", "2.386"]                       # the two copies agree
    assert read("runOpticalBand_1D.m").count(
        "{'2.386' '0' '0' '0' '2.386' '0' '0' '0' '2.386'}") == 2
    # lithium niobate and the other materials are untouched
    assert "noref_sel(lbd0)" in read("runOpticalBand_1D.m")
    assert "{'3.5' '0' '0' '0' '3.5'" in read("runOpticalBand_3D.m")


# ---- queue runner passes the options through ------------------------------------

def test_batch_queue_passes_xsect_and_tidy3d():
    import run_batch_rib as rb
    u = [round(x, 6) for x in _rib_u0()]
    cfg = dict(celltype="rib", optical="tidy3d", mech="surrogate_stub",
               xsect="wedge", tidy3d=dict(mode="dry_run", n_k=5),
               queue=[u, dict(u=u, xsect="rect", tidy3d=dict(n_k=7)),
                      dict(u=u, optical="none")])
    items = rb.build_queue(cfg, _ROOT)
    assert [it["xsect"] for it in items] == ["wedge", "rect", "wedge"]
    assert items[0]["tidy3d"] == dict(mode="dry_run", n_k=5)
    assert items[1]["tidy3d"] == dict(mode="dry_run", n_k=7)      # merged
    assert items[2]["tidy3d"] == {}                               # not tidy3d
    for bad in (dict(cfg, xsect="round"),
                dict(cfg, tidy3d=dict(alow_cloud=True)),
                dict(cfg, celltype="hole", optical="mpb",
                     queue=[[0.5] * 5])):
        _raises(ValueError, rb.build_queue, bad, _ROOT)
    # a run_one-style config brings its own xsect / tidy3d
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "one.yaml")
        with open(path, "w") as fh:
            fh.write("celltype: rib\nu: %s\noptical: tidy3d\nmech: surrogate_stub\n"
                     "xsect: wedge\ntidy3d: {mode: dry_run}\n" % u)
        items = rb.build_queue(dict(optical="none", mech="surrogate_stub",
                                    queue=[dict(config=path)]), tmp)
        assert items[0]["xsect"] == "wedge" and items[0]["tidy3d"] == {"mode": "dry_run"}
    # main(): evaluate_candidate receives both, as from run_one.py
    got = []

    def fake_eval(u, **kw):
        got.append((kw.get("xsect"), kw.get("tidy3d_options")))
        return dict(id="x%d" % len(got), status="success", score=0.0)

    real = (rb.evaluate_candidate, rb.save_result)
    rb.evaluate_candidate, rb.save_result = fake_eval, lambda rec, path: None
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "q.yaml")
            import yaml
            with open(path, "w") as fh:
                yaml.safe_dump(dict(cfg, quiet=True), fh)
            with mock.patch("sys.stdout"), mock.patch("sys.stderr"):
                assert rb.main(["--config", path]) == 0
    finally:
        rb.evaluate_candidate, rb.save_result = real
    assert got == [("wedge", dict(mode="dry_run", n_k=5)),
                   ("rect", dict(mode="dry_run", n_k=7)), ("wedge", None)]


# ---- rib air half-disk (replaces the air box) ---------------------------------

def test_rib_air_disk_template_check():
    """rib_air_is_disk / require_rib_air_disk read the .mph zip, no COMSOL."""
    import zipfile
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "omc_unitcell_rib.mph")
        with open(path, "w") as fh:
            fh.write("not a zip")
        assert comsol_client.rib_air_is_disk(path) is True       # unreadable
        assert comsol_client.rib_air_is_disk(os.path.join(tmp, "none.mph")) is True
        os.utime(path, (1, 1))
        with zipfile.ZipFile(path, "w") as z:                     # roof, old box
            z.writestr("dmodel.xml", 's("wdj") s("th") s("airblk")')
        assert comsol_client.rib_wedge_is_roof(path) is True
        assert comsol_client.rib_air_is_disk(path) is False
        err = _raises(RuntimeError, comsol_client.require_rib_air_disk, path)
        assert "air_r" in str(err) and "rebuild" in str(err)
        assert "Mechanical solves are unaffected" in str(err)
        os.utime(path, (2, 2))                                    # a new file
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("dmodel.xml", 's("wdj") s("air_gap") s("air_r")')
        assert comsol_client.rib_air_is_disk(path) is True
        assert comsol_client.require_rib_air_disk(path) == path


def test_builder_sets_the_comsol_kernel():
    # Without geomRep('comsol') the geometry takes the user's COMSOL preference;
    # with the CAD kernel every build and solve needs a CADIMPORT seat.
    with open(os.path.join(_ROOT, "comsol", "build_omc_unitcell_rib.m")) as fh:
        src = fh.read()
    code = "\n".join(line.split("%", 1)[0] for line in src.splitlines())
    assert "geom.geomRep('comsol');" in code
    assert "cadps" not in code
    assert code.index("comp.geom.create('geom1', 3)") < code.index("geom.geomRep('comsol')") \
        < code.index("geom.create('wp1', 'WorkPlane')")


def test_builder_has_the_air_half_disk():
    with open(os.path.join(_ROOT, "comsol", "build_omc_unitcell_rib.m")) as fh:
        src = fh.read()
    for needle in ("'air_r', 'sqrt((w/2)^2+(th/2-air_xc)^2)+air_gap'",
                   "model.param.set('air_gap', '1.55[um]'",
                   "airXc = '-wd/2';",
                   "airCyl.set('pos', {'air_xc', '0', '-a/2'});",
                   "airCyl.set('axistype', 'z');",
                   "air.selection('input').set({'airdsk'});",
                   "scat.set('subtract', {'zloA', 'zhiA', 'selYsymAll', 'beamA'});",
                   "xBeamLo = '-th/2-wd'", "xBeamLo = '-th/2'"):
        assert needle in src, needle
    assert "'airblk'" not in src and "2.5*th" not in src.split("%%", 1)[1]
    # every air size derives from air_gap / air_r: no number in the air pieces
    for line in src.splitlines():
        if re.match(r"\s*(airCyl|airCut)\.set\('(pos|size|r|h)'", line):
            assert "e-" not in line and "[nm]" not in line and "[um]" not in line, line


def test_batch_air_region_note_and_refusal():
    import comsol_batch
    disk = {"th": ("250[nm]", 250e-9), "air_gap": ("1.55[um]", 1.55e-6),
            "air_r": ("...", 2434.9e-9)}
    box = {"th": ("250[nm]", 250e-9)}
    note = comsol_batch.air_region_note
    assert note(box, ["std1", "std3"]) == "no optical study: air region not used"
    assert note(disk, ["std2"]).startswith("half-disk, air_r = 2434.9 nm")
    assert note({"hx": ("1", 1.0)}, ["std2"]) == \
        "as stored in the template (not a rib template)"
    assert note({}, ["std4"]) == "as stored in the template (not a rib template)"
    err = _raises(comsol_batch.BatchError, note, box, ["std2"])
    assert err.code == 2 and "--allow-old-air" in str(err) and "rebuild" in str(err)
    with mock.patch("sys.stderr"), mock.patch("sys.stdout"):
        assert note(box, ["std1", "std4"], allow_old=True).startswith("OLD AIR BOX")

    # main(): refused before any COMSOL process; --allow-old-air goes on
    import zipfile

    def blob(name, value, real):
        return ('{"description":"d","name":"%s","scalarImag":"0.0",'
                '"scalarReal":"%s","value":"%s"}' % (name, real, value))

    class Past(Exception):
        pass

    def stop(*a, **kw):
        raise Past()

    with tempfile.TemporaryDirectory() as tmp:
        mph = os.path.join(tmp, "old_box.mph")
        with zipfile.ZipFile(mph, "w") as z:
            z.writestr("smodel.json", "[" + ",".join([
                blob("a", "570[nm]", "5.7E-7"), blob("th", "250[nm]", "2.5E-7"),
                blob("kF", "pi/a", "5.5E6"),
                blob("n_diamond", "2.386", "2.386")]) + "]")
        args = ["--template", mph, "--study", "std2", "--outdir",
                os.path.join(tmp, "out")]
        with mock.patch.object(comsol_batch, "find_comsol_root", lambda a: tmp), \
                mock.patch.object(comsol_batch, "find_batch_command",
                                  lambda r: ["comsolbatch"]), \
                mock.patch.object(comsol_batch, "find_javac", stop), \
                mock.patch.object(comsol_batch, "run_logged", stop), \
                mock.patch("sys.stderr"), mock.patch("sys.stdout"):
            err = _raises(comsol_batch.BatchError, comsol_batch.main, args)
            assert err.code == 2 and "air_r" in str(err)
            assert not os.path.exists(os.path.join(tmp, "out"))   # nothing ran
            _raises(Past, comsol_batch.main, args + ["--allow-old-air"])
            # mechanics only: not refused
            _raises(Past, comsol_batch.main, ["--template", mph, "--study", "std1",
                                              "--outdir", os.path.join(tmp, "o2")])


# ---- crystal rotation rxt (rib, mechanics only) --------------------------------

def test_resolve_rxtal():
    r = geometry.resolve_rxtal
    assert geometry.RXTAL_ALLOWED_DEG == (0.0, 45.0)
    for x in (None, 0, 0.0, 90, -90, 180, "0", 360.0):
        assert r(x) == 0.0 and isinstance(r(x), float), x
    for x in (45, -45, 135, 45.0, 225):
        assert r(x) == 45.0, x
    for x in (22.5, 30, 44.9, float("nan"), float("inf"), "abc"):
        _raises(ValueError, r, x)
    assert "mirror" in str(_raises(ValueError, r, 30))


def test_rotated_tensor_closed_form_is_rotatextaltensor():
    """The DVo of the builder equals omc-comsol-chang/RotateXtalTensor.m
    (K*C*K', rotation about the slab normal) mapped to the template frame."""
    C11, C12, C44 = 1076.0, 125.0, 577.0
    C = np.zeros((6, 6))
    C[:3, :3] = C12
    np.fill_diagonal(C, [C11, C11, C11, C44, C44, C44])
    with open(os.path.join(_ROOT, "comsol", "build_omc_unitcell_rib.m")) as fh:
        src = fh.read()
    dvo_src = src[src.index("aniso.set('DVo', {"):]
    dvo_src = dvo_src[:dvo_src.index("});")]
    dvo = re.findall(r"'([^']*)'", dvo_src)[1:]                    # skip 'DVo'
    assert dvo == ['C11', 'C12', 'C11-dCrot', 'C12', 'C12+dCrot', 'C11-dCrot',
                   '0', 'gCrot', '-gCrot', 'C44+dCrot', '0', '0', '0', '0',
                   'C44', '0', '0', '0', '0', '0', 'C44']
    upper = [(i, j) for j in range(6) for i in range(j + 1)]       # COMSOL order
    want = {45.0: (1177.5, 23.5, 475.5, 0.0),
            22.5: (1126.75, 74.25, 526.25, -50.75)}
    for deg in (0.0, 22.5, 30.0, 45.0, 90.0):
        th = np.radians(deg)
        c, s = np.cos(th), np.sin(th)
        K = np.array([[c * c, s * s, 0, 0, 0, 2 * s * c],
                      [s * s, c * c, 0, 0, 0, -2 * s * c],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, -s, 0],
                      [0, 0, 0, s, c, 0],
                      [-s * c, s * c, 0, 0, 0, c * c - s * s]])
        P = [2, 1, 0, 5, 4, 3]                     # Voigt permutation [3 2 1 6 5 4]
        ref = (K @ C @ K.T)[np.ix_(P, P)]
        Cani = C11 - C12 - 2 * C44
        env = dict(C11=C11, C12=C12, C44=C44, Cani=Cani,
                   dCrot=Cani / 2 * np.sin(2 * th) ** 2,
                   gCrot=Cani / 4 * np.sin(4 * th))
        D = np.zeros((6, 6))
        for (i, j), e in zip(upper, dvo):
            D[i, j] = D[j, i] = eval(e, {}, env)
        assert np.abs(D - ref).max() < 1e-9, (deg, np.abs(D - ref).max())
        if deg in want:
            got = (D[1, 1], D[1, 2], D[3, 3], env["gCrot"])
            assert np.allclose(got, want[deg], atol=1e-9), (deg, got)
    # the four parameters and the label the index setters match on
    for needle in ("'rxt', '0[deg]'", "'Cani', 'C11-C12-2*C44'",
                   "'dCrot', 'Cani/2*sin(2*rxt)^2'", "'gCrot', 'Cani/4*sin(4*rxt)'",
                   "diamond.label('Diamond (cubic, rotated by rxt about x)');"):
        assert needle in src, needle
    label = re.search(r"diamond\.label\('([^']*)'\)", src).group(1)
    assert label.lower().startswith("diamond")      # set_diamond_index
    assert label.strip()[:7].lower() == "diamond"   # isDiamondLabel (strncmpi)


def test_set_crystal_rotation():
    m = FakeModel(params=("a", "w", "t", "s", "th", "r1", "r2", "kF", "rxt"))
    assert comsol_client.set_crystal_rotation(m, 45) == 45.0
    assert m._params["rxt"] == "45[deg]"
    assert comsol_client.set_crystal_rotation(m, 90) == 0.0
    assert m._params["rxt"] == "0[deg]"
    assert comsol_client.set_crystal_rotation(m) == 0.0
    _raises(ValueError, comsol_client.set_crystal_rotation, m, 30)
    old = mock.MagicMock()
    old.parameters.return_value = {"a": "1", "kF": "1"}
    assert comsol_client.set_crystal_rotation(old, 0) == 0.0
    assert comsol_client.set_crystal_rotation(old, 180) == 0.0
    old.parameter.assert_not_called()
    err = _raises(RuntimeError, comsol_client.set_crystal_rotation, old, 45)
    assert "rxt" in str(err)
    old.parameter.assert_not_called()
    # run_mechanical_comsol sets it right after the geometry, on every call
    import acoustic_comsol
    import inspect
    assert inspect.signature(acoustic_comsol.run_mechanical_comsol) \
        .parameters["rxtal_deg"].default == 0.0
    calls = []
    m2 = FakeModel(params=("a", "w", "t", "s", "th", "r1", "r2", "kF", "rxt"),
                   freqs=(7e9, 9e9))
    with mock.patch.object(acoustic_comsol, "_HAVE_MPH", True), \
            mock.patch.object(acoustic_comsol, "_get_model", lambda t: m2), \
            mock.patch.object(acoustic_comsol, "_set_geometry",
                              lambda model, g: calls.append("geom")), \
            mock.patch.object(acoustic_comsol, "_set_crystal_rotation",
                              lambda model, r: calls.append(("rxt", r))), \
            mock.patch("os.path.exists", lambda p: True):
        g = RibGeometry(**RIB_CELL)
        acoustic_comsol.run_mechanical_comsol(g, n_k=2, n_bands=2, template="x.mph",
                                              rxtal_deg=45.0)
        acoustic_comsol.run_mechanical_comsol(g, n_k=2, n_bands=2, template="x.mph")
    assert calls == ["geom", ("rxt", 45.0), "geom", ("rxt", 0.0)]


def test_rxtal_template_check():
    import zipfile
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "omc_unitcell_rib.mph")
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("dmodel.xml", 's("C11") s("C44") s("air_r")')
        assert comsol_client.template_has_rxtal(path) is False
        assert comsol_client.require_rxtal_template(path, 0) == path
        assert comsol_client.require_rxtal_template(path, 90) == path
        err = _raises(RuntimeError, comsol_client.require_rxtal_template, path, 45)
        assert "rebuild" in str(err) and "rxt" in str(err)
        _raises(ValueError, comsol_client.require_rxtal_template, path, 30)
        os.utime(path, (3, 3))
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("dmodel.xml", 's("rxt") s("Cani") s("dCrot") s("gCrot")')
        assert comsol_client.template_has_rxtal(path) is True
        assert comsol_client.require_rxtal_template(path, 45) == path
        with open(path, "w") as fh:
            fh.write("not a zip")
        os.utime(path, (4, 4))
        assert comsol_client.template_has_rxtal(path) is True     # unreadable


def test_rib_objective_rxtal_in_record_id_and_mechanics():
    import acoustic_comsol
    mech_kw, checked = [], []

    def fake_mech(g, *, n_k=15, n_bands=12, template=None, study="mech sym", **kw):
        mech_kw.append(kw)
        k = np.linspace(0, 1, 5)
        return dict(freqs_hz=np.column_stack([5e9 + 0 * k, 9e9 + 0 * k]),
                    k_frac=k, a=g.a)

    def fake_template(celltype="hole", xsect="rect", must_exist=False):
        return os.path.join("comsol", comsol_client.TEMPLATES[celltype, xsect])

    import yaml
    with open(os.path.join(_ROOT, "configs", "targets.yaml")) as fh:
        targets = yaml.safe_load(fh)
    bounds = geometry.load_bounds()

    def run(**kw):
        return objective_rib.evaluate_rib_candidate(
            _rib_u0(), bounds=bounds, targets=targets, optical_backend="none",
            mech_backend="comsol", mech_gap_mode="complete", **kw)

    with mock.patch.object(acoustic_comsol, "run_mechanical_comsol", fake_mech), \
            mock.patch.object(comsol_client, "template_path", fake_template), \
            mock.patch.object(objective_rib, "_check_rxtal_template",
                              lambda xs, r: checked.append((xs, r))):
        base = run()
        r0, r90, r45 = run(rxtal_deg=0), run(rxtal_deg=90), run(rxtal_deg=45)
        err = _raises(ValueError, run, rxtal_deg=30)
    assert "45" in str(err)
    assert base["status"] == r45["status"] == "success", (base, r45)
    assert base["rxtal_deg"] == r0["rxtal_deg"] == r90["rxtal_deg"] == 0.0
    assert r45["rxtal_deg"] == 45.0
    assert base["id"] == r0["id"] == r90["id"] != r45["id"]
    # unrotated ids carry no rotation tag (only the score version)
    legacy = objective_rib._rib_id(_rib_u0(), (
        "none", "comsol", "complete", 15, 20, "opt TE", 13, 10, 0.30, "score=v4"))
    assert base["id"] == legacy
    assert checked == [("rect", 45.0)]                 # before the solves, once
    assert [kw["rxtal_deg"] for kw in mech_kw] == [0.0] * 6 + [45.0] * 2
    # rib options flow through objective.evaluate_candidate
    with mock.patch.object(acoustic_comsol, "run_mechanical_comsol", fake_mech), \
            mock.patch.object(comsol_client, "template_path", fake_template), \
            mock.patch.object(objective_rib, "_check_rxtal_template",
                              lambda xs, r: None):
        rec = evaluate_candidate(_rib_u0(), celltype="rib", optical_backend="none",
                                 mech_backend="comsol",
                                 rib_options=dict(rxtal_deg=45))
    assert rec["rxtal_deg"] == 45.0 and rec["id"] == r45["id"]
    # the rib configs carry the key, and the loop accepts it
    import run_opt_comsol_rib as ror
    assert "rxtal_deg" in ror._RIB_OPTION_KEYS
    for name in ("run_one_rib", "run_loop_rib", "run_opt_comsol_rib",
                 "run_batch_rib", "run_one_rib_tidy3d", "run_opt_comsol_rib_tidy3d"):
        with open(os.path.join(_ROOT, "configs", name + ".yaml")) as fh:
            assert yaml.safe_load(fh)["rib"]["rxtal_deg"] == 0, name


def test_rib_warm_start_keeps_one_rotation():
    import run_opt_comsol_rib as ror
    added = []

    class Study:
        def add_trial(self, t):
            added.append(t)

    class Opt:
        study = Study()

    bounds = geometry.load_bounds()
    params = u_to_geometry(_rib_u0(), bounds, celltype="rib").as_dict()
    prior = [dict(status="success", params=params, score=0.1),
             dict(status="success", params=params, score=0.2, rxtal_deg=0.0),
             dict(status="success", params=params, score=0.3, rxtal_deg=45.0),
             dict(status="infeasible", params=params, score=-10.0)]
    try:
        import optuna  # noqa: F401
    except ImportError:
        return                                       # optuna not installed
    nd = geometry.n_dim("rib", bounds)
    assert ror._warm_start(Opt(), prior, bounds, nd) == 2
    assert ror._warm_start(Opt(), prior, bounds, nd, rxtal=45.0) == 1


TESTS = [v for k, v in sorted(globals().items(),
                              key=lambda kv: kv[1].__code__.co_firstlineno
                              if callable(kv[1]) and hasattr(kv[1], "__code__") else 0)
         if k.startswith("test_") and callable(v)
         and v.__module__ == __name__]

if __name__ == "__main__":
    for test in TESTS:
        print(test.__name__)
        test()
    print("all %d tests passed" % len(TESTS))
