"""Tests of scripts/plot_tidy3d_bands.py (offline analysis of Tidy3D results).

No network: the script blocks connections itself, and the tests that need
tidy3d (skipped where it is not installed; use the `tidy` conda env) build
SYNTHETIC SimulationData files in a temp folder. The real results under
results/tidy3d are only read (smoke test), never written.

Run:  python tests/test_plot_tidy3d_bands.py
"""
import csv
import io
import json
import os
import socket
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import numpy as np

import plot_tidy3d_bands as ptb

C0 = 299_792_458.0


def _have_tidy3d():
    import importlib.util
    try:
        return importlib.util.find_spec("tidy3d") is not None
    except (ImportError, ValueError):
        return False


def _raises(exc, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except exc as err:
        return err
    raise AssertionError(f"{fn.__name__} did not raise {exc.__name__}")


def _modules():
    """(diamond_bandStructure_sim on the real tidy3d, the repo scoring)."""
    ot = ptb.load_scoring()
    return ptb.ModuleLoader(ot).get("run this test"), ot


def _run(argv):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        rc = ptb.main(argv)
    return rc, out.getvalue(), err.getvalue()


# ---- pure parts (no tidy3d) ------------------------------------------------------

def test_script_imports_without_tidy3d_and_names_the_env():
    assert "tidy3d" not in ptb.__dict__                  # lazy import
    with open(os.path.join(_ROOT, "scripts", "plot_tidy3d_bands.py")) as fh:
        src = fh.read()
    code = "\n".join(line.split("#")[0] for line in src.splitlines())
    for banned in ("web.Batch", "web.run", ".upload(", ".download(", "run_cloud(",
                   "import tidy3d.web"):
        assert banned not in code, banned                # local files only
    assert "network_blocked()" in code
    if _have_tidy3d():
        print("  (tidy3d is installed: the missing-package message is not shown)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        folder = os.path.join(tmp, "rib", "cand_x")
        os.makedirs(folder)
        open(os.path.join(folder, "a.hdf5"), "w").close()
        rc, out, err = _run(["cand_x", "--results-dir", tmp])
    assert rc == 2 and "tidy" in err and "conda" in err and out == ""


def test_inputs_are_resolved_to_record_folders():
    with tempfile.TemporaryDirectory() as tmp:
        for rel in ("hole/cand_a/x.hdf5", "hole/cand_a/batch.hdf5",
                    "rib/cand_b/1.hdf5", "rib/cand_b/2.hdf5",
                    "rib/cand_empty/plan.json", "rib/cand_only_batch/batch.hdf5"):
            path = os.path.join(tmp, *rel.split("/"))
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").close()
        a, b = os.path.join(tmp, "hole", "cand_a"), os.path.join(tmp, "rib", "cand_b")
        assert ptb.record_folders(tmp) == [a, b]         # only folders with results
        res = ptb.resolve_inputs
        assert res([], tmp, all_records=True) == [dict(folder=a), dict(folder=b)]
        assert res(["cand_b"], tmp) == [dict(folder=b)]                 # record id
        assert res(["rib/cand_b"], tmp) == [dict(folder=b)]             # celltype/id
        assert res([b], tmp) == [dict(folder=b)]                        # folder
        assert res(["cand_a", a], tmp) == [dict(folder=a)]              # once only
        files = [os.path.join(b, "1.hdf5"), os.path.join(b, "2.hdf5")]
        assert res(files, tmp) == [dict(files=files, folder=b)]         # one job
        err = _raises(ptb.AnalysisError, res, ["cand_zzz"], tmp)
        assert "cand_a" in str(err) and "cand_b" in str(err)            # what exists
        _raises(ptb.AnalysisError, res, ["cand_empty"], tmp)            # no hdf5
        _raises(ptb.AnalysisError, res, ["cand_only_batch"], tmp)       # bookkeeping only
        _raises(ptb.AnalysisError, res, [], tmp)                        # nothing asked
        rc, out, err = _run(["cand_zzz", "--results-dir", tmp])
        assert rc == 2 and "cand_zzz" in err


def test_network_is_blocked_while_the_script_runs():
    with ptb.network_blocked():
        _raises(ptb.AnalysisError, socket.create_connection, ("127.0.0.1", 9))
        s = socket.socket()
        try:
            _raises(ptb.AnalysisError, s.connect, ("127.0.0.1", 9))
        finally:
            s.close()
    assert socket.create_connection.__module__ == "socket"              # restored


def test_comparison_with_the_recorded_gap_and_index():
    scored = dict(found=True, f_lower_THz=175.0, f_upper_THz=210.0,
                  f_center_THz=192.5, normalized_gap=35 / 192.5)
    same = ptb.compare_with_record(scored, dict(gap=dict(
        found=True, f_lower=175.0e12, f_upper=210.0e12, normalized_gap=35 / 192.5)))
    assert same["matches"] and abs(same["d_lower_THz"]) < 1e-9
    other = ptb.compare_with_record(scored, dict(gap=dict(
        found=True, f_lower=237.7e12, f_upper=262.6e12, normalized_gap=0.0998)))
    assert not other["matches"] and abs(other["d_lower_THz"] + 62.7) < 1e-6
    assert not ptb.compare_with_record(dict(found=False), dict(gap=dict(
        found=True, f_lower=1e14, f_upper=2e14, normalized_gap=0.5)))["matches"]
    assert ptb.compare_with_record(dict(found=False), dict(gap=dict(found=False)))["matches"]
    assert ptb.compare_with_record(scored, None) is None
    summary = dict(record_id="cand_x", celltype="rib", xsect="rect", a_um=0.65,
                   counts=dict(k_points=10, resonances=33, guided=22, guided_low_q=5),
                   scored_gap=dict(scored, freq_penalty=0.5, effective_gap=0.09),
                   x_point_gap=None, comparison_with_evaluation_json=same)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "sub", "index.csv")
        ptb.write_index(path, [(summary, "out")])
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    assert rows[0]["record_id"] == "cand_x" and rows[0]["matches_evaluation"] == "1"
    assert float(rows[0]["freq_penalty"]) == 0.5 and rows[0]["x_gap"] == ""


# ---- synthetic SimulationData (real tidy3d, offline) ----------------------------------

def _sim_data(bs, sim, modes):
    """SimulationData whose monitors ring at `modes` = [(f_Hz, Q), ...]."""
    td = bs.td
    rng = np.random.default_rng(int(round(bs.bloch_k(sim) * 1000)))
    datas = []
    for m in sim.monitors:
        t = np.linspace(m.start, sim.run_time, 5000)
        # low-Q modes get a larger amplitude: they have nearly decayed when
        # the monitors start, and the default filter drops amplitudes < 1
        sig = sum((a if q > 1000 else 200 * a) * np.exp(-np.pi * f * t / q)
                  * np.cos(2 * np.pi * f * t + ph)
                  for (f, q), a, ph in zip(modes, rng.uniform(5, 10, len(modes)),
                                           rng.uniform(0, 6, len(modes))))
        arr = td.ScalarFieldTimeDataArray(
            sig[None, None, None, :],
            coords=dict(x=[m.center[0]], y=[m.center[1]], z=[m.center[2]], t=t))
        datas.append(td.FieldTimeData(
            monitor=m, Ey=arr, symmetry=(0, 0, 0), symmetry_center=(0, 0, 0),
            grid_expanded=sim.discretize_monitor(m)))
    return td.SimulationData(simulation=sim, data=tuple(datas))


def _write_record(bs, folder, p, indices, modes_at, with_params=True, evaluation=None):
    """A record folder like the pipeline's: one .hdf5 per k (+ json files)."""
    os.makedirs(folder, exist_ok=True)
    sims = bs.build_simulations(p, x_only=(indices == [p.Nk]))
    for n, i in enumerate(indices):
        sim = sims["sim_%d" % i]
        _sim_data(bs, sim, modes_at(bs.bloch_k(sim))).to_file(
            os.path.join(folder, "fdve-%04d.hdf5" % n))
    with open(os.path.join(folder, "batch.hdf5"), "w") as fh:
        fh.write("not a simulation result")              # as the real folders have
    if with_params:
        with open(os.path.join(folder, "params.json"), "w") as fh:
            json.dump(bs.asdict(p), fh)
    if evaluation is not None:
        with open(os.path.join(folder, "evaluation.json"), "w") as fh:
            json.dump(evaluation, fh)


def _listing(root):
    return sorted((os.path.relpath(os.path.join(d, f), root),
                   os.path.getsize(os.path.join(d, f)), os.path.getmtime(os.path.join(d, f)))
                  for d, _, files in os.walk(root) for f in files)


def test_analysis_and_plots_on_synthetic_results():
    if not _have_tidy3d():
        print("  (skipped: tidy3d is not installed; use the `tidy` conda env)")
        return
    with ptb.network_blocked():
        bs, ot = _modules()
        hi_q = 5e3
        with tempfile.TemporaryDirectory() as tmp:
            results = os.path.join(tmp, "tidy3d")
            # ---- hole, X point only: a = 700 nm, light line at X 214.1 THz ----
            ph = bs.BandParams(geometry="holey", a=0.7, beam_width=0.66, thickness=0.325,
                               hole_x=0.46, hole_y=0.56, index=2.386, Nk=1,
                               target_freq_thz=193.414489)
            hole = os.path.join(results, "hole", "cand_hole")
            _write_record(bs, hole, ph, [1],
                          lambda k: [(174e12, hi_q), (210e12, hi_q), (250e12, 60.0)],
                          evaluation=dict(celltype="hole", gap=dict(
                              found=True, f_lower=174e12, f_upper=210e12,
                              f_center=192e12, normalized_gap=36 / 192)))
            # ---- rib, sweep of 6 k: k*a/pi = 0.5 .. 1, a = 570 nm ----
            pr = bs.BandParams(geometry="rib_block", a=0.57, t=0.391, w=1.747,
                               spine_width=0.254, thickness=0.25, fillet_r1=0.04,
                               fillet_r2=0.04, index=2.386, Nk=10,
                               target_freq_thz=193.414489)

            def rib_modes(k):
                ll = k * C0 / 0.57e-6
                # two flat bands, 140 and 200 THz (both above the 131.5 THz light
                # line of the first k, the lower end of the trusted range)
                modes = [(140e12, hi_q), (200e12, hi_q), (1.1 * ll, 50.0)]
                if abs(k - 0.5) < 1e-9:
                    modes.append((170e12, 40.0))         # spurious low-Q peak in the gap
                return modes

            rib = os.path.join(results, "rib", "cand_rib")
            _write_record(bs, rib, pr, [5, 6, 7, 8, 9, 10], rib_modes)
            # ---- the same rib files without params.json / evaluation.json ----
            bare = os.path.join(results, "rib", "cand_bare")
            _write_record(bs, bare, pr, [9, 10], rib_modes, with_params=False)
            before = _listing(results)

            out = os.path.join(tmp, "out")
            rc, text, err = _run(["--all", "--results-dir", results, "--outdir", out,
                                  "--min-guided-q", "1000", "--score-version", "2"])
            assert rc == 0, err
            assert _listing(results) == before           # the results are untouched
            for name in ("cand_hole", "cand_rib", "cand_bare"):
                for f in ("band_diagram.png", "resonances.csv", "summary.json"):
                    path = os.path.join(out, name, f)
                    assert os.path.getsize(path) > 0, path
            with open(os.path.join(out, "analysis_index.csv"), newline="") as fh:
                index = {r["record_id"]: r for r in csv.DictReader(fh)}
            assert set(index) == {"cand_hole", "cand_rib", "cand_bare"}

            # hole: the tracked gap at X = the module's two lowest guided
            with open(os.path.join(out, "cand_hole", "summary.json")) as fh:
                s = json.load(fh)
            g = s["scored_gap"]
            assert s["celltype"] == "hole" and s["counts"]["k_points"] == 1
            assert g["method"] == "optical_gap_tracked_x" and g["found"]
            assert abs(g["f_lower_THz"] - 174) < 0.05 and abs(g["f_upper_THz"] - 210) < 0.05
            assert "freq_penalty" not in g               # hole: no Gaussian weight
            assert abs(s["x_point_gap"]["normalized_gap"] - g["normalized_gap"]) < 1e-9
            assert s["comparison_with_evaluation_json"]["matches"]
            assert s["counts"]["guided"] == 2 and s["counts"]["resonances"] == 3
            assert s["symmetry"] == [0, -1, 1] and s["optical_tol_thz"] is None
            assert "MATCHES" in text

            # rib sweep: the pipeline's gap choice, with the Gaussian weight
            with open(os.path.join(out, "cand_rib", "summary.json")) as fh:
                s = json.load(fh)
            g = s["scored_gap"]
            assert s["celltype"] == "rib" and s["counts"]["k_points"] == 6
            assert np.allclose(s["k_2pi_over_a"], [0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
            assert g["method"] == "guided_interval_gap" and g["found"]
            # the spurious 170 THz peak at X splits 140-200 into 140-170 and
            # 170-200 (30 THz each > the 26.3 THz step): one of them is "the gap"
            assert round(g["f_upper_THz"] - g["f_lower_THz"]) == 30
            assert s["optical_tol_thz"] == 2.0 and g["freq_penalty"] < 1e-6
            assert len(g["candidates"]) == 2
            # ... and without guided |Q| < 1000 it is the real 140-200 THz gap
            c = s["scored_gap_min_guided_q"]
            assert abs(c["f_lower_THz"] - 140) < 0.05 and abs(c["f_upper_THz"] - 200) < 0.05
            assert s["counts"]["guided_low_q"] == 1
            assert s["x_point_gap"]["found"] and abs(s["x_point_gap"]["f_lower_THz"] - 140) < 0.05
            assert s["comparison_with_evaluation_json"] is None      # no evaluation.json
            with open(os.path.join(out, "cand_rib", "resonances.csv"), newline="") as fh:
                rows = list(csv.DictReader(fh))
            assert set(rows[0]) == {"k_2pi_over_a", "k_a_over_pi", "freq_THz", "Q",
                                    "amplitude", "error", "guided", "in_window"}
            at_x = [r for r in rows if float(r["k_2pi_over_a"]) == 0.5]
            assert len(at_x) == 4 and sum(r["guided"] == "1" for r in at_x) == 3
            leaky = [r for r in at_x if r["guided"] == "0"][0]
            assert abs(float(leaky["freq_THz"]) - 1.1 * 0.5 * C0 / 0.57e-6 * 1e-12) < 0.5
            assert abs(float(leaky["Q"])) < 100          # the low-Q leaky resonance
            assert float(rows[0]["k_a_over_pi"]) == 2 * float(rows[0]["k_2pi_over_a"])

            # no params.json: read from the simulation, and said so
            with open(os.path.join(out, "cand_bare", "summary.json")) as fh:
                s = json.load(fh)
            assert s["params_source"].startswith("simulation")
            assert abs(s["a_um"] - 0.57) < 1e-9 and abs(s["index"] - 2.386) < 1e-6
            assert s["celltype"] == "rib" and s["mode"] == "TE"     # from the folder
            assert s["counts"]["k_points"] == 2 and s["xsect"] == "rect"
            assert "no params.json" in text

            # filters: minQ = 1000 removes the leaky and the spurious peaks
            out2 = os.path.join(tmp, "out2")
            rc, text, err = _run(["rib/cand_rib", "--results-dir", results, "--outdir", out2,
                                  "--minQ", "1000", "--k-units", "ka/pi", "--no-geometry",
                                  "--score-version", "2"])
            assert rc == 0, err
            with open(os.path.join(out2, "summary.json")) as fh:
                s = json.load(fh)
            assert s["counts"]["resonances"] == 12 and s["counts"]["guided_low_q"] == 0
            assert abs(s["scored_gap"]["f_lower_THz"] - 140) < 0.05
            assert abs(s["scored_gap"]["f_upper_THz"] - 200) < 0.05
            assert s["filters"]["minQ"] == 1000.0
            # a wide tolerance: the same gap, a weight that is not ~ 0
            rc, text, err = _run(["cand_rib", "--results-dir", results, "--outdir", out2,
                                  "--minQ", "1000", "--optical-tol-thz", "60"])
            with open(os.path.join(out2, "summary.json")) as fh:
                s = json.load(fh)
            assert abs(s["scored_gap"]["freq_penalty"]
                       - np.exp(-((170 - 193.414489) / 60) ** 2)) < 1e-3

            # explicit files: grouped by the k stored in them; without the X
            # point there is no score, but resonances and plot are still made
            files = [os.path.join(rib, "fdve-%04d.hdf5" % n) for n in (0, 1, 1)]
            out3 = os.path.join(tmp, "out3")
            rc, text, err = _run(files + [os.path.join(rib, "batch.hdf5"),
                                          "--outdir", out3])
            assert rc == 0, err
            with open(os.path.join(out3, "summary.json")) as fh:
                s = json.load(fh)
            assert s["counts"]["k_points"] == 2 and s["x_point_gap"] is None
            assert not s["scored_gap"]["found"] and "X" in s["scored_gap"]["not_scored"]
            assert os.path.getsize(os.path.join(out3, "band_diagram.png")) > 0
            assert _listing(results) == before
            # default output: analysis/ next to the data, existing files kept
            rc, text, err = _run([hole])
            assert rc == 0, err
            assert sorted(os.listdir(os.path.join(hole, "analysis"))) == [
                "band_diagram.png", "resonances.csv", "summary.json"]
            kept = [e for e in _listing(results) if os.sep + "analysis" + os.sep not in e[0]]
            assert kept == before


def test_record_without_hdf5_is_plotted_from_the_saved_analysis():
    """A record pruned by tidy3d.keep_hdf5: bands.csv + params.json only."""
    if not _have_tidy3d():
        print("  (skipped: tidy3d is not installed; use the `tidy` conda env)")
        return
    with ptb.network_blocked():
        bs, ot = _modules()
        with tempfile.TemporaryDirectory() as tmp:
            results = os.path.join(tmp, "tidy3d")
            pr = bs.BandParams(geometry="rib_block", a=0.57, t=0.391, w=1.747,
                               spine_width=0.254, thickness=0.25, fillet_r1=0.04,
                               fillet_r2=0.04, index=2.386, Nk=10,
                               target_freq_thz=193.414489)
            rib = os.path.join(results, "rib", "cand_rib")
            _write_record(bs, rib, pr, [5, 6, 7, 8, 9, 10],
                          lambda k: [(140e12, 5e3), (200e12, 5e3)],
                          evaluation=dict(celltype="rib", gap=dict(
                              found=True, f_lower=140e12, f_upper=200e12,
                              f_center=170e12, normalized_gap=60 / 170)))
            out1 = os.path.join(tmp, "with_files")
            rc, text, err = _run([rib, "--outdir", out1, "--score-version", "2"])
            assert rc == 0, err
            # what the pipeline leaves behind: bands.csv in the module's format
            with open(os.path.join(out1, "resonances.csv"), newline="") as fh:
                rows = list(csv.DictReader(fh))
            with open(os.path.join(rib, "bands.csv"), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["k_2pi_over_a", "freq_THz", "Q", "below_light_line"])
                for r in rows:
                    w.writerow([r["k_2pi_over_a"], r["freq_THz"], r["Q"], r["guided"]])
            for name in os.listdir(rib):                 # ... and no .hdf5 files
                if name.endswith(".hdf5"):
                    os.remove(os.path.join(rib, name))
            assert ptb.record_folders(results) == [rib]  # still a record
            out2 = os.path.join(tmp, "pruned")
            rc, text, err = _run(["cand_rib", "--results-dir", results, "--outdir", out2,
                                  "--minQ", "1e6", "--min-guided-q", "1000",
                                  "--score-version", "2"])
            assert rc == 0, err
            assert "data: saved analysis" in text and "no .hdf5" in text
            assert "are ignored" in err                  # --minQ has no effect here
            with open(os.path.join(out1, "summary.json")) as fh:
                s1 = json.load(fh)
            with open(os.path.join(out2, "summary.json")) as fh:
                s2 = json.load(fh)
            assert s1["data_source"] == "hdf5" and "saved analysis" in s2["data_source"]
            assert s2["counts"] == s1["counts"] and s2["k_2pi_over_a"] == s1["k_2pi_over_a"]
            for key in ("f_lower_THz", "f_upper_THz", "normalized_gap", "freq_penalty"):
                assert abs(s2["scored_gap"][key] - s1["scored_gap"][key]) < 1e-9, key
            assert s2["x_point_gap"] == s1["x_point_gap"]
            assert s2["comparison_with_evaluation_json"]["matches"]
            assert s2["scored_gap_min_guided_q"]["found"]
            assert s2["symmetry"] == [0, -1, 1] and s2["filters"]["minQ"] == 1.0
            assert os.path.getsize(os.path.join(out2, "band_diagram.png")) > 10000
            with open(os.path.join(out2, "resonances.csv"), newline="") as fh:
                assert all(r["amplitude"] == "nan" for r in csv.DictReader(fh))
            # a folder with neither files nor a saved analysis is not a record
            os.remove(os.path.join(rib, "bands.csv"))
            rc, text, err = _run(["cand_rib", "--results-dir", results])
            assert rc == 2 and "cand_rib" in err


# ---- CSV / JSON records: the normal input, no tidy3d needed ---------------------

RIB_PARAMS = dict(geometry="rib_block", mode="TE", xsect="rect", index=2.386, a=0.57,
                  t=0.391, w=1.747, spine_width=0.254, thickness=0.25, fillet_r1=0.04,
                  fillet_r2=0.04, beam_width=0.9, hole_x=0.275, hole_y=0.5,
                  target_freq_thz=193.414489, freq_span_thz=100.0, Nk=10,
                  minQ=1.0, minamp=1.0, maxerr=1000.0)
HOLE_PARAMS = dict(RIB_PARAMS, geometry="holey", a=0.7, beam_width=0.66, thickness=0.325,
                   hole_x=0.46, hole_y=0.56, Nk=1)


def _csv_rows(a_um, k_list, bands, extra=()):
    """(k, f_THz, Q, guided) rows: flat `bands` (THz) at every k, plus `extra`
    rows (k, f_THz, Q)."""
    rows = [(k, f, 2e5) for k in k_list for f in bands] + list(extra)
    return [(k, f, q, int(f * 1e12 < k * C0 / (a_um * 1e-6))) for k, f, q in sorted(rows)]


def _write_csv_record(folder, params=None, rows=(), new_columns=False, evaluation=None,
                      plan=None, hdf5=False):
    """A record as the pipeline leaves it after pruning: bands.csv + JSON."""
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "bands.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        head = ["k_2pi_over_a", "freq_THz", "Q", "below_light_line"]
        w.writerow(head + (["amplitude", "error"] if new_columns else []))
        for k, f, q, g in rows:
            w.writerow([k, f, q, g] + ([12.5, 0.004] if new_columns else []))
    for name, obj in (("params.json", params), ("evaluation.json", evaluation),
                      ("plan.json", plan)):
        if obj is not None:
            with open(os.path.join(folder, name), "w") as fh:
                json.dump(obj, fh)
    if hdf5:
        for name in ("fdve-1.hdf5", "batch.hdf5"):
            with open(os.path.join(folder, name), "wb") as fh:
                fh.write(b"not a real result")


K6 = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]


def _rib_rows(extra=()):
    # bands at 140 and 200 THz, one leaky point per k, as in the synthetic sweep
    leaky = [(k, 1.1 * k * C0 / 0.57e-6 * 1e-12, 50.0) for k in K6]
    return _csv_rows(0.57, K6, [140.0, 200.0], extra=leaky + list(extra))


def _rib_evaluation(tol=None, found=True):
    detail = dict(window_THz=[93.414489, 293.414489], symmetry=[0, -1, 1], k_mode="sweep",
                  band_params={k: RIB_PARAMS[k] for k in (
                      "geometry", "mode", "xsect", "index", "a", "t", "w", "spine_width",
                      "thickness", "fillet_r1", "fillet_r2", "target_freq_thz")})
    if tol is not None:
        detail["optical_tol_thz"] = tol
    return dict(id="x", celltype="rib", detail=detail, gap=dict(
        found=found, f_lower=140e12, f_upper=200e12, f_center=170e12,
        normalized_gap=60 / 170))


def _summary(folder):
    with open(os.path.join(folder, "summary.json")) as fh:
        return json.load(fh)


def test_csv_record_is_plotted_without_tidy3d():
    """The saved analysis is enough: scoring, gaps, Q colours, outline, index."""
    had_tidy3d = "tidy3d" in sys.modules
    with tempfile.TemporaryDirectory() as tmp:
        results = os.path.join(tmp, "results", "tidy3d")
        rib = os.path.join(results, "rib", "cand_rib")
        spurious = [(0.5, 170.0, 40.0)]                  # low-Q guided peak at X
        _write_csv_record(rib, RIB_PARAMS, _rib_rows(spurious),
                          evaluation=_rib_evaluation(tol=5.0))
        hole = os.path.join(results, "hole", "cand_hole")
        _write_csv_record(hole, HOLE_PARAMS, _csv_rows(
            0.7, [0.5], [174.0, 210.0], extra=[(0.5, 250.0, 60.0)]), new_columns=True)
        os.makedirs(os.path.join(results, "rib", "cand_plan_only"))
        with open(os.path.join(results, "rib", "cand_plan_only", "plan.json"), "w") as fh:
            fh.write("{}")                               # a budget stop: no results
        before = _listing(results)
        assert ptb.record_folders(results) == [hole, rib]

        out = os.path.join(tmp, "out")
        rc, text, err = _run(["--all", "--results-dir", results, "--outdir", out,
                              "--min-guided-q", "1000", "--k-units", "ka/pi"])
        assert rc == 0, err
        assert _listing(results) == before               # read only
        if not _have_tidy3d():
            assert ("tidy3d" in sys.modules) == had_tidy3d   # never imported

        s = _summary(os.path.join(out, "cand_rib"))
        assert s["celltype"] == "rib" and s["params_source"] == "params.json"
        assert s["data_source"].startswith("saved analysis") and "no .hdf5" in s["data_source"]
        assert s["counts"]["k_points"] == 6 and s["counts"]["resonances"] == 19
        assert s["counts"]["guided_low_q"] == 1 and s["symmetry"] == [0, -1, 1]
        assert s["optical_tol_thz"] == 5.0 and "evaluation.json" in s["optical_tol_source"]
        assert s["score_version"] == 2                   # a Gaussian tolerance was recorded
        g = s["scored_gap"]                              # 140-170 / 170-200: split by the peak
        assert g["method"] == "guided_interval_gap" and g["found"]
        assert round(g["f_upper_THz"] - g["f_lower_THz"]) == 30 and len(g["candidates"]) == 2
        c = s["scored_gap_min_guided_q"]                 # the cleaning works on the CSV's Q
        assert abs(c["f_lower_THz"] - 140) < 1e-9 and abs(c["f_upper_THz"] - 200) < 1e-9
        assert abs(c["freq_penalty"] - np.exp(-((170 - 193.414489) / 5.0) ** 2)) < 1e-9
        assert s["x_point_gap"] == dict(found=True, f_lower_THz=140.0, f_upper_THz=170.0,
                                        f_center_THz=155.0, normalized_gap=30 / 155)
        assert s["comparison_with_evaluation_json"]["matches"] is False   # 140-200 recorded
        assert abs(s["index"] - 2.386) < 1e-12 and s["mode"] == "TE" and s["xsect"] == "rect"
        assert os.path.getsize(os.path.join(out, "cand_rib", "band_diagram.png")) > 10000
        with open(os.path.join(out, "cand_rib", "resonances.csv"), newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 19 and all(r["amplitude"] == "nan" for r in rows)   # old format
        # 140 THz is below the light line at 5 of the 6 k, 200 THz at 3, + the peak
        assert sum(r["guided"] == "1" for r in rows) == 9

        s = _summary(os.path.join(out, "cand_hole"))
        assert s["celltype"] == "hole" and s["optical_tol_thz"] is None
        assert s["scored_gap"]["method"] == "optical_gap_tracked_x"
        assert abs(s["scored_gap"]["normalized_gap"] - 36 / 192) < 1e-9
        assert "no amplitude" not in s["data_source"]    # the new columns are read
        with open(os.path.join(out, "cand_hole", "resonances.csv"), newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert {r["amplitude"] for r in rows} == {"12.5"} and {r["error"] for r in rows} == {"0.004"}

        with open(os.path.join(out, "analysis_index.csv"), newline="") as fh:
            index = {r["record_id"]: r for r in csv.DictReader(fh)}
        assert set(index) == {"cand_rib", "cand_hole"}   # the plan-only folder is no record
        assert index["cand_rib"]["data_source"].startswith("saved analysis")
        assert float(index["cand_rib"]["optical_tol_thz"]) == 5.0
        assert "optical_tol_thz = 5 (evaluation.json)" in text


def test_bands_csv_given_directly_and_parameter_fallbacks():
    with tempfile.TemporaryDirectory() as tmp:
        # 1) bands.csv path + params.json next to it
        a = os.path.join(tmp, "anywhere", "run_a")
        _write_csv_record(a, RIB_PARAMS, _rib_rows())
        rc, text, err = _run([os.path.join(a, "bands.csv"), "--outdir", os.path.join(tmp, "o1")])
        assert rc == 0, err
        s = _summary(os.path.join(tmp, "o1"))
        assert s["record_id"] == "run_a" and s["params_source"] == "params.json"
        assert s["celltype"] == "rib" and abs(s["scored_gap"]["f_lower_THz"] - 140) < 1e-9
        # nothing recorded about the scoring: the current version (4)
        assert s["score_version"] == 4 and s["optical_tol_thz"] is None
        assert "targets.yaml" in s["scoring"]
        # 2) no params.json: evaluation.json, then plan.json
        b = os.path.join(tmp, "run_b")
        _write_csv_record(b, None, _rib_rows(), evaluation=_rib_evaluation())
        rc, text, err = _run([b, "--outdir", os.path.join(tmp, "o2")])
        assert rc == 0, err
        s2 = _summary(os.path.join(tmp, "o2"))
        assert s2["params_source"].startswith("evaluation.json")
        assert s2["a_um"] == 0.57 and s2["window_THz"] == s["window_THz"]
        assert s2["scored_gap"]["f_upper_THz"] == s["scored_gap"]["f_upper_THz"]
        assert s2["comparison_with_evaluation_json"]["matches"]
        c = os.path.join(tmp, "run_c")
        _write_csv_record(c, None, _rib_rows(), plan=_rib_evaluation()["detail"])
        rc, text, err = _run([c, "--outdir", os.path.join(tmp, "o3")])
        assert rc == 0 and _summary(os.path.join(tmp, "o3"))["params_source"].startswith("plan.json")
        # 3) a bare bands.csv: the period must be given
        d = os.path.join(tmp, "run_d")
        _write_csv_record(d, None, _rib_rows())
        rc, text, err = _run([os.path.join(d, "bands.csv"), "--outdir", os.path.join(tmp, "o4")])
        assert rc == 2 and "--a-um" in err
        rc, text, err = _run([os.path.join(d, "bands.csv"), "--outdir", os.path.join(tmp, "o4"),
                              "--a-um", "0.57", "--celltype", "rib"])
        assert rc == 0, err
        s4 = _summary(os.path.join(tmp, "o4"))
        assert s4["params_source"].startswith("--a-um") and s4["celltype"] == "rib"
        assert s4["scored_gap"]["f_lower_THz"] == s["scored_gap"]["f_lower_THz"]
        assert "geometry not stored" in open(os.path.join(tmp, "o4", "summary.json")).read() \
            or s4["geometry"] is None
        # 4) not a bands.csv of the pipeline, and an empty one
        bad = os.path.join(tmp, "bad.csv")
        with open(bad, "w") as fh:
            fh.write("x,y\n1,2\n")
        rc, text, err = _run([bad, "--a-um", "0.57"])
        assert rc == 2 and "bands.csv" in err
        with open(bad, "w") as fh:
            fh.write("k_2pi_over_a,freq_THz,Q,below_light_line\n")
        rc, text, err = _run([bad, "--a-um", "0.57"])
        assert rc == 2 and "no resonances" in err
        # LiteParams: index from the material name, symmetry by mode and xsect
        lp = ptb.LiteParams(dict(a=0.6, material="diamond_telecom", mode="TM", xsect="wedge"))
        assert lp.n == 2.386 and lp.symmetry == (0, 1, 0)
        assert ptb.LiteParams(dict(a=0.6)).symmetry == (0, -1, 1)
        assert abs(ptb.LiteParams(dict(a=0.6)).freq_range[0] - 93.4e12) < 1e11


def test_csv_is_preferred_and_reanalyze_needs_the_files_and_tidy3d():
    with tempfile.TemporaryDirectory() as tmp:
        results = os.path.join(tmp, "tidy3d")
        both = os.path.join(results, "rib", "cand_both")       # CSV and (unreadable) .hdf5
        _write_csv_record(both, RIB_PARAMS, _rib_rows(), hdf5=True)
        pruned = os.path.join(results, "rib", "cand_pruned")   # CSV only
        _write_csv_record(pruned, RIB_PARAMS, _rib_rows())
        # default: the CSV, even with .hdf5 files next to it (they are not opened)
        rc, text, err = _run(["cand_both", "--results-dir", results, "--outdir",
                              os.path.join(tmp, "o1")])
        assert rc == 0, err
        s = _summary(os.path.join(tmp, "o1"))
        assert s["data_source"].startswith("saved analysis") and "--reanalyze" in s["data_source"]
        # --reanalyze / a filter option on a record without files: ignored, with a note
        for opt in (["--reanalyze"], ["--minQ", "1000"], ["--freq-max-thz", "250"]):
            rc, text, err = _run(["cand_pruned", "--results-dir", results, "--outdir",
                                  os.path.join(tmp, "o2")] + opt)
            assert rc == 0 and "are ignored" in err, opt
            s2 = _summary(os.path.join(tmp, "o2"))
            assert s2["data_source"].startswith("saved analysis")
            assert s2["filters"]["minQ"] == 1.0 and s2["window_THz"] == s["window_THz"]
        if _have_tidy3d():
            print("  (tidy3d is installed: the missing-package errors are not shown)")
            return
        # --reanalyze / filters with files but no tidy3d: a clear error, nothing written
        for opt in (["--reanalyze"], ["--minQ", "1000"]):
            rc, text, err = _run(["cand_both", "--results-dir", results, "--outdir",
                                  os.path.join(tmp, "o3")] + opt)
            assert rc == 2 and "tidy" in err and "conda" in err and "bands.csv" in err, opt
        assert not os.path.exists(os.path.join(tmp, "o3", "summary.json"))
        # explicit .hdf5 files need it too
        rc, text, err = _run([os.path.join(both, "fdve-1.hdf5")])
        assert rc == 2 and "tidy" in err


def test_tolerance_is_read_from_the_run_record_read_only():
    import sqlite3
    with tempfile.TemporaryDirectory() as tmp:
        results = os.path.join(tmp, "results")
        folder = os.path.join(results, "tidy3d", "rib", "cand_abc")
        _write_csv_record(folder, RIB_PARAMS, _rib_rows(), evaluation=_rib_evaluation())
        db = os.path.join(results, "runs_rib.sqlite3")
        con = sqlite3.connect(db)
        con.execute("CREATE TABLE candidates (id TEXT PRIMARY KEY, record TEXT)")
        con.execute("INSERT INTO candidates VALUES (?, ?)",
                    ("cand_abc", json.dumps(dict(id="cand_abc", optical_tol_thz=7.5))))
        con.commit()
        con.close()
        stamp = (os.path.getsize(db), os.path.getmtime(db))
        assert ptb.recorded_tolerance(folder, _rib_evaluation()) == (7.5, "runs_rib.sqlite3")
        assert ptb.recorded_tolerance(folder, _rib_evaluation(tol=3.0)) == (3.0, "evaluation.json")
        other = os.path.join(results, "tidy3d", "rib", "cand_unknown")
        assert ptb.recorded_tolerance(other, None) == (None, None)
        rc, text, err = _run([folder, "--outdir", os.path.join(tmp, "o")])
        assert rc == 0, err
        s = _summary(os.path.join(tmp, "o"))
        assert s["optical_tol_thz"] == 7.5 and "runs_rib.sqlite3" in s["optical_tol_source"]
        assert abs(s["scored_gap"]["freq_penalty"]
                   - np.exp(-((170 - 193.414489) / 7.5) ** 2)) < 1e-9
        # the command line wins
        rc, text, err = _run([folder, "--outdir", os.path.join(tmp, "o"),
                              "--optical-tol-thz", "60"])
        assert _summary(os.path.join(tmp, "o"))["optical_tol_thz"] == 60.0
        assert (os.path.getsize(db), os.path.getmtime(db)) == stamp   # not modified
        assert sorted(os.listdir(results)) == ["runs_rib.sqlite3", "tidy3d"]   # no journal


def test_rib_score_versions_are_shown_as_their_runs_scored_them():
    """Current rib scoring (version 4): the gap with the largest G_o -
    penalties with the flat-top detuning penalty, the terms in the summary;
    version 3 records keep the quadratic, version 2 records their Gaussian."""
    import sqlite3
    import yaml
    sys.path.insert(0, os.path.join(_ROOT, "src"))
    import objective_rib
    with open(os.path.join(_ROOT, "configs", "targets.yaml")) as fh:
        targets = yaml.safe_load(fh)
    w4 = objective_rib.rib_score_weights(targets)
    w_now = objective_rib.rib_score_weights(targets, version=3)   # as a v3 run had
    lof, tol = w4["lambda_optical_freq"], w4["optical_freq_tol_THz"]
    flat = lambda df, t=tol: lof * (1 - 1 / (1 + (abs(df) / t) ** 8))   # noqa: E731
    f_t = 193.414489
    with tempfile.TemporaryDirectory() as tmp:
        results = os.path.join(tmp, "results")
        # bands at 140 / 185 / 200 THz: gaps 140-185 (27.7 %, off target) and
        # 185-200 (7.8 %, near target). 6 k points: resolvable step 26.3 THz,
        # so only 140-185 qualifies; add a third record with 13 k for both.
        k13 = [0.5 * (0.478 + i * (1 - 0.478) / 12) for i in range(13)]
        rows = _csv_rows(0.57, k13, [140.0, 185.0, 200.0])
        new = os.path.join(results, "tidy3d", "rib", "cand_v3")
        _write_csv_record(new, RIB_PARAMS, rows, new_columns=True)
        old = os.path.join(results, "tidy3d", "rib", "cand_v2")
        _write_csv_record(old, RIB_PARAMS, rows, evaluation=_rib_evaluation(tol=5.0))
        db = os.path.join(results, "runs_rib.sqlite3")
        con = sqlite3.connect(db)
        con.execute("CREATE TABLE candidates (id TEXT PRIMARY KEY, record TEXT)")
        w_rec = dict(w_now, lambda_optical_freq=90.0)          # the run's own weights
        con.execute("INSERT INTO candidates VALUES (?, ?)", ("cand_v3", json.dumps(
            dict(id="cand_v3", score_version=3, score_weights=w_rec))))
        # a version-4 run with its own (very wide) flat top
        w4_rec = dict(w4, optical_freq_tol_THz=100.0)
        _write_csv_record(os.path.join(results, "tidy3d", "rib", "cand_v4"),
                          RIB_PARAMS, rows, new_columns=True)
        con.execute("INSERT INTO candidates VALUES (?, ?)", ("cand_v4", json.dumps(
            dict(id="cand_v4", score_version=4, score_weights=w4_rec))))
        con.commit()
        con.close()
        out = os.path.join(tmp, "out")
        rc, text, err = _run(["--all", "--results-dir", os.path.join(results, "tidy3d"),
                              "--outdir", out])
        assert rc == 0, err
        s3, s2 = _summary(os.path.join(out, "cand_v3")), _summary(os.path.join(out, "cand_v2"))
        # version 3, with the weights recorded by the run (strong detuning weight)
        assert s3["score_version"] == 3 and s3["score_weights"] == w_rec
        assert s3["optical_tol_thz"] is None and "runs_rib.sqlite3" in s3["scoring"]
        g = s3["scored_gap"]
        assert "freq_penalty" not in g and len(g["candidates"]) == 2
        assert abs(g["f_lower_THz"] - 185) < 1e-9 and abs(g["f_upper_THz"] - 200) < 1e-9
        gap = 15 / 192.5
        assert abs(g["pen_size_o"] - 10 * (0.2 - gap) ** 2) < 1e-12
        assert abs(g["pen_detune_o"] - 90 * ((192.5 - f_t) / f_t) ** 2) < 1e-12
        assert g["pen_lightcone_o"] == 0.0
        assert abs(g["value"] - (gap - g["pen_size_o"] - g["pen_detune_o"])) < 1e-12
        assert all("value" in c for c in g["candidates"]) and g["candidates"][0]["chosen"]
        assert "optical score term" in text and "light cone" in text
        assert "optical_freq_tol_THz" not in s3["score_weights"]
        assert "f_tol" not in s3["scoring"]
        # version 4 with the run's weights: flat top of 100 THz, so the wide
        # off-target gap 140-185 THz is the one that was scored
        s4 = _summary(os.path.join(out, "cand_v4"))
        assert s4["score_version"] == 4 and s4["score_weights"] == w4_rec
        assert "f_tol = 100 THz" in s4["scoring"] and "runs_rib.sqlite3" in s4["scoring"]
        g4 = s4["scored_gap"]
        assert abs(g4["f_lower_THz"] - 140) < 1e-9 and abs(g4["f_upper_THz"] - 185) < 1e-9
        assert abs(g4["pen_detune_o"] - flat(162.5 - f_t, 100.0)) < 1e-12
        assert "freq_penalty" not in g4 and "pen_center_m" not in g4
        # version 2 record: shown as it was scored, Gaussian weight and tolerance
        assert s2["score_version"] == 2 and s2["optical_tol_thz"] == 5.0
        assert abs(s2["scored_gap"]["freq_penalty"]
                   - np.exp(-((192.5 - f_t) / 5.0) ** 2)) < 1e-9
        assert "pen_detune_o" not in s2["scored_gap"]
        with open(os.path.join(out, "analysis_index.csv"), newline="") as fh:
            index = {r["record_id"]: r for r in csv.DictReader(fh)}
        assert index["cand_v3"]["score_version"] == "3" and index["cand_v2"]["score_version"] == "2"
        assert index["cand_v4"]["score_version"] == "4"
        assert abs(float(index["cand_v4"]["pen_detune_o"]) - g4["pen_detune_o"]) < 1e-12
        assert abs(float(index["cand_v3"]["pen_detune_o"]) - g["pen_detune_o"]) < 1e-12
        assert index["cand_v2"]["pen_detune_o"] == "" and index["cand_v3"]["freq_penalty"] == ""
        # nothing recorded: version 4 with today's targets.yaml. The on-target
        # gap 185-200 THz is scored: 140-185 THz is 31 THz off, far outside
        # the flat top
        other = os.path.join(results, "tidy3d", "rib", "cand_plain")
        _write_csv_record(other, RIB_PARAMS, rows)
        rc, text, err = _run([other, "--outdir", os.path.join(tmp, "o3")])
        assert rc == 0, err
        s = _summary(os.path.join(tmp, "o3"))
        assert s["score_version"] == 4 and s["score_weights"] == w4
        assert "f_tol = %g THz" % tol in s["scoring"] and "f_tol = %g THz" % tol in text
        if lof >= 0.5 and tol <= 15:
            assert abs(s["scored_gap"]["f_lower_THz"] - 185) < 1e-9      # 185-200 THz
            assert abs(s["scored_gap"]["pen_detune_o"] - flat(192.5 - f_t)) < 1e-12
            values = sorted(c["value"] for c in s["scored_gap"]["candidates"])
            assert abs(values[0] - (45 / 162.5 - flat(162.5 - f_t))) < 1e-9
        assert os.path.getsize(os.path.join(tmp, "o3", "band_diagram.png")) > 10000
        # overrides: version 3 forced (today's weights read as version 3:
        # quadratic detuning, the wide gap wins); version 4 forced on the
        # version-3 record (its recorded weights are not version 4's);
        # version 2 forced on the new record
        rc, text, err = _run([other, "--outdir", os.path.join(tmp, "o3b"),
                              "--score-version", "3"])
        s = _summary(os.path.join(tmp, "o3b"))
        assert s["score_version"] == 3 and s["score_weights"] == w_now
        assert abs(s["scored_gap"]["f_lower_THz"] - 140) < 1e-9          # 140-185 THz
        assert abs(s["scored_gap"]["pen_detune_o"]
                   - w_now["lambda_optical_freq"] * ((162.5 - f_t) / f_t) ** 2) < 1e-12
        rc, text, err = _run([new, "--outdir", os.path.join(tmp, "o3c"),
                              "--score-version", "4"])
        s = _summary(os.path.join(tmp, "o3c"))
        assert s["score_version"] == 4 and s["score_weights"] == w4
        assert "configs/targets.yaml" in s["scoring"]
        rc, text, err = _run([new, "--outdir", os.path.join(tmp, "o4"),
                              "--score-version", "2", "--min-guided-q", "1000"])
        s = _summary(os.path.join(tmp, "o4"))
        assert s["score_version"] == 2 and s["optical_tol_thz"] == 2.0
        assert "freq_penalty" in s["scored_gap"] and "freq_penalty" in s["scored_gap_min_guided_q"]
        # a rib record with no gap at all: the size penalty, centre = target
        empty = os.path.join(results, "tidy3d", "rib", "cand_none")
        _write_csv_record(empty, RIB_PARAMS, _csv_rows(0.57, k13, [140.0]))
        rc, text, err = _run([empty, "--outdir", os.path.join(tmp, "o5")])
        g = _summary(os.path.join(tmp, "o5"))["scored_gap"]
        assert not g["found"] and abs(g["pen_size_o"] - 0.4) < 1e-12 and g["pen_detune_o"] == 0.0
        assert os.path.getsize(os.path.join(tmp, "o5", "band_diagram.png")) > 10000


def test_reanalyze_with_real_tidy3d_and_new_bands_csv_columns():
    """With tidy3d: the CSV is still the default, --reanalyze reads the files,
    and the module now saves amplitude and fit error in bands.csv."""
    if not _have_tidy3d():
        print("  (skipped: tidy3d is not installed; use the `tidy` conda env)")
        return
    with ptb.network_blocked():
        bs, ot = _modules()
        with tempfile.TemporaryDirectory() as tmp:
            results = os.path.join(tmp, "tidy3d")
            pr = bs.BandParams(**{k: v for k, v in RIB_PARAMS.items()})
            rib = os.path.join(results, "rib", "cand_rib")
            _write_record(bs, rib, pr, [5, 6, 7, 8, 9, 10],
                          lambda k: [(140e12, 5e3), (200e12, 5e3)])
            # bands.csv exactly as the pipeline writes it now
            result = bs.analyze_bands(bs.load_run(rib), pr)
            assert set(result.bands[0]) == {"k", "freq", "Q", "amplitude", "error"}
            bs.save_bands_csv(os.path.join(rib, "bands.csv"), result, pr)
            with open(os.path.join(rib, "bands.csv"), newline="") as fh:
                head = next(csv.reader(fh))
            assert head == ["k_2pi_over_a", "freq_THz", "Q", "below_light_line",
                            "amplitude", "error"]
            f, q = bs.find_resonances(next(iter(bs.load_run(rib).values())).data, pr)
            assert len(f) == len(q) == 2                 # the old function is unchanged
            # a BandResult without the new arrays (older callers) still saves
            old = bs.BandResult(bands=[dict(k=0.5, freq=np.array([1.4e14]), Q=np.array([9.0]))],
                                f_light_X=pr.f_light_X, has_X=True)
            bs.save_bands_csv(os.path.join(tmp, "old.csv"), old, pr)
            assert ptb.rows_from_bands_csv(os.path.join(tmp, "old.csv"))[0]["freq"][0] == 1.4e14

            rc, text, err = _run([rib, "--outdir", os.path.join(tmp, "csv")])
            assert rc == 0, err
            rc, text, err = _run([rib, "--outdir", os.path.join(tmp, "h5"), "--reanalyze"])
            assert rc == 0, err
            s_csv, s_h5 = _summary(os.path.join(tmp, "csv")), _summary(os.path.join(tmp, "h5"))
            assert s_csv["data_source"].startswith("saved analysis")
            assert "no amplitude" not in s_csv["data_source"] and s_h5["data_source"] == "hdf5"
            assert s_csv["counts"] == s_h5["counts"]
            for key in ("f_lower_THz", "f_upper_THz", "normalized_gap"):
                assert abs(s_csv["scored_gap"][key] - s_h5["scored_gap"][key]) < 1e-9
            with open(os.path.join(tmp, "csv", "resonances.csv"), newline="") as fh:
                a_csv = [float(r["amplitude"]) for r in csv.DictReader(fh)]
            with open(os.path.join(tmp, "h5", "resonances.csv"), newline="") as fh:
                a_h5 = [float(r["amplitude"]) for r in csv.DictReader(fh)]
            assert np.allclose(a_csv, a_h5) and min(a_csv) > 0
            # a filter option re-reads the files without --reanalyze
            rc, text, err = _run([rib, "--outdir", os.path.join(tmp, "f"), "--minQ", "1e6"])
            s_f = _summary(os.path.join(tmp, "f"))
            assert s_f["data_source"] == "hdf5" and s_f["counts"]["resonances"] == 0


def test_real_results_smoke():
    """The records in results/tidy3d, read only, output to a temp folder.
    Without tidy3d only records with a saved bands.csv can be read."""
    records = [f for f in ptb.record_folders(ptb.DEFAULT_RESULTS)
               if _have_tidy3d() or ptb._has_cache(f)]
    if not records:
        print("  (skipped: no readable record under results/tidy3d)")
        return
    smallest = min(records, key=lambda f: sum(
        os.path.getsize(os.path.join(f, n)) for n in os.listdir(f)))
    before = _listing(smallest)
    with tempfile.TemporaryDirectory() as tmp:
        rc, text, err = _run([smallest, "--outdir", tmp])
        assert rc == 0, err
        with open(os.path.join(tmp, "summary.json")) as fh:
            s = json.load(fh)
        assert s["counts"]["k_points"] >= 1 and s["counts"]["resonances"] >= 1
        assert s["record_id"] == os.path.basename(smallest)
        assert os.path.getsize(os.path.join(tmp, "band_diagram.png")) > 10000
    assert _listing(smallest) == before
    print("  %s: %d k, %d resonances, gap found: %s"
          % (s["record_id"], s["counts"]["k_points"], s["counts"]["resonances"],
             s["scored_gap"]["found"]))


TESTS = [test_script_imports_without_tidy3d_and_names_the_env,
         test_inputs_are_resolved_to_record_folders,
         test_network_is_blocked_while_the_script_runs,
         test_comparison_with_the_recorded_gap_and_index,
         test_analysis_and_plots_on_synthetic_results,
         test_record_without_hdf5_is_plotted_from_the_saved_analysis,
         test_csv_record_is_plotted_without_tidy3d,
         test_bands_csv_given_directly_and_parameter_fallbacks,
         test_csv_is_preferred_and_reanalyze_needs_the_files_and_tidy3d,
         test_tolerance_is_read_from_the_run_record_read_only,
         test_rib_score_versions_are_shown_as_their_runs_scored_them,
         test_reanalyze_with_real_tidy3d_and_new_bands_csv_columns,
         test_real_results_smoke]

if __name__ == "__main__":
    for test in TESTS:
        print(test.__name__)
        test()
    print("all %d tests passed" % len(TESTS))
