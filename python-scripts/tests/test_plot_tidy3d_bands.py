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
        bs, ot = ptb.load_module()
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
                                  "--min-guided-q", "1000"])
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
                                  "--minQ", "1000", "--k-units", "ka/pi", "--no-geometry"])
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
        bs, ot = ptb.load_module()
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
            rc, text, err = _run([rib, "--outdir", out1])
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
                                  "--minQ", "1e6", "--min-guided-q", "1000"])
            assert rc == 0, err
            assert "ResonanceFinder cannot be re-run" in text
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


def test_real_results_smoke():
    """The records in results/tidy3d, read only, output to a temp folder."""
    if not _have_tidy3d():
        print("  (skipped: tidy3d is not installed; use the `tidy` conda env)")
        return
    records = ptb.record_folders(ptb.DEFAULT_RESULTS)
    if not records:
        print("  (skipped: no result files under results/tidy3d)")
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
         test_real_results_smoke]

if __name__ == "__main__":
    for test in TESTS:
        print(test.__name__)
        test()
    print("all %d tests passed" % len(TESTS))
