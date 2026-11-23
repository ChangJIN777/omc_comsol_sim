"""Tests for scripts/run_batch_rib.py and comsol_client.release_client.

No COMSOL, no MPh, no licence: the solver-free backends (mech: surrogate_stub,
optical: none) and a fake client stand in. Run:  python tests/test_batch_queue.py
"""
import contextlib
import io
import json
import os
import sys
import tempfile

_ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import comsol_client          # noqa: E402
import database               # noqa: E402
import run_batch_rib as rb    # noqa: E402

U0 = [0.385714, 0.575254, 0.358947, 0.214737]      # configs/run_one_rib.yaml


class FakeServer:
    def __init__(self):
        self.stopped = False

    def running(self):
        return not self.stopped

    def stop(self):
        self.stopped = True


class FakeClient:
    """Stands in for mph.Client: records what release_client does to it."""
    def __init__(self):
        self.port, self.cleared, self.disconnected = 2036, False, False

    def clear(self):
        self.cleared = True

    def disconnect(self):
        self.disconnected, self.port = True, None

    def remove(self, model):
        pass


@contextlib.contextmanager
def _sandbox(queue_yaml):
    """Temp dir with a queue config; rib records go to a temp database."""
    old = (database._DB_RIB, comsol_client._client, comsol_client._released,
           rb.evaluate_candidate)
    with tempfile.TemporaryDirectory() as tmp:
        database._DB_RIB = os.path.join(tmp, "runs_rib.sqlite3")
        comsol_client._client, comsol_client._released = None, False
        cfg = os.path.join(tmp, "queue.yaml")
        with open(cfg, "w") as fh:
            fh.write(queue_yaml)
        try:
            yield tmp, cfg
        finally:
            (database._DB_RIB, comsol_client._client, comsol_client._released,
             rb.evaluate_candidate) = old


def _run(cfg, *extra):
    """main() with stdout/stderr captured: (exit code, stdout JSON, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = rb.main(["--config", cfg, *extra])
    return rc, json.loads(out.getvalue()), err.getvalue()


_STUB = "celltype: rib\noptical: none\nmech: surrogate_stub\nquiet: true\n"


def _fake_eval(script):
    """evaluate_candidate stand-in. script[i] is what item i does:
    'ok', 'licence', 'interrupt', 'boom'. Opens a fake session like a real
    COMSOL solve would."""
    calls = []

    def fake(u, **kw):
        i = len(calls)
        calls.append(dict(u=list(u), session_before=comsol_client.session_active(), **kw))
        if kw["mech_backend"] == "comsol" or kw["optical_backend"] == "comsol":
            if comsol_client._client is None:
                comsol_client._client = FakeClient()
        what = script[i]
        if what == "interrupt":
            raise KeyboardInterrupt
        if what == "boom":
            raise RuntimeError("boom")
        rec = dict(id=f"cand_fake{i:04d}", u=list(u), params={}, celltype="rib",
                   optical_backend=kw["optical_backend"],
                   mech_backend=kw["mech_backend"], status="success",
                   optical_gap=None, mechanical_gap=0.1, score=0.1,
                   optical_center_frequency=None, mechanical_center_frequency=8e9)
        if what == "licence":
            rec.update(status="mech_failed", score=-2.0, mechanical_gap=0.0,
                       error="License error: -4. Licensed number of users "
                             "already reached. Feature: STRUCTURALMECHANICS")
        return rec
    fake.calls = calls
    return fake


def test_release_client_is_safe_and_final():
    old = (comsol_client._client, comsol_client._released)
    try:
        comsol_client._client, comsol_client._released = None, False
        assert comsol_client.release_client() is False      # nothing open
        assert comsol_client._released is False              # ... so not "released"
        c = comsol_client._client = FakeClient()
        assert comsol_client.session_active()
        assert comsol_client.release_client() is True
        assert c.cleared and c.disconnected
        assert not comsol_client.session_active()
        assert comsol_client.release_client() is False       # twice is harmless
        if comsol_client.have_mph():
            try:
                comsol_client.get_client()                   # must not restart
            except RuntimeError as e:
                assert "released" in str(e)
            else:
                raise AssertionError("get_client restarted a released session")
    finally:
        comsol_client._client, comsol_client._released = old


def test_queue_forms_and_validation():
    with _sandbox(_STUB) as (tmp, cfg):
        with open(os.path.join(tmp, "designs.csv"), "w") as fh:
            fh.write("ua,uw,ut,us\n# comment\n0.4,0.55,0.36,0.22\n\n0.5,0.5,0.5,0.5\n")
        with open(os.path.join(tmp, "one.yaml"), "w") as fh:
            fh.write(f"celltype: rib\nu: {U0}\noptical: none\nmech: surrogate_stub\n"
                     "rib: {n_k: 7}\nquiet: true\n")
        cfg_d = dict(celltype="rib", optical="none", mech="surrogate_stub",
                     rib=dict(n_k=9, mech_gap_mode="sym"),
                     queue=[U0, dict(u=U0, rib=dict(n_k=5), label="x"),
                            dict(config="one.yaml")],
                     queue_csv="designs.csv")
        items = rb.build_queue(cfg_d, tmp)
        assert len(items) == 5
        assert items[0]["rib"] == dict(n_k=9, mech_gap_mode="sym")
        assert items[1]["rib"] == dict(n_k=5, mech_gap_mode="sym")   # merged
        assert items[2]["rib"] == dict(n_k=7)                        # config's own
        assert items[3]["u"] == [0.4, 0.55, 0.36, 0.22]
        for bad in (dict(cfg_d, queue=[[0.1, 0.2]], queue_csv=None),       # length
                    dict(cfg_d, queue=[[0.1, 0.2, 0.3, 1.5]], queue_csv=None),
                    dict(cfg_d, queue=[dict(u=U0, typo=1)], queue_csv=None),
                    dict(cfg_d, queue=[], queue_csv=None),
                    dict(celltype="rib", queue=[U0])):                    # no backends
            try:
                rb.build_queue(bad, tmp)
            except ValueError:
                pass
            else:
                raise AssertionError(f"accepted a bad queue: {bad}")


def test_queue_saves_one_record_per_item_like_run_one():
    infeasible = [0.0, 0.0, 1.0, 1.0]          # t > a: fails the geometry gate
    with _sandbox(_STUB + f"queue:\n  - {U0}\n  - {infeasible}\n"
                          "  - {u: [0.40, 0.55, 0.36, 0.22], label: b}\n") as (tmp, cfg):
        rc, out, err = _run(cfg)
        assert rc == 0 and out["n_done"] == 3 and out["stopped"] is None
        assert out["not_run"] == []
        saved = database.load_completed_results(database.db_path("rib"))
        assert sorted(r["id"] for r in saved) == sorted(r["id"] for r in out["results"])
        assert len({r["id"] for r in saved}) == 3
        by_id = {r["id"]: r for r in saved}
        statuses = [by_id[r["id"]]["status"] for r in out["results"]]
        assert statuses[1] in ("infeasible", "unfabricable")     # failure saved too
        assert all(r["celltype"] == "rib" for r in saved)
        # identical record to what run_one.py would save for the same item
        from objective import evaluate_candidate
        ref = evaluate_candidate(U0, optical_backend="none",
                                 mech_backend="surrogate_stub", celltype="rib",
                                 rib_options=None)
        assert by_id[ref["id"]] == ref
        assert not comsol_client.session_active()                # never opened one
        assert "COMSOL session released" not in err


def test_dry_run_solves_nothing():
    with _sandbox(_STUB + f"queue:\n  - {U0}\n") as (tmp, cfg):
        rb.evaluate_candidate = _fake_eval(["boom"])
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            assert rb.main(["--config", cfg, "--dry-run"]) == 0
        assert json.loads(out.getvalue())[0]["u"] == U0
        assert rb.evaluate_candidate.calls == []
        assert not os.path.exists(database.db_path("rib"))


_COMSOL = "celltype: rib\noptical: none\nmech: comsol\nquiet: true\n"


def test_session_is_reused_then_released_when_queue_is_empty():
    with _sandbox(_COMSOL + f"queue:\n  - {U0}\n  - {U0}\n  - {U0}\n") as (tmp, cfg):
        fake = rb.evaluate_candidate = _fake_eval(["ok", "ok", "ok"])
        rc, out, err = _run(cfg)
        assert rc == 0 and out["n_done"] == 3
        # one session for all three: open before items 2 and 3, not before 1
        assert [c["session_before"] for c in fake.calls] == [False, True, True]
        assert not comsol_client.session_active()
        assert "COMSOL session released (queue empty)" in err


def test_session_released_early_when_rest_needs_no_comsol():
    q = (f"queue:\n  - {U0}\n  - {{u: {U0}, mech: surrogate_stub}}\n"
         f"  - {{u: {U0}, mech: surrogate_stub}}\n")
    with _sandbox(_COMSOL + q) as (tmp, cfg):
        fake = rb.evaluate_candidate = _fake_eval(["ok", "ok", "ok"])
        rc, out, err = _run(cfg)
        assert rc == 0 and out["n_done"] == 3
        assert [c["session_before"] for c in fake.calls] == [False, False, False]
        assert "no remaining item needs COMSOL" in err


def test_licence_error_stops_queue_and_releases():
    with _sandbox(_COMSOL + f"queue:\n  - {U0}\n  - {U0}\n  - {U0}\n") as (tmp, cfg):
        fake = rb.evaluate_candidate = _fake_eval(["ok", "licence", "ok"])
        rc, out, err = _run(cfg)
        assert rc == 3 and out["n_done"] == 2 and len(out["not_run"]) == 1
        assert len(fake.calls) == 2
        saved = database.load_completed_results(database.db_path("rib"))
        assert sorted(r["status"] for r in saved) == ["mech_failed", "success"]
        assert not comsol_client.session_active()
        assert "License error" in err            # for wait_for_comsol_seat.sh -r
    with _sandbox(_COMSOL + "stop_on_license_error: false\n"
                  + f"queue:\n  - {U0}\n  - {U0}\n  - {U0}\n") as (tmp, cfg):
        rb.evaluate_candidate = _fake_eval(["ok", "licence", "ok"])
        rc, out, err = _run(cfg)
        assert rc == 0 and out["n_done"] == 3


def test_ctrl_c_and_errors_release_the_session():
    with _sandbox(_COMSOL + f"queue:\n  - {U0}\n  - {U0}\n  - {U0}\n") as (tmp, cfg):
        rb.evaluate_candidate = _fake_eval(["ok", "interrupt", "ok"])
        rc, out, err = _run(cfg)
        assert rc == 130 and out["n_done"] == 1 and len(out["not_run"]) == 2
        assert not comsol_client.session_active()
        assert len(database.load_completed_results(database.db_path("rib"))) == 1
    with _sandbox(_COMSOL + f"queue:\n  - {U0}\n  - {U0}\n") as (tmp, cfg):
        rb.evaluate_candidate = _fake_eval(["ok", "boom"])
        client = None
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                rb.main(["--config", cfg])
        except RuntimeError as e:
            assert str(e) == "boom"
        else:
            raise AssertionError("an unexpected error must propagate")
        assert not comsol_client.session_active()
        assert len(database.load_completed_results(database.db_path("rib"))) == 1


if __name__ == "__main__":
    test_release_client_is_safe_and_final()
    test_queue_forms_and_validation()
    test_queue_saves_one_record_per_item_like_run_one()
    test_dry_run_solves_nothing()
    test_session_is_reused_then_released_when_queue_is_empty()
    test_session_released_early_when_rest_needs_no_comsol()
    test_licence_error_stops_queue_and_releases()
    test_ctrl_c_and_errors_release_the_session()
    print("all tests passed")
