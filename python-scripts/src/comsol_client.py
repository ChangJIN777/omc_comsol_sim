"""Shared COMSOL/MPh client singleton.

Both acoustic_comsol and optical_comsol import from here so the whole Python
process shares ONE COMSOL instance and ONE loaded model.
"""
from __future__ import annotations

import os

try:
    import mph
    _HAVE_MPH = True
except Exception:
    _HAVE_MPH = False

_client = None
_model = None
_model_path = None
_released = False


def have_mph() -> bool:
    return _HAVE_MPH


def get_client():
    global _client
    if not _HAVE_MPH:
        raise RuntimeError("MPh not installed. Run: pip install MPh")
    if _client is None:
        if _released:
            raise RuntimeError(
                "the COMSOL session of this process was released "
                "(comsol_client.release_client); start a new Python process "
                "for more COMSOL work")
        _client = mph.start()
    return _client


def session_active() -> bool:
    """True while this process holds a COMSOL session (and so licences)."""
    return _client is not None


def get_model(template: str):
    """Return a loaded mph.Model, reloading only when the path changes."""
    global _model, _model_path
    client = get_client()
    if _model is None or _model_path != template:
        if _model is not None:
            try:
                client.remove(_model)
            except Exception:
                pass
        _model = client.load(template)
        _model_path = template
    return _model


def release_model():
    """Unload the current model (e.g. to force a reload after template rebuild)."""
    global _model, _model_path
    if _model is not None and _client is not None:
        try:
            _client.remove(_model)
        except Exception:
            pass
    _model = None
    _model_path = None


def release_client() -> bool:
    """End the COMSOL session of this process so its licences are returned.

    Unloads every model, disconnects the client and stops the local COMSOL
    server that mph.start() launched. Call it when no more COMSOL work is
    queued (scripts/run_batch_rib.py does, in a finally block): MPh would
    otherwise keep the session, and the seats, until the interpreter exits.

    Returns True if a session was open. Never raises, and calling it without
    a session (or twice) does nothing. A released session is not restarted:
    get_client() raises afterwards.

    With MPh's stand-alone client (mph.option('session', 'stand-alone')) there
    is no server to stop; the module seats go back when the models are
    cleared, the COMSOL seat when the process exits.
    """
    global _client, _released
    release_model()
    if _client is None:
        return False
    client, _client = _client, None
    _released = True
    try:
        client.clear()
    except Exception:
        pass
    try:
        if getattr(client, "port", None):
            client.disconnect()
    except Exception:
        pass
    try:
        import mph.session as _session
        server = getattr(_session, "server", None)
        if server is not None and server.running():
            server.stop()
    except Exception:
        pass
    return True


# One template per cell type (see geometry.CELLTYPES).
_COMSOL_DIR = os.path.join(os.path.dirname(__file__), "..", "comsol")
TEMPLATES = {
    "hole": "omc_unitcell_iso.mph",   # isotropic diamond, elliptical hole
    "rib": "omc_unitcell_rib.mph",    # cubic diamond, build_omc_unitcell_rib.m
}


def template_path(celltype: str = "hole") -> str:
    """Path of the .mph template for a cell type."""
    if celltype not in TEMPLATES:
        raise ValueError(f"unknown celltype {celltype!r}; expected one of "
                         f"{sorted(TEMPLATES)}")
    return os.path.join(_COMSOL_DIR, TEMPLATES[celltype])


def set_geometry(model, g) -> None:
    """Set the geometry parameters of `g` (meters) on a loaded template.

    The names in g.as_dict() must all be global parameters of the model
    already. COMSOL would otherwise create the missing ones silently, and a
    hole geometry sent to the rib template (or the reverse) would solve with
    the wrong meaning of `t` (hole: thickness, rib: rib width).
    """
    params = g.as_dict()
    have = set(model.parameters())
    missing = sorted(name for name in params if name not in have)
    if missing:
        raise RuntimeError(
            f"template has no parameter(s) {missing}: "
            f"{getattr(g, 'celltype', '?')!r} geometry sent to the wrong "
            f"template? Template parameters: {sorted(have)}")
    for name, val in params.items():
        model.parameter(name, f"{val}[m]")
