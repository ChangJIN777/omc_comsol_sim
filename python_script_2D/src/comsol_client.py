"""Shared COMSOL/MPh client singleton.

Ported unchanged from python-scripts/src/comsol_client.py -- acoustic_comsol_2d
(and, later, optical_comsol_2d) import from here so the whole Python process
shares ONE COMSOL instance and ONE loaded model.
"""
from __future__ import annotations

try:
    import mph
    _HAVE_MPH = True
except Exception:
    _HAVE_MPH = False

_client = None
_model = None
_model_path = None


def have_mph() -> bool:
    return _HAVE_MPH


def get_client():
    global _client
    if not _HAVE_MPH:
        raise RuntimeError("MPh not installed. Run: pip install MPh")
    if _client is None:
        _client = mph.start()
    return _client


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
