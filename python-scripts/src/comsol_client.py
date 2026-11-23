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
