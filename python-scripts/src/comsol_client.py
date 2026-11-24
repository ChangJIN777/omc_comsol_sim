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


# One template per cell type (geometry.CELLTYPES) and cross-section
# (geometry.XSECTS). Checked against the .mph contents (dmodel.xml):
#   omc_unitcell.mph      Block(t, w, a) minus the scaled hole cylinder: a
#                         clean RECTANGULAR beam.
#   omc_unitcell_iso.mph  the same plus Polygon 1 / Extrude 1, the triangular
#                         WEDGE under the slab. ("iso" = isotropic elasticity;
#                         omc_unitcell.mph is isotropic too.)
#   omc_unitcell_rib.mph        build_omc_unitcell_rib()        (rectangular)
#   omc_unitcell_rib_wedge.mph  build_omc_unitcell_rib('wedge') (see BUILD)
_COMSOL_DIR = os.path.join(os.path.dirname(__file__), "..", "comsol")
TEMPLATES = {
    ("hole", "rect"): "omc_unitcell.mph",
    ("hole", "wedge"): "omc_unitcell_iso.mph",
    ("rib", "rect"): "omc_unitcell_rib.mph",
    ("rib", "wedge"): "omc_unitcell_rib_wedge.mph",
}
# The hole loop scripts/run_opt_comsol.py has always used the wedge template.
LEGACY_XSECT = {"hole": "wedge", "rib": "rect"}
# How to get a missing template.
BUILD = {
    ("hole", "rect"): "build it in the COMSOL GUI (comsol/README_template.md)",
    ("hole", "wedge"): "it is the rectangular template plus the wedge polygon "
                       "(comsol/README_template.md)",
    ("rib", "rect"): "in MATLAB with LiveLink: cd python-scripts/comsol; "
                     "build_omc_unitcell_rib();",
    ("rib", "wedge"): "in MATLAB with LiveLink: cd python-scripts/comsol; "
                      "build_omc_unitcell_rib('wedge');",
}
# Name of the slab-thickness parameter in the templates of each cell type.
THICKNESS_PARAM = {"hole": "t", "rib": "th"}


def template_path(celltype: str = "hole", xsect: str = "rect",
                  must_exist: bool = False) -> str:
    """Path of the .mph template for a cell type and cross-section.

    must_exist=True raises FileNotFoundError, with the build step, when the
    template has not been built (the rib wedge template is not shipped).
    """
    key = (celltype, "rect" if xsect is None else xsect)
    if key not in TEMPLATES:
        raise ValueError(f"no template for celltype {celltype!r}, xsect "
                         f"{xsect!r}; expected one of {sorted(TEMPLATES)}")
    path = os.path.join(_COMSOL_DIR, TEMPLATES[key])
    if must_exist and not os.path.exists(path):
        raise FileNotFoundError(
            f"COMSOL template for celltype {key[0]!r}, xsect {key[1]!r} not "
            f"found: {path}\nTo get it: {BUILD[key]}")
    if must_exist and key == ("rib", "wedge") and not rib_wedge_is_roof(path):
        raise RuntimeError(
            f"{path} was built with the first wedge definition (straight "
            "prisms with a flat end, parameter wd = max(s,t)/2*tan(wedge_angle)). "
            "The wedge is now the roof of the filleted outline "
            "(geometry.rib_roof_distance). Delete the file and rebuild it: "
            f"{BUILD[key]}")
    return path


_mentions_checked = {}


def _template_mentions(path: str, word: str) -> bool:
    """Does the model description inside this .mph mention `word`?

    Read from the .mph (a zip: dmodel.xml / smodel.json) without COMSOL, to
    tell which builder version made a template from the global parameters it
    defines. Cached by (path, mtime, word). A file that cannot be read, or
    has neither entry, counts as True (current), so an unusual .mph layout
    never blocks a run.
    """
    try:
        st = os.stat(path)       # size too: a rewrite within one mtime tick
        stamp = (path, st.st_mtime_ns, st.st_size, word)
    except OSError:
        return True
    if stamp not in _mentions_checked:
        import zipfile
        try:
            with zipfile.ZipFile(path) as z:
                names = [n for n in ("dmodel.xml", "smodel.json")
                         if n in z.namelist()]
                text = "".join(z.read(n).decode("utf-8", "replace") for n in names)
            _mentions_checked[stamp] = (not names) or (word in text)
        except Exception:
            _mentions_checked[stamp] = True
    return _mentions_checked[stamp]


def rib_wedge_is_roof(path: str) -> bool:
    """Was this rib wedge template built with the roof construction?

    The roof builder defines the global parameter `wdj` (wedge depth under
    the junction); the first, straight-prism builder did not.
    """
    return _template_mentions(path, "wdj")


def rib_air_is_disk(path: str) -> bool:
    """Was this rib template built with the air half-disk?

    The half-disk builder defines the global parameter `air_r` (and air_gap,
    air_xc); the templates built before it have an air box of 5*th by 1.5*w
    and no air_r.
    """
    return _template_mentions(path, "air_r")


def require_rib_air_disk(path: str) -> str:
    """Return `path`, or raise if this rib template still has the air box.

    For the COMSOL optical solves only: the old box put the scattering
    boundary inside the evanescent tail of the guided modes. Mechanics does
    not use the air, so it does not call this.
    """
    if not rib_air_is_disk(path):
        raise RuntimeError(
            f"{path} was built with the old air box (no global parameter "
            "air_r): its outer faces sit inside the evanescent tail of the "
            "guided modes. The rib templates now have an air half-disk of "
            "radius air_r = farthest beam corner + air_gap. Delete the file "
            "and rebuild it in MATLAB with LiveLink: "
            "build_omc_unitcell_rib() (rect) or build_omc_unitcell_rib('wedge'). "
            "Mechanical solves are unaffected.")
    return path


# Crystal rotation of the rib templates (rib.rxtal_deg; mechanics only). The
# builder defines the global parameter `rxt` and the derived dCrot / gCrot;
# templates built before it have the unrotated tensor and no rxt.
RXTAL_PARAM = "rxt"
_RXTAL_MARKER = "dCrot"


def template_has_rxtal(path: str) -> bool:
    """Was this rib template built with the crystal rotation `rxt`?"""
    return _template_mentions(path, _RXTAL_MARKER)


def require_rxtal_template(path: str, rxtal_deg=0.0) -> str:
    """Return `path`, or raise if a nonzero rotation needs a rebuilt template."""
    from geometry import resolve_rxtal
    r = resolve_rxtal(rxtal_deg)
    if r != 0.0 and not template_has_rxtal(path):
        raise RuntimeError(
            f"rxtal_deg = {r:g} needs a rib template with the crystal "
            f"rotation parameter {RXTAL_PARAM!r}, but {path} was built "
            "before it (unrotated tensor, beam along <100>). Delete the file "
            "and rebuild it in MATLAB with LiveLink: build_omc_unitcell_rib() "
            "(rect) or build_omc_unitcell_rib('wedge').")
    return path


def set_crystal_rotation(model, rxtal_deg=0.0) -> float:
    """Set the crystal rotation of a loaded rib template; returns it in deg.

    The rotation is reduced by geometry.resolve_rxtal (0 or 45 only). A
    template with the global parameter `rxt` gets it set on every call (the
    model is a shared singleton, so a previous rotation must not stay). A
    template without it (hole templates, rib templates built before rxt)
    has the unrotated tensor: 0 is accepted, anything else raises.
    """
    from geometry import resolve_rxtal
    r = resolve_rxtal(rxtal_deg)
    if RXTAL_PARAM in set(model.parameters()):
        model.parameter(RXTAL_PARAM, f"{r:g}[deg]")
    elif r != 0.0:
        raise RuntimeError(
            f"rxtal_deg = {r:g} but the loaded template has no parameter "
            f"{RXTAL_PARAM!r} (hole template, or a rib template built before "
            "the crystal rotation): refusing to solve the unrotated crystal. "
            "Rebuild the rib template with build_omc_unitcell_rib.")
    return r


def set_diamond_index(model, n=None) -> float:
    """Set the diamond refractive index of a loaded template; returns n.

    The templates store relpermittivity 5.76 (n = 2.40) and the hole ones are
    binary files that are not rebuilt, so the value of configs/materials.yaml
    (2.386) is written at run time, before every optical solve. Mechanics does
    not use it. A template with a global parameter `n_diamond` (rib templates
    built after 2026-11-23, hole templates patched with
    comsol/patch_hole_templates_index.m, both with relpermittivity =
    n_diamond^2) gets the parameter set; any other gets the relpermittivity
    of its diamond material (label starting with "Diamond", else mat1) set
    to n^2, so unpatched templates keep working. Raises if neither works: solving with another index
    than the one recorded must not happen silently.
    """
    if n is None:
        from geometry import diamond_index
        n = diamond_index()
    n = float(n)
    if "n_diamond" in set(model.parameters()):
        model.parameter("n_diamond", repr(n))
        return n
    try:
        comp = model.java.component("comp1")
        tags = [str(t) for t in comp.material().tags()]
        diamond = [t for t in tags
                   if str(comp.material(t).label()).lower().startswith("diamond")]
        tag = diamond[0] if diamond else "mat1"
        eps, zero = repr(n * n), "0"
        tensor = [eps, zero, zero, zero, eps, zero, zero, zero, eps]
        try:                            # as MPh does: an explicit String[]
            from jpype import JArray, JString
            tensor = JArray(JString)(tensor)
        except Exception:               # no jpype / no JVM: plain sequence
            pass
        comp.material(tag).propertyGroup("def").set("relpermittivity", tensor)
    except Exception as e:
        raise RuntimeError(
            f"could not set the diamond refractive index {n} on the COMSOL "
            f"template ({type(e).__name__}: {e}); refusing to solve optics "
            "with the template's own value") from e
    return n


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
