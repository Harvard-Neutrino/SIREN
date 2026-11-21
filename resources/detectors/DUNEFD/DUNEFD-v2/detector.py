"""
DUNE far-detector loader.

Loads the official dunecore GDML (FD1-HD horizontal drift and/or FD2-VD vertical
drift, nowires variants) into a single SIREN coordinate frame, optionally with the
Homestake overburden (real topography, geological formations, PREM Earth).

SITE FRAME (shared by all modules and the overburden):
    +y = local up (vertical at SURF)
    +z = LBNF neutrino-beam direction (downstream, traveling away from Fermilab);
         horizontal projection, geographic azimuth ~277.15 deg at the far detector
    +x = right-handed (horizontal transverse, ~ geographic south)
The real beam travels +z tilted UP by ~5.71 deg (it dives under from FNAL and
arrives going upward at the FD). This is the LArSoft DUNE convention
(z = neutrino travel, y = up, x = right-handed).

The HD GDML is already in this frame (y=up, z=beam). The VD GDML is built x=up
(vertical drift), so it is rotated x->y (GDML rotation z=-90) to match.

Usage:
    from siren._util import load_detector
    model = load_detector("DUNEFD", detector="VD")                 # VD only
    model = load_detector("DUNEFD", detector="HD")                 # HD only
    model = load_detector("DUNEFD", detector="both")               # FD1-HD + FD2-VD
    model = load_detector("DUNEFD", detector="both", earth_model=True)  # + overburden

The module GDMLs are the official DUNE/dunecore far-detector geometries, mirrored
(sha256-verified) in SIREN-data under detectors/DUNEFD/v2/{HD,VD}/ with a README
recording the exact upstream commit. They are fetched on first use, or via
``siren-download --fetch DUNEFD``.
"""
import importlib.util
import os
import sys

from siren.download import writable_data_dir, ensure_files, resolve_data_path

_THIS_DIR = os.path.dirname(os.path.realpath(__file__))
# Module GDMLs are mirrored in SIREN-data (each with a provenance README) and
# fetched, sha256-verified, on first use. _INSTALL_GDML holds any copy shipped
# with the package; _GDML_DIR is the writable download/cache location and is
# also where the composite stub is written.
_ABS_DIR = writable_data_dir(_THIS_DIR)
_INSTALL_GDML = os.path.join(_THIS_DIR, "gdml")
_GDML_DIR = os.path.join(_ABS_DIR, "gdml")

# SIREN-data mirror of the official DUNE/dunecore far-detector GDMLs. Each module
# lives under detectors/DUNEFD/v2/<HD|VD>/ alongside a README documenting the
# exact upstream dunecore commit it was taken from.
_DATA_BASE = (
    "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/"
    "main/detectors/DUNEFD/v2"
)

# ---- site-frame / beam geometry (computed from the FNAL->SURF chord) ----
BEAM_AZIMUTH_DEG = 277.15   # +z: beam downstream azimuth at the FD (from N toward E)
BEAM_DIP_DEG = 5.71         # beam travels UP by this much at the FD (info only)
# HD (north cavern) <-> VD (south cavern) transverse center-to-center spacing.
# The N-S spacing is not published; ~75 m from cavern widths 19.8 m + CUC 19.3 m
# + rock pillars. Transverse = +x in the site frame.
CAVERN_SPACING_M = 75.0


def _load_sibling(name, filename):
    fqn = f"siren._dunefd.{name}"
    if fqn in sys.modules:
        return sys.modules[fqn]
    spec = importlib.util.spec_from_file_location(fqn, os.path.join(_THIS_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fqn] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        del sys.modules[fqn]
        raise
    return mod


# Per-module spec. `active_center` is the active-LAr center in the module's OWN
# GDML frame; `rotation_deg` is the GDML <rotation> (passive) that brings the
# module into the site frame; HD is identity, VD is x->up rotated to y->up.
_MODULES = {
    "HD": {
        "name": "dune10kt_v6_refactored_1x2x6_nowires",
        "url": f"{_DATA_BASE}/HD/dune10kt_v6_refactored_1x2x6_nowires.gdml",
        "sha256": "b801f3eb28c5af58205cbe393039fc6f414e384fe6ee66a81ebe8c3b578b7e68",
        "active_center": (0.0, 0.21, 5.70),    # HD already y=up, z=beam
        "rotation_deg": (0.0, 0.0, 0.0),
    },
    "VD": {
        "name": "dunevd10kt_3view_30deg_v7_refactored_1x8x14_nowires",
        "url": f"{_DATA_BASE}/VD/dunevd10kt_3view_30deg_v7_refactored_1x8x14_nowires.gdml",
        "sha256": "24fffafeeb17ceea3760ab378bf4ed76803a072cc9060704d7810ea9d1d85b56",
        "active_center": (0.0, 0.0, 10.5),     # VD x=up -> rotate to y=up
        "rotation_deg": (0.0, 0.0, -90.0),
    },
}


def _ensure_gdml(spec):
    """Resolve a module GDML, fetching it from SIREN-data if absent.

    Prefers a copy shipped in the package's gdml/ directory; otherwise downloads
    the sha256-verified file from SIREN-data into a writable mirror. The files
    originate from DUNE/dunecore -- see the per-module README under
    detectors/DUNEFD/v2/<HD|VD>/ in SIREN-data for the exact upstream commit.
    Returns the absolute path to the resolved file.
    """
    fn = spec["name"] + ".gdml"
    path = resolve_data_path(_INSTALL_GDML, _GDML_DIR, fn)
    if not os.path.isfile(path):
        ensure_files([{"path": path, "url": spec["url"], "sha256": spec["sha256"]}])
    return path


def fetch_data():
    """Pre-fetch the FD module GDMLs from SIREN-data.

    Entry point for ``siren-download --fetch DUNEFD``; downloads both modules so
    later ``load_detector`` calls run offline.
    """
    for spec in _MODULES.values():
        _ensure_gdml(spec)


def _write_composite(placements, out_path, world_half=400.0):
    """Write a stub GDML referencing each module's GDML with a position+rotation.

    placements: list of (gdml_path, (rx,ry,rz) deg, (tx,ty,tz) m). gdml_path is
    an absolute path so the composite need not be co-located with the modules
    (they may be resolved from the package dir or a writable download mirror).
    """
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    physvols = []
    for i, (fn, rot, off) in enumerate(placements):
        # as_assembly="true" unwraps the sub-geometry into this frame so the
        # position/rotation propagates to containment (not just the bounding box).
        pv = [f'      <physvol name="pv_module_{i}">',
              f'        <file name="{fn}" as_assembly="true"/>',
              f'        <position unit="m" x="{off[0]:.6f}" y="{off[1]:.6f}" z="{off[2]:.6f}"/>']
        if any(abs(a) > 1e-9 for a in rot):
            pv.append(f'        <rotation unit="deg" x="{rot[0]:.6f}" y="{rot[1]:.6f}" z="{rot[2]:.6f}"/>')
        pv.append("      </physvol>")
        physvols.append("\n".join(pv))
    stub = (
        '<?xml version="1.0"?>\n<gdml>\n  <define/>\n  <materials>\n'
        '    <material name="WorldVacuum" Z="1"><D value="1e-25" unit="g/cm3"/>'
        '<atom value="1.008"/></material>\n  </materials>\n  <solids>\n'
        f'    <box name="sol_world" lunit="m" x="{2*world_half}" y="{2*world_half}" '
        f'z="{2*world_half}"/>\n  </solids>\n  <structure>\n'
        '    <volume name="World">\n      <materialref ref="WorldVacuum"/>\n'
        '      <solidref ref="sol_world"/>\n'
        + "\n".join(physvols)
        + '\n    </volume>\n  </structure>\n'
        '  <setup name="Default" version="1.0"><world ref="World"/></setup>\n</gdml>\n')
    with open(out_path, "w") as f:
        f.write(stub)
    return out_path


def load_detector(detector=None, earth_model=False):
    """Load a DUNE far-detector model in the site frame (+y up, +z beam).

    detector : "HD", "VD", or "both" (FD1-HD + FD2-VD).
    earth_model : if True, overlay the Homestake overburden (real topography,
        formations, PREM). DEM tiles download on first use. Default False.
    """
    if detector is None:
        raise TypeError('"detector" is required: "HD", "VD", or "both".')
    key = str(detector).lower()
    if key in ("hd", "vd"):
        modules = [key.upper()]
    elif key == "both":
        modules = ["HD", "VD"]
    else:
        raise ValueError(f'Unknown detector "{detector}". Choose "HD", "VD", or "both".')

    # transverse offsets: single module at origin; for "both", HD at -spacing/2,
    # VD at +spacing/2 along +x (north cavern / south cavern).
    offsets = {m: (0.0, 0.0, 0.0) for m in modules}
    if modules == ["HD", "VD"]:
        offsets["HD"] = (-CAVERN_SPACING_M / 2.0, 0.0, 0.0)
        offsets["VD"] = (+CAVERN_SPACING_M / 2.0, 0.0, 0.0)

    placements = []
    for m in modules:
        spec = _MODULES[m]
        gdml = _ensure_gdml(spec)
        placements.append((gdml, spec["rotation_deg"], offsets[m]))

    from siren.detector import DetectorModel, GeometryPosition
    from siren.math import Vector3D

    composite = os.path.join(_GDML_DIR, f"composite_{key}.gdml")
    _write_composite(placements, composite)
    model = DetectorModel()
    model.LoadGDML(composite)

    if earth_model:
        earth = _load_sibling("earth_model", "earth_model.py")
        earth.add_earth_model(model)

    # DetectorOrigin = active-LAr center of the primary module (HD for "both"),
    # mapped into the site frame (rotation of HD/VD primary is about the long
    # axis, so the on-axis active center only shifts by the transverse offset).
    primary = modules[0]
    ac = _MODULES[primary]["active_center"]
    off = offsets[primary]
    model.DetectorOrigin = GeometryPosition(Vector3D(ac[0] + off[0], ac[1] + off[1], ac[2] + off[2]))
    return model
