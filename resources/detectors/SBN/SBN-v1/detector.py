"""
SBN detector loader.

Composes BNB beamline + NuMI beamline + detector GDML files into a single
composite geometry in BNB-frame coordinates, then loads via the SIREN GDML
parser. Missing GDML files are downloaded automatically.

Usage:
    from siren._util import load_detector
    model = load_detector("SBN", detector="ICARUS")
    model = load_detector("SBN", detector="SBND")
"""

import importlib.util
import os

_THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def _load_sibling(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(_THIS_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


geo = _load_sibling("sbn_geometry", "sbn_geometry.py")
sbn_loader = _load_sibling("sbn_loader", "sbn_loader.py")

from siren.download import writable_data_dir
_ABS_DIR = writable_data_dir(_THIS_DIR)

_T_numi = geo.transform("NuMI", "BNB")
_numi_origin_bnb = _T_numi.apply([0.0, 0.0, 0.0])
_numi_rx, _numi_ry, _numi_rz = geo.gdml_rotation_angles(_T_numi.R)

_DATA_BASE = (
    "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/"
    "main/detectors/SBN/v1"
)

_BEAMLINE_SOURCES = [
    {
        "file": "gdml/BooNE_50m.gdml",
        "prefix": "bnb",
        "position": (0.0, 0.0, -geo.BNB_TARGET_Z_SAND_M),
        "rotation": None,
        "unwrap": False,
        "url": f"{_DATA_BASE}/BNB/BooNE_50m.gdml",
        "sha256": "70d3a5ef55062b8bfc3c94cc7ddc559dce09dec03d670dddd5286518d80f12da",
    },
    {
        "file": "gdml/numi_g4export_2026-05-19.gdml",
        "prefix": "numi",
        "position": tuple(_numi_origin_bnb),
        "rotation": (_numi_rx, _numi_ry, _numi_rz),
        "unwrap": False,
        "url": f"{_DATA_BASE}/NuMI/numi_g4export_2026-05-19.gdml",
        "sha256": "39670d52a6181352a8ae7c798387a9c58de950462c634e57da7d39fb23abe30a",
    },
]

_DETECTOR_SPECS = {
    "ICARUS": {
        "file": "gdml/icarus_refactored_nounderscore_20230918_nowires.gdml",
        "prefix": "icarus",
        "unwrap": True,
        "url": f"{_DATA_BASE}/ICARUS/icarus_refactored_nounderscore_20230918_nowires.gdml",
        "sha256": "3f68961a21fa037d3278c3235e48920cfc85a306adc5437f296dd7c09f095cd1",
    },
    "SBND": {
        "file": "gdml/sbnd_v02_06.gdml",
        "prefix": "sbnd",
        "unwrap": True,
        "url": f"{_DATA_BASE}/SBND/sbnd_v02_06.gdml",
        "sha256": "224a0efa55e66b1fb3b527937e4a899854c2bbab367db3659e55466c4e7f013a",
    },
}


def fetch_data():
    """Download GDML files for all detectors (called by siren-download --fetch)."""
    all_sources = list(_BEAMLINE_SOURCES)
    for spec in _DETECTOR_SPECS.values():
        if spec.get("file"):
            all_sources.append(spec)
    sbn_loader._ensure_gdml_files(_ABS_DIR, all_sources)


def load_detector(detector=None):
    if detector is None:
        raise TypeError(
            '"detector" is a required argument. '
            'Choose from: ' + ", ".join(_DETECTOR_SPECS.keys()))

    detector = str(detector)
    if detector not in _DETECTOR_SPECS:
        raise ValueError(
            f'Unknown detector "{detector}". '
            f'Choose from: {", ".join(_DETECTOR_SPECS.keys())}')

    spec = _DETECTOR_SPECS[detector]

    # Full rigid transform from the detector's native (LArSoft) frame to
    # BNB frame: r_BNB = R @ r_det + t.  Used for both placing the detector
    # GDML in the composite and for telling SIREN how to convert between
    # detector-local and geometry-global coordinates.
    T_det_to_bnb = geo.detector_transform(detector, "BNB")
    origin_bnb = T_det_to_bnb.t
    qx, qy, qz, qw = geo.quaternion_from_matrix(T_det_to_bnb.R)

    # GDML Euler angles for the <physvol> placement rotation.  None when
    # the detector frame is axis-aligned with BNB (pure translation).
    det_rotation = None
    if abs(qx) > 1e-12 or abs(qy) > 1e-12 or abs(qz) > 1e-12:
        rx, ry, rz = geo.gdml_rotation_angles(T_det_to_bnb.R)
        det_rotation = (rx, ry, rz)

    sources = list(_BEAMLINE_SOURCES)
    if spec["file"] is not None:
        sources.append({
            "file": spec["file"],
            "prefix": spec["prefix"],
            "position": tuple(origin_bnb),
            "rotation": det_rotation,
            "unwrap": spec["unwrap"],
            "url": spec.get("url"),
            "sha256": spec.get("sha256", ""),
        })

    cache_name = f"composite_{detector.lower()}.gdml"
    cache_path = sbn_loader.build_composite(_ABS_DIR, sources, cache_name)

    from siren.detector import DetectorModel
    from siren.math import Vector3D, Quaternion
    from siren.detector import GeometryPosition

    model = DetectorModel()
    model.LoadGDML(cache_path)

    # DetectorOrigin + DetectorRotation define the rigid transform between
    # detector-local coordinates and geometry (BNB frame) coordinates:
    #   ToGeo(det_pos) = R @ det_pos + origin
    #   ToDet(geo_pos) = R^T @ (geo_pos - origin)
    # Without these, SIREN assumes the detector is at the BNB target
    # with axes aligned to BNB.
    model.DetectorOrigin = GeometryPosition(
        Vector3D(origin_bnb[0], origin_bnb[1], origin_bnb[2]))
    model.DetectorRotation = Quaternion(qx, qy, qz, qw)

    return model
