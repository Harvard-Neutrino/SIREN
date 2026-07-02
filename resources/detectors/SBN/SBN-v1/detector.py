"""
SBN detector loader.

Composes BNB beamline + NuMI beamline + detector GDML files into a single
composite geometry in BNB-frame coordinates, then loads via the SIREN GDML
parser. Missing GDML files are downloaded automatically.

Usage:
    from siren._util import load_detector

    # Local-only (default): just the GDML site geology (~500 m around
    # detectors). Use for beam neutrino studies where all interactions
    # are near the beamline and detector.
    model = load_detector("SBN", detector="ICARUS")

    # Full Earth model (PREM + atmosphere + local Moho correction).
    # Use for atmospheric neutrinos or BSM searches with interactions
    # far from the detector.
    model = load_detector("SBN", detector="ICARUS", earth_model=True)
"""

import importlib.util
import os
import sys

_THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def _load_sibling(name, filename):
    fqn = f"siren._sbn.{name}"
    if fqn in sys.modules:
        return sys.modules[fqn]
    spec = importlib.util.spec_from_file_location(
        fqn, os.path.join(_THIS_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fqn] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        del sys.modules[fqn]
        raise
    return mod


geo = _load_sibling("sbn_geometry", "sbn_geometry.py")
sbn_loader = _load_sibling("sbn_loader", "sbn_loader.py")

from siren.download import writable_data_dir
_ABS_DIR = writable_data_dir(_THIS_DIR)

_DATA_BASE = (
    "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/"
    "main/detectors/SBN/v1"
)


def _beamline_sources():
    T_numi = geo.transform("NuMI", "BNB")
    numi_origin_bnb = T_numi.apply([0.0, 0.0, 0.0])
    numi_rx, numi_ry, numi_rz = geo.gdml_rotation_angles(T_numi.R.T)
    return [
        {
            "file": "gdml/BooNE_50m.gdml",
            "prefix": "bnb",
            "position": (0.0, 0.0, geo.SAND_CENTER_Z_BNB_M),
            "rotation": None,
            "unwrap": False,
            "url": f"{_DATA_BASE}/BNB/BooNE_50m.gdml",
            "sha256": "70d3a5ef55062b8bfc3c94cc7ddc559dce09dec03d670dddd5286518d80f12da",
        },
        {
            "file": "gdml/numi_g4export_2026-05-19.gdml",
            "prefix": "numi",
            "position": tuple(numi_origin_bnb),
            "rotation": (numi_rx, numi_ry, numi_rz),
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
    # MiniBooNE: spherical mineral-oil tank (inner oil + veto shell + steel
    # shell, in an air vault), generated locally by
    # sbn_loader.ensure_miniboone_gdml (no remote URL). Placed at the G4BNB
    # bsim::Location (0, 1.896, 541.34) m via the MiniBooNE_local frame.
    # unwrap=True: the MB_AIR bounding world box is dropped (as_assembly) so
    # the enclosure volumes sit directly in the site geology; keeping it as a
    # real sector would carve an air box out of the surrounding glacial till.
    "MiniBooNE": {
        "file": "gdml/miniboone_tank.gdml",
        "prefix": "miniboone",
        "unwrap": True,
        "url": None,
        "sha256": "",
    },
}


def fetch_data():
    """Download GDML files for all detectors (called by siren-download --fetch)."""
    # MiniBooNE's tank GDML has no remote URL; generate it locally first so
    # _ensure_gdml_files finds it present rather than failing to download.
    sbn_loader.ensure_miniboone_gdml(_ABS_DIR)
    all_sources = list(_beamline_sources())
    for spec in _DETECTOR_SPECS.values():
        if spec.get("file"):
            all_sources.append(spec)
    sbn_loader._ensure_gdml_files(_ABS_DIR, all_sources)


def load_detector(detector=None, earth_model=False):
    """Load an SBN detector model.

    Parameters
    ----------
    detector : str
        Which detector to load ("ICARUS" or "SBND").
    earth_model : bool, optional
        If *False* (default), only load the GDML site-geology volume
        (beamlines, detector, local stratigraphy within ~500 m).
        Use this for beam-neutrino studies where all interactions
        are near the beamline.

        If *True*, embed the GDML geometry in a full PREM Earth
        model with a multi-shell atmosphere and a local Moho
        correction for the Midwest crust.  Use this for atmospheric
        neutrino studies or BSM searches where interactions can occur
        far from the detector.

    Returns
    -------
    DetectorModel
    """
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

    # MiniBooNE's placeholder tank GDML is generated locally (no download).
    if detector == "MiniBooNE":
        sbn_loader.ensure_miniboone_gdml(_ABS_DIR)

    from siren.detector import DetectorModel, GeometryPosition
    from siren.math import Vector3D, Quaternion, Matrix3D

    T_det_to_bnb = geo.detector_transform(detector, "BNB")
    origin_bnb = T_det_to_bnb.t

    # DetectorRotation is an active rotation (r_geo = R @ r_det + origin),
    # so it uses the det-to-bnb rotation directly.
    R = T_det_to_bnb.R
    det_quat = Quaternion()
    det_quat.SetMatrix(Matrix3D(
        R[0, 0], R[0, 1], R[0, 2],
        R[1, 0], R[1, 1], R[1, 2],
        R[2, 0], R[2, 1], R[2, 2]))

    # GDML <physvol> rotation is passive (SIREN applies M^T), so we need
    # M = R_det_to_bnb^T (the inverse rotation). None when identity.
    det_rotation = None
    is_rotated = (abs(det_quat.X) > 1e-12
                  or abs(det_quat.Y) > 1e-12
                  or abs(det_quat.Z) > 1e-12)
    if is_rotated:
        rx, ry, rz = geo.gdml_rotation_angles(T_det_to_bnb.R.T)
        det_rotation = (rx, ry, rz)

    sources = list(_beamline_sources())
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

    model = DetectorModel()
    model.LoadGDML(cache_path)

    if earth_model:
        earth = _load_sibling("earth_model", "earth_model.py")
        earth.add_earth_model(model, sbn_loader._FNAL_SITE_GRADE_Y + geo.T_MiniBooNE_local[1])

    # DetectorOrigin is the point in BNB (geometry) coordinates that
    # corresponds to (0,0,0) in DetectorCoordinates. We place it at
    # the geometric center of the LAr volume (from the Detector entry
    # in sbn_geometry) so that injectors aiming at the detector origin
    # hit the middle of the active volume, not e.g. the cathode plane.
    #
    # DetectorRotation is the active rotation from detector-local axes
    # to geometry (BNB) axes: r_geo = R @ r_det + DetectorOrigin.
    det_info = geo.DETECTORS[detector]
    det_center_bnb = T_det_to_bnb.apply(det_info.center_native)
    model.DetectorOrigin = GeometryPosition(
        Vector3D(det_center_bnb[0], det_center_bnb[1], det_center_bnb[2]))
    model.DetectorRotation = det_quat

    return model
