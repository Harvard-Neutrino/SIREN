"""End-to-end test with real detector GDML: gold nugget in LAr.

Uses the actual SBN detector loader infrastructure (sbn_loader,
sbn_geometry, detector.py constants) to build a composite GDML with
the real detector geometry plus a 1cm gold cube, then verifies:
  - Gold found at the expected position via DetectorCoordinates
  - Gold found via BNB (GeometryCoordinates)
  - Gold found via NuMI coordinates
  - LAr (not air/rock/steel) found 2cm away in every direction
  - DetectorPosition(0,0,0) maps to the LAr volume center
  - All coordinate transforms are consistent across frames

Requires --run-network to download GDML files (~5 MB total).
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
import pytest

from siren.detector import DetectorModel, DetectorPosition, GeometryPosition
from siren.math import Quaternion, Vector3D, Matrix3D

pytestmark = pytest.mark.network

GOLD_DENSITY = 19.3
LAR_DENSITY = 1.39
GOLD_SIZE_M = 0.01

_SBN_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "resources", "detectors",
    "SBN", "SBN-v1")


def _load_sbn_module(name, filename):
    fqn = f"siren._sbn.{name}"
    path = os.path.join(_SBN_DIR, filename)
    spec = importlib.util.spec_from_file_location(fqn, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fqn] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def sbn():
    """Load the three SBN modules (geo, sbn_loader, detector) as a bundle."""
    module_names = [
        "siren._sbn.sbn_geometry",
        "siren._sbn.sbn_loader",
        "siren._sbn.sbn_detector",
    ]
    previous = {n: sys.modules.get(n) for n in module_names}
    geo = _load_sbn_module("sbn_geometry", "sbn_geometry.py")
    loader = _load_sbn_module("sbn_loader", "sbn_loader.py")
    det = _load_sbn_module("sbn_detector", "detector.py")
    yield geo, loader, det
    for n in module_names:
        sys.modules.pop(n, None)
        if previous[n] is not None:
            sys.modules[n] = previous[n]


def _gold_gdml():
    return f"""\
<?xml version="1.0"?>
<gdml><define/><materials>
<isotope N="197" Z="79" name="Au197"><atom unit="g/mole" value="196.967"/></isotope>
<element name="Gold"><fraction n="1.0" ref="Au197"/></element>
<material name="GoldMetal" state="solid">
  <D value="{GOLD_DENSITY}" unit="g/cm3"/><fraction n="1.0" ref="Gold"/>
</material>
</materials>
<solids><box name="gold_box" lunit="m" x="{GOLD_SIZE_M}" y="{GOLD_SIZE_M}" z="{GOLD_SIZE_M}"/></solids>
<structure>
<volume name="vol_gold"><materialref ref="GoldMetal"/><solidref ref="gold_box"/></volume>
</structure>
<setup name="Default" version="1.0"><world ref="vol_gold"/></setup>
</gdml>"""


def _build_model_with_gold(sbn, detector_name, gold_pos_larsoft, tmpdir):
    """Build a detector model with a gold nugget, using the same integration
    path as detector.py:load_detector but with an extra gold physvol.

    Uses the real _beamline_sources(), _DETECTOR_SPECS, and sbn_loader from
    the detector module, plus the same DetectorOrigin/Rotation logic.
    """
    geo, loader, det = sbn

    spec = det._DETECTOR_SPECS[detector_name]
    T_det_to_bnb = geo.detector_transform(detector_name, "BNB")
    origin_bnb = T_det_to_bnb.t
    det_quat = Quaternion()
    det_quat.SetMatrix(Matrix3D(*T_det_to_bnb.R.flatten()))

    det_rotation = None
    is_rotated = (abs(det_quat.X) > 1e-12
                  or abs(det_quat.Y) > 1e-12
                  or abs(det_quat.Z) > 1e-12)
    if is_rotated:
        rx, ry, rz = geo.gdml_rotation_angles(T_det_to_bnb.R.T)
        det_rotation = (rx, ry, rz)

    # Same sources as load_detector
    sources = list(det._beamline_sources())
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

    # Add a gold nugget physvol at the gold position in BNB coords
    gold_pos_bnb = T_det_to_bnb.apply(gold_pos_larsoft)
    gold_gdml_path = os.path.join(tmpdir, "gold.gdml")
    with open(gold_gdml_path, "w") as f:
        f.write(_gold_gdml())
    sources.append({
        "file": "gold.gdml",
        "prefix": "gold",
        "position": tuple(gold_pos_bnb),
        "rotation": None,
        "unwrap": False,
    })

    cache_name = f"composite_{detector_name.lower()}_gold.gdml"
    cache_path = loader.build_composite(tmpdir, sources, cache_name)

    dm = DetectorModel()
    dm.LoadGDML(cache_path)

    # Same DetectorOrigin/Rotation as load_detector
    det_info = geo.DETECTORS[detector_name]
    det_center_bnb = T_det_to_bnb.apply(det_info.center_native)
    dm.DetectorOrigin = GeometryPosition(
        Vector3D(det_center_bnb[0], det_center_bnb[1], det_center_bnb[2]))
    dm.DetectorRotation = det_quat

    return dm, det_info.center_native


def _vec(v):
    if hasattr(v, 'get'):
        v = v.get()
    return np.array([v.GetX(), v.GetY(), v.GetZ()])


# ======================================================================
# ICARUS
# ======================================================================

# ICARUS center is all LAr (verified by probing the real GDML)
ICARUS_GOLD_POS = np.array([0.0, -0.202, 0.0])


@pytest.fixture(scope="module")
def icarus_model(sbn, tmp_path_factory):
    tmpdir = str(tmp_path_factory.mktemp("icarus"))
    return _build_model_with_gold(sbn, "ICARUS", ICARUS_GOLD_POS, tmpdir)


class TestICARUSGold:

    def test_gold_in_det_coords(self, icarus_model):
        dm, center = icarus_model
        gold_det = ICARUS_GOLD_POS - center
        rho = dm.GetMassDensity(DetectorPosition(Vector3D(*gold_det)))
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_lar_2cm_away_all_directions(self, icarus_model):
        dm, center = icarus_model
        gold_det = ICARUS_GOLD_POS - center
        for axis in range(3):
            for sign in [-1, 1]:
                offset = np.zeros(3)
                offset[axis] = sign * 0.02
                rho = dm.GetMassDensity(
                    DetectorPosition(Vector3D(*(gold_det + offset))))
                label = "xyz"[axis]
                assert abs(rho - LAR_DENSITY) < 0.01, (
                    f"Expected LAr at {label}{'+' if sign > 0 else '-'}2cm, "
                    f"got {rho:.4f}")

    def test_gold_via_bnb(self, icarus_model, sbn):
        dm, _ = icarus_model
        geo = sbn[0]
        gold_bnb = geo.transform("ICARUS_LArSoft", "BNB").apply(ICARUS_GOLD_POS)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*gold_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_gold_via_numi(self, icarus_model, sbn):
        dm, _ = icarus_model
        geo = sbn[0]
        gold_numi = geo.transform("ICARUS_LArSoft", "NuMI").apply(ICARUS_GOLD_POS)
        gold_bnb = geo.transform("NuMI", "BNB").apply(gold_numi)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*gold_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_det_to_geo_matches_frame_graph(self, icarus_model, sbn):
        dm, center = icarus_model
        geo = sbn[0]
        gold_det = ICARUS_GOLD_POS - center
        expected_bnb = geo.transform("ICARUS_LArSoft", "BNB").apply(ICARUS_GOLD_POS)
        actual_bnb = _vec(dm.DetPositionToGeoPosition(
            DetectorPosition(Vector3D(*gold_det))))
        np.testing.assert_allclose(actual_bnb, expected_bnb, atol=1e-10)

    def test_detector_origin_is_in_lar(self, icarus_model):
        """A point near DetectorPosition(0,0,0) should be in LAr."""
        dm, _ = icarus_model
        for dx in [-0.1, 0.1]:
            rho = dm.GetMassDensity(DetectorPosition(Vector3D(dx, 0, 0)))
            assert abs(rho - LAR_DENSITY) < 0.01, (
                f"Expected LAr near detector origin (dx={dx}), got {rho:.4f}")


# ======================================================================
# SBND
# ======================================================================

# SBND LArSoft origin is at the cathode (x=0). Gold placed at x=0.5
# to be inside the drift volume, away from the cathode.
SBND_GOLD_POS = np.array([0.5, 0.0, 0.0])


@pytest.fixture(scope="module")
def sbnd_model(sbn, tmp_path_factory):
    tmpdir = str(tmp_path_factory.mktemp("sbnd"))
    return _build_model_with_gold(sbn, "SBND", SBND_GOLD_POS, tmpdir)


class TestSBNDGold:

    def test_gold_in_det_coords(self, sbnd_model):
        dm, center = sbnd_model
        gold_det = SBND_GOLD_POS - center
        rho = dm.GetMassDensity(DetectorPosition(Vector3D(*gold_det)))
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_lar_2cm_away_all_directions(self, sbnd_model):
        dm, center = sbnd_model
        gold_det = SBND_GOLD_POS - center
        for axis in range(3):
            for sign in [-1, 1]:
                offset = np.zeros(3)
                offset[axis] = sign * 0.02
                rho = dm.GetMassDensity(
                    DetectorPosition(Vector3D(*(gold_det + offset))))
                label = "xyz"[axis]
                assert abs(rho - LAR_DENSITY) < 0.01, (
                    f"Expected LAr at {label}{'+' if sign > 0 else '-'}2cm, "
                    f"got {rho:.4f}")

    def test_gold_via_bnb(self, sbnd_model, sbn):
        dm, _ = sbnd_model
        geo = sbn[0]
        gold_bnb = geo.transform("SBND_LArSoft", "BNB").apply(SBND_GOLD_POS)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*gold_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_gold_via_numi(self, sbnd_model, sbn):
        dm, _ = sbnd_model
        geo = sbn[0]
        gold_numi = geo.transform("SBND_LArSoft", "NuMI").apply(SBND_GOLD_POS)
        gold_bnb = geo.transform("NuMI", "BNB").apply(gold_numi)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*gold_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_det_to_geo_matches_frame_graph(self, sbnd_model, sbn):
        dm, center = sbnd_model
        geo = sbn[0]
        gold_det = SBND_GOLD_POS - center
        expected_bnb = geo.transform("SBND_LArSoft", "BNB").apply(SBND_GOLD_POS)
        actual_bnb = _vec(dm.DetPositionToGeoPosition(
            DetectorPosition(Vector3D(*gold_det))))
        np.testing.assert_allclose(actual_bnb, expected_bnb, atol=1e-10)

    def test_detector_origin_is_in_lar(self, sbnd_model):
        """Points near DetectorPosition(0,0,0) should be in LAr.

        The detector origin is at the LAr volume center. The LAr center
        sits at x=0 (the cathode) so we offset into one drift volume
        (x=+0.5m) and check that nearby points are all LAr.
        """
        dm, _ = sbnd_model
        base = np.array([0.5, 0.0, 0.0])
        for axis in range(3):
            for sign in [-1, 1]:
                offset = np.zeros(3)
                offset[axis] = sign * 0.1
                rho = dm.GetMassDensity(
                    DetectorPosition(Vector3D(*(base + offset))))
                label = "xyz"[axis]
                assert abs(rho - LAR_DENSITY) < 0.01, (
                    f"Expected LAr near detector center "
                    f"({label}{'+' if sign > 0 else '-'}10cm), "
                    f"got {rho:.4f}")
