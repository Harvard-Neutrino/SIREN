"""End-to-end test: place a gold nugget in ICARUS and find it via coordinate transforms.

Creates a minimal GDML with a 10cm gold cube at the ICARUS detector center,
places it in a composite world in BNB coordinates using the SBN geometry
transforms, then queries material density at the predicted gold location
using BNB coordinates, NuMI coordinates, and both DetectorPosition and
GeometryPosition interfaces.

This validates the full rotation/translation pipeline:
  sbn_geometry frame graph -> gdml_rotation_angles (Euler decomposition)
  -> GDML <rotation> (passive convention) -> SIREN GDML parser (QFromXYZs)
  -> DetectorModel coordinate transforms (active convention)
"""
from __future__ import annotations

import importlib.util
import math
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

from siren.detector import (
    DetectorModel, DetectorPosition, GeometryPosition,
)
from siren.math import Quaternion, Vector3D


@pytest.fixture(scope="module")
def geo():
    sbn_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", "resources", "detectors",
        "SBN", "SBN-v1")
    spec = importlib.util.spec_from_file_location(
        "sbn_geometry", os.path.join(sbn_dir, "sbn_geometry.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["sbn_geometry"] = mod
    spec.loader.exec_module(mod)
    return mod


# Gold at the center of ICARUS in ICARUS_LArSoft coordinates
GOLD_POS_ICARUS = np.array([0.0, -0.202, 0.0])
GOLD_SIZE_M = 0.1
GOLD_DENSITY = 19.3
AIR_DENSITY = 0.001225


def _make_detector_gdml(gold_pos):
    """Minimal GDML: a gold cube at gold_pos inside an air world."""
    return f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <isotope N="197" Z="79" name="Au197">
      <atom unit="g/mole" value="196.967"/>
    </isotope>
    <element name="Gold"><fraction n="1.0" ref="Au197"/></element>
    <material name="GoldMetal" state="solid">
      <D value="{GOLD_DENSITY}" unit="g/cm3"/>
      <fraction n="1.0" ref="Gold"/>
    </material>
    <isotope N="14" Z="7" name="N14">
      <atom unit="g/mole" value="14"/>
    </isotope>
    <element name="Nitrogen"><fraction n="1.0" ref="N14"/></element>
    <material name="TestAir" state="gas">
      <D value="{AIR_DENSITY}" unit="g/cm3"/>
      <fraction n="1.0" ref="Nitrogen"/>
    </material>
  </materials>
  <solids>
    <box name="world_box" lunit="m" x="20" y="20" z="20"/>
    <box name="gold_box" lunit="m"
         x="{GOLD_SIZE_M}" y="{GOLD_SIZE_M}" z="{GOLD_SIZE_M}"/>
  </solids>
  <structure>
    <volume name="vol_gold">
      <materialref ref="GoldMetal"/>
      <solidref ref="gold_box"/>
    </volume>
    <volume name="volDetWorld">
      <materialref ref="TestAir"/>
      <solidref ref="world_box"/>
      <physvol name="pv_gold">
        <volumeref ref="vol_gold"/>
        <position unit="m"
          x="{gold_pos[0]:.10f}" y="{gold_pos[1]:.10f}" z="{gold_pos[2]:.10f}"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="volDetWorld"/>
  </setup>
</gdml>
"""


def _make_composite_gdml(det_filename, origin, rotation_angles):
    """Composite GDML: place the detector file in a BNB-frame world."""
    rot_line = ""
    if rotation_angles is not None:
        rx, ry, rz = rotation_angles
        rot_line = (
            f'        <rotation unit="rad" '
            f'x="{rx:.10f}" y="{ry:.10f}" z="{rz:.10f}"/>')

    return f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <isotope N="14" Z="7" name="env_N14">
      <atom unit="g/mole" value="14"/>
    </isotope>
    <element name="env_Nitrogen">
      <fraction n="1.0" ref="env_N14"/>
    </element>
    <material name="EnvAir" state="gas">
      <D value="{AIR_DENSITY}" unit="g/cm3"/>
      <fraction n="1.0" ref="env_Nitrogen"/>
    </material>
  </materials>
  <solids>
    <box name="sol_world" lunit="m" x="1000" y="1000" z="2000"/>
  </solids>
  <structure>
    <volume name="volCompositeWorld">
      <materialref ref="EnvAir"/>
      <solidref ref="sol_world"/>
      <physvol name="pv_detector">
        <file name="{det_filename}" as_assembly="true"/>
        <position unit="m"
          x="{origin[0]:.10f}" y="{origin[1]:.10f}" z="{origin[2]:.10f}"/>
{rot_line}
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="volCompositeWorld"/>
  </setup>
</gdml>
"""


@pytest.fixture(scope="module")
def model_and_transforms(geo):
    """Build a composite GDML with gold in ICARUS, load it, return
    (DetectorModel, transform_dict)."""

    T_ic_bnb = geo.transform("ICARUS_LArSoft", "BNB")
    origin_bnb = T_ic_bnb.t

    # GDML physvol is passive: SIREN applies M^T.
    # For ICARUS->BNB (identity rotation), angles are zero.
    rx, ry, rz = geo.gdml_rotation_angles(T_ic_bnb.R.T)
    rotation_angles = None
    if abs(rx) > 1e-12 or abs(ry) > 1e-12 or abs(rz) > 1e-12:
        rotation_angles = (rx, ry, rz)

    tmpdir = tempfile.mkdtemp(prefix="siren_gold_test_")

    det_path = os.path.join(tmpdir, "detector.gdml")
    with open(det_path, "w") as f:
        f.write(_make_detector_gdml(GOLD_POS_ICARUS))

    comp_path = os.path.join(tmpdir, "composite.gdml")
    with open(comp_path, "w") as f:
        f.write(_make_composite_gdml("detector.gdml", origin_bnb, rotation_angles))

    dm = DetectorModel()
    dm.LoadGDML(comp_path)

    # DetectorOrigin is at the LAr center in BNB coords (matching detector.py).
    # DetectorRotation is the active rotation from det-local to BNB.
    det_info = geo.DETECTORS["ICARUS"]
    det_center_bnb = T_ic_bnb.apply(det_info.center_native)
    qx, qy, qz, qw = geo.quaternion_from_matrix(T_ic_bnb.R)
    dm.DetectorOrigin = GeometryPosition(
        Vector3D(det_center_bnb[0], det_center_bnb[1], det_center_bnb[2]))
    dm.DetectorRotation = Quaternion(qx, qy, qz, qw)

    transforms = {
        "T_ic_bnb": T_ic_bnb,
        "T_ic_numi": geo.transform("ICARUS_LArSoft", "NuMI"),
        "T_bnb_numi": geo.transform("BNB", "NuMI"),
        "T_numi_bnb": geo.transform("NuMI", "BNB"),
        "T_numi_ic": geo.transform("NuMI", "ICARUS_LArSoft"),
    }

    yield dm, transforms

    shutil.rmtree(tmpdir)


def _vec(v3):
    """Extract (x,y,z) tuple from a Vector3D or position wrapper."""
    if hasattr(v3, 'get'):
        v3 = v3.get()
    return np.array([v3.GetX(), v3.GetY(), v3.GetZ()])


# ======================================================================
# Tests
# ======================================================================

def _to_det_coords(pos_larsoft, geo, det_name="ICARUS"):
    """Convert LArSoft-frame position to DetectorCoordinates (offset by center)."""
    return pos_larsoft - geo.DETECTORS[det_name].center_native


class TestGoldInDetectorCoordinates:
    """Query the gold using DetectorCoordinates (centered on LAr volume)."""

    def test_gold_found(self, model_and_transforms, geo):
        dm, _ = model_and_transforms
        gold_det = _to_det_coords(GOLD_POS_ICARUS, geo)
        p = DetectorPosition(Vector3D(*gold_det))
        rho = dm.GetMassDensity(p)
        assert abs(rho - GOLD_DENSITY) < 0.5, f"Expected gold ({GOLD_DENSITY}), got {rho}"

    def test_air_nearby(self, model_and_transforms, geo):
        dm, _ = model_and_transforms
        gold_det = _to_det_coords(GOLD_POS_ICARUS, geo)
        off = gold_det + np.array([1.0, 0.0, 0.0])
        rho = dm.GetMassDensity(DetectorPosition(Vector3D(*off)))
        assert rho < 0.01, f"Expected air, got {rho}"


class TestGoldInBNBCoordinates:
    """Find the gold by computing its BNB position from the frame graph,
    then converting BNB (geometry) -> detector coords for the query."""

    def test_gold_found_via_bnb(self, model_and_transforms):
        dm, xforms = model_and_transforms
        gold_bnb = xforms["T_ic_bnb"].apply(GOLD_POS_ICARUS)

        p_geo = GeometryPosition(Vector3D(*gold_bnb))
        p_det = dm.GeoPositionToDetPosition(p_geo)
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5, f"Expected gold via BNB, got {rho}"

    def test_det_to_geo_matches_frame_graph(self, model_and_transforms, geo):
        """DetPositionToGeoPosition must map det coords to the correct BNB position."""
        dm, xforms = model_and_transforms
        gold_det = _to_det_coords(GOLD_POS_ICARUS, geo)
        gold_bnb_expected = xforms["T_ic_bnb"].apply(GOLD_POS_ICARUS)

        p_det = DetectorPosition(Vector3D(*gold_det))
        p_geo = dm.DetPositionToGeoPosition(p_det)
        gold_bnb_actual = _vec(p_geo)

        np.testing.assert_allclose(gold_bnb_actual, gold_bnb_expected, atol=1e-10)

    def test_geo_to_det_roundtrip(self, model_and_transforms, geo):
        dm, xforms = model_and_transforms
        gold_det = _to_det_coords(GOLD_POS_ICARUS, geo)
        gold_bnb = xforms["T_ic_bnb"].apply(GOLD_POS_ICARUS)

        p_geo = GeometryPosition(Vector3D(*gold_bnb))
        p_det = dm.GeoPositionToDetPosition(p_geo)
        gold_det_roundtrip = _vec(p_det)

        np.testing.assert_allclose(gold_det_roundtrip, gold_det, atol=1e-10)

    def test_air_in_bnb_coords(self, model_and_transforms):
        dm, xforms = model_and_transforms
        off_ic = GOLD_POS_ICARUS + np.array([1.0, 0.0, 0.0])
        off_bnb = xforms["T_ic_bnb"].apply(off_ic)

        p_geo = GeometryPosition(Vector3D(*off_bnb))
        p_det = dm.GeoPositionToDetPosition(p_geo)
        rho = dm.GetMassDensity(p_det)
        assert rho < 0.01, f"Expected air at offset BNB position, got {rho}"


class TestGoldInNuMICoordinates:
    """Find the gold by computing its NuMI position, transforming
    NuMI -> BNB -> detector coords."""

    def test_gold_found_via_numi(self, model_and_transforms):
        dm, xforms = model_and_transforms
        gold_numi = xforms["T_ic_numi"].apply(GOLD_POS_ICARUS)

        # NuMI -> BNB -> detector
        gold_bnb = xforms["T_numi_bnb"].apply(gold_numi)
        p_geo = GeometryPosition(Vector3D(*gold_bnb))
        p_det = dm.GeoPositionToDetPosition(p_geo)
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5, f"Expected gold via NuMI, got {rho}"

    def test_numi_position_consistency(self, model_and_transforms):
        """NuMI position computed directly vs via BNB must match."""
        _, xforms = model_and_transforms
        gold_numi_direct = xforms["T_ic_numi"].apply(GOLD_POS_ICARUS)

        gold_bnb = xforms["T_ic_bnb"].apply(GOLD_POS_ICARUS)
        gold_numi_via_bnb = xforms["T_bnb_numi"].apply(gold_bnb)

        np.testing.assert_allclose(gold_numi_direct, gold_numi_via_bnb, atol=1e-9)

    def test_numi_roundtrip(self, model_and_transforms):
        """ICARUS -> NuMI -> ICARUS must return the original position."""
        _, xforms = model_and_transforms
        gold_numi = xforms["T_ic_numi"].apply(GOLD_POS_ICARUS)
        gold_ic_roundtrip = xforms["T_numi_ic"].apply(gold_numi)
        np.testing.assert_allclose(gold_ic_roundtrip, GOLD_POS_ICARUS, atol=1e-9)

    def test_air_in_numi_coords(self, model_and_transforms):
        dm, xforms = model_and_transforms
        off_numi = xforms["T_ic_numi"].apply(GOLD_POS_ICARUS + np.array([1.0, 0.0, 0.0]))
        off_bnb = xforms["T_numi_bnb"].apply(off_numi)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*off_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert rho < 0.01, f"Expected air at offset NuMI position, got {rho}"


class TestTransformConsistency:
    """Verify the transforms themselves are internally consistent."""

    def test_icarus_bnb_is_pure_translation(self, model_and_transforms):
        _, xforms = model_and_transforms
        np.testing.assert_allclose(xforms["T_ic_bnb"].R, np.eye(3), atol=1e-15)
        np.testing.assert_allclose(xforms["T_ic_bnb"].t, [0, 0, 600], atol=1e-10)

    def test_numi_has_rotation(self, model_and_transforms):
        _, xforms = model_and_transforms
        assert not np.allclose(xforms["T_ic_numi"].R, np.eye(3), atol=1e-3)

    def test_detector_origin_at_lar_center(self, model_and_transforms, geo):
        """DetectorOrigin should be at the LAr center in BNB coords."""
        dm, xforms = model_and_transforms
        origin = _vec(dm.GetDetectorOrigin())
        det_info = geo.DETECTORS["ICARUS"]
        expected = xforms["T_ic_bnb"].apply(det_info.center_native)
        np.testing.assert_allclose(origin, expected, atol=1e-10)

    def test_compose_ic_bnb_numi_equals_ic_numi(self, model_and_transforms, geo):
        """ICARUS->BNB composed with BNB->NuMI must equal ICARUS->NuMI."""
        _, xforms = model_and_transforms
        T_composed = xforms["T_ic_bnb"].compose(xforms["T_bnb_numi"])
        T_direct = xforms["T_ic_numi"]
        np.testing.assert_allclose(T_composed.R, T_direct.R, atol=1e-9)
        np.testing.assert_allclose(T_composed.t, T_direct.t, atol=1e-9)
