"""End-to-end test with real detector GDML: gold nugget in LAr.

Downloads the real ICARUS and SBND GDML files, places a 1cm gold cube
at a known position inside the active LAr volume, then verifies:
  - Gold found at the expected position via DetectorCoordinates
  - Gold found via BNB (GeometryCoordinates)
  - Gold found via NuMI coordinates
  - LAr (not air/rock/steel) found 2cm away in every direction
  - All coordinate transforms are consistent across frames

Requires network access to download GDML files (~5 MB total).
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

from siren.detector import DetectorModel, DetectorPosition, GeometryPosition
from siren.download import ensure_files
from siren.math import Quaternion, Vector3D

_DATA_BASE = (
    "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/"
    "main/detectors/SBN/v1"
)

GOLD_DENSITY = 19.3
LAR_DENSITY = 1.39
GOLD_SIZE_M = 0.01


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


def _composite_gdml(det_filename, det_origin, gold_origin):
    """Composite GDML with a real detector + gold nugget in BNB frame."""
    return f"""\
<?xml version="1.0"?>
<gdml><define/><materials>
<isotope N="14" Z="7" name="env_N14"><atom unit="g/mole" value="14"/></isotope>
<element name="env_N"><fraction n="1.0" ref="env_N14"/></element>
<material name="EnvAir" state="gas">
  <D value="0.001225" unit="g/cm3"/><fraction n="1.0" ref="env_N"/>
</material>
</materials>
<solids><box name="sol_world" lunit="m" x="1000" y="1000" z="2000"/></solids>
<structure>
<volume name="volCompositeTestWorld">
  <materialref ref="EnvAir"/><solidref ref="sol_world"/>
  <physvol name="pv_det">
    <file name="{det_filename}" as_assembly="true"/>
    <position unit="m" x="{det_origin[0]:.10f}" y="{det_origin[1]:.10f}" z="{det_origin[2]:.10f}"/>
  </physvol>
  <physvol name="pv_gold">
    <file name="gold.gdml"/>
    <position unit="m" x="{gold_origin[0]:.10f}" y="{gold_origin[1]:.10f}" z="{gold_origin[2]:.10f}"/>
  </physvol>
</volume>
</structure>
<setup name="Default" version="1.0"><world ref="volCompositeTestWorld"/></setup>
</gdml>"""


def _vec(v):
    if hasattr(v, 'get'):
        v = v.get()
    return np.array([v.GetX(), v.GetY(), v.GetZ()])


def _build_model(geo, det_name, det_gdml_url, det_sha256, gold_pos_det, tmpdir):
    """Download detector GDML, add gold, load model, return (model, transforms)."""
    det_filename = f"{det_name.lower()}.gdml"
    det_path = os.path.join(tmpdir, det_filename)
    ensure_files([{"path": det_path, "url": det_gdml_url, "sha256": det_sha256}])

    gold_path = os.path.join(tmpdir, "gold.gdml")
    with open(gold_path, "w") as f:
        f.write(_gold_gdml())

    frame_name = f"{det_name}_LArSoft"
    T_det_bnb = geo.transform(frame_name, "BNB")
    det_origin_bnb = T_det_bnb.t
    gold_pos_bnb = T_det_bnb.apply(gold_pos_det)

    comp_path = os.path.join(tmpdir, f"composite_{det_name.lower()}.gdml")
    with open(comp_path, "w") as f:
        f.write(_composite_gdml(det_filename, det_origin_bnb, gold_pos_bnb))

    dm = DetectorModel()
    dm.LoadGDML(comp_path)
    qx, qy, qz, qw = geo.quaternion_from_matrix(T_det_bnb.R)
    dm.DetectorOrigin = GeometryPosition(Vector3D(*det_origin_bnb))
    dm.DetectorRotation = Quaternion(qx, qy, qz, qw)

    return dm


# ======================================================================
# ICARUS
# ======================================================================

# ICARUS center is all LAr (verified by probing)
ICARUS_GOLD_POS = np.array([0.0, -0.202, 0.0])


@pytest.fixture(scope="module")
def icarus_model(geo, tmp_path_factory):
    tmpdir = str(tmp_path_factory.mktemp("icarus"))
    return _build_model(
        geo, "ICARUS",
        f"{_DATA_BASE}/ICARUS/icarus_refactored_nounderscore_20230918_nowires.gdml",
        "3f68961a21fa037d3278c3235e48920cfc85a306adc5437f296dd7c09f095cd1",
        ICARUS_GOLD_POS, tmpdir)


class TestICARUSGold:

    def test_gold_in_det_coords(self, icarus_model):
        rho = icarus_model.GetMassDensity(
            DetectorPosition(Vector3D(*ICARUS_GOLD_POS)))
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_lar_2cm_away_all_directions(self, icarus_model):
        for axis in range(3):
            for sign in [-1, 1]:
                offset = np.zeros(3)
                offset[axis] = sign * 0.02
                pos = ICARUS_GOLD_POS + offset
                rho = icarus_model.GetMassDensity(
                    DetectorPosition(Vector3D(*pos)))
                label = "xyz"[axis]
                assert abs(rho - LAR_DENSITY) < 0.01, (
                    f"Expected LAr at {label}{'+' if sign > 0 else '-'}2cm, "
                    f"got {rho:.4f}")

    def test_gold_via_bnb(self, icarus_model, geo):
        T = geo.transform("ICARUS_LArSoft", "BNB")
        gold_bnb = T.apply(ICARUS_GOLD_POS)
        p_det = icarus_model.GeoPositionToDetPosition(
            GeometryPosition(Vector3D(*gold_bnb)))
        rho = icarus_model.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_gold_via_numi(self, icarus_model, geo):
        gold_numi = geo.transform("ICARUS_LArSoft", "NuMI").apply(ICARUS_GOLD_POS)
        gold_bnb = geo.transform("NuMI", "BNB").apply(gold_numi)
        p_det = icarus_model.GeoPositionToDetPosition(
            GeometryPosition(Vector3D(*gold_bnb)))
        rho = icarus_model.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_det_to_geo_matches_frame_graph(self, icarus_model, geo):
        T = geo.transform("ICARUS_LArSoft", "BNB")
        expected_bnb = T.apply(ICARUS_GOLD_POS)
        actual_bnb = _vec(icarus_model.DetPositionToGeoPosition(
            DetectorPosition(Vector3D(*ICARUS_GOLD_POS))))
        np.testing.assert_allclose(actual_bnb, expected_bnb, atol=1e-10)


# ======================================================================
# SBND
# ======================================================================

# SBND origin has a cathode at x=0; offset into the LAr drift region
SBND_GOLD_POS = np.array([0.5, 0.0, 0.0])


@pytest.fixture(scope="module")
def sbnd_model(geo, tmp_path_factory):
    tmpdir = str(tmp_path_factory.mktemp("sbnd"))
    return _build_model(
        geo, "SBND",
        f"{_DATA_BASE}/SBND/sbnd_v02_06.gdml",
        "224a0efa55e66b1fb3b527937e4a899854c2bbab367db3659e55466c4e7f013a",
        SBND_GOLD_POS, tmpdir)


class TestSBNDGold:

    def test_gold_in_det_coords(self, sbnd_model):
        rho = sbnd_model.GetMassDensity(
            DetectorPosition(Vector3D(*SBND_GOLD_POS)))
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_lar_2cm_away_all_directions(self, sbnd_model):
        for axis in range(3):
            for sign in [-1, 1]:
                offset = np.zeros(3)
                offset[axis] = sign * 0.02
                pos = SBND_GOLD_POS + offset
                rho = sbnd_model.GetMassDensity(
                    DetectorPosition(Vector3D(*pos)))
                label = "xyz"[axis]
                assert abs(rho - LAR_DENSITY) < 0.01, (
                    f"Expected LAr at {label}{'+' if sign > 0 else '-'}2cm, "
                    f"got {rho:.4f}")

    def test_gold_via_bnb(self, sbnd_model, geo):
        T = geo.transform("SBND_LArSoft", "BNB")
        gold_bnb = T.apply(SBND_GOLD_POS)
        p_det = sbnd_model.GeoPositionToDetPosition(
            GeometryPosition(Vector3D(*gold_bnb)))
        rho = sbnd_model.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_gold_via_numi(self, sbnd_model, geo):
        gold_numi = geo.transform("SBND_LArSoft", "NuMI").apply(SBND_GOLD_POS)
        gold_bnb = geo.transform("NuMI", "BNB").apply(gold_numi)
        p_det = sbnd_model.GeoPositionToDetPosition(
            GeometryPosition(Vector3D(*gold_bnb)))
        rho = sbnd_model.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_det_to_geo_matches_frame_graph(self, sbnd_model, geo):
        T = geo.transform("SBND_LArSoft", "BNB")
        expected_bnb = T.apply(SBND_GOLD_POS)
        actual_bnb = _vec(sbnd_model.DetPositionToGeoPosition(
            DetectorPosition(Vector3D(*SBND_GOLD_POS))))
        np.testing.assert_allclose(actual_bnb, expected_bnb, atol=1e-10)
