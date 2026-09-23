"""End-to-end test with real detector GDML: gold nugget in LAr.

Uses the actual SBN detector loader infrastructure (sbn_loader,
sbn_geometry, detector.py constants) to build a composite GDML with
the real detector geometry plus a 1cm gold cube, then verifies:
  - Gold found at the expected position via DetectorCoordinates
  - Gold found via BNB (GeometryCoordinates)
  - Gold found via NuMI coordinates
  - LAr (not air/rock/steel) found 2cm away in every direction
  - DetectorPosition(0,0,0) maps to the LAr volume center (for ICARUS,
    to the air gap between its two cold vessels)
  - All coordinate transforms are consistent across frames

Requires --run-network to download GDML files (~5 MB total).
"""
from __future__ import annotations

import hashlib
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
AIR_DENSITY = 0.001205   # LArSoft "Air"
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

# The T600 is two separate cold vessels, so ICARUS' centre is the air gap
# between them: argon starts at |x| = 0.29 m and volTPCActive runs
# 0.62 - 3.585 m either side of a cathode at |x| = 2.10215 m. Put the
# nugget mid-drift in the +x module (ICARUS_C1).
ICARUS_GOLD_POS = np.array([1.5, -0.202, 0.0])


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

    def test_detector_origin_is_between_the_cryostats(self, icarus_model):
        """DetectorPosition(0,0,0) is the T600 centre: the inter-module gap.

        Pin air at the origin and LAr inside either module, so a shift in
        the detector placement fails here rather than passing quietly.
        """
        dm, _ = icarus_model
        rho = dm.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
        assert abs(rho - AIR_DENSITY) < 1e-5, (
            f"Expected the inter-cryostat air gap at the detector origin, "
            f"got {rho:.4f}")
        # z=1m avoids the nugget at z=0 in the +x module.
        for dx in [-1.5, 1.5]:
            rho = dm.GetMassDensity(DetectorPosition(Vector3D(dx, 0, 1.0)))
            assert abs(rho - LAR_DENSITY) < 0.01, (
                f"Expected LAr inside the module at dx={dx}, got {rho:.4f}")


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


# ======================================================================
# MicroBooNE
# ======================================================================

# Active-volume centre in the LArSoft world: the TPC box is centred at
# (1.28175, 0, 5.185) m and volTPCActive is offset (-1.55, 0.97, 0) cm from
# it. Gold goes half a metre away, well inside volTPCActive.
MICROBOONE_CENTER = np.array([1.28175 - 0.0155, 0.0097, 5.185])
MICROBOONE_GOLD_POS = MICROBOONE_CENTER + np.array([0.5, 0.0, 0.0])

# Detector coordinates and the density the production geometry puts there:
# steel vessel (an 11 mm wall, so the probe sits in the middle of it), foam,
# LArTF concrete, ground ring.
MICROBOONE_MATERIALS = [
    ((0.0, 1.90, 0.0), 7.93),
    ((0.0, 2.10, 0.0), 0.0384),
    ((7.6, 0.0, 0.0), 2.3),
    ((10.0, 0.0, 0.0), 1.7),
]

# The composite's atmosphere, as opposed to the 0.001205 LArSoft Air inside
# volDetEnclosure. uboonecode's LAr is 1.4, not the 1.39 ICARUS and SBND use.
ATMOSPHERE_DENSITY = 0.001225
MICROBOONE_LAR_DENSITY = 1.4


@pytest.fixture(scope="module")
def microboone_dir(sbn, tmp_path_factory):
    """Derive the SIREN copy of the production GDML, as load_detector does."""
    tmpdir = str(tmp_path_factory.mktemp("microboone"))
    _, _, det = sbn
    det._GENERATED_GDML["MicroBooNE"](tmpdir)
    return tmpdir


@pytest.fixture(scope="module")
def microboone_model(sbn, microboone_dir):
    return _build_model_with_gold(
        sbn, "MicroBooNE", MICROBOONE_GOLD_POS, microboone_dir)


class TestMicroBooNEGold:

    def test_gold_in_det_coords(self, microboone_model):
        dm, center = microboone_model
        gold_det = MICROBOONE_GOLD_POS - center
        rho = dm.GetMassDensity(DetectorPosition(Vector3D(*gold_det)))
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_lar_2cm_away_all_directions(self, microboone_model):
        dm, center = microboone_model
        gold_det = MICROBOONE_GOLD_POS - center
        for axis in range(3):
            for sign in [-1, 1]:
                offset = np.zeros(3)
                offset[axis] = sign * 0.02
                rho = dm.GetMassDensity(
                    DetectorPosition(Vector3D(*(gold_det + offset))))
                label = "xyz"[axis]
                assert abs(rho - MICROBOONE_LAR_DENSITY) < 0.01, (
                    f"Expected LAr at {label}{'+' if sign > 0 else '-'}2cm, "
                    f"got {rho:.4f}")

    def test_gold_via_bnb(self, microboone_model, sbn):
        dm, _ = microboone_model
        geo = sbn[0]
        gold_bnb = geo.transform("MicroBooNE_LArSoft", "BNB").apply(
            MICROBOONE_GOLD_POS)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*gold_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_gold_via_numi(self, microboone_model, sbn):
        dm, _ = microboone_model
        geo = sbn[0]
        gold_numi = geo.transform("MicroBooNE_LArSoft", "NuMI").apply(
            MICROBOONE_GOLD_POS)
        gold_bnb = geo.transform("NuMI", "BNB").apply(gold_numi)
        p_det = dm.GeoPositionToDetPosition(GeometryPosition(Vector3D(*gold_bnb)))
        rho = dm.GetMassDensity(p_det)
        assert abs(rho - GOLD_DENSITY) < 0.5

    def test_det_to_geo_matches_frame_graph(self, microboone_model, sbn):
        dm, center = microboone_model
        geo = sbn[0]
        gold_det = MICROBOONE_GOLD_POS - center
        expected_bnb = geo.transform("MicroBooNE_LArSoft", "BNB").apply(
            MICROBOONE_GOLD_POS)
        actual_bnb = _vec(dm.DetPositionToGeoPosition(
            DetectorPosition(Vector3D(*gold_det))))
        np.testing.assert_allclose(actual_bnb, expected_bnb, atol=1e-10)

    def test_detector_origin_is_in_lar(self, microboone_model):
        """The detector origin is the active-volume centre, so it and the
        points around it are liquid argon."""
        dm, _ = microboone_model
        for dx in [-0.1, 0.1]:
            rho = dm.GetMassDensity(DetectorPosition(Vector3D(dx, 0, 0)))
            assert abs(rho - MICROBOONE_LAR_DENSITY) < 0.01, (
                f"Expected LAr near detector origin (dx={dx}), got {rho:.4f}")

    def test_baseline_to_the_bnb_target(self, microboone_model):
        """The active volume is the published 468.5 m from the BNB target."""
        dm, _ = microboone_model
        origin = _vec(dm.GetDetectorOrigin())
        np.testing.assert_allclose(
            origin, [0.023, 0.0190, 468.548525], atol=1e-6)

    @pytest.mark.parametrize("point,density", MICROBOONE_MATERIALS,
                             ids=["steel", "foam", "concrete", "ground"])
    def test_building_materials(self, microboone_model, point, density):
        """Boolean-solid volumes of the production file parse and place."""
        dm, _ = microboone_model
        rho = dm.GetMassDensity(DetectorPosition(Vector3D(*point)))
        assert rho == pytest.approx(density, rel=2e-2)

    def test_vacuum_box_is_replaced_by_the_site_atmosphere(self,
                                                           microboone_model):
        """Left in place, LArSoft's vacuum box would own everything above
        grade, including the air over SBND 358.5 m upstream."""
        dm, _ = microboone_model
        for point in [(0.0, 20.0, 0.0), (0.0, 30.0, -358.5)]:
            rho = dm.GetMassDensity(DetectorPosition(Vector3D(*point)))
            assert rho == pytest.approx(ATMOSPHERE_DENSITY, rel=1e-3), (
                f"Expected the composite atmosphere at {point}, got {rho:.6f}")


def test_microboone_derives_from_the_pinned_production_asset(sbn,
                                                             microboone_dir):
    """The parsed file is the pinned asset, stripped of one placement."""
    _, loader, _ = sbn
    raw = os.path.join(microboone_dir, loader._MICROBOONE_SOURCE["file"])
    with open(raw, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    assert digest == loader._MICROBOONE_SOURCE["sha256"]

    derived = os.path.join(
        microboone_dir, "gdml", "microboonev12_nowires_siren.gdml")
    with open(derived, encoding="utf-8") as f:
        text = f.read()
    assert "volVacuumSpace" in text, "the volume definition is kept"
    assert "Derived by SIREN" in text.splitlines()[1]

    # Count as XML: a text search finds one more, because the geometry also
    # carries a commented-out volOverburden placement.
    def placements(path):
        import xml.etree.ElementTree as ET
        found = ET.parse(path).getroot().findall(".//physvol")
        refs = [p.find("volumeref").get("ref") for p in found
                if p.find("volumeref") is not None]
        return len(found), refs.count("volVacuumSpace")

    raw_total, raw_vacuum = placements(raw)
    derived_total, derived_vacuum = placements(derived)
    assert (raw_total, raw_vacuum) == (3945, 1)
    assert (derived_total, derived_vacuum) == (raw_total - 1, 0)


def test_numi_me_asset_target_placement(sbn, tmp_path):
    """The pinned ME asset places all 48 fins upstream of horn 1."""
    import xml.etree.ElementTree as ET

    _, loader, det = sbn
    source = next(s for s in det._beamline_sources(numi_config="ME")
                  if s["prefix"] == "numi")
    loader._ensure_gdml_files(str(tmp_path), [source])
    path = tmp_path / source["file"]
    root = ET.parse(path).getroot()
    volumes = {v.get("name"): v for v in root.find("structure")}
    positions = []

    def walk(volume, offset):
        for pv in volume.findall("physvol"):
            position = pv.find("position")
            local = [float(position.get(axis, "0")) if position is not None else 0.
                     for axis in "xyz"]
            world = offset + np.array(local)
            child = volumes[pv.find("volumeref").get("ref")]
            if pv.get("name", "").startswith("TGT10x"):
                positions.append(world / 1000.)
            walk(child, world)

    # This pinned export's target and ancestors have identity rotations and mm units.
    walk(volumes[root.find("setup/world").get("ref")], np.zeros(3))
    assert len(positions) == 48
    np.testing.assert_allclose(sorted(p[2] for p in positions),
                               (-1363.5 + 24.5*np.arange(48))/1000., atol=1e-10)
    model = DetectorModel()
    model.LoadGDML(str(path), True)
    for point in positions:
        assert model.GetMassDensity(DetectorPosition(Vector3D(*point))) == pytest.approx(1.78)
