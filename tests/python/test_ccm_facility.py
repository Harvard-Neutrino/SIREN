"""CCM-v3 resource contracts, using the installed native geometry engine."""

import hashlib
import importlib.util
import math
from pathlib import Path
import shutil

import numpy as np
import pytest


RESOURCE_ROOT = Path(__file__).resolve().parents[2] / "resources"
RESOURCE = RESOURCE_ROOT / "detectors" / "CCM" / "CCM-v3"


def module_at(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def bundle():
    return module_at(RESOURCE / "detector.py", "ccm_v3_test")


@pytest.fixture(scope="module")
def model(bundle):
    return bundle.load_detector()


def position(xyz):
    from siren.detector import DetectorPosition
    from siren.math import Vector3D
    return DetectorPosition(Vector3D(*xyz))


def xyz(v):
    return np.array([v.GetX(), v.GetY(), v.GetZ()])


def test_resource_selection_and_fiducial(monkeypatch, bundle):
    from siren import _util
    monkeypatch.setattr(_util, "resource_package_dir", lambda: str(RESOURCE_ROOT))
    assert Path(_util.get_detector_model_path("CCM")) == RESOURCE
    model = _util.load_detector("CCM-v3")
    assert model.GetMassDensity(position([0, 0, 0])) == pytest.approx(1.3954)
    fiducial = _util.get_fiducial_volume("CCM-v3")
    assert fiducial.Radius == pytest.approx(1.03385)
    assert fiducial.Z == pytest.approx(1.2396)
    assert xyz(fiducial.placement.Position) == pytest.approx([0, 0, 0])
    assert bundle.fiducial_volume().Radius == fiducial.Radius
    # Explicit v2 keeps its historical facility and detector coordinates.
    old = _util.load_detector("CCM-v2")
    assert xyz(old.DetectorOrigin.get()) == pytest.approx([23, 0, -.65])
    assert len(old.Sectors) < 100


def test_legacy_detector_layers_preserved(model, bundle):
    from siren.detector import DetectorModel
    old = DetectorModel()
    old.LoadMaterialModel(str(RESOURCE.parent / "CCM-v2" / "materials.dat"))
    old.LoadDetectorModel(str(RESOURCE.parent / "CCM-v2" / "densities.dat"))
    # Probe each cylindrical interface on both sides, radially and axially.
    points = []
    for radius, height in [(1.38, 2.62), (1.35, 2.52), (1.25, 2.40),
                           (1.20, 2.30), (1.130076, 1.24261), (1.12776, 1.2396)]:
        for delta in [-1e-6, 1e-6]:
            if radius + delta < 1.38:
                points.append([radius + delta, 0, 0])
            if height / 2 + delta < 1.31:
                points.extend([[0, 0, height / 2 + delta],
                               [0, 0, -height / 2 - delta]])
    for point in points:
        expected = old.GetMassDensity(position(point))
        if expected == 0:
            expected = 1e-25
        assert model.GetMassDensity(position(point)) == pytest.approx(expected, rel=1e-12, abs=0)
    sectors = {s.name: s for s in model.Sectors}
    for layer in bundle.geometry_description()["layers"]:
        sector = sectors[layer["name"]]
        old_sector = next(s for s in old.Sectors if s.name == layer["name"])
        assert sector.geo.Radius == old_sector.geo.Radius
        assert sector.geo.Z == old_sector.geo.Z
        for target in old.Materials.GetMaterialTargets(old_sector.material_id):
            expected = old.Materials.GetTargetMassFraction(old_sector.material_id, target)
            if layer["material"] == "STEEL":
                expected /= 1.0003
            assert model.Materials.GetTargetMassFraction(sector.material_id, target) == pytest.approx(expected)


def test_analytic_detector_column(model):
    # Independent chord calculation through the complete central radial stack.
    expected = 200 * ((1.37999 - 1.35) * 7.83 + (1.35 - 1.25) * 1e-25
                      + (1.25 - 1.20) * 7.83 + (1.20 - 1.130076) * 1.3954
                      + (1.130076 - 1.12776) * 2.70 + 1.12776 * 1.3954)
    actual = model.GetColumnDepthInCGS(position([-1.37999, 0, 0]),
                                     position([1.37999, 0, 0]))
    assert actual == pytest.approx(expected, rel=1e-11)


def test_placement_and_frame(model, bundle):
    from siren.detector import GeometryPosition
    from siren.math import Vector3D
    spec = bundle.geometry_description()["placement"]
    center = xyz(model.DetPositionToGeoPosition(position([0, 0, 0])).get())
    assert center == pytest.approx(spec["center_facility_m"])
    assert center[2] - 1.31 == pytest.approx(spec["floor_z_m"])
    target = xyz(model.GeoPositionToDetPosition(GeometryPosition(Vector3D(0, 0, 0))).get())
    baseline = np.linalg.norm(center[:2])
    theta = .654498  # Darcy's horizontal target-relative angle.
    assert target == pytest.approx([-baseline * math.cos(theta),
                                    -baseline * math.sin(theta), -center[2]])
    point = position([.123, -.456, .789])
    assert xyz(model.GeoPositionToDetPosition(model.DetPositionToGeoPosition(point)).get()) == pytest.approx(xyz(point.get()))


def test_regeneration_and_legacy_inputs(tmp_path, bundle):
    for name in ["facility.gdml", "geometry.json"]:
        shutil.copy2(RESOURCE / name, tmp_path / name)
    builder = module_at(RESOURCE / "build_geometry.py", "ccm_v3_builder_test")
    builder.build(tmp_path)
    for name in ["ccm-facility.gdml", "densities.dat"]:
        assert (tmp_path / name).read_bytes() == (RESOURCE / name).read_bytes()
    spec = bundle.geometry_description()
    for name, key in [("densities.dat", "legacy_density_sha256"),
                      ("materials.dat", "legacy_materials_sha256")]:
        old_bytes = (RESOURCE.parent / "CCM-v2" / name).read_bytes()
        assert hashlib.sha256(old_bytes).hexdigest() == spec[key]


def test_combination_retains_facility(model):
    names = {s.name for s in model.Sectors}
    assert len(model.Sectors) == 557
    assert "ccm_inner_argon" in names
    assert not any("CCMFloorLead" in name for name in names)
    assert "Fe_around_TMRS" not in names  # The old approximate target is absent.
    assert any("Facility_floorplan" in name for name in names)
