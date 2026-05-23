"""Offline fixtures for the SBN detector loader.

These tests exercise detector.py and sbn_loader.py with tiny local GDML files
pre-seeded at the same relative paths as the downloadable SBN data.
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
import pytest

from siren.detector import DetectorPosition, GeometryPosition
from siren.math import Vector3D

LAR_DENSITY = 1.39

_SBN_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "resources", "detectors",
    "SBN", "SBN-v1")


@pytest.fixture
def sbn_detector_module():
    module_names = ("sbn_detector_offline",
                     "siren._sbn.sbn_geometry", "siren._sbn.sbn_loader",
                     "siren._sbn.earth_model")
    previous_modules = {name: sys.modules.pop(name, None) for name in module_names}
    try:
        spec = importlib.util.spec_from_file_location(
            "sbn_detector_offline", os.path.join(_SBN_DIR, "detector.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        yield mod
    finally:
        for name in module_names:
            sys.modules.pop(name, None)
        for name, module in previous_modules.items():
            if module is not None:
                sys.modules[name] = module


def _material_block(name, density):
    return f"""\
    <isotope N="14" Z="7" name="{name}_N14">
      <atom unit="g/mole" value="14"/>
    </isotope>
    <element name="{name}_N">
      <fraction n="1.0" ref="{name}_N14"/>
    </element>
    <material name="{name}" state="gas">
      <D value="{density}" unit="g/cm3"/>
      <fraction n="1.0" ref="{name}_N"/>
    </material>"""


def _beamline_fixture_gdml(prefix, material):
    return f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
{_material_block(material, 0.001225)}
  </materials>
  <solids>
    <box name="{prefix}_world_box" lunit="m" x="10" y="10" z="10"/>
  </solids>
  <structure>
    <volume name="{prefix}_world">
      <materialref ref="{material}"/>
      <solidref ref="{prefix}_world_box"/>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="{prefix}_world"/>
  </setup>
</gdml>
"""


def _detector_fixture_gdml(prefix):
    return f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <isotope N="40" Z="18" name="{prefix}_Ar40">
      <atom unit="g/mole" value="39.95"/>
    </isotope>
    <element name="{prefix}_Ar">
      <fraction n="1.0" ref="{prefix}_Ar40"/>
    </element>
    <material name="{prefix}_LAr" state="liquid">
      <D value="{LAR_DENSITY}" unit="g/cm3"/>
      <fraction n="1.0" ref="{prefix}_Ar"/>
    </material>
    <isotope N="14" Z="7" name="{prefix}_N14">
      <atom unit="g/mole" value="14"/>
    </isotope>
    <element name="{prefix}_N">
      <fraction n="1.0" ref="{prefix}_N14"/>
    </element>
    <material name="{prefix}_Air" state="gas">
      <D value="0.001225" unit="g/cm3"/>
      <fraction n="1.0" ref="{prefix}_N"/>
    </material>
  </materials>
  <solids>
    <box name="{prefix}_world_box" lunit="m" x="40" y="40" z="40"/>
    <box name="{prefix}_lar_box" lunit="m" x="10" y="10" z="10"/>
  </solids>
  <structure>
    <volume name="{prefix}_lar">
      <materialref ref="{prefix}_LAr"/>
      <solidref ref="{prefix}_lar_box"/>
    </volume>
    <volume name="{prefix}_world">
      <materialref ref="{prefix}_Air"/>
      <solidref ref="{prefix}_world_box"/>
      <physvol name="pv_{prefix}_lar">
        <volumeref ref="{prefix}_lar"/>
        <position unit="m" x="0" y="0" z="0"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="{prefix}_world"/>
  </setup>
</gdml>
"""


def _write_fixture_file(root, rel_path, content):
    path = root / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


@pytest.fixture
def offline_sbn_cache(tmp_path):
    _write_fixture_file(
        tmp_path, "gdml/BooNE_50m.gdml",
        _beamline_fixture_gdml("bnb_fixture", "BNBFixtureAir"))
    _write_fixture_file(
        tmp_path, "gdml/numi_g4export_2026-05-19.gdml",
        _beamline_fixture_gdml("numi_fixture", "NuMIFixtureAir"))
    _write_fixture_file(
        tmp_path, "gdml/icarus_refactored_nounderscore_20230918_nowires.gdml",
        _detector_fixture_gdml("icarus_fixture"))
    _write_fixture_file(
        tmp_path, "gdml/sbnd_v02_06.gdml",
        _detector_fixture_gdml("sbnd_fixture"))
    return tmp_path


def _forbid_download(*args, **kwargs):
    raise AssertionError("offline SBN fixture test attempted a network download")


@pytest.mark.parametrize("detector_name", ["ICARUS", "SBND"])
def test_load_detector_with_preseeded_gdml_offline(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name)
    rho = model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
    assert abs(rho - LAR_DENSITY) < 1e-12

    geo = sbn_detector_module.geo
    transform = geo.detector_transform(detector_name, "BNB")
    expected_origin = transform.apply(geo.DETECTORS[detector_name].center_native)
    actual_origin = model.GetDetectorOrigin().get()
    np.testing.assert_allclose(
        [actual_origin.GetX(), actual_origin.GetY(), actual_origin.GetZ()],
        expected_origin,
        atol=1e-12)

    composite = offline_sbn_cache / f"composite_{detector_name.lower()}.gdml"
    assert composite.is_file()


def test_fetch_data_uses_preseeded_gdml_offline(
        sbn_detector_module, offline_sbn_cache, monkeypatch):
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    sbn_detector_module.fetch_data()


def test_public_load_detector_with_preseeded_gdml_offline(
        offline_sbn_cache, monkeypatch):
    import siren.download as download
    from siren import _util

    resources_root = os.path.abspath(os.path.join(_SBN_DIR, "..", "..", ".."))
    module_names = ("siren-detector-SBN",
                     "siren._sbn.sbn_geometry", "siren._sbn.sbn_loader",
                     "siren._sbn.earth_model")
    previous_modules = {name: sys.modules.pop(name, None) for name in module_names}

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(download, "writable_data_dir", lambda module_dir: str(offline_sbn_cache))
    monkeypatch.setattr(_util, "resource_package_dir", lambda: resources_root)

    try:
        model = _util.load_detector("SBN", detector="ICARUS")
        rho = model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
        assert abs(rho - LAR_DENSITY) < 1e-12
    finally:
        for name in module_names:
            sys.modules.pop(name, None)
        for name, module in previous_modules.items():
            if module is not None:
                sys.modules[name] = module


# ---------------------------------------------------------------------------
# Earth model integration tests
# ---------------------------------------------------------------------------

def _geo_density(model, x, y, z):
    """Get mass density at a point given in BNB (geometry) coordinates."""
    gp = GeometryPosition(Vector3D(x, y, z))
    dp = model.GeoPositionToDetPosition(gp)
    return model.GetMassDensity(dp)


def _geo_sector_name(model, x, y, z):
    """Get containing sector name at a point in BNB (geometry) coordinates."""
    gp = GeometryPosition(Vector3D(x, y, z))
    dp = model.GeoPositionToDetPosition(gp)
    return model.GetContainingSector(dp).name


def _load_earth_constants():
    """Import earth_model module constants without full load."""
    import importlib.util as ilu
    spec = ilu.spec_from_file_location(
        "earth_model_consts", os.path.join(_SBN_DIR, "earth_model.py"))
    mod = ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_earth_model_sector_count(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """The loaded model should have PREM + atmosphere sectors alongside GDML."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    sectors = model.Sectors
    names = [s.name for s in sectors]

    earth = _load_earth_constants()
    expected_earth_count = len(earth._PREM_LAYERS) + len(earth._ATMO_LAYERS)
    earth_names = [name for name, *_ in earth._all_layers()]
    for en in earth_names:
        assert en in names, f"Missing earth sector: {en}"

    assert len(sectors) > expected_earth_count


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_detector_center_is_lar(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Detector origin (0,0,0 in detector coords) should still be LAr."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    rho = model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
    assert abs(rho - LAR_DENSITY) < 1e-12


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_atmosphere_column_depth(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Vertical column depth through atmosphere should be ~1030 g/cm2."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    earth = _load_earth_constants()

    column_depth = 0.0
    for h1, h2 in earth._ATMO_SHELL_ALTS:
        avg_rho = earth._atmo_avg_density(h1, h2)
        thickness_cm = (h2 - h1) * 100.0
        column_depth += avg_rho * thickness_cm

    assert 980 < column_depth < 1080, f"Column depth {column_depth:.1f} g/cm2 out of range"


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_prem_density_upper_crust(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Upper crust density should be 2.6 g/cm3 (PREM constant layer)."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    earth = _load_earth_constants()

    r_crust = earth.R_PREM - 5000.0
    y_bnb = r_crust - (earth.R_PREM - earth._GRADE_Y_BNB)

    rho = _geo_density(model, 0, y_bnb, 0)
    assert abs(rho - 2.6) < 0.01, f"Upper crust density {rho} != 2.6"


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_prem_density_upper_mantle(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Upper mantle (just below local Moho) should follow PREM polynomial."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    earth = _load_earth_constants()

    r_mantle = earth.R_PREM - earth._LOCAL_MOHO_DEPTH - 1000.0
    y_bnb = r_mantle - (earth.R_PREM - earth._GRADE_Y_BNB)

    rho = _geo_density(model, 0, y_bnb, 0)
    coeffs = [2.691, 1.08679956050855438e-07]
    expected = coeffs[0] + coeffs[1] * r_mantle
    assert abs(rho - expected) < 0.01, f"Mantle density {rho} vs expected {expected}"


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_local_moho_at_45km(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Near Fermilab, the local Moho correction should place the
    crust/mantle boundary at 45 km (CRUST1.0), not PREM's 24.4 km."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    earth = _load_earth_constants()

    # 30 km depth: between PREM Moho (24.4 km) and local Moho (45 km).
    # Should be ROCK thanks to the local correction.
    r_30 = earth.R_PREM - 30000.0
    y_30 = r_30 - (earth.R_PREM - earth._GRADE_Y_BNB)
    sector_30 = _geo_sector_name(model, 0, y_30, 0)
    assert sector_30 == "local_thick_crust", \
        f"30 km depth near Fermilab should be local_thick_crust, got {sector_30}"

    # 44 km depth: still above local Moho -> local_thick_crust
    r_44 = earth.R_PREM - 44000.0
    y_44 = r_44 - (earth.R_PREM - earth._GRADE_Y_BNB)
    sector_44 = _geo_sector_name(model, 0, y_44, 0)
    assert sector_44 == "local_thick_crust", \
        f"44 km depth near Fermilab should be local_thick_crust, got {sector_44}"

    # 46 km depth: below local Moho -> PREM mantle
    r_46 = earth.R_PREM - 46000.0
    y_46 = r_46 - (earth.R_PREM - earth._GRADE_Y_BNB)
    sector_46 = _geo_sector_name(model, 0, y_46, 0)
    assert sector_46 == "moho_boundary", \
        f"46 km depth near Fermilab should be moho_boundary, got {sector_46}"


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_prem_moho_far_from_fermilab(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Far from Fermilab (beyond the theta cap), the standard PREM Moho
    at 24.4 km should apply -- 30 km depth should be mantle, not crust."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    earth = _load_earth_constants()

    # A point at 30 km depth, 2500 km arc distance from Fermilab.
    # This is outside the local correction theta cap (0.3 rad ~ 1900 km).
    r = earth.R_PREM - 30000.0
    import math
    alpha = 2500000.0 / earth.R_PREM
    x = r * math.sin(alpha)
    y = r * math.cos(alpha) - (earth.R_PREM - earth._GRADE_Y_BNB)
    sector = _geo_sector_name(model, x, y, 0)
    assert sector == "moho_boundary", \
        f"30 km depth at 2500 km from Fermilab should be PREM mantle, got {sector}"


@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_innercore_density(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """Inner core (Earth center) density should be ~13 g/cm3."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=True)
    earth = _load_earth_constants()

    y_center = -(earth.R_PREM - earth._GRADE_Y_BNB)
    rho = _geo_density(model, 0, y_center, 0)
    assert 12.5 < rho < 14.0, f"Inner core density {rho} not in [12.5, 14.0]"


# ---------------------------------------------------------------------------
# earth_model=False (beam-only) tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("detector_name", ["ICARUS"])
def test_load_without_earth_model(
        sbn_detector_module, offline_sbn_cache, monkeypatch, detector_name):
    """With earth_model=False, the model should load without PREM sectors."""
    import siren.download as download
    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector(detector_name, earth_model=False)

    # Detector center should still be LAr
    rho = model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
    assert abs(rho - LAR_DENSITY) < 1e-12

    # There should be no PREM sectors
    sector_names = [s.name for s in model.Sectors]
    earth = _load_earth_constants()
    for layer_name, *_ in earth._all_layers():
        assert layer_name not in sector_names, \
            f"earth_model=False should not have sector {layer_name}"
    assert "local_thick_crust" not in sector_names
