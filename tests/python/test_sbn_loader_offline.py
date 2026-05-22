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

from siren.detector import DetectorPosition
from siren.math import Vector3D

LAR_DENSITY = 1.39

_SBN_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "resources", "detectors",
    "SBN", "SBN-v1")


@pytest.fixture
def sbn_detector_module():
    module_names = ("sbn_detector_offline", "sbn_geometry", "sbn_loader")
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
    module_names = ("siren-detector-SBN", "sbn_geometry", "sbn_loader")
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
