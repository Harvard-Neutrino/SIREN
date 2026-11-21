"""Offline fixtures for the SBN detector loader.

These tests exercise detector.py and sbn_loader.py with tiny local GDML files
pre-seeded at the same relative paths as the downloadable SBN data.
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from threading import Barrier, Event

import numpy as np
import pytest

from siren.detector import DetectorPosition, GeometryPosition
from siren.math import Vector3D

LAR_DENSITY = 1.39

_SBN_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "resources", "detectors",
    "SBN", "SBN-v1")


def _bnb_ground_level(sbn_detector_module):
    """BNB-frame ground level used by the earth model.

    Mirrors the value detector.add_earth_model() is called with
    (sbn_loader._FNAL_SITE_GRADE_Y + geo.T_MiniBooNE_local[1]) so density
    tests can convert PREM radii to BNB-frame y-coordinates.
    """
    return (sbn_detector_module.sbn_loader._FNAL_SITE_GRADE_Y
            + sbn_detector_module.geo.T_MiniBooNE_local[1])


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


MICROBOONE_LAR_DENSITY = 1.40


def _microboone_fixture_gdml():
    """Stand-in for the uboonecode export: a LAr box at the TPC-box centre
    and, as in the real file, a volVacuumSpace placement to drop."""
    return f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
    <isotope N="40" Z="18" name="ub_fixture_Ar40">
      <atom unit="g/mole" value="39.95"/>
    </isotope>
    <element name="ub_fixture_Ar">
      <fraction n="1.0" ref="ub_fixture_Ar40"/>
    </element>
    <material name="ub_fixture_LAr" state="liquid">
      <D value="{MICROBOONE_LAR_DENSITY}" unit="g/cm3"/>
      <fraction n="1.0" ref="ub_fixture_Ar"/>
    </material>
    <material name="ub_fixture_Vacuum" state="gas">
      <D value="1e-25" unit="g/cm3"/>
      <fraction n="1.0" ref="ub_fixture_Ar"/>
    </material>
{_material_block("ub_fixture_Air", 0.001205)}
  </materials>
  <solids>
    <box name="ub_fixture_world_box" lunit="m" x="80" y="80" z="80"/>
    <box name="ub_fixture_lar_box" lunit="m" x="3" y="3" z="12"/>
    <box name="ub_fixture_vacuum_box" lunit="m" x="60" y="20" z="60"/>
  </solids>
  <structure>
    <volume name="volTPCActive">
      <materialref ref="ub_fixture_LAr"/>
      <solidref ref="ub_fixture_lar_box"/>
    </volume>
    <volume name="volVacuumSpace">
      <materialref ref="ub_fixture_Vacuum"/>
      <solidref ref="ub_fixture_vacuum_box"/>
    </volume>
    <volume name="volWorld">
      <materialref ref="ub_fixture_Air"/>
      <solidref ref="ub_fixture_world_box"/>
      <physvol name="pv_ub_fixture_vacuum">
        <volumeref ref="volVacuumSpace"/>
        <position unit="m" x="0" y="16.25" z="0"/>
      </physvol>
      <physvol name="pv_ub_fixture_lar">
        <volumeref ref="volTPCActive"/>
        <position unit="m" x="1.28175" y="0" z="5.185"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="volWorld"/>
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
        tmp_path, "gdml/numi_ME_g4export_2026-09-17.gdml",
        _beamline_fixture_gdml("numi_fixture", "NuMIFixtureAir"))
    _write_fixture_file(
        tmp_path, "gdml/g4lbnf.gdml",
        _beamline_fixture_gdml("lbnf_fixture", "LBNFFixtureAir"))
    _write_fixture_file(
        tmp_path, "gdml/icarus_refactored_nounderscore_20230918_nowires.gdml",
        _detector_fixture_gdml("icarus_fixture"))
    _write_fixture_file(
        tmp_path, "gdml/sbnd_v02_06.gdml",
        _detector_fixture_gdml("sbnd_fixture"))
    _write_fixture_file(
        tmp_path, "gdml/nd_hall_with_lar_tms_sand.gdml",
        _detector_fixture_gdml("dune_nd_fixture"))
    # The raw uboonecode download; the loader derives its SIREN copy from it.
    _write_fixture_file(
        tmp_path, "gdml/microboonev12_nowires.gdml",
        _microboone_fixture_gdml())
    return tmp_path


def _forbid_download(*args, **kwargs):
    raise AssertionError("offline SBN fixture test attempted a network download")


@pytest.fixture
def microboone_pin(sbn_detector_module, monkeypatch):
    """Pin the loader's digest to the stand-in the offline fixture seeds."""
    digest = hashlib.sha256(
        _microboone_fixture_gdml().encode("utf-8")).hexdigest()
    monkeypatch.setitem(
        sbn_detector_module.sbn_loader._MICROBOONE_SOURCE, "sha256", digest)
    return digest


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

    composites = list(offline_sbn_cache.glob(f"composite_{detector_name.lower()}_*.gdml"))
    assert len(composites) == 1


def test_load_microboone_with_preseeded_gdml_offline(
        sbn_detector_module, offline_sbn_cache, monkeypatch, microboone_pin):
    """MicroBooNE derives its SIREN GDML from the pre-seeded uboonecode file
    without downloading, drops the LArSoft vacuum box, and sits at the
    surveyed baseline."""
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    model = sbn_detector_module.load_detector("MicroBooNE")

    derived = offline_sbn_cache / "gdml" / "microboonev12_nowires_siren.gdml"
    assert derived.is_file()
    text = derived.read_text()
    assert 'volumeref ref="volVacuumSpace"' not in text
    assert 'volumeref ref="volTPCActive"' in text
    assert '<volume name="volVacuumSpace">' in text  # definition kept, unplaced
    assert text.startswith("<?xml")
    assert "<!-- Derived by SIREN from microboonev12_nowires.gdml" in text.splitlines()[1]

    rho = model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
    assert rho == pytest.approx(MICROBOONE_LAR_DENSITY, rel=1e-12, abs=0)
    assert model.GetContainingSector(DetectorPosition(Vector3D(0, 0, 0))).name == "volTPCActive"

    origin = model.GetDetectorOrigin().get()
    geo = sbn_detector_module.geo
    expected = geo.detector_center("MicroBooNE", "BNB")
    np.testing.assert_allclose(
        [origin.GetX(), origin.GetY(), origin.GetZ()], expected, atol=1e-12)
    assert abs(expected[2] - 468.5) < 0.1  # the published MicroBooNE baseline

    # Where the fixture's vacuum box would sit, the composite atmosphere remains.
    assert _geo_density(model, 0.0, 12.0, expected[2]) == pytest.approx(0.001225, rel=1e-12, abs=0)
    assert _geo_sector_name(model, 0.0, 12.0, expected[2]) == "vol_atmosphere"

    # A second load reuses the derived file without touching the raw one.
    raw = offline_sbn_cache / "gdml" / "microboonev12_nowires.gdml"
    stamp = (raw.stat().st_mtime_ns, derived.stat().st_mtime_ns)
    sbn_detector_module.load_detector("MicroBooNE")
    assert (raw.stat().st_mtime_ns, derived.stat().st_mtime_ns) == stamp


def test_microboone_rejects_a_cached_file_with_the_wrong_digest(
        sbn_detector_module, offline_sbn_cache, monkeypatch, microboone_pin):
    """The loader checks the pin itself: ensure_files trusts any file
    already on disk."""
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))
    raw = offline_sbn_cache / "gdml" / "microboonev12_nowires.gdml"
    raw.write_text(
        _microboone_fixture_gdml().replace("volTPCActive", "volTampered"))

    with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
        sbn_detector_module.load_detector("MicroBooNE")


def test_microboone_rebuilds_a_derived_copy_from_another_source(
        sbn_detector_module, offline_sbn_cache, monkeypatch, microboone_pin):
    """A derived copy made from a different source file is replaced."""
    loader = sbn_detector_module.sbn_loader
    loader.ensure_microboone_gdml(str(offline_sbn_cache))
    derived = offline_sbn_cache / "gdml" / "microboonev12_nowires_siren.gdml"
    assert microboone_pin in derived.read_text().splitlines()[1]

    derived.write_text(
        '<?xml version="1.0"?>\n'
        "<!-- Derived by SIREN from microboonev12_nowires.gdml "
        "(uboone/ubcore microboonev12, sha256 " + "0" * 64 + "): stale. -->\n"
        "<gdml/>\n")
    loader.ensure_microboone_gdml(str(offline_sbn_cache))
    text = derived.read_text()
    assert microboone_pin in text.splitlines()[1]
    assert 'volumeref ref="volTPCActive"' in text


def test_microboone_rebuilds_an_edited_copy_with_a_current_marker(
        sbn_detector_module, offline_sbn_cache, microboone_pin):
    """A copy that still places the vacuum box is replaced even though its
    marker line is current."""
    loader = sbn_detector_module.sbn_loader
    loader.ensure_microboone_gdml(str(offline_sbn_cache))
    derived = offline_sbn_cache / "gdml" / "microboonev12_nowires_siren.gdml"
    good = derived.read_bytes()
    marker = good.decode("utf-8").splitlines(keepends=True)[1]
    assert microboone_pin in marker

    head, sep, rest = _microboone_fixture_gdml().partition("?>\n")
    derived.write_text(head + sep + marker + rest)  # vacuum placement back
    assert 'volumeref ref="volVacuumSpace"' in derived.read_text()

    loader.ensure_microboone_gdml(str(offline_sbn_cache))
    assert derived.read_bytes() == good


def test_microboone_concurrent_rebuild_of_a_stale_copy(
        sbn_detector_module, offline_sbn_cache, microboone_pin):
    """Concurrent loads that find a stale copy all succeed; none finds the
    file missing while another replaces it."""
    loader = sbn_detector_module.sbn_loader
    derived = offline_sbn_cache / "gdml" / "microboonev12_nowires_siren.gdml"
    stale = ('<?xml version="1.0"?>\n<!-- Derived by SIREN (sha256 '
             + "0" * 64 + "): stale. -->\n<gdml/>\n")
    workers = 8
    for _ in range(20):
        derived.write_text(stale)
        barrier = Barrier(workers)

        def load(_):
            barrier.wait()
            return loader.ensure_microboone_gdml(str(offline_sbn_cache))

        with ThreadPoolExecutor(workers) as pool:
            results = list(pool.map(load, range(workers)))
        assert set(results) == {"gdml/microboonev12_nowires_siren.gdml"}
        assert microboone_pin in derived.read_text().splitlines()[1]


def test_strip_physvols_removes_only_named_placements(sbn_detector_module):
    loader = sbn_detector_module.sbn_loader
    text = _microboone_fixture_gdml()
    stripped, removed = loader._strip_physvols(text, ("volVacuumSpace",))
    assert removed == 1
    assert 'volumeref ref="volVacuumSpace"' not in stripped
    assert stripped.count("<physvol") == text.count("<physvol") - 1
    same, none = loader._strip_physvols(text, ("volNotPlaced",))
    assert none == 0 and same == text


def test_fetch_data_uses_preseeded_gdml_offline(
        sbn_detector_module, offline_sbn_cache, monkeypatch, microboone_pin):
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))

    sbn_detector_module.fetch_data()


@pytest.mark.parametrize("detector_name", ["ICARUS", "MiniBooNE"])
@pytest.mark.parametrize("lbnf_flags", [(False, True), (False, False), (True, True)],
                         ids=["distinct", "identical-sbn", "identical-lbnf"])
def test_concurrent_compositions_load_requested_geometry(
        sbn_detector_module, offline_sbn_cache, monkeypatch,
        detector_name, lbnf_flags):
    """Both writes finish before either native reader opens its composition."""
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(offline_sbn_cache))
    loader = sbn_detector_module.sbn_loader
    build = loader.build_composite
    written = Barrier(2)

    def build_before_read(*args, **kwargs):
        path = build(*args, **kwargs)
        written.wait(timeout=15)
        return path

    monkeypatch.setattr(loader, "build_composite", build_before_read)
    if detector_name == "MiniBooNE":
        # Force both callers past the missing-file check on a cold cache.
        assert not (offline_sbn_cache / "gdml/miniboone_tank.gdml").exists()
        atomic_output = download.atomic_output_path
        tank_writers_met = Event()
        tank_writers = Barrier(2, action=tank_writers_met.set)

        @contextmanager
        def simultaneous_tank_writes(path):
            with atomic_output(path) as tmp:
                if Path(path).name == "miniboone_tank.gdml":
                    tank_writers.wait(timeout=15)
                yield tmp

        monkeypatch.setattr(download, "atomic_output_path", simultaneous_tank_writes)

    with ThreadPoolExecutor(max_workers=2) as pool:
        pending = [pool.submit(sbn_detector_module.load_detector, detector_name,
                               lbnf=flag) for flag in lbnf_flags]
        models = [future.result(timeout=30) for future in pending]

    if detector_name == "MiniBooNE":
        assert tank_writers_met.is_set(), "Cold-cache tank writers did not synchronize"

    for lbnf, model in zip(lbnf_flags, models):
        assert any("lbnf_fixture_world" in s.name for s in model.Sectors) == lbnf
        assert model.Materials.HasMaterial("LBNFFixtureAir") == lbnf
        expected_density = 0.845 if detector_name == "MiniBooNE" else LAR_DENSITY
        assert model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0))) == pytest.approx(
            expected_density, rel=1e-12, abs=0)
        if lbnf:
            origin = sbn_detector_module.geo.transform("LBNF", "BNB").t
            assert "lbnf_fixture_world" in _geo_sector_name(model, *origin)
            assert _geo_density(model, *origin) == pytest.approx(0.001225, rel=1e-12, abs=0)


@pytest.mark.parametrize("field,value", [
    ("file", "gdml/g4lbnf.gdml"),
    ("prefix", "other_beam"),
    ("position", (15.0, 20.0, 25.0)),
    ("rotation", (0.1, 0.2, 0.3)),
    ("unwrap", True),
])
def test_composition_variants_preserve_previous_file(
        sbn_detector_module, offline_sbn_cache, monkeypatch, field, value):
    import siren.download as download

    monkeypatch.setattr(download, "download_file", _forbid_download)
    loader = sbn_detector_module.sbn_loader
    sources = sbn_detector_module._beamline_sources()
    first = Path(loader.build_composite(str(offline_sbn_cache), sources))
    contents = first.read_bytes()
    changed = [dict(source) for source in sources]
    changed[0][field] = value
    second = Path(loader.build_composite(str(offline_sbn_cache), changed))

    assert first != second
    assert first.read_bytes() == contents
    assert second.read_bytes() != contents
    assert Path(loader.build_composite(str(offline_sbn_cache), sources)) == first


def test_composition_identity_includes_site_geometry(
        sbn_detector_module, offline_sbn_cache, monkeypatch):
    loader = sbn_detector_module.sbn_loader
    first = Path(loader.build_composite(str(offline_sbn_cache), []))
    contents = first.read_bytes()
    monkeypatch.setattr(loader, "_TILL_THICKNESS", loader._TILL_THICKNESS + 1)
    second = Path(loader.build_composite(str(offline_sbn_cache), []))
    assert second != first
    assert first.read_bytes() == contents
    assert second.read_bytes() != contents


@pytest.mark.parametrize("numi_options", [{}, {"numi_config": "ME"}, {"numi_config": "me"}])
def test_public_load_detector_with_preseeded_gdml_offline(
        offline_sbn_cache, monkeypatch, numi_options):
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
        model = _util.load_detector("SBN", detector="ICARUS", **numi_options)
        rho = model.GetMassDensity(DetectorPosition(Vector3D(0, 0, 0)))
        assert abs(rho - LAR_DENSITY) < 1e-12
        composite, = offline_sbn_cache.glob("composite_icarus_*.gdml")
        assert 'name="gdml/numi_ME_g4export_2026-09-17.gdml"' in composite.read_text()
        assert "numi_g4export_2026-05-19.gdml" not in composite.read_text()
    finally:
        for name in module_names:
            sys.modules.pop(name, None)
        for name, module in previous_modules.items():
            if module is not None:
                sys.modules[name] = module


# ---------------------------------------------------------------------------
# NuMI configuration selection
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("config", ["LE", "HE", "me000z200i", "", None])
def test_unsupported_numi_configuration_has_no_geometry_side_effects(
        sbn_detector_module, monkeypatch, tmp_path, config):
    monkeypatch.setattr(sbn_detector_module, "_ABS_DIR", str(tmp_path))
    monkeypatch.setattr(sbn_detector_module.sbn_loader, "_ensure_gdml_files", _forbid_download)
    monkeypatch.setattr(sbn_detector_module.sbn_loader, "ensure_miniboone_gdml", _forbid_download)
    with pytest.raises(ValueError, match="Unsupported NuMI configuration.*Supported configurations: ME"):
        sbn_detector_module.load_detector("MiniBooNE", numi_config=config)
    assert list(tmp_path.iterdir()) == []


def test_numi_me_selection_preserves_other_beamlines(sbn_detector_module):
    ordinary = sbn_detector_module._beamline_sources()
    explicit = sbn_detector_module._beamline_sources(numi_config="me")
    assert ordinary == explicit
    extended = sbn_detector_module._beamline_sources(lbnf=True, numi_config="ME")
    assert extended[:2] == ordinary
    assert [s["prefix"] for s in extended] == ["bnb", "numi", "lbnf"]
    assert ordinary[1]["sha256"] == "730466f287196d65a7fee074203014471faee6be0fbfa3da4769046d92355ed7"


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
    y_bnb = r_crust - (earth.R_PREM - _bnb_ground_level(sbn_detector_module))

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
    y_bnb = r_mantle - (earth.R_PREM - _bnb_ground_level(sbn_detector_module))

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
    y_30 = r_30 - (earth.R_PREM - _bnb_ground_level(sbn_detector_module))
    sector_30 = _geo_sector_name(model, 0, y_30, 0)
    assert sector_30 == "local_thick_crust", \
        f"30 km depth near Fermilab should be local_thick_crust, got {sector_30}"

    # 44 km depth: still above local Moho -> local_thick_crust
    r_44 = earth.R_PREM - 44000.0
    y_44 = r_44 - (earth.R_PREM - _bnb_ground_level(sbn_detector_module))
    sector_44 = _geo_sector_name(model, 0, y_44, 0)
    assert sector_44 == "local_thick_crust", \
        f"44 km depth near Fermilab should be local_thick_crust, got {sector_44}"

    # 46 km depth: below local Moho -> PREM mantle
    r_46 = earth.R_PREM - 46000.0
    y_46 = r_46 - (earth.R_PREM - _bnb_ground_level(sbn_detector_module))
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
    y = r * math.cos(alpha) - (earth.R_PREM - _bnb_ground_level(sbn_detector_module))
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

    y_center = -(earth.R_PREM - _bnb_ground_level(sbn_detector_module))
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
