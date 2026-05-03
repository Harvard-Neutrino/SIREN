"""Tests for the high-level resource-loading API in siren.utilities."""
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def utilities():
    from siren import utilities
    return utilities


@pytest.fixture(scope="module")
def installed_resources_root(utilities):
    return Path(utilities.get_resource_package_dir())


def _has_resource(installed_resources_root: Path, kind: str, name: str) -> bool:
    return (installed_resources_root / kind / name).exists()


# Resources that have shipped with siren long enough that their absence
# from a built install is a packaging regression, not a stale checkout.
@pytest.mark.parametrize(
    "kind,name",
    [
        ("detectors", "IceCube"),
        ("detectors", "MINERvA"),
        ("detectors", "MiniBooNE"),
        ("fluxes", "BNB"),
        ("fluxes", "NUMI"),
    ],
)
def test_canonical_resource_present(installed_resources_root, kind, name):
    assert (installed_resources_root / kind / name).exists(), (
        f"{kind}/{name} missing from installed siren - packaging regression?"
    )


@pytest.mark.parametrize(
    "kind,name",
    [
        ("detectors", "IceCube"),
        ("detectors", "ND280"),
        ("detectors", "KM3NeTORCA"),
        ("detectors", "SINE"),
        ("detectors", "UNDINE"),
        ("detectors", "MINERvA"),
        ("detectors", "MiniBooNE"),
    ],
)
def test_detector_model_path_resolution(utilities, installed_resources_root, kind, name):
    if not _has_resource(installed_resources_root, kind, name):
        pytest.skip(f"installed siren has no {kind}/{name}")
    path = utilities.get_detector_model_path(name + "-v1")
    assert Path(path).exists(), f"resolved path {path!r} does not exist on disk"


def test_unknown_detector_raises(utilities):
    with pytest.raises(ValueError):
        utilities.get_detector_model_path("NoSuchDetector-v1")


def test_version_selection_picks_latest(utilities, tmp_path, monkeypatch):
    fake_root = tmp_path / "resources"
    for v in ("Fake-v1", "Fake-v2"):
        d = fake_root / "detectors" / "Fake" / v
        d.mkdir(parents=True)
        (d / "densities.dat").write_text("")
        (d / "materials.dat").write_text("")

    from siren import _util
    monkeypatch.setattr(_util, "resource_package_dir", lambda: str(fake_root))

    path = utilities.get_detector_model_path("Fake")
    assert path.endswith("Fake-v2"), f"latest version not selected: {path!r}"


def test_load_detector_returns_built_model(utilities, installed_resources_root):
    if not _has_resource(installed_resources_root, "detectors", "IceCube"):
        pytest.skip("installed siren has no detectors/IceCube")
    dm = utilities.load_detector("IceCube-v1")
    assert hasattr(dm, "Sectors")
    assert len(dm.Sectors) > 0


def test_load_flux_t2k_kaons(utilities, installed_resources_root, tmp_path):
    if not _has_resource(installed_resources_root, "fluxes", "T2K_Kaons"):
        pytest.skip("installed siren has no fluxes/T2K_Kaons")
    # Mirror the source data into tmp so flux.py's output stays out of site-packages.
    src = installed_resources_root / "fluxes" / "T2K_Kaons" / "T2K_Kaons-v1.0"
    for fname in ("kaon-flux-data.dat", "ratio.dat"):
        (tmp_path / fname).write_bytes((src / fname).read_bytes())
    out = utilities.load_flux("T2K_Kaons", "numu_PLUS", abs_flux_dir=str(tmp_path))
    assert Path(out).exists()
    assert Path(out).stat().st_size > 0
