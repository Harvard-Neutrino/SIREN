"""End-to-end parsing tests for shipped detector geometries via DetectorModel."""
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def DetectorModel():
    from siren.detector import DetectorModel as DM
    return DM


def _shipped_detector_versions(detectors_dir: Path):
    pairs = []
    for det_root in sorted(p for p in detectors_dir.iterdir()
                           if p.is_dir() and not p.name.startswith(".") and p.name != "visuals"):
        for ver_dir in sorted(p for p in det_root.iterdir() if p.is_dir()):
            if (ver_dir / "densities.dat").exists() and (ver_dir / "materials.dat").exists():
                pairs.append((det_root.name, ver_dir.name))
    return pairs


def pytest_generate_tests(metafunc):
    if "detector" not in metafunc.fixturenames:
        return
    # Resolve detectors_dir at collection time (fixtures aren't available here).
    here = Path(__file__).resolve()
    for parent in (here, *here.parents):
        if (parent / "resources").exists():
            detectors_dir = parent / "resources" / "detectors"
            break
    else:
        metafunc.parametrize("detector", [], ids=[])
        return
    pairs = _shipped_detector_versions(detectors_dir)
    metafunc.parametrize("detector", pairs, ids=[f"{n}/{v}" for n, v in pairs])


def test_detector_loads(DetectorModel, detectors_dir, detector):
    name, version = detector
    materials = detectors_dir / name / version / "materials.dat"
    densities = detectors_dir / name / version / "densities.dat"

    dm = DetectorModel()
    dm.LoadMaterialModel(str(materials))
    dm.LoadDetectorModel(str(densities))

    sectors = dm.Sectors
    assert len(sectors) > 0, f"{name}/{version}: parsed but produced 0 sectors"


def test_icecube_v1_has_deepcore_sectors(DetectorModel, detectors_dir):
    """Regression test for the DeepCore extension added in #115."""
    base = detectors_dir / "IceCube" / "IceCube-v1"
    dm = DetectorModel()
    dm.LoadMaterialModel(str(base / "materials.dat"))
    dm.LoadDetectorModel(str(base / "densities.dat"))

    sector_names = {s.name for s in dm.Sectors}
    for expected in ("icecube", "icecube_shell", "deepcore", "deepcore_shell"):
        assert expected in sector_names, f"missing sector {expected!r}: have {sorted(sector_names)}"
