"""Smoke tests for siren.visualization's optional-dependency-free surface."""
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest


def test_import_does_not_require_optional_deps():
    """siren.visualization must import even when every optional dep is absent.

    Simulated by blocking the optional imports outright; anything else in the
    siren import chain that merely prefers them (e.g. DarkNews) must already
    tolerate their absence.
    """
    code = (
        "import sys\n"
        "BLOCKED = ('matplotlib', 'pyg4ometry', 'pyvista', 'pygltflib', 'jinja2')\n"
        "class Blocker:\n"
        "    def find_spec(self, name, path=None, target=None):\n"
        "        if name.split('.')[0] in BLOCKED:\n"
        "            raise ImportError(name + ' is blocked for this test')\n"
        "        return None\n"
        "sys.meta_path.insert(0, Blocker())\n"
        "import siren.visualization\n"
        "leaked = [m for m in BLOCKED if m in sys.modules]\n"
        "assert not leaked, 'optional deps imported eagerly: %s' % leaked\n"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


def test_solid_xml_box_emits_full_widths():
    """A Box exports its full widths (GDML <box> x == Box.X)."""
    from siren import visualization
    from siren.geometry import Box

    box = Box(2.0, 4.0, 6.0)
    el = ET.fromstring(visualization._solid_xml(box, "TESTBOX"))
    assert el.tag == "box"
    assert float(el.get("x")) == pytest.approx(2.0)
    assert float(el.get("y")) == pytest.approx(4.0)
    assert float(el.get("z")) == pytest.approx(6.0)


@pytest.fixture(scope="module")
def ccm_model(detectors_dir):
    from siren.detector import DetectorModel
    base = detectors_dir / "CCM" / "CCM-v1"
    dm = DetectorModel()
    dm.LoadMaterialModel(str(base / "materials.dat"))
    dm.LoadDetectorModel(str(base / "densities.dat"))
    return dm


def test_to_gdml_writes_readable_gdml(ccm_model, tmp_path):
    from siren import visualization

    out = tmp_path / "ccm.gdml"
    stats = visualization.to_gdml(ccm_model, str(out))

    assert stats["path"] == str(out)
    assert stats["n_exact"] + stats["n_bbox"] > 0, "no sectors exported"
    assert stats["n_solids"] > 0 and stats["n_materials"] > 0

    root = ET.parse(out).getroot()
    assert root.tag == "gdml"

    solids = {s.get("name") for s in root.find("solids")}
    materials = {m.get("name") for m in root.find("materials")}
    volumes = {v.get("name") for v in root.find("structure")}
    assert "WORLD" in solids and "World" in volumes
    assert len(volumes) - 1 == stats["n_exact"] + stats["n_bbox"]

    # every reference in the structure must resolve to a defined name
    for vol in root.find("structure"):
        for ref in vol.iter("solidref"):
            assert ref.get("ref") in solids
        for ref in vol.iter("materialref"):
            assert ref.get("ref") in materials
        for ref in vol.iter("volumeref"):
            assert ref.get("ref") in volumes
    world_ref = root.find("setup/world").get("ref")
    assert world_ref in volumes


def test_describe_returns_material_table(ccm_model):
    from siren import visualization

    table = visualization.describe(ccm_model, printout=False)
    assert table, "describe() returned an empty material table"
    for name, density, count in table:
        assert isinstance(name, str) and name
        assert density >= 0.0
        assert count > 0
    assert sum(r[2] for r in table) == len(list(ccm_model.Sectors))
