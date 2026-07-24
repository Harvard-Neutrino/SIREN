"""Module-only resources: modules that expose model classes directly and set
_SIREN_RESOURCE_MODULE_ONLY instead of providing a load_<type> factory."""

import pytest

siren = pytest.importorskip("siren")
from siren import _util


def test_module_only_resource_loads_as_module():
    module = siren.resources.processes.BeamDecays
    assert module._SIREN_RESOURCE_MODULE_ONLY
    assert hasattr(module, "MesonTwoBodyLeptonicDecay")
    assert hasattr(module, "MuonThreeBodyDecay")


def test_module_only_resource_rejects_factory_loading():
    with pytest.raises(
            TypeError,
            match=r"siren\.resources\.processes\.BeamDecays"):
        siren.load_processes("BeamDecays")


def test_module_only_resource_documents_itself():
    doc = _util.process_docs("BeamDecays")
    assert "decay" in doc.lower()
