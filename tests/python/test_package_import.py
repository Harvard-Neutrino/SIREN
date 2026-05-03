"""Smoke tests for the top-level siren package surface."""
import importlib

import pytest


def test_top_level_import():
    siren = importlib.import_module("siren")
    assert hasattr(siren, "__version__")
    assert isinstance(siren.__version__, str)


@pytest.mark.parametrize(
    "submodule",
    [
        "utilities",
        "math",
        "dataclasses",
        "geometry",
        "detector",
        "interactions",
        "distributions",
        "injection",
        "resources",
        "_util",
    ],
)
def test_submodule_imports(submodule):
    importlib.import_module(f"siren.{submodule}")


@pytest.mark.parametrize(
    "name",
    ["load_flux", "load_detector", "load_processes",
     "get_flux_model_path", "get_detector_model_path", "get_processes_model_path",
     "get_resource_package_dir", "get_fiducial_volume"],
)
def test_utilities_public_helpers(name):
    from siren import utilities
    assert hasattr(utilities, name), f"siren.utilities is missing public helper {name!r}"
    assert callable(getattr(utilities, name))


def test_resources_public_helpers():
    from siren import resources
    for name in ("load_flux", "load_detector", "load_processes"):
        assert callable(getattr(resources, name)), f"siren.resources.{name} not callable"
    for name in ("fluxes", "detectors", "processes"):
        assert hasattr(resources, name), f"siren.resources missing {name}"
