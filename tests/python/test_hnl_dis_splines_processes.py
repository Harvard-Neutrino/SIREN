"""Validation tests for the HNL DIS splines processes.py loader scripts.

The full success path (constructing HNLDISFromSpline / HNLDipoleDISFromSpline)
needs C++ classes that aren't on main until #122 lands. Here we just regression-
protect the review fixes that landed alongside these files.
"""
import importlib.util

import pytest


HNL_REL_DIRS = [
    "HNLDISSplines/HNLDISSplines-v1.0",
    "HNLDISSplines/HNLDISSplines-v2.0",
    "DipoleHNLDISSplines/DipoleHNLDISSplines-v1.0",
    "DipoleHNLDISSplines/DipoleHNLDISSplines-v2.0",
]


@pytest.fixture(params=HNL_REL_DIRS)
def processes_mod(request, processes_dir):
    rel = request.param
    path = processes_dir / rel / "processes.py"
    if not path.exists():
        raise FileNotFoundError(f"Could not find processes.py at {path}")

    spec = importlib.util.spec_from_file_location(
        f"_test_processes_{rel.replace('/', '_').replace('.', '_')}",
        str(path),
    )
    if spec is None:
        raise ImportError(f"Could not create import spec for {path}")
    if spec.loader is None:
        raise ImportError(f"Import spec for {path} has no loader")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_get_primary_types_rejects_unsupported(processes_mod):
    import siren
    with pytest.raises(ValueError, match="not supported"):
        processes_mod._get_primary_types([siren.dataclasses.Particle.ParticleType.EPlus])


def test_get_isoscalar_rejects_false(processes_mod):
    with pytest.raises(ValueError, match="Non-isoscalar splines"):
        processes_mod._get_isoscalar(False)


def test_isoscalar_error_mentions_correct_file_name(processes_mod):
    """Regression: the error said CSMSDISSplines-v1.0 in all four files
    before #116's review fixes."""
    with pytest.raises(ValueError) as exc:
        processes_mod._get_isoscalar(False)
    assert "CSMSDISSplines" not in str(exc.value), f"stale CSMSDISSplines name in {exc.value!r}"


def test_load_processes_requires_m4_mev(processes_mod):
    """Regression: int(m4_MeV) blew up if m4_MeV was None before #116's fix."""
    with pytest.raises(ValueError, match="m4_MeV"):
        processes_mod.load_processes()
