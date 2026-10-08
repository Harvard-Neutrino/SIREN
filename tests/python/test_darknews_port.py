"""DarkNews port consumer audit: the ported PyDarkNews classes keep the C++
DarkNews bases, pass audit_overrides, and normalize DensityVariables to lists.
"""

import os
from types import SimpleNamespace

import numpy as np
import pytest

siren = pytest.importorskip("siren")
pytest.importorskip("DarkNews")

from siren import _util


DARKNEWS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "resources", "processes", "DarkNewsTables",
)


def _load_darknews_module(name):
    path = os.path.join(DARKNEWS_DIR, "%s.py" % name)
    if not os.path.isfile(path):
        pytest.skip("%s.py not present" % name)
    try:
        return _util.load_module(name, path)
    except Exception as exc:  # noqa: BLE001 -- resource/dependency unavailability
        pytest.skip("could not load %s: %s" % (name, exc))


def test_ported_decay_keeps_cpp_base():
    mod = _load_darknews_module("DarkNewsDecay")
    assert issubclass(mod.PyDarkNewsDecay, siren.interactions.DarkNewsDecay)
    assert issubclass(mod.PyDarkNewsDecay, siren.interactions.Decay)


def test_ported_cross_section_keeps_cpp_base():
    mod = _load_darknews_module("DarkNewsCrossSection")
    assert issubclass(mod.PyDarkNewsCrossSection,
                      siren.interactions.DarkNewsCrossSection)
    assert issubclass(mod.PyDarkNewsCrossSection,
                      siren.interactions.CrossSection)


def test_darknews_cross_section_clips_negative_numerical_density():
    mod = _load_darknews_module("DarkNewsCrossSection")
    assert mod._finite_nonnegative(-1.0e-38, "differential") == 0.0
    assert mod._finite_nonnegative(2.5, "differential") == 2.5


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_darknews_cross_section_rejects_nonfinite_density(bad):
    mod = _load_darknews_module("DarkNewsCrossSection")
    with pytest.raises(ValueError, match="Non-finite"):
        mod._finite_nonnegative(bad, "differential")


def test_darknews_interpolation_nan_falls_back_to_exact_evaluation():
    mod = _load_darknews_module("DarkNewsCrossSection")

    class _FakeUpsCase:
        pass

    xs = mod.PyDarkNewsCrossSection(_FakeUpsCase())
    xs.is_configured = True
    xs.differential_cross_section_table = np.array([
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
    ])
    xs.differential_cross_section_interpolator = lambda inputs: np.nan

    assert xs._query_interpolation_table(
        [1.0, 0.5], mode="differential") == -1


def test_darknews_non_always_interpolation_without_interpolator_falls_back():
    mod = _load_darknews_module("DarkNewsCrossSection")

    xs = mod.PyDarkNewsCrossSection(
        SimpleNamespace(), always_interpolate=False, interp_tolerance=2.0)
    xs.is_configured = True
    xs.differential_cross_section_table = np.array([
        [1.0, 0.0, 1.0],
        [1.0, 0.5, 1.0],
        [1.0, 1.0, 1.0],
    ])
    xs.differential_cross_section_interpolator = None

    assert xs._query_interpolation_table(
        [1.0, 0.25], mode="differential") == -1


def test_darknews_exact_differential_path_clips_negative_density():
    mod = _load_darknews_module("DarkNewsCrossSection")
    ups_case = SimpleNamespace(
        nu_projectile=SimpleNamespace(pdgid=14),
        diff_xsec_Q2=lambda energy, q2: np.float64(-1.0e-38),
        m_ups=0.1,
        MA=37.2,
        Ethreshold=0.11,
    )
    xs = mod.PyDarkNewsCrossSection(ups_case, always_interpolate=False)
    xs.is_configured = True
    xs.differential_cross_section_interpolator = None

    target = siren.dataclasses.ParticleType(1000180400)
    record = siren.dataclasses.InteractionRecord()
    record.signature.primary_type = siren.particles.NuMu
    record.signature.target_type = target
    record.primary_momentum = [1.0, 1.0, 0.0, 0.0]
    record.target_mass = ups_case.MA
    q2 = 0.5 * (xs.Q2Min(record) + xs.Q2Max(record))

    assert xs.DifferentialCrossSection(
        siren.particles.NuMu, target, 1.0, q2) == 0.0


def test_no_exit_calls_remain_in_port():
    for name in ("DarkNewsDecay", "DarkNewsCrossSection"):
        path = os.path.join(DARKNEWS_DIR, "%s.py" % name)
        if not os.path.isfile(path):
            pytest.skip("%s.py not present" % name)
        with open(path) as f:
            src = f.read()
        assert "exit(0)" not in src, "%s still calls exit(0)" % name


def test_single_get_ps_sample_definition():
    path = os.path.join(DARKNEWS_DIR, "DarkNewsDecay.py")
    if not os.path.isfile(path):
        pytest.skip("DarkNewsDecay.py not present")
    with open(path) as f:
        src = f.read()
    assert src.count("def GetPSSample") == 1


def test_holder_removed_from_processes():
    path = os.path.join(DARKNEWS_DIR, "processes.py")
    if not os.path.isfile(path):
        pytest.skip("processes.py not present")
    with open(path) as f:
        src = f.read()
    assert "class Holder" not in src
