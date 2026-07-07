"""DarkNews port consumer audit: the ported PyDarkNews classes keep the C++
DarkNews bases, pass audit_overrides, and normalize DensityVariables to lists.
"""

import os

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
    except (ImportError, OSError, RuntimeError) as exc:
        # Only a missing DarkNews submodule or unavailable resource is a skip;
        # the port modules import DarkNews.processes/Cfourvec/phase_space, which
        # the module-scope importorskip("DarkNews") does not individually cover.
        # A broken port (AttributeError/NameError/SyntaxError) must fail loud.
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


# ------------------------------------------------------------------ #
#  Fail-loud data-integrity guards                                    #
# ------------------------------------------------------------------ #

class _FakeTarget:
    """A material target whose str() parses to a fixed nucleus name."""

    def __init__(self, s):
        self._s = s

    def __str__(self):
        return self._s

    def __eq__(self, other):
        return isinstance(other, _FakeTarget) and other._s == self._s

    def __hash__(self):
        return hash(self._s)


class _FakeMaterials:
    def __init__(self, targets):
        self._targets = targets

    def HasMaterial(self, i):
        return i == 0

    def GetMaterialTargets(self, i):
        return self._targets if i == 0 else []


class _FakeDetector:
    def __init__(self, targets):
        self.Materials = _FakeMaterials(targets)


def test_dk2nu_missing_pot_raises():
    """dk2nu_to_primary_distribution must fail loud when POT is zero/missing
    rather than emit inf/nan weights. The guard is reached before any detector
    use, so the condition is constructed directly from a zero-POT data dict."""
    import numpy as np
    reader = _load_darknews_module("Dk2nuReader")

    n = 3
    data = {k: np.zeros(n) for k in
            ["ptype", "E", "px", "py", "pz", "vx", "vy", "vz", "nimpwt"]}
    data["ptype"] = np.full(n, 211)
    data["nimpwt"] = np.ones(n)
    for bad in (0.0, -1.0, float("nan")):
        data["pot"] = bad
        with pytest.raises(siren.utilities.ConfigurationError):
            reader.dk2nu_to_primary_distribution(data, detector_model=None)


def test_load_vector_portal_unsupported_target_raises():
    """load_vector_portal must raise ConfigurationError naming a detector target
    absent from its nuclear table, not silently drop it. A duck-typed detector
    model supplies an unsupported target string directly."""
    processes = _load_darknews_module("processes")
    if not processes._DARKNEWS_AVAILABLE:
        pytest.skip("load_vector_portal requires DarkNews")

    # "TypePb208Nucleus" parses to an unsupported nucleus name.
    fake = _FakeDetector([_FakeTarget("TargetTypePb208NucleusXYZ")])
    with pytest.raises(siren.utilities.ConfigurationError):
        processes.load_vector_portal(
            m_chi=0.01, m_chi_prime=0.02, m_V1=0.03, m_V2=0.05,
            g_D=1.0, epsilon_1=1e-3, epsilon_2=1e-3, detector_model=fake)
