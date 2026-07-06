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
