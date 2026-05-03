"""Tests for siren._util model-name and version parsing."""
import pytest


@pytest.fixture(scope="module")
def util():
    from siren import _util
    return _util


@pytest.mark.parametrize(
    "spec,expected_name,expected_release",
    [
        ("ND280-v1",                  "ND280",                   "1"),
        ("ND280UPGRD-v1.0",           "ND280UPGRD",              "1.0"),
        ("KM3NeTORCA-v1",             "KM3NeTORCA",              "1"),
        ("SINE-v1",                   "SINE",                    "1"),
        ("UNDINE-v1",                 "UNDINE",                  "1"),
        ("IceCube-v1",                "IceCube",                 "1"),
        ("MINERvA-v1",                "MINERvA",                 "1"),
        ("T2K_Kaons-v1.0",            "T2K_Kaons",               "1.0"),
        ("Lake_Geneva-v2",            "Lake_Geneva",             "2"),
        ("HNL_Decays-v1.0",           "HNL_Decays",              "1.0"),
        ("foo-bar-v1",                "foo-bar",                 "1"),
        ("a.b.c-v3",                  "a.b.c",                   "3"),
        ("just_name",                 "just_name",               None),
        ("just-name",                 "just-name",               None),
    ],
)
def test_model_regex_parses(util, spec, expected_name, expected_release):
    m = util._model_regex.match(spec)
    assert m is not None, f"regex did not match {spec!r}"
    assert m.group("model_name") == expected_name
    assert m.group("release") == expected_release


@pytest.mark.parametrize(
    "version,expected",
    [
        ("1",     {"release": "1",     "dev": None, "dev_num": None}),
        ("1.0",   {"release": "1.0",   "dev": None, "dev_num": None}),
        ("v1",    {"release": "1",     "dev": None, "dev_num": None}),
        ("v2.3",  {"release": "2.3",   "dev": None, "dev_num": None}),
    ],
)
def test_decompose_version(util, version, expected):
    parts = util.decompose_version(version)
    assert parts is not None, f"could not decompose {version!r}"
    for key, val in expected.items():
        assert parts.get(key) == val, f"{version!r}: {key}={parts.get(key)!r}, expected {val!r}"


@pytest.mark.parametrize(
    "v1,v2,greater",
    [
        ("1.0",  "0.9",  True),
        ("v1.1", "v1.0", True),
        ("2",    "10",   False),  # numeric, not lexicographic
    ],
)
def test_version_tokenize_orderable(util, v1, v2, greater):
    t1 = util.tokenize_version(v1)
    t2 = util.tokenize_version(v2)
    if greater:
        assert t1 > t2, f"{v1!r} should rank above {v2!r}"
    else:
        assert not (t1 > t2), f"{v1!r} should not rank above {v2!r}"
