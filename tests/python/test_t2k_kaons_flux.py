"""Tests for the T2K_Kaons flux module."""
import importlib.util
import os

import pytest


T2K_RELDIR = "T2K_Kaons/T2K_Kaons-v1.0"


@pytest.fixture(scope="module")
def t2k_dir(fluxes_dir):
    return fluxes_dir / T2K_RELDIR


@pytest.fixture(scope="module")
def flux_mod(t2k_dir):
    spec = importlib.util.spec_from_file_location("t2k_flux_test", str(t2k_dir / "flux.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.fetch_data()
    return mod


@pytest.mark.parametrize("name", ["MakeFluxFile", "load_flux", "bar_scaling"])
def test_module_exposes_expected_callables(flux_mod, name):
    assert callable(getattr(flux_mod, name))


@pytest.mark.parametrize(
    "tag,err_substr",
    [
        ("",                "Tag must be"),
        ("noseparator",     "Tag must be"),
        ("a_b_c",           "Tag must be"),
        ("numu_FOO",        "Enhance tag"),
        ("numubar_BAR",     "Enhance tag"),
        ("foo_PLUS",        "is not valid"),
    ],
)
def test_invalid_tag_raises_value_error(flux_mod, tmp_path, tag, err_substr):
    with pytest.raises(ValueError, match=err_substr):
        flux_mod.MakeFluxFile(tag, output_dir=str(tmp_path))


@pytest.mark.parametrize(
    "tag,expected_basename",
    [
        ("numu_PLUS",      "kaon-flux-numu.dat"),
        ("numubar_MINUS",  "kaon-flux-numubar_bar.dat"),
        ("numubar_PLUS",   "kaon-flux-numubar.dat"),
        ("numu_MINUS",     "kaon-flux-numu_bar.dat"),
    ],
)
def test_valid_tag_produces_flux_file(flux_mod, tmp_path, tag, expected_basename):
    out = flux_mod.MakeFluxFile(tag, output_dir=str(tmp_path))
    assert os.path.basename(out) == expected_basename
    assert os.path.getsize(out) > 0
    with open(out) as f:
        first = f.readline().strip().split()
    assert len(first) == 2
    energy, flux = float(first[0]), float(first[1])
    assert energy > 0 and flux > 0


def test_load_flux_with_output_dir(flux_mod, tmp_path):
    out = flux_mod.load_flux("numu_PLUS", abs_flux_dir=str(tmp_path))
    assert os.path.exists(out)
    assert os.path.basename(out) == "kaon-flux-numu.dat"
    assert str(tmp_path) in out
