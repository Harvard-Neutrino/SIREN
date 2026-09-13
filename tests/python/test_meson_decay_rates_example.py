"""Regression tests for the MesonDecay_rates.py example script.

Checks that its SM two-body denominator uses the same meson decay
constants MesonProduction uses internally, and that its ratio
calculation runs to a finite, positive result for every meson and
mediator-type configuration.
"""

import math

import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def example_script(resources_dir):
    from siren import _util

    script_path = resources_dir / "examples" / "example4" / "MesonDecay_rates.py"
    return _util.load_module("example4_meson_decay_rates", str(script_path))


def test_configs_use_meson_production_constants(example_script):
    """Each config's f_M/V_Mq must equal MesonProduction's own
    _meson_params(m_M), so the denominator cannot drift out of sync
    with the value MesonThreeBodyDecay uses."""
    meson_module = example_script._meson_module
    for label, m_M, m_l, f_M, V_Mq, mtype in example_script.configs:
        expected_f_M, expected_V_Mq = meson_module._meson_params(m_M)
        assert f_M == expected_f_M, label
        assert V_Mq == expected_V_Mq, label


def test_kaon_decay_constant_matches_pdg_value(example_script):
    """The script's kaon decay constant matches the PDG value
    MesonProduction uses internally."""
    assert example_script.FK == pytest.approx(0.1557)


def test_ratio_calculation_runs(example_script):
    """Run the script's ratio calculation for one mediator mass per
    config and check it returns a finite, positive number."""
    m_phi = 0.01  # GeV, below threshold for every config in this script
    for label, m_M, m_l, f_M, V_Mq, mtype in example_script.configs:
        decay = example_script.MesonThreeBodyDecay(
            m_M, m_l, m_phi, example_script.g_mu, mtype)
        w3 = decay.total_width()
        w2 = example_script.sm_width_2body(m_M, m_l, f_M, V_Mq)
        assert w2 > 0.0, label
        ratio = w3 / w2
        assert math.isfinite(ratio) and ratio > 0.0, label
