"""C++ MultiChannelPhaseSpace.UpdateWeights == Python _kp_update.

Phase B moved the Kleiss-Pittau accumulation + weight update into C++.  These
tests pin that the C++ update reproduces the reference Python `_kp_update`
(plus the caller's damping blend) bit-for-bit on identical inputs, and that the
nested inner/outer point counts stay in lockstep (so each level divides its
accumulator by the same count, exactly as the old single n_nonzero did).
"""
import math

import pytest

import siren
from siren.tune import _kp_update

PT = siren.dataclasses.ParticleType
M_V1 = 0.017


def _iso_mixture(weights):
    """A minimal mixture; UpdateWeights only reads kp_accumulator/kp_count and
    weights, so the channel kinematics are irrelevant here."""
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [siren.injection.Isotropic2BodyChannel(0) for _ in weights]
    mc.weights = list(weights)
    return mc


@pytest.mark.parametrize("rule", ["sqrt_W", "alpha_sqrt_W"])
@pytest.mark.parametrize("min_weight", [0.0, 0.1])
@pytest.mark.parametrize("damping", [1.0, 0.5])
def test_cpp_updateweights_matches_python_kp_update(rule, min_weight, damping):
    weights = [0.3, 0.7]
    W = [4.0, 1.0]
    count = 1000

    mc = _iso_mixture(weights)
    mc.kp_accumulator = [W[i] * count for i in range(len(W))]
    mc.kp_count = count
    mc.UpdateWeights(rule, damping, min_weight)

    new_alpha = _kp_update(list(weights), W, rule, min_weight)
    expected = [damping * new_alpha[i] + (1.0 - damping) * weights[i]
                for i in range(len(weights))]

    assert mc.weights == pytest.approx(expected, abs=1e-12)
    # UpdateWeights resets its accumulators at the end of the call.
    assert mc.kp_count == 0


def test_cpp_updateweights_degenerate_keeps_weights():
    weights = [0.4, 0.6]
    mc = _iso_mixture(weights)
    mc.kp_accumulator = [0.0, 0.0]
    mc.kp_count = 500
    mc.UpdateWeights("sqrt_W", 1.0, 0.0)
    assert mc.weights == pytest.approx(weights, abs=1e-15)


def test_cpp_updateweights_zero_count_keeps_weights():
    weights = [0.4, 0.6]
    mc = _iso_mixture(weights)
    mc.UpdateWeights("alpha_sqrt_W", 0.5, 0.01)
    assert mc.weights == pytest.approx(weights, abs=1e-15)


def test_cpp_updateweights_unknown_rule_raises():
    mc = _iso_mixture([0.5, 0.5])
    with pytest.raises(ValueError):
        mc.UpdateWeights("not_a_rule", 1.0, 0.0)


def _make_record(m_chi=0.008, E_V1=0.030):
    pz = math.sqrt(max(E_V1 * E_V1 - M_V1 * M_V1, 0.0))
    rec = siren.dataclasses.InteractionRecord()
    rec.signature.primary_type = PT.N4
    rec.signature.target_type = PT.Decay
    rec.signature.secondary_types = [PT.NuLight, PT.Gamma]
    rec.primary_mass = M_V1
    rec.primary_momentum = [E_V1, 0.0, 0.0, pz]
    rec.secondary_masses = [m_chi, m_chi]
    rec.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
    rec.secondary_helicities = [0, 0]
    rec.interaction_vertex = [0.0, 0.0, 0.0]
    rec.primary_initial_position = [0.0, 0.0, 0.0]
    return rec


def test_inner_outer_count_lockstep():
    """Accumulate must bump a nested mixture's count once per outer point, so
    inner and outer divide their accumulators by the same count."""
    inner = siren.injection.MultiChannelPhaseSpace()
    inner.channels = [siren.injection.Isotropic2BodyChannel(0),
                      siren.injection.Isotropic2BodyChannel(0)]
    inner.weights = [0.5, 0.5]
    nested = siren.injection.NestedMixtureChannel(inner)

    outer = siren.injection.MultiChannelPhaseSpace()
    outer.channels = [siren.injection.Isotropic2BodyChannel(0), nested]
    outer.weights = [0.5, 0.5]

    outer.ResetAccumulators()
    rng = siren.utilities.SIREN_random(3)
    n_counted = 0
    for _ in range(500):
        r = _make_record()
        outer.Sample(rng, None, r)
        g = outer.Density(None, r)
        if g <= 0 or not math.isfinite(g):
            continue
        outer.Accumulate(None, r, 1.0)
        n_counted += 1

    assert outer.kp_count == n_counted
    assert inner.kp_count == outer.kp_count
