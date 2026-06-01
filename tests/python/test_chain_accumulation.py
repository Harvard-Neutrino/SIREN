"""Phase C: generation-time accumulation pieces.

The chain optimizer now feeds the mixtures' accumulators in C++ instead of
re-walking trees in Python.  These tests pin the two genuinely-new C++ pieces:

1. the chain failure penalty folded into UpdateWeights (W_i inflated by
   1/(1-f_i) from per-channel success/failure selection sums), and
2. AccumulateSelection accumulating the per-channel selection probability
   p_i = alpha_i*g_i/g into the success or failure vector.

The reference Python `_kp_update` is reused as the oracle for (1).
"""
import math

import pytest

import siren
from siren.optimize import _kp_update

PT = siren.dataclasses.ParticleType
M_V1 = 0.017


def _iso_mixture(weights):
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [siren.injection.Isotropic2BodyChannel(0) for _ in weights]
    mc.weights = list(weights)
    return mc


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


def test_updateweights_failure_penalty_inflates_W():
    """With failure-selection data present, UpdateWeights inflates W_i by
    1/(1-f_i) before the KP step -- exactly the chain failure penalty."""
    weights = [0.5, 0.5]
    W_bare = [1.0, 1.0]
    count = 100
    succ = [0.8, 0.2]
    fail = [0.2, 0.8]  # channel 1 fails much more: f = [0.2, 0.8]

    mc = _iso_mixture(weights)
    mc.kp_accumulator = [W_bare[i] * count for i in range(2)]
    mc.kp_count = count
    mc.kp_succ_select = succ
    mc.kp_fail_select = fail
    mc.UpdateWeights("alpha_sqrt_W", 1.0, 0.0)  # damping=1 -> result == _kp_update

    W_inflated = []
    for i in range(2):
        f_i = fail[i] / (succ[i] + fail[i])
        W_inflated.append(W_bare[i] / (1.0 - f_i) if f_i < 1.0 else W_bare[i])
    expected = _kp_update(list(weights), W_inflated, "alpha_sqrt_W", 0.0)

    assert mc.weights == pytest.approx(expected, abs=1e-12)
    # Inflating W_i raises the channel's KP weight (higher variance contribution
    # draws more samples), so the more-failing channel gets the larger weight --
    # this is the existing behavior the C++ reproduces, not new behavior.
    assert mc.weights[1] > mc.weights[0]


def test_updateweights_no_penalty_without_failure_data():
    """No failure-selection data -> no inflation (single-vertex optimizer path)."""
    weights = [0.5, 0.5]
    W = [4.0, 1.0]
    count = 100
    mc = _iso_mixture(weights)
    mc.kp_accumulator = [W[i] * count for i in range(2)]
    mc.kp_count = count
    mc.UpdateWeights("alpha_sqrt_W", 1.0, 0.0)
    expected = _kp_update(list(weights), W, "alpha_sqrt_W", 0.0)
    assert mc.weights == pytest.approx(expected, abs=1e-12)


def test_accumulate_selection_sums_selection_probability():
    """AccumulateSelection adds p_i = alpha_i*g_i/g (= DensityBreakdown[i]/g) to
    the success or failure vector."""
    target = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 2.0)),
        1.0, 1.0, 1.0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [siren.injection.Isotropic2BodyChannel(0),
                   siren.injection.DetectorDirected2BodyChannel(target, 0)]
    mc.weights = [0.4, 0.6]
    mc.ResetAccumulators()

    rng = siren.utilities.SIREN_random(5)
    r = _make_record()
    mc.Sample(rng, None, r)

    contribs = mc.DensityBreakdown(None, r)
    g = sum(contribs)
    expected_p = [contribs[i] / g for i in range(2)]

    mc.AccumulateSelection(None, r, False)
    assert list(mc.kp_succ_select) == pytest.approx(expected_p, abs=1e-12)
    assert all(v == 0.0 for v in mc.kp_fail_select)

    mc.AccumulateSelection(None, r, True)
    assert list(mc.kp_fail_select) == pytest.approx(expected_p, abs=1e-12)
    # A second success sample doubles the success vector.
    mc.AccumulateSelection(None, r, False)
    assert list(mc.kp_succ_select) == pytest.approx(
        [2.0 * p for p in expected_p], abs=1e-12)


def _make_routing_injector():
    """A minimal C++ injector whose primary process has one >=2-channel mixture
    registered for a single signature (no GenerateEvent needed -- we hand-build
    trees).  Uses siren.injection._Injector, the raw C++ class (the public
    siren.Injector is a higher-level wrapper)."""
    dm = siren.detector.DetectorModel()
    proc = siren.injection.PrimaryInjectionProcess()
    proc.primary_type = PT.N4
    # A vertex distribution is required by SetPrimaryProcess.
    proc.distributions = [siren.distributions.PrimaryPhysicalVertexDistribution()]

    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = PT.N4
    sig.target_type = PT.Decay
    sig.secondary_types = [PT.NuLight, PT.Gamma]

    target = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 2.0)), 1.0, 1.0, 1.0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [siren.injection.Isotropic2BodyChannel(0),
                   siren.injection.DetectorDirected2BodyChannel(target, 0)]
    mc.weights = [0.5, 0.5]
    proc.SetPhaseSpace(sig, mc)

    rng = siren.utilities.SIREN_random(1)
    inj = siren.injection._Injector(10, dm, proc, rng)
    return inj, mc, rng


def test_injector_routing_event_credits_all_datums_selection_once():
    """AccumulateEventToMixtures routes each datum to its mixture and credits
    EVERY matching datum (no break); AccumulateSelectionToMixtures credits at
    most once per tree per mixture (break).  Two same-signature datums in one
    tree expose the asymmetry."""
    inj, mc, rng = _make_routing_injector()

    assert len(inj.GetPhaseSpaces()) == 1

    def sampled_record():
        r = _make_record()
        mc.channels[0].Sample(rng, None, r)  # set on-shell kinematics
        return r

    r1 = sampled_record()
    r2 = sampled_record()
    tree = siren.dataclasses.InteractionTree()
    tree.add_entry(r1, None)
    tree.add_entry(r2, None)
    assert [d.depth() for d in tree.tree] == [0, 0]

    # Weighted accumulation: both datums credited.
    mc.ResetAccumulators()
    w = 2.0
    inj.AccumulateEventToMixtures(tree, w, discount_fallback=False)
    assert mc.kp_count == 2

    def credited(r):
        c = mc.DensityBreakdown(None, r)
        g = sum(c)
        return [w * w * (c[i] / mc.weights[i]) / g for i in range(2)]
    expected = [credited(r1)[i] + credited(r2)[i] for i in range(2)]
    assert list(mc.kp_accumulator) == pytest.approx(expected, abs=1e-12)

    # Selection accumulation: one sample per tree per mixture.
    mc.ResetAccumulators()
    inj.AccumulateSelectionToMixtures(tree, False)
    c = mc.DensityBreakdown(None, r1)
    g = sum(c)
    assert list(mc.kp_succ_select) == pytest.approx(
        [c[i] / g for i in range(2)], abs=1e-12)
    assert all(v == 0.0 for v in mc.kp_fail_select)


def test_injector_routing_recurses_into_nested_group():
    """recurse=True makes AccumulateEventToMixtures descend into a
    NestedMixtureChannel, so the chain optimizer can tune the inner per-target
    weights; recurse=False leaves the group's inner accumulators untouched."""
    dm = siren.detector.DetectorModel()
    proc = siren.injection.PrimaryInjectionProcess()
    proc.primary_type = PT.N4
    proc.distributions = [siren.distributions.PrimaryPhysicalVertexDistribution()]

    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = PT.N4
    sig.target_type = PT.Decay
    sig.secondary_types = [PT.NuLight, PT.Gamma]

    box_a = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 2.0)), 1.0, 1.0, 1.0)
    box_b = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, -2.0)), 1.0, 1.0, 1.0)
    inner = siren.injection.MultiChannelPhaseSpace()
    inner.channels = [siren.injection.DetectorDirected2BodyChannel(box_a, 0),
                      siren.injection.DetectorDirected2BodyChannel(box_b, 0)]
    inner.weights = [0.5, 0.5]
    group = siren.injection.NestedMixtureChannel(inner)

    outer = siren.injection.MultiChannelPhaseSpace()
    outer.channels = [siren.injection.Isotropic2BodyChannel(0), group]
    outer.weights = [0.5, 0.5]
    proc.SetPhaseSpace(sig, outer)

    rng = siren.utilities.SIREN_random(2)
    inj = siren.injection._Injector(10, dm, proc, rng)

    r = _make_record()
    outer.channels[0].Sample(rng, None, r)
    tree = siren.dataclasses.InteractionTree()
    tree.add_entry(r, None)

    # recurse=True credits the outer mixture AND the nested inner group (in
    # lockstep), so the inner per-target weights can be tuned.
    outer.ResetAccumulators(True)
    inj.AccumulateEventToMixtures(tree, 1.5, discount_fallback=False, recurse=True)
    assert outer.kp_count == 1
    assert inner.kp_count == 1
    assert any(v != 0.0 for v in inner.kp_accumulator)

    # recurse=False leaves the nested group untouched (pre-existing behavior).
    outer.ResetAccumulators(True)
    inj.AccumulateEventToMixtures(tree, 1.5, discount_fallback=False, recurse=False)
    assert outer.kp_count == 1
    assert inner.kp_count == 0
