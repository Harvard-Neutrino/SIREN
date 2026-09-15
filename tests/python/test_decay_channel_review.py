"""Review witnesses for validation, stationary parents, and configuration changes."""

import math

import numpy as np
import pytest
import siren
from test_decay_channels import (
    HBARC,
    MOMENTUM,
    PT,
    PartialDecay,
    d,
    generate,
    generation_probability,
    injector,
    models,
    primary_distributions,
    record,
)


@pytest.mark.parametrize("fixed", [False, True])
@pytest.mark.parametrize("width", [-0.8, float("nan"), float("inf")])
@pytest.mark.parametrize("separate_total", [False, True])
def test_invalid_unselected_widths_rejected_before_generation(
    fixed, width, separate_total
):
    class SeparateTotal(PartialDecay):
        def TotalDecayWidthAllFinalStates(self, record):
            return MOMENTUM * HBARC

    ms = models()
    if separate_total:
        ms[1] = SeparateTotal(0.8, (PT.MuMinus, PT.MuPlus))
    ms[1].width = width * MOMENTUM * HBARC
    inj = injector(ms, fixed=fixed)
    with pytest.raises(siren.errors.ConfigurationError, match="width"):
        generate(inj, 1)
    assert inj.injected_events == 0
    with pytest.raises(siren.errors.ConfigurationError, match="width"):
        generation_probability(inj, record(ms[0]))


def test_invalid_all_final_state_width_is_checked_separately():
    class BadTotal(PartialDecay):
        def TotalDecayWidthAllFinalStates(self, record):
            return float("inf")

    ms = [PartialDecay(0.2), BadTotal(0.8, (PT.MuMinus, PT.MuPlus))]
    inj = injector(ms)
    with pytest.raises(siren.errors.ConfigurationError, match="width"):
        generate(inj, 1)


def rest_injector(ms, selection, kinematics=None):
    source = d.PrimaryExternalDistribution(
        ["E", "px", "py", "pz", "x", "y", "z", "m"], [[1, 0, 0, 0, 0, 0, 0, 1]]
    )
    vertex = siren.Vertex(
        PT.N4,
        ms,
        decay_channels=selection,
        distributions=[source],
        weighting=siren.Fixed(),
        kinematics=kinematics,
    )
    return siren.injection.Injector(
        detector=siren.detector.DetectorModel(), primary=vertex, events=32, seed=71
    )


@pytest.mark.parametrize(
    "model_count,biased", [(1, False), (1, True), (2, False), (2, True), (3, False)]
)
@pytest.mark.parametrize("selection", ["none", "first", "all"])
def test_fixed_decay_at_rest_has_correct_branching_and_shape(
    model_count, selection, biased
):
    ms = (models() + [PartialDecay(0.5)])[:model_count]
    selected = None if selection == "none" else [ms[0]]
    if selection == "all":
        # The third model shares the first model's observable signature.
        selected = [m.GetPossibleSignatures()[0] for m in ms[:2]]
    target = siren.geometry.Box(widths=[1, 1, 1], center=[0, 0, 4])
    kin = siren.channels.toward(0, target, fraction=0.6) if biased else None
    inj = rest_injector(ms, selected, kin)
    trees = generate(inj, 32)
    weighter = siren.injection.Weighter(inj)
    full_width = sum(m.width for m in ms)
    allowed_width = (
        sum(m.width for m in ms if m.daughters == ms[0].daughters)
        if selection == "first"
        else full_width
    )
    process = inj.engine.GetPrimaryProcess()
    for tree, weight in zip(trees, weighter.weight_all(trees)):
        r = tree.tree[0].record
        ps = process.GetPhaseSpace(r.signature)
        proposal = ps.Density(inj.detector_model, r) if ps else 1 / (4 * math.pi)
        expected = allowed_width / full_width / (4 * math.pi * proposal)
        assert weight * inj.injection_attempts == pytest.approx(
            expected, rel=2e-12, abs=0
        )
        assert weighter.explain(tree).total == pytest.approx(weight, rel=2e-12, abs=0)
    if model_count == 2 and selection == "first":
        physical = [PartialDecay(0.6), PartialDecay(1.4, (PT.MuMinus, PT.MuPlus))]
        reweighter = siren.injection.Weighter(
            inj, overrides={"primary_interactions": physical}
        )
        for tree, weight in zip(trees, reweighter.weight_all(trees)):
            r = tree.tree[0].record
            ps = process.GetPhaseSpace(r.signature)
            proposal = ps.Density(inj.detector_model, r) if ps else 1 / (4 * math.pi)
            assert weight * inj.injection_attempts == pytest.approx(
                0.3 / (4 * math.pi * proposal), rel=2e-12, abs=0
            )


@pytest.mark.parametrize("other_width", [0.0, 0.8])
def test_closed_selected_rest_decay_retains_no_support_failure(other_width):
    ms = [PartialDecay(0), PartialDecay(other_width, (PT.MuMinus, PT.MuPlus))]
    inj = rest_injector(ms, ms[0])
    with pytest.raises(siren.errors.GenerationFailure):
        generate(inj, 1)
    assert inj.injected_events == 0
    r = record(ms[0])
    r.primary_momentum = [1, 0, 0, 0]
    assert generation_probability(inj, r) == 0


class N5Decay(PartialDecay):
    parent = PT.N5


def test_secondary_replacement_validates_all_before_updating_any():
    a = models()
    b = [N5Decay(0.2), N5Decay(0.8, (PT.MuMinus, PT.MuPlus))]
    sigs = {
        PT.N4: [a[0].GetPossibleSignatures()[0]],
        PT.N5: [b[0].GetPossibleSignatures()[0]],
    }
    inj = siren.injection.Injector(
        detector=siren.detector.DetectorModel(),
        events=4,
        seed=71,
        primary_type=PT.N4,
        primary_interactions=a,
        primary_injection_distributions=primary_distributions(True),
        primary_weighting_mode=siren.Fixed(),
        secondary_interactions={PT.N4: a, PT.N5: b},
        secondary_injection_distributions={
            pt: [d.SecondaryBoundedVertexDistribution(2)] for pt in sigs
        },
        secondary_decay_channels=sigs,
        stopping_condition=lambda *args: True,
    )
    processes = inj.engine.GetSecondaryProcessMap()
    original = {pt: p.interactions for pt, p in processes.items()}
    replacements = {PT.N4: models(scale=2), PT.N5: [b[1]]}
    with pytest.raises(siren.errors.ConfigurationError, match="absent"):
        inj.secondary_interactions = replacements
    for pt, p in processes.items():
        assert p.interactions is original[pt]
    assert inj.secondary_interactions == {PT.N4: a, PT.N5: b}
    replacements[PT.N5] = [N5Decay(0.6), N5Decay(2.4, (PT.MuMinus, PT.MuPlus))]
    inj.secondary_interactions = replacements
    for pt, factor in [(PT.N4, 2), (PT.N5, 3)]:
        r = record(replacements[pt][0])
        assert processes[pt].interactions.TotalDecayWidthAllFinalStates(
            r
        ) == pytest.approx(factor * MOMENTUM * HBARC, rel=1e-14, abs=0)
        assert processes[pt].interactions.GetDecayChannels() == sigs[pt]


def test_selection_order_preserves_equality_samples_weights_and_merge():
    ms = models()
    sigs = [m.GetPossibleSignatures()[0] for m in ms]
    collections, results = [], []
    detector = siren.detector.DetectorModel()
    for selected in [sigs, sigs[::-1]]:
        c = siren.interactions.InteractionCollection(PT.N4, ms)
        c.SetDecayChannels(selected)
        collections.append(c)
        vertex = siren.Vertex(
            PT.N4,
            ms,
            decay_channels=selected,
            distributions=primary_distributions(True),
            weighting=siren.Fixed(),
        )
        results.append(
            siren.Simulation(detector=detector, primary=vertex, events=8, seed=71).run()
        )
    assert collections[0] == collections[1]
    assert collections[0].GetDecayChannels() == collections[1].GetDecayChannels()
    np.testing.assert_array_equal(results[0].weights, results[1].weights)
    for a, b in zip(results[0].events, results[1].events):
        assert a.tree[0].record.signature == b.tree[0].record.signature
        np.testing.assert_array_equal(
            a.tree[0].record.secondary_momenta, b.tree[0].record.secondary_momenta
        )
    merged = siren.Results.merge(results)
    assert len(merged.events) == 16
    assert sum(merged.weights) == pytest.approx(1.0, rel=1e-12)


def test_vertex_selection_assignment_uses_constructor_validation():
    ms = models()
    vertex = siren.Vertex(PT.N4, ms, decay_channels=ms[0])
    vertex.decay_channels = ms[1]
    expected = ms[1].GetPossibleSignatures()
    assert vertex.decay_channels == expected
    for invalid in ([], "unknown", object()):
        with pytest.raises(siren.errors.ConfigurationError):
            vertex.decay_channels = invalid
        assert vertex.decay_channels == expected
    assert vertex.compile(is_primary=True).interactions.GetDecayChannels() == expected


def test_finite_model_widths_cannot_overflow_the_propagation_total():
    ms = models()
    for model in ms:
        model.width = 0.75 * np.finfo(float).max
    inj = injector(ms)
    with pytest.raises(siren.errors.ConfigurationError, match="width"):
        generate(inj, 1)
    assert inj.injected_events == 0


def test_native_multisignature_decays_at_rest_keep_partial_branching():
    decay = siren.interactions.HNLDipoleDecay(
        1.0, [1e-7, 2e-7, 3e-7], siren.interactions.HNLDipoleDecay.Majorana
    )
    for signature in decay.GetPossibleSignaturesFromParent(PT.N4):
        inj = rest_injector([decay], signature)
        trees = generate(inj, 32)
        weighter = siren.injection.Weighter(inj)
        weights = weighter.weight_all(trees)
        r = trees[0].tree[0].record
        expected = decay.TotalDecayWidth(r) / decay.TotalDecayWidthAllFinalStates(r)
        assert all(t.tree[0].record.signature == signature for t in trees)
        assert sum(weights) == pytest.approx(expected, rel=2e-12, abs=0)
