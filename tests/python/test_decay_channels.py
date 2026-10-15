"""Selection changes the proposal; inclusive widths still set physical rates."""

import math
import pickle

import pytest
import siren
from siren import dataclasses as dc, distributions as dist, injection as inj
from siren import interactions as xs, utilities
from siren.math import Vector3D

P = dc.Particle.ParticleType


def channel(daughters, partial, total):
    signature = dc.InteractionSignature()
    signature.primary_type = P.N4
    signature.target_type = P.Decay
    signature.secondary_types = daughters
    return inj.PhaseSpaceDecay(
        signature, [0.0, 0.0], partial, total, inj.Isotropic2BodyChannel()
    )


def setup(models, selected=None, n=32, energy=2.0, mixture=False):
    detector = siren.detector.DetectorModel()
    collection = xs.InteractionCollection(P.N4, models)
    collection.SetDecayChannels(selected)
    process = inj.PrimaryInjectionProcess(P.N4, collection)
    process.weighting_mode = inj.VertexWeightingMode.Fixed()
    process.distributions = [
        dist.PrimaryMass(1.0),
        dist.Monoenergetic(energy),
        dist.PrimaryNeutrinoHelicityDistribution(),
        dist.FixedDirection(Vector3D(0, 0, 1)),
        dist.SphereVolumePositionDistribution(siren.geometry.Sphere(1.0, 0.0)),
    ]
    if mixture:
        for model in models:
            signature = model.GetPossibleSignatures()[0]
            process.SetPhaseSpace(
                signature,
                inj.MultiChannelPhaseSpace(
                    [
                        inj.PhysicalDecayChannel(model, signature),
                        inj.Isotropic2BodyChannel(),
                    ],
                    [0.3, 0.7],
                ),
            )
    generator = inj._Injector(n, detector, process, utilities.SIREN_random(752))
    physical = inj.PhysicalProcess(P.N4, collection)
    physical.weighting_mode = inj.VertexWeightingMode.Fixed()
    physical.distributions = process.distributions
    return generator, inj._Weighter([generator], detector, physical), collection


@pytest.mark.parametrize("energy", [1.0, 2.0])
@pytest.mark.parametrize("width", [1e-16, 1e-280])
@pytest.mark.parametrize("mixture", [False, True])
def test_selected_channel_keeps_parent_lifetime_and_branch(energy, width, mixture):
    a = channel([P.NuMu, P.Gamma], 0.2 * width, width)
    b = channel([P.NuE, P.Gamma], 0.3 * width, width)
    signature = a.GetPossibleSignatures()[0]
    generator, weighter, collection = setup(
        [a, b], [signature], energy=energy, mixture=mixture
    )
    events = [generator.GenerateEvent() for _ in range(32)]
    assert generator.InjectionAttempts() == 32
    for event in events:
        assert len(event.tree) == 1
        record = event.tree[0].record
        assert record.signature == signature
        assert collection.TotalDecayWidthAllFinalStates(record) == width
        assert weighter.EventWeight(event) == pytest.approx(0.2 / 32, rel=2e-12, abs=0)
        assert weighter.EventWeightWithBreakdown(event).total == weighter.EventWeight(
            event
        )
        record.signature = b.GetPossibleSignatures()[0]
        process = generator.GetPrimaryProcess()
        per_process = inj.PrimaryProcessWeighter(
            weighter.GetPrimaryPhysicalProcess(), process, generator.GetDetectorModel()
        )
        datum = dc.InteractionTreeDatum(record)
        assert per_process.GenerationProbability(datum) == 0.0


def test_unrestricted_generation_conditions_on_all_modeled_channels():
    a = channel([P.NuMu, P.Gamma], 0.2e-16, 1e-16)
    b = channel([P.NuE, P.Gamma], 0.3e-16, 1e-16)
    generator, weighter, _ = setup([a, b], n=256)
    events = [generator.GenerateEvent() for _ in range(256)]
    seen = {tuple(event.tree[0].record.signature.secondary_types) for event in events}
    assert len(seen) == 2
    for event in events:
        assert weighter.EventWeight(event) == pytest.approx(0.5 / 256, rel=2e-12, abs=0)


def test_selection_updates_validate_before_replacing_live_state():
    model = channel([P.NuMu, P.Gamma], 1e-16, 1e-16)
    signature = model.GetPossibleSignatures()[0]
    _, _, collection = setup([model], [signature])
    foreign = channel([P.NuE, P.Gamma], 1e-16, 1e-16).GetPossibleSignatures()[0]
    for invalid in ([], [signature, signature], [foreign]):
        with pytest.raises(RuntimeError, match="decay_channels"):
            collection.SetDecayChannels(invalid)
        assert collection.GetDecayChannels() == [signature]
    collection.SetDecayChannels(None)
    assert not collection.HasDecayChannels()


def test_selected_native_archive_and_pickle_continue(tmp_path):
    from test_native_phase_space import record_state

    a = channel([P.NuMu, P.Gamma], 0.2e-16, 1e-16)
    b = channel([P.NuE, P.Gamma], 0.3e-16, 1e-16)
    signature = a.GetPossibleSignatures()[0]
    generator, weighter, _ = setup([a, b], [signature], mixture=True)
    for _ in range(3):
        generator.GenerateEvent()
    path = str(tmp_path / "selected.injector")
    generator.SaveInjector(path)
    restored = inj._Injector(32, path, utilities.SIREN_random(1))
    copied = pickle.loads(pickle.dumps(generator))
    weights_path = str(tmp_path / "selected")
    weighter.SaveWeighter(weights_path)
    loaded_weighter = inj._Weighter([], weights_path)
    for _ in range(4):
        original = generator.GenerateEvent()
        assert (
            record_state(original)
            == record_state(restored.GenerateEvent())
            == record_state(copied.GenerateEvent())
        )
        assert loaded_weighter.EventWeight(original) == pytest.approx(
            0.2 / 3, rel=2e-12, abs=0
        )
    assert restored.GetPrimaryProcess().interactions.GetDecayChannels() == [signature]


class WidthOverride(xs.Decay):
    def __init__(self, model, partial=None, owned=None, parent=0.0):
        super().__init__()
        self.model, self.partial, self.owned, self.parent = (
            model,
            partial,
            owned,
            parent,
        )

    def TotalDecayWidth(self, record):
        return (
            self.model.TotalDecayWidth(record) if self.partial is None else self.partial
        )

    def TotalDecayWidthAllFinalStates(self, record):
        return (
            self.model.TotalDecayWidthAllFinalStates(record)
            if self.owned is None
            else self.owned
        )

    def ParentDecayWidth(self, record):
        return self.parent

    def GetPossibleSignatures(self):
        return self.model.GetPossibleSignatures()

    def GetPossibleSignaturesFromParent(self, parent):
        return self.model.GetPossibleSignaturesFromParent(parent)

    def SecondaryMasses(self, types):
        return self.model.SecondaryMasses(types)

    def Measure(self):
        return self.model.Measure()

    def Topology(self):
        return self.model.Topology()

    def FinalStateProbability(self, record):
        return self.model.FinalStateProbability(record)

    def DifferentialDecayWidth(self, record):
        return self.TotalDecayWidth(record) * self.FinalStateProbability(record)

    def SampleFinalState(self, record, random):
        return self.model.SampleFinalState(record, random)

    def DensityVariables(self):
        return self.model.DensityVariables()

    def equal(self, other):
        return self is other


@pytest.mark.parametrize("bad", [-1.0, math.inf, math.nan])
@pytest.mark.parametrize("field", ["partial", "owned"])
@pytest.mark.parametrize("selected", [False, True])
def test_invalid_width_is_visible_even_for_unselected_channel(bad, field, selected):
    a = channel([P.NuMu, P.Gamma], 0.2e-16, 1e-16)
    b = WidthOverride(channel([P.NuE, P.Gamma], 0.3e-16, 1e-16), **{field: bad})
    selected_model = b if selected else a
    generator, _, _ = setup([a, b], selected_model.GetPossibleSignatures())
    with pytest.raises(RuntimeError, match="finite|nonnegative|width"):
        generator.GenerateEvent()


def test_python_parent_width_and_legacy_owned_width_share_one_lifetime():
    a = WidthOverride(channel([P.NuMu, P.Gamma], 0.1e-16, 1e-16), parent=1e-16)
    b = WidthOverride(channel([P.NuE, P.Gamma], 0.9e-16, 1e-16))
    generator, weighter, collection = setup([a, b], a.GetPossibleSignatures())
    events = [generator.GenerateEvent() for _ in range(32)]
    assert collection.ParentDecayWidth(events[0].tree[0].record) == 1e-16
    for event in events:
        assert weighter.EventWeight(event) == pytest.approx(0.1 / 32, rel=2e-12, abs=0)


def test_darknews_trampoline_preserves_explicit_parent_width():
    class RareDecay(xs.DarkNewsDecay):
        def ParentDecayWidth(self, record):
            return 1e-16

        def TotalDecayWidthAllFinalStates(self, record):
            return 0.1e-16

    model = RareDecay()
    collection = xs.InteractionCollection(P.N4, [model])
    record = dc.InteractionRecord()
    record.signature.primary_type = P.N4
    assert collection.ParentDecayWidth(record) == 1e-16
    assert collection.TotalDecayWidthAllFinalStates(record) == 1e-16


def test_legacy_moving_fixed_process_does_not_query_unrequested_lifetime():
    class LegacyDecay(WidthOverride):
        def TotalDecayWidthAllFinalStates(self, record):
            raise AssertionError("No lifetime was requested")

    model = LegacyDecay(channel([P.NuMu, P.Gamma], 1e-16, 1e-16))
    generator, weighter, _ = setup([model])
    events = [generator.GenerateEvent() for _ in range(32)]
    for event in events:
        assert weighter.EventWeight(event) == pytest.approx(1 / 32, rel=2e-12, abs=0)


def test_native_wrapper_cannot_silently_archive_live_python_physics():
    model = WidthOverride(channel([P.NuMu, P.Gamma], 1e-16, 1e-16))
    native = inj.PhaseSpaceDecay(
        model.GetPossibleSignatures()[0],
        [0.0, 0.0],
        1e-16,
        1e-16,
        inj.PhysicalDecayChannel(model),
    )
    with pytest.raises(RuntimeError):
        pickle.dumps(native)


def test_rejected_parent_change_preserves_process_and_collection():
    model = channel([P.NuMu, P.Gamma], 1e-16, 1e-16)
    generator, _, collection = setup([model], model.GetPossibleSignatures())
    process = generator.GetPrimaryProcess()
    with pytest.raises(RuntimeError, match="decay_channels"):
        process.primary_type = P.Pi0
    assert process.primary_type == collection.GetPrimaryType() == P.N4


def test_native_decay_equality_is_conservative_for_live_python_channels():
    model = WidthOverride(channel([P.NuMu, P.Gamma], 1e-16, 1e-16))
    signature = model.GetPossibleSignatures()[0]
    physical = inj.PhysicalDecayChannel(model)
    first = inj.PhaseSpaceDecay(signature, [0.0, 0.0], 1e-16, 1e-16, physical)
    shared = inj.PhaseSpaceDecay(signature, [0.0, 0.0], 1e-16, 1e-16, physical)
    separate = inj.PhaseSpaceDecay(
        signature, [0.0, 0.0], 1e-16, 1e-16, inj.PhysicalDecayChannel(model))
    assert first == first
    assert first == shared
    assert first != separate
    with pytest.raises((RuntimeError, TypeError)):
        pickle.dumps(first)
