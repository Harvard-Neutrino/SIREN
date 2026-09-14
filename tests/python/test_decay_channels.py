"""Forced decay signatures retain the complete model's lifetime and branching."""

import math
import pickle

import numpy as np
import pytest

import siren
from siren import distributions as d

PT = siren.dataclasses.ParticleType
HBARC = siren.utilities.Constants.hbarc
MASS = 1.0
MOMENTUM = 3.0
LENGTH = 2.0


class PartialDecay(siren.DecayModel):
    parent = PT.N4
    measure = siren.Measure.SolidAngleRest()

    def __init__(self, fraction, daughters=(PT.EMinus, PT.EPlus), scale=1.0):
        super().__init__()
        self.width = fraction * MOMENTUM * HBARC * scale
        self.daughters = daughters

    def total_width(self):
        return self.width

    def differential_width(self, record):
        return self.width / (4 * math.pi)

    def sample(self, record, random):
        self.sample_isotropic(record, random)


class UniformFlight(d.VertexPositionDistribution):
    """Independent position proposal with a known density and finite bounds."""

    def SamplePosition(self, random, detector, interactions, record):
        return siren.math.Vector3D(), siren.math.Vector3D(0, 0, LENGTH * random.Uniform(0, 1))

    def InjectionBounds(self, detector, interactions, record):
        return siren.math.Vector3D(), siren.math.Vector3D(0, 0, LENGTH)

    def GenerationProbability(self, detector, interactions, record):
        return 1 / LENGTH if 0 <= record.interaction_vertex[2] <= LENGTH else 0

    def Name(self):
        return "UniformFlight"

    def clone(self):
        return UniformFlight()

    def equal(self, other):
        return isinstance(other, UniformFlight)

    def less(self, other):
        return False


def models(scale=1.0):
    return [PartialDecay(0.2, scale=scale),
            PartialDecay(0.8, (PT.MuMinus, PT.MuPlus), scale=scale)]


def record(model):
    r = siren.dataclasses.InteractionRecord()
    r.signature = model.GetPossibleSignatures()[0]
    r.primary_mass = MASS
    r.primary_momentum = [math.sqrt(MASS**2 + MOMENTUM**2), 0, 0, MOMENTUM]
    r.interaction_vertex = [0, 0, 0.7]
    return r


def primary_distributions(fixed):
    if fixed:
        return [d.PrimaryExternalDistribution(
            ["E", "px", "py", "pz", "x", "y", "z", "m"],
            [[math.sqrt(10), 0, 0, 3, 0, 0, 0, 1]])]
    return [d.PrimaryMass(MASS), d.Monoenergetic(math.sqrt(10)),
            d.FixedDirection(siren.math.Vector3D(0, 0, 1)),
            d.PrimaryNeutrinoHelicityDistribution(), UniformFlight()]


def injector(ms, *, selected=True, fixed=True, kinematics=None, events=64, seed=71):
    kwargs = {"decay_channels": [ms[0]]} if selected else {}
    v = siren.Vertex(PT.N4, ms, distributions=primary_distributions(fixed),
                     weighting=siren.Fixed() if fixed else siren.Propagated(),
                     kinematics=kinematics, **kwargs)
    return siren.injection.Injector(detector=siren.detector.DetectorModel(),
                                    primary=v, events=events, seed=seed)


def generate(inj, n=64):
    return inj.generate(n, on_failure="raise", on_shortfall="raise")


def generation_probability(inj, r):
    process = inj.engine.GetPrimaryProcess()
    physical = siren.injection.PhysicalProcess(process.primary_type, process.interactions)
    datum = siren.dataclasses.InteractionTreeDatum(r)
    return siren.injection.PrimaryProcessWeighter(
        physical, process, inj.detector_model).GenerationProbability(datum)


def test_collection_separates_channel_width_from_lifetime():
    ms = models()
    collection = siren.interactions.InteractionCollection(PT.N4, ms)
    r = record(ms[0])
    width = sum(m.total_width() for m in ms)
    assert collection.TotalDecayWidthAllFinalStates(r) == pytest.approx(width, rel=1e-14)
    assert collection.TotalDecayLengthAllFinalStates(r) == pytest.approx(1, rel=1e-14)
    collection.SetDecayChannels([r.signature])
    assert collection.TotalDecayWidthAllFinalStates(r) == pytest.approx(width, rel=1e-14)
    assert collection.TotalDecayLengthAllFinalStates(r) == pytest.approx(1, rel=1e-14)
    assert collection.GetDecays() == ms
    assert not collection.AllowsDecay(ms[1].GetPossibleSignatures()[0])
    collection.SetDecayChannels(None)
    assert not collection.HasDecayChannels()


@pytest.mark.parametrize("fixed", [False, True])
@pytest.mark.parametrize("kinematics", [False, True])
def test_forced_channel_weights_and_lifetime_reweighting(fixed, kinematics):
    ms = models()
    inj = injector(ms, fixed=fixed,
                   kinematics=siren.channels.isotropic(0) if kinematics else None)
    trees = generate(inj)
    weighter = siren.injection.Weighter(inj)
    weights = weighter.weight_all(trees)
    for tree, weight in zip(trees, weights):
        r = tree.tree[0].record
        assert r.signature == ms[0].GetPossibleSignatures()[0]
        expected = 0.2 if fixed else LENGTH * 0.2 * math.exp(-r.interaction_vertex[2])
        assert weight * inj.injection_attempts == pytest.approx(expected, rel=2e-12, abs=0)
        assert weighter.explain(tree).total == pytest.approx(weight, rel=2e-12, abs=0)

    # Change both the selected partial width and the full physical lifetime.
    target = [PartialDecay(0.6), PartialDecay(1.4, (PT.MuMinus, PT.MuPlus))]
    reweighter = siren.injection.Weighter(inj, overrides={"primary_interactions": target})
    for tree, weight in zip(trees, reweighter.weight_all(trees)):
        z = tree.tree[0].record.interaction_vertex[2]
        expected = 0.3 if fixed else LENGTH * 0.6 * math.exp(-2 * z)
        assert weight * inj.injection_attempts == pytest.approx(expected, rel=2e-12, abs=0)


def test_unrestricted_competition_and_forced_channel_have_same_absolute_yield():
    ms = models()
    n = 4096
    inj = injector(ms, selected=False, events=n)
    trees = generate(inj, n)
    weights = siren.injection.Weighter(inj).weight_all(trees)
    np.testing.assert_allclose(weights, np.full(n, 1 / n), rtol=2e-12, atol=0)
    count = sum(t.tree[0].record.signature == ms[0].GetPossibleSignatures()[0] for t in trees)
    assert abs(count - n * 0.2) < 5 * math.sqrt(n * 0.2 * 0.8)
    forced = injector(ms)
    forced_trees = generate(forced)
    assert sum(siren.injection.Weighter(forced).weight_all(forced_trees)) == pytest.approx(0.2, rel=2e-12)


def test_same_signature_contributions_are_selected_together():
    ms = models() + [PartialDecay(0.5)]
    inj = injector(ms)
    trees = generate(inj)
    assert sum(siren.injection.Weighter(inj).weight_all(trees)) == pytest.approx(0.7 / 1.5, rel=2e-12)


def test_forced_channel_with_biased_kinematics_includes_both_proposal_factors():
    target = siren.geometry.Box(widths=[1, 1, 1], center=[0, 0, 4])
    ms = models()
    inj = injector(ms, kinematics=siren.channels.toward(0, target, fraction=0.6))
    trees = generate(inj)
    weighter = siren.injection.Weighter(inj)
    process = inj.engine.GetPrimaryProcess()
    for tree, weight in zip(trees, weighter.weight_all(trees)):
        r = tree.tree[0].record
        proposal = process.GetPhaseSpace(r.signature).Density(inj.detector_model, r)
        # The analytic physical density is isotropic with branching 0.2.
        expected = 0.2 / (4 * math.pi * proposal)
        assert weight * inj.injection_attempts == pytest.approx(expected, rel=2e-12, abs=0)
        assert weighter.explain(tree).total == pytest.approx(weight, rel=2e-12, abs=0)


def test_lower_level_signature_selection_and_model_replacement():
    ms = models()
    sig = ms[0].GetPossibleSignatures()[0]
    inj = siren.injection.Injector(
        detector=siren.detector.DetectorModel(), events=8, seed=51,
        primary_type=PT.N4, primary_interactions=ms,
        primary_injection_distributions=primary_distributions(True),
        primary_weighting_mode=siren.Fixed(), primary_decay_channels=[sig])
    trees = generate(inj, 8)
    assert all(t.tree[0].record.signature == sig for t in trees)
    inj.primary_interactions = models(scale=2)
    assert inj.engine.GetPrimaryProcess().interactions.GetDecayChannels() == [sig]
    with pytest.raises(siren.errors.ConfigurationError, match="absent"):
        inj.primary_interactions = [ms[1]]
    assert len(inj.engine.GetPrimaryProcess().interactions.GetDecays()) == 2
    with pytest.raises(siren.errors.ConfigurationError, match="absent"):
        inj.primary_type = PT.N5
    assert inj.primary_type == PT.N4
    assert inj.engine.GetPrimaryProcess().primary_type == PT.N4


@pytest.mark.parametrize("selection", [[], "unknown", "foreign", "duplicate"])
def test_invalid_selections_fail_without_changing_existing_configuration(selection):
    ms = models()
    c = siren.interactions.InteractionCollection(PT.N4, ms)
    sig = ms[0].GetPossibleSignatures()[0]
    c.SetDecayChannels([sig])
    if selection == "duplicate":
        bad = [sig, sig]
    elif selection == "unknown":
        bad_sig = ms[1].GetPossibleSignatures()[0]
        bad_sig.secondary_types = [PT.Gamma, PT.Gamma]
        bad = [bad_sig]
    elif selection == "foreign":
        bad_sig = ms[1].GetPossibleSignatures()[0]
        bad_sig.primary_type = PT.N5
        bad = [bad_sig]
    else:
        bad = selection
    with pytest.raises(siren.errors.ConfigurationError, match="decay_channels"):
        c.SetDecayChannels(bad)
    assert c.GetDecayChannels() == [sig]
    with pytest.raises(siren.errors.ConfigurationError):
        c.SetPrimaryType(PT.N5)
    assert c.GetPrimaryType() == PT.N4


def test_closed_selected_channel_has_no_support_or_substitute_events():
    ms = models()
    ms[0].width = 0
    inj = injector(ms)
    with pytest.raises(siren.errors.GenerationFailure):
        generate(inj, 1)
    assert inj.injection_attempts == 1
    assert inj.injected_events == 0
    assert generation_probability(inj, record(ms[0])) == 0


def test_records_outside_selection_have_zero_generation_density():
    ms = models()
    inj = injector(ms)
    assert generation_probability(inj, record(ms[1])) == 0


@pytest.mark.parametrize("width", [-1.0, float("nan"), float("inf")])
def test_invalid_selected_rates_fail_loudly(width):
    ms = models()
    ms[0].width = width
    inj = injector(ms)
    with pytest.raises(siren.errors.ConfigurationError, match="rates|rate"):
        generate(inj, 1)
    assert inj.injected_events == 0


class FeedDown(PartialDecay):
    parent = PT.N5

    def SecondaryMasses(self, types):
        return [1.0, 0.0]


@pytest.mark.parametrize("facade", [False, True])
def test_secondary_selection_preserves_full_lifetime_and_branching(facade):
    source = d.PrimaryExternalDistribution(
        ["E", "px", "py", "pz", "x", "y", "z", "m"],
        [[math.sqrt(13), 0, 0, 3, 0, 0, 0, 2]])
    feed = FeedDown(1.0, (PT.N4, PT.Gamma))
    ms = models()
    primary = siren.Vertex(PT.N5, feed, distributions=[source],
                           weighting=siren.Fixed(), expand=[siren.expand.child(PT.N4, index=0)])
    secondary = siren.Vertex(PT.N4, ms, decay_channels=ms[0],
                             position=d.SecondaryBoundedVertexDistribution(LENGTH),
                             expand=[siren.expand.depth_below(0)],
                             kinematics=siren.channels.isotropic(0))
    kwargs = dict(detector=siren.detector.DetectorModel(), primary=primary,
                  secondaries=[secondary], events=32, seed=122)
    if facade:
        result = siren.Simulation(**kwargs).run(on_failure="raise", on_shortfall="raise")
    else:
        inj = siren.injection.Injector(**kwargs)
        result = siren.generate(inj, siren.injection.Weighter(inj), events=32,
                                on_failure="raise", on_shortfall="raise")
    for tree, weight in zip(result.events, result.weights):
        assert len(tree.tree) == 2
        r = tree.tree[1].record
        assert r.signature == ms[0].GetPossibleSignatures()[0]
        lab_length = np.linalg.norm(r.primary_momentum[1:]) / MOMENTUM
        expected = 0.2 * -math.expm1(-LENGTH / lab_length)
        assert weight * result.attempts == pytest.approx(expected, rel=3e-12, abs=0)
    if not facade:
        inj.secondary_interactions = {PT.N4: models(scale=2)}
        collection = inj.engine.GetSecondaryProcessMap()[PT.N4].interactions
        assert collection.GetDecayChannels() == ms[0].GetPossibleSignatures()
        assert collection.TotalDecayLengthAllFinalStates(record(ms[0])) == pytest.approx(0.5, rel=1e-14)


def test_vertex_compiles_and_expands_only_selected_decay_signatures():
    ms = models()
    v = siren.Vertex(PT.N4, ms, decay_channels=ms[0], kinematics=siren.channels.isotropic(0))
    for primary in (False, True):
        process = v.compile(is_primary=primary)
        assert list(process.GetPhaseSpaceMap()) == ms[0].GetPossibleSignatures()
        assert process.interactions.GetDecays() == ms
    assert set(v.as_vertex_spec().secondary_types) == set(ms[0].GetPossibleSignatures()[0].secondary_types)


def test_simulation_preserves_selection_and_results_reject_different_selections():
    ms = models()
    def run(selection):
        v = siren.Vertex(PT.N4, ms, decay_channels=selection,
                         distributions=primary_distributions(True), weighting=siren.Fixed(),
                         kinematics=siren.channels.isotropic(0))
        sim = siren.Simulation(detector=siren.detector.DetectorModel(), primary=v, events=16, seed=15)
        return sim.run()
    selected = run(ms[0])
    assert sum(selected.weights) == pytest.approx(0.2, rel=2e-12)
    full = run(None)
    with pytest.raises(siren.errors.ConfigurationError, match="configuration|config|compatible"):
        siren.Results.merge([selected, full])


def native_injector(selected=True):
    decay = siren.interactions.HNLDipoleDecay(1.0, [1e-7, 2e-7, 3e-7],
                                            siren.interactions.HNLDipoleDecay.Majorana)
    signatures = decay.GetPossibleSignaturesFromParent(PT.N4)
    kwargs = {"decay_channels": [signatures[0]]} if selected else {}
    dists = primary_distributions(False)[:-1] + [d.SphereVolumePositionDistribution(
        siren.geometry.Sphere(siren.geometry.Placement(), 1.0, 0.0))]
    v = siren.Vertex(PT.N4, decay, distributions=dists,
                     weighting=siren.Fixed(), kinematics=siren.channels.isotropic(0), **kwargs)
    return siren.injection.Injector(detector=siren.detector.DetectorModel(),
                                    primary=v, events=64, seed=91)


def test_save_keeps_selection_models_and_resumable_rng(tmp_path):
    inj = native_injector()
    generate(inj, 3)
    path = str(tmp_path / "selected")
    inj.save(path)
    restored = siren.injection.Injector.load(path)
    selected = inj.engine.GetPrimaryProcess().interactions.GetDecayChannels()
    assert restored.engine.GetPrimaryProcess().interactions.GetDecayChannels() == selected
    expected, actual = generate(inj, 6), generate(restored, 6)
    for a, b in zip(expected, actual):
        ar, br = a.tree[0].record, b.tree[0].record
        assert ar.signature == br.signature
        np.testing.assert_array_equal(ar.primary_momentum, br.primary_momentum)
        np.testing.assert_array_equal(ar.secondary_momenta, br.secondary_momenta)
        np.testing.assert_array_equal(ar.interaction_vertex, br.interaction_vertex)
    w = siren.injection.Weighter(restored)
    before = w.weight_all(actual)
    path = str(tmp_path / "weighter")
    w.save(path)
    loaded = siren.injection.Weighter()
    loaded.load(path)
    np.testing.assert_allclose(loaded.weight_all(actual), before, rtol=2e-12, atol=0)
    restored.reset(events=4)
    assert all(t.tree[0].record.signature == selected[0] for t in generate(restored, 4))


def test_pickle_keeps_selection_with_registered_native_distributions():
    inj = native_injector()
    # PointSourcePositionDistribution has the existing portable-binary archive
    # registration. Keep this a configuration round trip, without sampling its
    # material-bounded position distribution as a vacuum decay proposal.
    inj.primary.distributions = primary_distributions(False)[:-1] + [
        d.PointSourcePositionDistribution(siren.math.Vector3D(), LENGTH)]
    original = inj.engine.GetPrimaryProcess().interactions
    restored = pickle.loads(pickle.dumps(inj))
    collection = restored.engine.GetPrimaryProcess().interactions
    assert collection.GetDecayChannels() == original.GetDecayChannels()
    r = record(models()[0])
    assert collection.TotalDecayWidthAllFinalStates(r) == original.TotalDecayWidthAllFinalStates(r)
