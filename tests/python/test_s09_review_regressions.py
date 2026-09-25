"""Regression coverage for the S09 review's merge and facade failures."""

import copy
import pickle

import numpy as np
import pytest

import siren
from siren import channels, distributions as d, injection
from siren.errors import ConfigurationError, GenerationFailure, NotSerializableError
from siren.tune import tune
from simulation_fixtures import offline_detector
from test_generation_policy import _native_sampler
from test_injector_api import _depth_ge_1, _distributions
from test_simulation_physical_contracts import Factor, simulation
from test_variance_report_integration import _IsoDecay, M_N4


def _box(z=8):
    return siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, z)), 6, 6, 6)


def _vertex(model, kinematics=None):
    return siren.Vertex(
        model.parent, model,
        distributions=[d.PrimaryMass(M_N4), d.PowerLaw(2, .03, .05),
                       d.IsotropicDirection(),
                       d.PointSourcePositionDistribution(siren.math.Vector3D(), 25)],
        kinematics=kinematics)


def _sim(vertex, seed=71, **kwargs):
    return siren.Simulation(events=5, detector=offline_detector('CCM'),
                            primary=vertex, seed=seed, **kwargs)


@pytest.mark.parametrize('secondary', [False, True])
def test_merge_rejects_different_flux_normalizations(secondary):
    a, b = simulation(chain=secondary), simulation(chain=secondary)
    a.run()
    b.run()
    low, high = d.PowerLaw(2, .5, 5), d.PowerLaw(2, .5, 5)
    low.normalization, high.normalization = 1, 10
    assert low == high  # Native equality alone cannot protect pooling.
    if secondary:
        ra = a.reweight(secondary_physical_distributions={siren.particles.NuMu: [low]})
        rb = b.reweight(secondary_physical_distributions={siren.particles.NuMu: [high]})
    else:
        ra = a.reweight(physical_distributions=[low])
        rb = b.reweight(physical_distributions=[high])
    np.testing.assert_allclose(rb.weights, ra.weights * 10, rtol=1e-12, atol=0)
    with pytest.raises(ConfigurationError, match='differing config'):
        siren.Results.merge([ra, rb])


def test_merge_rejects_different_physical_detectors():
    a, b = simulation(chain=False), simulation(chain=False, seed=72)
    ra, original = a.run(), b.run()
    weighter = injection.Weighter(
        b.injector, primary_physical=b.physical_distributions,
        overrides={'detector_model': offline_detector('IceCube')})
    rb = siren.Results(original.events, weighter.weight_all(original.events),
                       original.gen_times, weighter, b.injector,
                       requested=original.requested)
    assert np.any(rb.weights != original.weights)
    with pytest.raises(ConfigurationError, match='differing config'):
        siren.Results.merge([ra, rb])


@pytest.mark.parametrize('nested', [False, True])
def test_merge_compiled_channels_agrees_with_pooled_weighter(nested):
    model, target = _IsoDecay(), _box()
    if nested:
        def group(signature, **kwargs):
            return injection.NestedMixtureChannel(injection.MultiChannelPhaseSpace(
                [injection.Isotropic2BodyChannel(),
                 injection.DetectorDirected2BodyChannel(target)], [.25, .75]))
        directed = channels.Channel(group)
    else:
        directed = channels.toward('NuLight', target)
    vertex = _vertex(model, .1 * channels.physical() + .9 * directed)
    a, b = _sim(vertex), _sim(vertex, seed=72)
    ra, rb = a.run(), b.run()
    if nested:
        group_b = next(iter(b.phase_spaces.values())).channels[1]
        group_b.label = 'different diagnostic label'
        group_b.mixture.kp_count = 77
        group_b.mixture.kp_accumulator = [3, 4]
        rb = siren.Results(rb.events, rb.weights, rb.gen_times, b.weighter, b.injector,
                           requested=rb.requested)
    merged = siren.Results.merge([ra, rb])
    np.testing.assert_array_equal(merged.weights, np.r_[ra.weights, rb.weights] / 2)
    np.testing.assert_allclose(merged.weights, merged._weighter.weight_all(merged.events),
                               rtol=1e-12, atol=0)
    assert merged.explain(len(ra)).total == pytest.approx(merged.weights[len(ra)])


@pytest.mark.parametrize('change', ['geometry', 'weights', 'nested_weights'])
def test_merge_rejects_different_channel_proposals(change):
    model = _IsoDecay()

    def vertex(second):
        target = _box(18 if second and change == 'geometry' else 8)
        if change == 'nested_weights':
            def group(signature, **kwargs):
                return injection.NestedMixtureChannel(injection.MultiChannelPhaseSpace(
                    [injection.Isotropic2BodyChannel(),
                     injection.DetectorDirected2BodyChannel(target)],
                    [.75, .25] if second else [.25, .75]))
            directed = channels.Channel(group)
        else:
            directed = channels.toward('NuLight', target)
        fraction = .2 if second and change == 'weights' else .9
        return _vertex(model, (1 - fraction) * channels.physical() + fraction * directed)

    a, b = _sim(vertex(False)), _sim(vertex(True), seed=72)
    with pytest.raises(ConfigurationError, match='differing config'):
        siren.Results.merge([a.run(), b.run()])


def test_shared_opaque_channel_still_merges():
    class CustomChannel(injection.Isotropic2BodyChannel):
        pass

    raw = CustomChannel()
    vertex = _vertex(_IsoDecay(), channels.Channel(lambda signature, **kwargs: raw))
    a, b = _sim(vertex), _sim(vertex, seed=72)
    assert len(siren.Results.merge([a.run(), b.run()])) == 10


def test_merge_honors_model_equality_across_subclasses_of_one_native_base():
    class EquivalentDecay(_IsoDecay):
        def equal(self, other):
            return isinstance(other, EquivalentDecay)

    class First(EquivalentDecay):
        pass

    class Second(EquivalentDecay):
        pass

    a = _sim(_vertex(First(), channels.physical()))
    b = _sim(_vertex(Second(), channels.physical()), seed=72)
    merged = siren.Results.merge([a.run(), b.run()])
    np.testing.assert_allclose(merged.weights, merged._weighter.weight_all(merged.events),
                               rtol=1e-12, atol=0)


def test_merge_mismatched_native_model_bases_raises_configuration_error():
    class NuMuDecay(_IsoDecay):
        parent = 'NuMu'

    decay = _sim(_vertex(NuMuDecay())).run()
    scatter = simulation(chain=False).run()
    with pytest.raises(ConfigurationError, match='differing config'):
        siren.Results.merge([scatter, decay])


@pytest.mark.parametrize('distinct', [False, True])
def test_simulation_and_vertex_reject_colliding_signatures(distinct):
    model = _IsoDecay()
    vertex = _vertex(model, channels.isotropic())
    vertex.interactions = [model, _IsoDecay() if distinct else model]
    with pytest.raises(ConfigurationError, match='both produce signature'):
        vertex.compile(is_primary=True, detector=offline_detector('CCM'))
    with pytest.raises(ConfigurationError, match='both produce signature'):
        _sim(vertex)


def test_primary_biasing_cannot_replace_vertex_kinematics():
    vertex = _vertex(_IsoDecay(), channels.isotropic())
    with pytest.raises(ConfigurationError, match='Multiple kinematics/biasing'):
        _sim(vertex, biasing=channels.Mixture([(1, channels.physical())]))


@pytest.mark.parametrize('route', ['biasing', 'bias_targets'])
def test_secondary_biasing_cannot_replace_vertex_kinematics(route):
    class ChildDecay(_IsoDecay):
        parent = 'NuLight'
        daughters = ('NuMu', 'Gamma')

    model = ChildDecay()
    secondary = siren.Vertex('NuLight', model, position=d.SecondaryPhysicalVertexDistribution(),
                             kinematics=channels.isotropic())
    extra = ({'biasing': channels.Mixture([(1, channels.physical())])}
             if route == 'biasing' else
             {'bias_targets': {model.GetPossibleSignatures()[0]:
                               injection.MultiChannelPhaseSpace([injection.Isotropic2BodyChannel()])},
              'bias_daughter': 'NuMu'})
    with pytest.raises(ConfigurationError, match='Multiple kinematics/biasing'):
        _sim(_vertex(_IsoDecay()), secondaries=[secondary],
             stopping_condition=_depth_ge_1, **extra)


def test_distinct_primary_and_secondary_biasing_remain_supported():
    class ChildDecay(_IsoDecay):
        parent = 'NuLight'
        daughters = ('NuMu', 'Gamma')

    secondary = siren.Vertex('NuLight', ChildDecay(),
                             position=d.SecondaryPhysicalVertexDistribution())
    sim = _sim(_vertex(_IsoDecay(), channels.isotropic()), secondaries=[secondary],
               stopping_condition=_depth_ge_1,
               biasing=channels.Mixture([(1, channels.physical())]))
    assert len(sim.phase_spaces) == 2
    assert len(sim.injector.engine.GetPrimaryProcess().GetPhaseSpaceMap()) == 1
    assert len(sim.injector.engine.GetSecondaryProcesses()[0].GetPhaseSpaceMap()) == 1


@pytest.mark.parametrize('shared', [False, True])
@pytest.mark.parametrize('route', ['biasing', 'vertices', 'geometry', 'geometry_by_type'])
def test_secondary_biasing_uses_each_types_own_model(shared, route):
    first = siren.interactions.DummyCrossSection()
    models = {siren.particles.NuMu: first,
              siren.particles.NuE: first if shared else siren.interactions.DummyCrossSection()}
    position = d.SecondaryPhysicalVertexDistribution()
    bias = channels.Mixture([(1, channels.physical())])
    extra = {'secondary_interactions': {pt: [model] for pt, model in models.items()},
             'secondary_position': position}
    if route == 'vertices':
        extra = {'secondaries': [siren.Vertex(pt, model, position=position, kinematics=bias)
                                 for pt, model in models.items()]}
    elif route == 'biasing':
        extra['biasing'] = bias
    else:
        target = _box()
        extra['bias_targets'] = ({pt: target for pt in models}
                                 if route == 'geometry_by_type' else target)
        extra['bias_daughter'] = first.GetPossibleSignatures()[0].target_type
    sim = siren.Simulation(
        events=5, detector=offline_detector('CCM'), primary='NuMu', interactions=[first],
        injection_distributions=_distributions(25), stopping_condition=_depth_ge_1,
        seed=71, **extra)
    assert len(first.GetPossibleSignatures()) > len(models)
    assert {sig.primary_type for sig in sim.phase_spaces} == set(models)
    for sig, mixture in sim.phase_spaces.items():
        assert mixture.channels[0].GetCrossSection() is models[sig.primary_type]
    if route in ('biasing', 'vertices'):
        processes = sim.injector.engine.GetSecondaryProcesses()
        assert len(processes) == len(models)
        for process in processes:
            registered = process.GetPhaseSpaceMap()
            assert registered
            assert all(sig.primary_type == process.primary_type for sig in registered)


def test_phase_spaces_prefers_secondary_when_primary_signature_matches():
    model = siren.interactions.DummyCrossSection()
    primary = siren.Vertex(
        'NuMu', model, distributions=_distributions(25),
        kinematics=.2 * channels.physical() + .8 * channels.physical())
    secondary = siren.Vertex(
        'NuMu', model, position=d.SecondaryPhysicalVertexDistribution(),
        kinematics=channels.physical())
    sim = _sim(primary, secondaries=[secondary], stopping_condition=_depth_ge_1)
    engine = sim.injector.engine
    primary_spaces = engine.GetPrimaryProcess().GetPhaseSpaceMap()
    secondary_spaces = engine.GetSecondaryProcesses()[0].GetPhaseSpaceMap()
    assert secondary_spaces
    for sig, secondary_mixture in secondary_spaces.items():
        assert sim.phase_spaces[sig] is secondary_mixture
        assert primary_spaces[sig] is not secondary_mixture
        assert list(primary_spaces[sig].weights) == pytest.approx([.2, .8])
        assert list(secondary_mixture.weights) == [1.]


def test_phase_spaces_includes_primary_and_returns_a_copy():
    sim = _sim(_vertex(_IsoDecay(), channels.isotropic()))
    registered = sim.injector.engine.GetPrimaryProcess().GetPhaseSpaceMap()
    assert sim.phase_spaces == registered
    sim.phase_spaces.clear()
    assert len(sim.phase_spaces) == 1


def test_shallow_copies_preserve_factor_without_invoking_it():
    factor = Factor()
    sim = simulation(factor=factor, chain=False)
    result = sim.run()
    calls = factor.calls
    sim_copy, result_copy = copy.copy(sim), copy.copy(result)
    assert sim_copy is not sim and sim_copy.__dict__ is not sim.__dict__
    assert result_copy is not result and result_copy.__dict__ is not result.__dict__
    assert sim_copy._event_factor is factor
    assert result_copy._weighter.event_factor is factor
    assert result_copy.events is result.events
    assert result_copy.weights is result.weights and not result_copy.weights.flags.writeable
    assert factor.calls == calls
    for obj in (sim_copy, result_copy):
        with pytest.raises(NotSerializableError, match='event_factor'):
            pickle.dumps(obj)


def test_shallow_copy_preserves_subclass_slots():
    class SlottedResults(siren.Results):
        __slots__ = ('marker',)

    result = SlottedResults([], [], [], None, None)
    result.marker = object()
    clone = copy.copy(result)
    assert type(clone) is SlottedResults and clone.marker is result.marker


def test_drift_exception_retains_bad_factor_diagnostics():
    factor = Factor()
    result = simulation(factor=factor, chain=False).run()
    factor.value = float('nan')
    with pytest.raises(ConfigurationError, match='event factor must be finite') as exc:
        result.explain(0)
    assert exc.value.breakdown.flags
    assert np.isnan(exc.value.breakdown.total)
    assert result.weights[0] > 0


def test_failed_tuning_restores_quota_and_keeps_failure_report():
    inj, _ = _native_sampler(0)
    mixture = next(iter(inj.primary_phase_spaces.values()))
    mixture.channels = list(mixture.channels) * 2
    mixture.weights = [.5, .5]
    inj.engine.ResetInjectedEvents(99)
    with pytest.raises(GenerationFailure, match='KinematicallyForbidden') as exc:
        tune(inj, lambda tree: 1, rounds=1, events=7, on_failure='raise')
    assert inj.engine.EventsToInject() == 99
    assert inj.engine.InjectionAttempts() == 0
    assert exc.value.report.attempts == 1
    assert exc.value.report.by_vertex[0].reason == injection.FailureReason.KinematicallyForbidden
    assert len(exc.value.report.last_failed_tree.tree) == 1


def test_tuning_restores_quota_when_weighting_raises():
    sim = _sim(_vertex(_IsoDecay(), .1 * channels.physical() + .9 * channels.toward('NuLight', _box())))
    inj = sim.injector
    inj.engine.ResetInjectedEvents(99)

    def fail(tree):
        raise RuntimeError('weighting failed')

    with pytest.raises(RuntimeError, match='weighting failed'):
        tune(inj, fail, rounds=1, events=7)
    assert inj.engine.EventsToInject() == 99
    assert inj.engine.InjectionAttempts() == 0
