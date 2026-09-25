"""Facade weights, failures, and exported snapshots agree with Layer 2."""
import pickle

import awkward as ak
import h5py
import numpy as np
import pytest

import siren
from siren import distributions as d, injection
from siren.errors import ConfigurationError, GenerationFailure, NotSerializableError, WeightCalculationError
from siren.tune import Plan
from simulation_fixtures import offline_detector
from test_injector_api import _distributions, _depth_ge_1
from test_generation_policy import _native_sampler

PT = siren.particles


class Factor:
    def __init__(self, value=3.0):
        self.value = value
        self.calls = 0

    def __call__(self, tree):
        self.calls += 1
        return self.value


def simulation(*, factor=None, chain=True, physical=True, seed=71, events=5):
    primary = siren.Vertex(
        PT.NuMu, siren.interactions.DummyCrossSection(),
        distributions=_distributions(25.0),
        physical=[d.NormalizationConstant(0.2)] if physical else [],
        physical_interactions=siren.interactions.DummyCrossSection())
    secondary = siren.Vertex(
        PT.NuMu, siren.interactions.DummyCrossSection(),
        position=d.SecondaryPhysicalVertexDistribution(),
        physical=[d.NormalizationConstant(0.7)] if physical else [],
        physical_interactions=siren.interactions.DummyCrossSection())
    return siren.Simulation(
        detector=offline_detector('CCM'), primary=primary,
        secondaries=[secondary] if chain else [],
        stopping_condition=_depth_ge_1 if chain else None,
        seed=seed, events=events, event_factor=factor)


def read_hdf5(path):
    with h5py.File(str(path) + '.hdf5') as f:
        group = f['Events']
        return ak.from_buffers(group.attrs['form'], int(group.attrs['length']),
                               {k: np.asarray(v) for k, v in group.items()}), dict(group.attrs)


def test_physical_declarations_and_event_factor_match_explicit_weighter():
    factor = Factor()
    sim = simulation(factor=factor)
    results = sim.run(on_failure='raise', on_shortfall='raise')
    assert factor.calls == len(results)
    assert all(len(tree.tree) == 2 for tree in results.events)
    direct = injection.Weighter(
        sim.injector,
        primary_physical=sim.physical_distributions,
        secondary_physical=sim._secondary_physical_distributions,
        overrides={'primary_interactions': sim._physical_interactions,
                   'secondary_interactions': sim._secondary_physical_interactions},
        event_factor=factor)
    np.testing.assert_array_equal(results.weights, direct.weight_all(results.events))
    assert sim.weighter.primary_interactions[0] is sim._physical_interactions[0]
    assert sim.weighter.secondary_interactions[PT.NuMu][0] is sim._secondary_physical_interactions[PT.NuMu][0]
    assert sim.weighter.primary_interactions[0] is not sim.injector.primary_interactions[0]
    conditional = injection.Weighter(sim.injector)
    np.testing.assert_allclose(results.weights / conditional.weight_all(results.events), 0.2 * 0.7 * 3)
    assert results.explain(0).event_factor == 3


def test_vertex_empty_physical_list_matches_direct_injector():
    sim = simulation(chain=False, physical=False)
    results = sim.run(on_failure='raise')
    assert sim.physical_distributions == []
    np.testing.assert_array_equal(results.weights, injection.Weighter(sim.injector).weight_all(results.events))


def test_vertex_physical_declarations_with_named_injection_slots():
    primary = siren.Vertex(PT.NuMu, siren.interactions.DummyCrossSection(),
                           physical=[d.NormalizationConstant(7)])
    sim = siren.Simulation(events=2, detector=offline_detector('CCM'), primary=primary,
                           injection_energy=d.Monoenergetic(2),
                           injection_direction=d.IsotropicDirection(),
                           position=d.PointSourcePositionDistribution(siren.math.Vector3D(), 25))
    assert sim.physical_distributions == primary.physical
    results = sim.run(on_failure='raise')
    direct = injection.Weighter(sim.injector, primary_physical=primary.physical)
    np.testing.assert_array_equal(results.weights, direct.weight_all(results.events))


def test_reweight_inherits_replaces_and_removes_factor_independently():
    factor = Factor()
    sim = simulation(factor=factor)
    results = sim.run()
    np.testing.assert_array_equal(sim.reweight().weights, results.weights)
    np.testing.assert_allclose(sim.reweight(event_factor=None).weights, results.weights / 3)
    np.testing.assert_allclose(sim.reweight(event_factor=Factor(6)).weights, results.weights * 2)
    np.testing.assert_allclose(sim.reweight(secondary_physical_distributions={}).weights, results.weights / 0.7)
    np.testing.assert_allclose(sim.reweight(physical_distributions=[]).weights, results.weights / 0.2)
    np.testing.assert_array_equal(sim.reweight().weights, results.weights)
    sim.injector.reset()
    with pytest.raises(ConfigurationError, match='changed since run'):
        sim.reweight()


@pytest.mark.parametrize('bad', [-1.0, float('nan'), float('inf')])
def test_bad_event_factor_raises_with_event_context(bad):
    sim = simulation(factor=Factor(bad), chain=False)
    with pytest.raises(WeightCalculationError, match='Event 0') as exc:
        sim.run()
    assert exc.value.event_index == 0
    assert exc.value.event is not None


@pytest.mark.parametrize('optimize', [False, Plan(rounds=1, events=2)])
def test_strict_policy_reaches_native_generation_and_tuning(optimize):
    inj, _ = _native_sampler(0)
    # Wrap the failing native channel in the public channel vocabulary.
    mixture = next(iter(inj.primary_phase_spaces.values()))
    channel = mixture.channels[0]
    primary = siren.Vertex(
        inj.primary_type, inj.primary_interactions,
        distributions=inj.primary_injection_distributions,
        kinematics=siren.channels.Channel(lambda signature, **kwargs: channel))
    sim = siren.Simulation(events=1, detector=inj.detector_model,
                           primary=primary, seed=331663)
    with pytest.raises(GenerationFailure, match='KinematicallyForbidden') as exc:
        sim.run(optimize=optimize, on_failure='raise')
    assert exc.value.report.attempts == 1


def test_tuning_uses_corrected_event_weights():
    from test_variance_report_integration import _build_results
    existing = _build_results(offline_detector('CCM'), n_events=5)
    factor = Factor(2)
    existing._weighter.event_factor = factor
    from siren.tune import tune
    report = tune(existing._injector, existing._weighter, rounds=1, events=8,
                  on_failure='raise')
    assert len(report.ess_trajectory) == 1
    assert factor.calls > 0


@pytest.mark.parametrize('selected', [False, True])
@pytest.mark.parametrize('native', [False, True])
def test_save_uses_weight_and_count_snapshots(tmp_path, selected, native):
    factor = Factor()
    sim = simulation(factor=factor)
    original = sim.run()
    results = original[1:4].where(lambda tree, weight: True) if selected else original
    weights = results.weights.copy()
    headers = [list(tree.header.weights) for tree in results.events]
    factor.value = 17
    sim.injector.reset(123)
    calls = factor.calls
    out = tmp_path / 'snapshot'
    results.save(out, save_parquet=True, save_siren_events=native, pot=9e20)
    assert factor.calls == calls
    table, attrs = read_hdf5(out)
    np.testing.assert_array_equal(ak.to_numpy(table.event_weight), weights)
    np.testing.assert_array_equal(ak.to_numpy(ak.from_parquet(str(out) + '.parquet').event_weight), weights)
    assert attrs['attempted_events'] == results.attempts
    assert attrs['accepted_events'] == results.injected
    assert attrs['events_to_inject'] == results.requested
    assert attrs['pot'] == 9e20
    if native:
        trees = siren.dataclasses.LoadInteractionTrees(str(out))
        np.testing.assert_array_equal([tree.header.weights[0] for tree in trees], weights)
    assert [list(tree.header.weights) for tree in results.events] == headers
    with pytest.raises(ConfigurationError, match='snapshot'):
        results.explain(0)


def test_merge_matches_pooled_weighter_and_saves_scaled_weights(tmp_path):
    factor = Factor()
    sim = simulation(factor=factor, chain=False)
    a = sim.run()
    b = sim.run()
    merged = siren.Results.merge([a, b])
    pooled = injection.Weighter(a._injector, b._injector,
                               primary_physical=sim.physical_distributions,
                               overrides={'primary_interactions': sim._physical_interactions},
                               event_factor=factor)
    np.testing.assert_allclose(merged.weights, pooled.weight_all(merged.events), rtol=1e-12)
    np.testing.assert_array_equal(merged.weights, np.concatenate([a.weights, b.weights]) / 2)
    out = tmp_path / 'pooled'
    merged[1:].save(out, save_parquet=False, save_siren_events=False)
    table, attrs = read_hdf5(out)
    np.testing.assert_array_equal(ak.to_numpy(table.event_weight), merged.weights[1:])
    assert attrs['accepted_events'] == a.injected + b.injected
    assert merged.explain(0).total == pytest.approx(merged.weights[0], rel=1e-12)
    assert merged.explain(0).event_factor == 3
    with pytest.raises(ConfigurationError, match='independent'):
        siren.Results.merge([a, a[:2]])


@pytest.mark.parametrize('change', ['physical', 'factor', 'energy', 'detector'])
def test_merge_rejects_different_physics_or_generation(change):
    sim = simulation(factor=Factor(), chain=False)
    a = sim.run()
    if change == 'physical':
        sim._physical_distributions = [d.NormalizationConstant(8)]
    elif change == 'factor':
        sim._event_factor = Factor(8)
    elif change == 'energy':
        sim._injection_distributions[1] = d.Monoenergetic(1)
    else:
        sim._detector_model = offline_detector()
    b = sim.run()
    with pytest.raises(ConfigurationError, match='differing config'):
        siren.Results.merge([a, b])


def test_result_array_does_not_alias_input_and_rejects_mismatched_lengths():
    weights = np.array([1.0, 2.0])
    result = siren.Results([1, 2], weights, [0, 0], None, None)
    weights[0] = 9
    assert result.weights[0] == 1
    with pytest.raises(ValueError, match='one weight'):
        siren.Results([1, 2], [1], [0, 0], None, None)
    with pytest.raises(ConfigurationError, match='unknown configuration'):
        siren.Results.merge([result])


def test_archiving_with_event_factor_is_explicitly_rejected(tmp_path):
    sim = simulation(factor=Factor(), chain=False)
    result = sim.run()
    for obj in (sim, result, sim.weighter):
        with pytest.raises(NotSerializableError, match='event_factor'):
            pickle.dumps(obj)
    with pytest.raises(NotSerializableError, match='event_factor'):
        sim.weighter.save(str(tmp_path / 'weighter'))


def test_fixed_weighting_uses_the_same_compatibility_check_as_layer2():
    from test_generation_policy import _Source, _Decay

    class PointSource(_Source):
        def DensityVariables(self):
            return []  # The prescribed point has no differential position measure.

    source = PointSource()
    primary = siren.Vertex(PT.NuMu, _Decay(PT.NuMu, PT.NuE, source),
                           distributions=[source], weighting=siren.Fixed())
    detector = siren.detector.DetectorModel()
    sim = siren.Simulation(events=2, detector=detector, primary=primary)
    assert not sim.weighter.engine.GetPrimaryPhysicalProcess().GetWeightingMode().compute_position_probability
    direct = injection.Weighter(sim.injector)
    assert not direct.engine.GetPrimaryPhysicalProcess().GetWeightingMode().compute_position_probability
    primary.weighting = siren.Propagated()
    with pytest.raises(ValueError, match='not covered by injection'):
        siren.Simulation(events=2, detector=detector, primary=primary).weighter.engine


def test_secondary_weighting_and_expansion_declarations_reach_engine():
    from test_generation_policy import _Source, _Decay
    from test_failure_report import _FixedSecondary
    from siren.expand import child, depth_below
    source = _Source()
    primary = siren.Vertex(PT.NuMu, _Decay(PT.NuMu, PT.NuE, source),
                           distributions=[source], expand=[child('NuE')])
    secondary = siren.Vertex(PT.NuE, _Decay(PT.NuE, PT.NuTau, source),
                             position=_FixedSecondary(), weighting=siren.Fixed(),
                             expand=[depth_below(0)])
    sim = siren.Simulation(detector=siren.detector.DetectorModel(), primary=primary,
                           secondaries=[secondary], events=2)
    engine = sim.injector.engine
    assert not engine.GetSecondaryProcessMap()[PT.NuE].GetWeightingMode().compute_interaction_probability
    trees = sim.injector.generate(on_failure='raise')
    assert all(len(tree.tree) == 2 for tree in trees)
    assert not sim.weighter.engine.GetSecondaryPhysicalProcesses()[0].GetWeightingMode().compute_interaction_probability


def test_zero_factor_retains_events_and_exposure():
    result = simulation(factor=Factor(0), chain=False).run(on_failure='raise')
    assert len(result) == result.injected == 5
    np.testing.assert_array_equal(result.weights, np.zeros(5))


def test_save_rejects_weight_override_before_writing(tmp_path):
    result = simulation(chain=False).run()
    with pytest.raises(ValueError, match='hepmc3_weights'):
        result.save(tmp_path / 'bad', hepmc3_weights='none')
    assert not list(tmp_path.iterdir())


def test_native_snapshot_reexport_trusts_corrected_headers(tmp_path):
    from siren._util import resolve_hepmc3_weight_policy
    result = simulation(chain=False, factor=Factor()).run()
    for tree in result.events:
        tree.header.provenance = {'siren.weights_state': 'unweighted'}
    result.save(tmp_path / 'trusted', save_hdf5=False, save_parquet=False)
    trees = siren.dataclasses.LoadInteractionTrees(str(tmp_path / 'trusted'))
    weights, state = resolve_hepmc3_weight_policy(trees, 'auto', None)
    assert state == 'header'
    np.testing.assert_array_equal(weights, result.weights)
    assert all(tree.header.provenance['siren.weights_state'] == 'unweighted' for tree in result.events)


def test_reweight_requires_a_run_even_after_inspecting_pipeline():
    sim = simulation(chain=False)
    assert sim.injector is not None and sim.weighter is not None
    with pytest.raises(RuntimeError, match='Must call run'):
        sim.reweight()
