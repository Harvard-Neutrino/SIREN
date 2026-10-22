"""Attempt-normalized snapshots, public compiler parity, persistence and pooling."""
import pickle

import awkward as ak
import h5py
import numpy as np
import pytest
import siren

from simulation_fixtures import offline_detector
from test_injector_api import _distributions, _working_injector, _forced_failure_injector


def simulation(events=5, seed=71, factor=None):
    primary = siren.Vertex('NuMu', siren.interactions.DummyCrossSection(),
                           distributions=_distributions(25),
                           physical=[siren.dist.NormalizationConstant(0.2)])
    return siren.Simulation(events=events, seed=seed, detector=offline_detector('CCM'),
                            primary=primary, event_factor=factor)


def test_facade_equals_direct_compiler_and_seeded_stream():
    sim = simulation()
    result = sim.run(on_failure='raise', on_shortfall='raise')
    direct = siren.injection.Injector(events=5, seed=71, detector=sim.detector_model,
                                      primary=sim._primary)
    trees = direct.generate(on_failure='raise', on_shortfall='raise')
    weights = siren.injection.Weighter(direct).weight_all(trees)
    assert result.attempts == direct.injection_attempts
    assert result.injected == len(trees)
    np.testing.assert_array_equal(result.weights, weights)
    for left, right in zip(result.events, trees):
        for a, b in zip(left.tree, right.tree):
            np.testing.assert_array_equal(a.record.primary_momentum, b.record.primary_momentum)
            np.testing.assert_array_equal(a.record.secondary_momenta, b.record.secondary_momenta)
            np.testing.assert_array_equal(a.record.interaction_vertex, b.record.interaction_vertex)


def test_weights_and_counts_survive_later_generation_and_reset():
    sim = simulation()
    first = sim.run()
    stored = first.weights.copy()
    counts = (first.attempts, first.injected, first.requested)
    sim.injector.generate(3)
    sim.injector.reset(91)
    sim.injector.generate(5)
    np.testing.assert_array_equal(first.weights, stored)
    np.testing.assert_array_equal(first.reweight().weights, stored)
    assert (first.attempts, first.injected, first.requested) == counts
    with pytest.raises(ValueError):
        first.weights.setflags(write=True)


def test_physical_source_mutation_cannot_change_snapshot():
    sim = simulation()
    result = sim.run()
    sim._primary.physical[0].normalization = 90
    np.testing.assert_array_equal(result.reweight().weights, result.weights)
    changed = result.reweight(primary_physical=[siren.dist.NormalizationConstant(0.4)])
    np.testing.assert_allclose(changed.weights, 2 * result.weights, rtol=1e-14)


def test_pooling_uses_attempts_and_matches_native_pooled_weighter():
    a, b = simulation(3, 71).run(), simulation(8, 72).run()
    assert a.attempts != b.attempts
    merged = siren.Results.merge([a, b])
    total = a.attempts + b.attempts
    np.testing.assert_allclose(merged.weights,
        np.concatenate([a.weights * a.attempts / total, b.weights * b.attempts / total]), rtol=1e-14)
    pooled = siren.injection.Weighter(*a._weighter.injectors, *b._weighter.injectors,
        primary_physical=a._weighter.primary_physical_distributions)
    np.testing.assert_allclose(merged.weights, pooled.weight_all(merged.events), rtol=1e-13)
    np.testing.assert_allclose(merged.reweight().weights, merged.weights, rtol=1e-13)
    assert merged.attempts == total


def test_merge_rejects_reset_epochs_from_same_source():
    gen = _working_injector(events=3)
    a = siren.generate(gen, siren.injection.Weighter(gen), events=3)
    gen.reset()
    b = siren.generate(gen, siren.injection.Weighter(gen), events=3)
    with pytest.raises(siren.errors.ConfigurationError, match='independent'):
        siren.Results.merge([a, b])
    with pytest.raises(siren.errors.ConfigurationError, match='independent'):
        siren.Results.merge([a, pickle.loads(pickle.dumps(a))])


def test_single_call_helper_refuses_unrepresented_prior_batches():
    gen = _working_injector(events=3)
    gen.generate(2)
    with pytest.raises(siren.errors.ConfigurationError, match='fresh counters'):
        siren.generate(gen, siren.injection.Weighter(gen), events=3)


@pytest.mark.parametrize('save_native', [False, True])
def test_export_and_snapshot_round_trip_preserve_weights_and_counts(tmp_path, save_native):
    sim = simulation()
    result = sim.run()[1:4].where(lambda tree, weight: weight > 0)
    sim.injector.reset()
    path = tmp_path / 'snapshot'
    result.save(path, snapshot=True, save_siren_events=save_native, pot=1e20)
    restored = siren.Results.load(path)
    np.testing.assert_array_equal(restored.weights, result.weights)
    np.testing.assert_array_equal(restored.reweight().weights, result.weights)
    assert restored.attempts == result.attempts
    with h5py.File(str(path)+'.hdf5') as stream:
        group = stream['Events']
        values = ak.from_buffers(group.attrs['form'], int(group.attrs['length']),
                                 {k:np.asarray(v) for k,v in group.items()})
        np.testing.assert_array_equal(ak.to_numpy(values.event_weight), result.weights)
        assert group.attrs['attempted_events'] == result.attempts
        assert group.attrs['accepted_events'] == result.injected
        assert group.attrs['events_to_inject'] == result.requested
        assert group.attrs['pot'] == 1e20
    np.testing.assert_array_equal(ak.to_numpy(ak.from_parquet(str(path)+'.parquet').event_weight), result.weights)
    if save_native:
        trees = siren.dataclasses.LoadInteractionTrees(str(path))
        np.testing.assert_array_equal([tree.header.weights[0] for tree in trees], result.weights)


def test_zero_success_run_keeps_attempts():
    gen = _forced_failure_injector(events=2, max_attempts=5)
    result = siren.generate(gen, siren.injection.Weighter(gen), events=2, on_shortfall='ignore')
    assert len(result) == 0
    assert result.attempts == 5 and result.injected == 0 and result.requested == 2


def test_unknown_keywords_fail_before_loading_detector():
    with pytest.raises(TypeError, match='event_factor'):
        siren.Simulation(events=2, detector='no-such-detector', primary='NuMu', event_fator=lambda x:1)


def test_legacy_slots_use_same_compiler():
    dists = _distributions(25)
    with pytest.warns(DeprecationWarning):
        legacy = siren.Simulation(events=3, seed=4, detector=offline_detector('CCM'),
            primary='NuMu', interactions=[siren.interactions.DummyCrossSection()],
            injection_distributions=dists)
    direct = siren.Simulation(events=3, seed=4, detector=legacy.detector_model,
        primary=siren.Vertex('NuMu', legacy._primary.interactions, distributions=dists))
    np.testing.assert_array_equal(legacy.run().weights, direct.run().weights)


def test_event_factor_snapshot_export_does_not_reevaluate(tmp_path):
    calls=[]
    sim=simulation(factor=lambda tree: calls.append(tree) or 2.0)
    result=sim.run()
    before=len(calls)
    result.save(tmp_path/'factor', save_parquet=False)
    assert len(calls)==before
    with pytest.raises(siren.errors.NotSerializableError):
        result.save(tmp_path/'factor-config', snapshot=True)
    assert not (tmp_path/'factor-config.siren_events').exists()
    np.testing.assert_allclose(result.reweight(event_factor=None).weights,result.weights/2)


def test_reweighted_target_is_frozen_and_can_be_pooled():
    norm = siren.dist.NormalizationConstant(0.4)
    a, b = simulation(3, 13).run(), simulation(3, 14).run()
    ra, rb = a.reweight(primary_physical=[norm]), b.reweight(primary_physical=[norm])
    norm.normalization = 12
    np.testing.assert_allclose(ra.reweight().weights, ra.weights, rtol=1e-14)
    merged = siren.Results.merge([ra, rb])
    np.testing.assert_allclose(merged.reweight().weights, merged.weights, rtol=1e-14)
    with pytest.raises(siren.errors.ConfigurationError, match='differing config'):
        siren.Results.merge([a, rb])


def test_generate_rejects_unrelated_weighter_before_attempting():
    a, b = _working_injector(events=2), _working_injector(events=2)
    with pytest.raises(siren.errors.ConfigurationError, match='same injector'):
        siren.generate(a, siren.injection.Weighter(b), events=2)
    assert a.injection_attempts == 0


def test_keyword_box_geometry_and_copy():
    box = siren.geometry.Box(widths=[2, 4, 6], center=[1, 2, 3])
    assert box.X == 2 and box.Y == 4 and box.Z == 6
    assert siren.geometry.Box(box).X == 2
    with pytest.raises(TypeError):
        siren.geometry.Box([2, 4, 6], [1, 2, 3])


def missed_ray_simulation(seed, attempts=30):
    detector = siren.detector.DetectorModel()
    materials = detector.Materials
    materials.AddMaterial('ARGON', {1000180400: 1.0})
    detector.Materials = materials
    sector = siren.detector.DetectorSector()
    sector.name, sector.level = 'argon', 0
    sector.material_id = materials.GetMaterialId('ARGON')
    sector.geo = siren.geometry.Sphere(siren.geometry.Placement(), 1.0, 0.0)
    sector.density = siren.detector.ConstantDensityDistribution(1.4)
    detector.AddSector(sector)
    dists = _distributions(10)
    dists[-1] = siren.distributions.PointSourcePositionDistribution(siren.math.Vector3D(2,0,0), 10)
    return siren.Simulation(events=1, seed=seed, max_attempts=attempts, detector=detector,
        primary=siren.Vertex('NuMu', siren.interactions.DummyCrossSection(), distributions=dists))


def test_pooling_includes_failed_and_zero_success_attempts(tmp_path):
    path = tmp_path/'rays.csv'
    path.write_text('E,px,py,pz,x0,y0,z0,m\n1,0,0,1,0,0,-2,0\n1,0,0,1,10,0,-2,0\n')
    runs=[]
    for seed in (1,2,3):
        sim = missed_ray_simulation(seed,2)
        sim._primary.distributions = [siren.dist.External(str(path)),
            siren.dist.BoundedPrimaryVertex(siren.geometry.Sphere(siren.geometry.Placement(),1.,0.),10)]
        runs.append(sim.run(on_shortfall='ignore'))
    assert [(r.attempts,r.injected) for r in runs] == [(1,1),(2,0),(2,1)]
    merged=siren.Results.merge(runs)
    assert (merged.attempts,merged.injected,merged.requested) == (5,2,3)
    np.testing.assert_allclose(merged.weights, [runs[0].weights[0]/5, runs[2].weights[0]*2/5],rtol=1e-14)
    np.testing.assert_allclose(merged.reweight().weights,merged.weights,rtol=1e-14)
    np.testing.assert_allclose(pickle.loads(pickle.dumps(merged)).reweight().weights,merged.weights,rtol=1e-14)


def test_legacy_loader_options_are_explicit_and_metadata_survives(monkeypatch):
    from siren import _simulation_legacy as legacy
    calls=[]
    def loader(name, **options):
        calls.append((name,options))
        return siren.ProcessBundle({siren.particles.NuMu:[siren.interactions.DummyCrossSection()]},{},'metadata')
    monkeypatch.setattr(legacy,'load_processes',loader)
    with pytest.warns(DeprecationWarning):
        sim=siren.Simulation(events=1, detector=offline_detector(), primary='NuMu',
            interactions='mock', targets='Nucleon', process='CC',
            process_options={'coupling':0.2}, injection_distributions=_distributions(25))
    assert calls == [('mock',{'primary_types':[siren.particles.NuMu],
        'target_types':[siren.particles.Nucleon], 'process_types':['CC'],'coupling':0.2})]
    assert sim.process_metadata == ('metadata',)


def test_merge_rejects_identical_seeds_from_fresh_injectors():
    sim = simulation(4, 71)
    a, b = sim.run(), sim.run()
    np.testing.assert_array_equal(a.weights, b.weights)
    with pytest.raises(siren.errors.ConfigurationError, match='independent seeds'):
        siren.Results.merge([a, b])


@pytest.mark.parametrize('change', ['length', 'fiducial'])
def test_merge_rejects_different_fixed_target_bounds(change):
    cylinder=siren.geometry.Cylinder(siren.geometry.Placement(),0.05,0.,0.1)
    fid=siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(0,0,5)),1.,0.)
    positions=[siren.distributions.FixedTargetPositionDistribution(cylinder,fid,25.)]
    positions.append(siren.distributions.FixedTargetPositionDistribution(cylinder,fid,5.)
                     if change=='length' else
                     siren.distributions.FixedTargetPositionDistribution(cylinder,25.))
    results=[]
    for seed,pos in zip((81,82),positions):
        sim=simulation(3,seed)
        sim._primary.distributions[-1]=pos
        results.append(sim.run())
    with pytest.raises(siren.errors.ConfigurationError,match='differing config'):
        siren.Results.merge(results)


@pytest.mark.parametrize('option', ['targets', 'process', 'isoscalar'])
@pytest.mark.parametrize('use_vertex', [False, True])
def test_legacy_loader_selectors_are_not_ignored(option,use_vertex):
    original=simulation(2)
    kwargs={option: {'targets':'Proton','process':'NC','isoscalar':True}[option]}
    if use_vertex:
        kwargs['primary']=original._primary
    else:
        kwargs.update(primary='NuMu',interactions=original._primary.interactions,
                      injection_distributions=original._primary.distributions)
    with pytest.warns(DeprecationWarning), pytest.raises(TypeError,match='named interaction loader'):
        siren.Simulation(events=2,detector=original.detector_model,**kwargs)


@pytest.mark.parametrize('pooled', [False, True])
def test_results_constructor_checks_weighter_source(pooled):
    a,b=simulation(3,11),simulation(3,12)
    ra,rb=a.run(),b.run()
    wrong=siren.injection.Weighter(a.injector,b.injector) if pooled else b.weighter
    with pytest.raises(siren.errors.ConfigurationError,match='same single injector'):
        siren.Results(ra.events,wrong.weight_all(ra.events),ra.gen_times,wrong,a.injector)


@pytest.mark.parametrize('zero_success', [False, True])
def test_empty_results_export_preserves_run_accounting(tmp_path,zero_success):
    if zero_success:
        gen=_forced_failure_injector(events=2,max_attempts=5)
        result=siren.generate(gen,siren.injection.Weighter(gen),events=2,on_shortfall='ignore')
    else:
        result=simulation(3).run().where(lambda event,weight:False)
    path=tmp_path/'empty'
    result.save(path,snapshot=True)
    restored=siren.Results.load(path)
    assert len(restored)==0 and restored.attempts==result.attempts
    assert restored.injected==result.injected
    assert len(restored.reweight())==0
    with h5py.File(str(path)+'.hdf5') as stream:
        assert stream['Events'].attrs['length']==0
        assert stream['Events'].attrs['attempted_events']==result.attempts
        assert stream['Events'].attrs['accepted_events']==result.injected


@pytest.mark.parametrize('scale',[1.,1e-280,1e280])
def test_effective_sample_size_is_scale_invariant(scale):
    result=siren.Results([None,None],[scale,scale],[0.,0.],None,None)
    assert result.ess()==2.


def _constant_factor(tree):
    return 2.


def test_pickleable_event_factor_still_obeys_archive_refusal(tmp_path):
    result=simulation(factor=_constant_factor).run()
    with pytest.raises(siren.errors.NotSerializableError):
        result.save(tmp_path/'factor',snapshot=True)
    assert not list(tmp_path.iterdir())


def test_process_bundle_unpacking_exposes_metadata_explicitly(monkeypatch):
    from siren import _util
    monkeypatch.setattr(_util,'load_resource',lambda *args,**kwargs:({'primary':[]},{},'ups','dec'))
    bundle=siren.load_processes('mock')
    primary,secondary=bundle
    assert primary=={'primary':[]} and secondary=={}
    assert bundle.metadata==('ups','dec')
    with pytest.raises(ValueError):
        primary,secondary,ups,dec=bundle
