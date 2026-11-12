"""Pilot isolation, support, exact weight explanations and frozen Results."""
import math
import pickle
from fractions import Fraction

import numpy as np
import pytest
import siren
from siren.tune import Plan, tune, _ess
from test_authoring_creation_survival import facade
from test_results_simulation import simulation
from test_weighter_inheritance import _injector


def pair(tmp_path, seed=31):
    gen = facade(tmp_path)
    gen.seed = seed
    gen.primary.kinematics = .5*siren.channels.physical() + .5*siren.channels.isotropic()
    return gen, siren.injection.Weighter(gen)


def test_pilot_isolation_and_production_snapshot(tmp_path):
    gen, weight = pair(tmp_path)
    gen.engine
    before = pickle.dumps(gen)
    report = tune(gen, weight, events=64, rounds=3, seed=91)
    assert pickle.dumps(gen) == before
    assert report.attempts == [64]*3
    assert report.successes == [64]*3
    assert all(v <= 64 for v in report.ess_trajectory)
    assert report.injector.seed == 31
    assert report.injector.injection_attempts == 0
    assert report.injector.number_of_events == 8
    assert all(w >= .001 for level in report.alpha_history[-1] for w in level)
    a = siren.generate(report.injector, report.weighter, events=8, on_failure='raise')
    expected = .2 * -math.expm1(-1.) * math.exp(-5.)
    assert sum(a.weights) == pytest.approx(expected, rel=3e-12, abs=0)
    assert a.attempts == a.requested == a.injected == 8
    old = a.weights.copy()
    for i,w in enumerate(old):
        bd = a.explain(i)
        assert bd.total == pytest.approx(w, rel=1e-12, abs=0)
        assert bd.vertices[0].survival == pytest.approx(math.exp(-5.),rel=1e-12)
    # Subsequent tuning uses another private pair and cannot affect this snapshot.
    tune(gen, weight, events=17, rounds=2, seed=92)
    np.testing.assert_array_equal(a.reweight().weights, old)
    b = pickle.loads(pickle.dumps(a))
    np.testing.assert_array_equal(b.weights, old)
    assert b.explain(0).total == pytest.approx(old[0], rel=1e-12, abs=0)


def test_simulation_runs_tuned_pair_and_keeps_pilot_counts_separate(tmp_path):
    gen, _ = pair(tmp_path)
    sim = siren.Simulation(events=8, detector=gen.detector_model, primary=gen.primary, seed=31)
    result = sim.run(optimize=Plan(events=32, rounds=2, seed=93),on_failure='raise')
    assert result.attempts == 8
    assert sim.optimization_report.attempts == [32,32]
    first = result.weights.copy()
    sim.run(optimize=Plan(events=40, rounds=1, seed=94),on_failure='raise')
    np.testing.assert_array_equal(result.reweight().weights, first)


@pytest.mark.parametrize('kwargs',[{'events':0},{'events':2.5},{'rounds':False},
    {'seed':-1},{'rule':'bad'},{'damping':float('nan')},{'min_weight':0},
    {'min_weight':float('inf')},{'metric':3},{'recurse':1}])
def test_plan_rejects_invalid_configuration(kwargs):
    with pytest.raises((ValueError,TypeError)):
        Plan(**kwargs)


def test_pilot_errors_leave_original_state_unchanged(tmp_path):
    gen, weight = pair(tmp_path)
    before=pickle.dumps(gen)
    with pytest.raises(siren.errors.WeightCalculationError, match='Metric'):
        tune(gen, weight, events=4, rounds=1, metric=lambda ev,w: math.nan)
    assert pickle.dumps(gen)==before
    with pytest.raises(siren.errors.ConfigurationError, match='distinct'):
        tune(gen,weight,seed=31)
    with pytest.raises(ValueError, match='simplex'):
        tune(gen,weight,min_weight=.9)
    gen.generate(1)
    with pytest.raises(siren.errors.ConfigurationError, match='fresh'):
        tune(gen,weight)


def test_explanations_event_factors_and_pooled_injectors():
    a,b=simulation(seed=71).run(),simulation(seed=72).run()
    pooled=siren.Results.merge([a,b])
    bd=pooled.explain(0)
    assert {v.injector_index for v in bd.vertices}=={0,1}
    assert bd.total==pytest.approx(pooled.weights[0],rel=1e-12,abs=0)
    gen=_injector()
    tree=gen.generate(1)[0]
    for factor in [0.,1.,3.,-1.,math.inf,math.nan,Fraction(1,10**400)]:
        w=siren.injection.Weighter(gen,event_factor=lambda ev:factor)
        bd=w.explain(tree)
        if factor in (0.,1.,3.):
            assert bd.total==w(tree)
        else:
            assert math.isnan(bd.total)
            assert bd.flags
            with pytest.raises(siren.errors.WeightCalculationError):
                w(tree)


def test_positive_physical_products_do_not_become_zero_support():
    gen=_injector()
    tree=gen.generate(1)[0]
    w=siren.injection.Weighter(gen,
        primary_physical=[siren.dist.NormalizationConstant(1e-200)],
        secondary_physical={gen.primary_type:[siren.dist.NormalizationConstant(1e-200)]})
    with pytest.raises(siren.errors.WeightCalculationError, match='underflow'):
        w(tree)
    bd=w.explain(tree)
    assert math.isnan(bd.total)
    assert any('underflow' in f for v in bd.vertices for f in v.flags)


def test_scaled_ess_and_signed_metrics():
    assert _ess([1e300,1e300,0])==2
    assert _ess([1e-300,1e-300])==2
    assert _ess([-1e300,1e300])==0
    assert _ess([])==0


def iso_mixture(weights=(.5,.5)):
    return siren.injection.MultiChannelPhaseSpace(
        [siren.injection.Isotropic2BodyChannel(0) for _ in weights], list(weights))


def record():
    r=siren.dataclasses.InteractionRecord()
    r.primary_mass=1
    r.primary_momentum=[1,0,0,0]
    r.secondary_masses=[0,0]
    r.secondary_momenta=[[.5,0,0,.5],[.5,0,0,-.5]]
    r.signature.primary_type=siren.particles.N4
    r.signature.target_type=siren.particles.Decay
    r.signature.secondary_types=[siren.particles.NuMu,siren.particles.Gamma]
    return r


@pytest.mark.parametrize('w',[1.,1e300,1e-300,-1e300])
def test_native_tuning_extreme_weights_and_pending_archive(w):
    mc=iso_mixture((.2,.8))
    mc.Accumulate(None,record(),w)
    assert mc.tuning_count==1
    with pytest.raises(RuntimeError,match='tuning'):
        pickle.dumps(mc)
    mc.UpdateWeights(min_weight=.01)
    assert mc.tuning_count==0
    assert all(a>=.01 for a in mc.weights)
    restored=pickle.loads(pickle.dumps(mc))
    assert restored.weights==mc.weights


def test_native_nested_tuning_and_transactional_validation():
    inner=iso_mixture((.2,.8))
    outer=siren.injection.MultiChannelPhaseSpace(
        [siren.injection.Isotropic2BodyChannel(0),siren.injection.NestedMixtureChannel(inner)],[.5,.5])
    outer.Accumulate(None,record(),2.)
    assert inner.tuning_count==outer.tuning_count==1
    before=outer.weights
    with pytest.raises(ValueError):
        outer.UpdateWeights(damping=math.nan)
    assert outer.weights==before and inner.tuning_count==1
    outer.UpdateWeights(min_weight=.02)
    assert inner.tuning_count==outer.tuning_count==0
    assert min(inner.weights)>=.02
    outer.Accumulate(None,record(),2.)
    inner.weights=[.3,.7]
    with pytest.raises(ValueError,match='changed'):
        outer.UpdateWeights()
    assert outer.weights==before
    outer.ResetAccumulators()


def test_terminal_counts_remain_separate_from_failures():
    sim = simulation(events=2)
    results = sim.run(on_failure='raise', on_shortfall='raise')
    report = sim.injector.report()
    assert report.unregistered_secondaries > 0
    assert report.successes == 2
    assert 'terminal daughters' in str(report)


def test_common_enum_values_remain_compatible():
    assert int(siren.injection.FailureReason.UnregisteredSecondaryType) == 5


@pytest.mark.parametrize('value',[Fraction(1,10**400),Fraction(-1,10**400),[1.],complex(1,1)])
def test_invalid_metric_values_are_visible(tmp_path,value):
    gen,weight=pair(tmp_path)
    with pytest.raises(siren.errors.WeightCalculationError,match='Metric'):
        tune(gen,weight,events=2,rounds=1,metric=lambda ev,w:value)
    assert gen.injection_attempts==0


def test_signed_metric_has_same_squared_statistic(tmp_path):
    gen,weight=pair(tmp_path)
    a=tune(gen,weight,events=16,rounds=2,seed=201)
    b=tune(gen,weight,events=16,rounds=2,seed=201,metric=lambda ev,w:-w)
    assert a.alpha_history==b.alpha_history
    assert a.ess_trajectory==b.ess_trajectory


def test_all_miss_pilot_retains_owned_failure_reports(tmp_path):
    gen,weight=pair(tmp_path)
    far=siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(1e6,0,11)),1.,0.)
    gen.primary.distributions[-1]=siren.distributions.PrimaryBoundedVertexDistribution(far,100.)
    report=tune(gen,weight,events=9,rounds=2,seed=810,on_failure='raise')
    assert report.attempts==[9,9]
    assert report.successes==[0,0]
    assert report.ess_trajectory==[0.,0.]
    assert all(r.failures==9 for r in report.failures)
    assert all(len(r.last_failed_tree.tree)==1 for r in report.failures)
    assert report.injector.injection_attempts==0
    report.injector.reset()
    assert all(r.attempts==9 for r in report.failures)
    assert all(len(r.last_failed_tree.tree)==1 for r in report.failures)
