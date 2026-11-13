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
    inner.SetWeights([.3,.7])
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
    gen=facade(tmp_path)
    gen.primary.kinematics = .5*siren.channels.physical() + .5*siren.channels.isotropic()
    far=siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(1e6,0,11)),1.,0.)
    gen.primary.distributions[-1]=siren.distributions.PrimaryBoundedVertexDistribution(far,100.)
    weight=siren.injection.Weighter(gen)
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


def _constant_factor(event):
    return 2.0


@pytest.mark.parametrize('chain',[False,True])
def test_tuning_retains_declared_physical_models_and_stopping_condition(chain):
    gen=_injector(chain=chain)
    weight=siren.injection.Weighter(gen)
    report=tune(gen,weight,events=4,rounds=1,seed=103)
    trees=report.injector.generate(5,on_failure='raise',on_shortfall='raise')
    assert all(len(t.tree)==(2 if chain else 1) for t in trees)
    expected=weight.weight_all(trees)*gen.number_of_events/report.injector.injection_attempts
    np.testing.assert_allclose(report.weighter.weight_all(trees),expected,rtol=2e-12,atol=0)
    assert report.injector is report.weighter.injectors[0]
    assert gen.injection_attempts==0


def anisotropic_pair(n=4096,seed=601,selected=False,condition=None):
    from test_native_phase_space import fixed_decay_setup,P
    model=siren.interactions.HNLDipoleDecay(1.,[1e-6]*3,
        siren.interactions.HNLDipoleDecay.ChiralNature.Dirac)
    raw,_=fixed_decay_setup(model,P.N4,1.,n)
    process=raw.GetPrimaryProcess()
    model=process.interactions.GetDecays()[0]
    choices=model.GetPossibleSignaturesFromParent(siren.particles.N4)
    vertex=siren.Vertex('N4',model,distributions=list(process.distributions),
        physical=list(process.distributions),weighting=siren.Fixed(),
        decay_channels=[choices[0]] if selected else None,
        kinematics=.3*siren.channels.physical()+.7*siren.channels.isotropic())
    secondary=[] if condition is None else [siren.Vertex(p,siren.interactions.DummyCrossSection(),
        position=siren.dist.SecondaryPhysicalVertexDistribution()) for p in ('NuE','NuMu','NuTau')]
    gen=siren.injection.Injector(events=n,seed=seed,detector=raw.GetDetectorModel(),
        primary=vertex,secondaries=secondary,stopping_condition=condition)
    return gen,siren.injection.Weighter(gen)


@pytest.mark.parametrize('selected',[False,True])
def test_changed_proposal_survives_production_and_native_archives(tmp_path,selected):
    gen,weight=anisotropic_pair(selected=selected)
    report=tune(gen,weight,events=1024,rounds=6,seed=104,damping=1.)
    active=[a for a in report.alpha_history[-1] if abs(a[0]-.3)>.05]
    assert active
    assert all(a[0]>.3 for a in active)
    assert all(list(m.weights)==[.3,.7] for m in gen.engine.GetPrimaryProcess().GetPhaseSpaceMap().values())
    result=siren.generate(report.injector,report.weighter,events=4096,on_failure='raise',on_shortfall='raise')
    expected=1/3 if selected else 1.
    se=np.std(result.weights*4096,ddof=1)/64
    assert abs(sum(result.weights)-expected)<6*se
    assert result.explain(1).total==pytest.approx(result.weights[1],rel=2e-12)
    np.testing.assert_allclose(result.reweight().weights,result.weights,rtol=2e-12,atol=0)
    restored=pickle.loads(pickle.dumps(result))
    np.testing.assert_allclose(restored.reweight().weights,result.weights,rtol=2e-12,atol=0)
    name=str(tmp_path/'tuned')
    with pytest.raises(siren.errors.NotSerializableError,match='tuning weighter'):
        report.injector.save(name)
    report.weighter.save(name)
    loaded_weight=siren.injection.Weighter();loaded_weight.load(name)
    np.testing.assert_allclose(loaded_weight.weight_all(result.events),result.weights,rtol=2e-12,atol=0)
    copied=pickle.loads(pickle.dumps(report.weighter))
    assert [list(m.weights) for m in copied.injectors[0].GetPrimaryProcess().GetPhaseSpaceMap().values()]==report.alpha_history[-1]


def test_event_factor_is_cloned_or_refused_explicitly(tmp_path):
    gen,weight=pair(tmp_path)
    weight.event_factor=_constant_factor
    report=tune(gen,weight,events=16,rounds=2,seed=107)
    event=report.injector.generate(1,on_failure='raise')[0]
    plain=siren.injection.Weighter(report.injector)
    assert report.weighter(event)==pytest.approx(2*plain(event))
    weight.event_factor=lambda event:1
    with pytest.raises(siren.errors.NotSerializableError,match='event_factor'):
        tune(gen,weight,events=4,rounds=1,seed=108)


def test_skipped_positive_vertex_is_not_a_culprit():
    from siren.report import VertexWeightLine,WeightBreakdown
    skipped=VertexWeightLine(0,14,1.,math.nan,{},[],['physical factors not evaluated (no usable proposal)'])
    broken=VertexWeightLine(1,14,1.,-1.,{},[],['physical density negative'])
    assert skipped.is_ok
    assert WeightBreakdown(math.nan,[skipped,broken]).culprit() is broken


def test_internal_positive_decay_product_does_not_become_physical_zero():
    from test_decay_channels import WidthOverride,channel,setup,P
    class RareDensity(WidthOverride):
        def FinalStateProbability(self,record):return 1e-200
    model=channel([P.NuMu,P.Gamma],1e-16,1e-16)
    gen,_,_=setup([model]);tree=gen.GenerateEvent()
    rare=RareDensity(model,partial=1e-216,owned=1e-216,parent=1e-16)
    physical=siren.injection.PhysicalProcess(P.N4,siren.interactions.InteractionCollection(P.N4,[rare]))
    physical.weighting_mode=siren.Fixed()
    physical.distributions=gen.GetPrimaryProcess().distributions
    weight=siren.injection._Weighter([gen],gen.GetDetectorModel(),physical)
    with pytest.raises(siren.errors.WeightCalculationError,match='underflow'):
        weight.EventWeight(tree)
    bd=weight.EventWeightWithBreakdown(tree)
    assert math.isnan(bd.total)
    assert any('underflow' in f for v in bd.vertices for f in v.flags)


class _AbortLaterPilot:
    def __init__(self):
        self.calls=0
    def __call__(self,tree,datum,index):
        self.calls+=1
        if self.calls>8:
            raise siren.errors.GenerationFailure('deliberate second-round pilot failure')
        return True


def test_later_pilot_failure_does_not_mutate_input_callback_or_proposal():
    condition=_AbortLaterPilot()
    gen,weight=anisotropic_pair(n=8,condition=condition)
    before=[list(m.weights) for m in gen.engine.GetPrimaryProcess().GetPhaseSpaceMap().values()]
    with pytest.raises(siren.errors.GenerationFailure):
        tune(gen,weight,events=8,rounds=2,seed=111)
    assert condition.calls==0
    assert gen.injection_attempts==0
    assert before==[list(m.weights) for m in gen.engine.GetPrimaryProcess().GetPhaseSpaceMap().values()]


def test_production_callback_starts_from_pristine_state():
    condition=_AbortLaterPilot()
    gen,weight=anisotropic_pair(n=2,condition=condition)
    report=tune(gen,weight,events=4,rounds=2,seed=112)
    assert condition.calls==0
    assert report.injector.stopping_condition.calls==0
    report.injector.generate(2,on_failure='raise',on_shortfall='raise')
    assert report.injector.stopping_condition.calls==2


def test_tuned_physical_declarations_cannot_be_lost_through_injector(tmp_path):
    gen=_injector(chain=False)
    report=tune(gen,siren.injection.Weighter(gen),events=4,rounds=1,seed=113)
    with pytest.raises(siren.errors.NotSerializableError,match='tuning weighter'):
        pickle.dumps(report.injector)
    with pytest.raises(siren.errors.NotSerializableError,match='tuning weighter'):
        report.injector.save(str(tmp_path/'bare'))
    with pytest.raises(siren.errors.ConfigurationError,match='tuning weighter'):
        siren.injection.Weighter(report.injector)
    event=report.injector.generate(1,on_failure='raise')[0]
    name=str(tmp_path/'physical')
    report.weighter.save(name)
    restored=siren.injection.Weighter();restored.load(name)
    assert restored(event)==pytest.approx(report.weighter(event),rel=2e-12)


def test_tuned_injector_without_physical_declarations_still_archives(tmp_path):
    gen,weight=pair(tmp_path)
    report=tune(gen,weight,events=8,rounds=1,seed=114)
    name=str(tmp_path/'sampling')
    event=report.injector.generate(1,on_failure='raise')[0]
    report.injector.save(name)
    loaded=siren.injection.Injector.load(name)
    assert siren.injection.Weighter(loaded)(event)==pytest.approx(report.weighter(event),rel=2e-12)


@pytest.mark.parametrize('density',[-1.,math.nan,math.inf])
def test_invalid_native_generation_retains_diagnostic_culprit(tmp_path,density):
    from test_python_distributions import InitialTimeShift
    class BadDensity(InitialTimeShift):
        def GenerationProbability(self,detector,interactions,record):return density
    gen=facade(tmp_path)
    gen.primary.distributions.append(BadDensity())
    weight=siren.injection.Weighter(gen)
    event=gen.generate(1,on_failure='raise')[0]
    with pytest.raises(siren.errors.WeightCalculationError):weight(event)
    breakdown=weight.explain(event)
    assert math.isnan(breakdown.total)
    assert breakdown.culprit() is not None
    assert not breakdown.culprit().is_ok


def test_legacy_unbuilt_input_is_refused_without_compiling_it(tmp_path):
    from test_creation_survival import setup,P
    raw,_=setup(tmp_path,5.,True)
    gen=facade(tmp_path)
    weight=siren.injection.Weighter(injectors=[gen],detector_model=raw.GetDetectorModel(),
        primary_type=P.N4,primary_interactions=list(raw.GetPrimaryProcess().interactions.GetDecays()),
        primary_physical_distributions=[])
    with pytest.raises(siren.errors.ConfigurationError,match='Build the injector'):
        tune(gen,weight,events=4,rounds=1,seed=115)
    gen.number_of_events=9


def test_tuned_python_stopping_callback_has_typed_archive_refusal(tmp_path):
    gen=_injector(chain=True)
    report=tune(gen,siren.injection.Weighter(gen),events=4,rounds=1,seed=116)
    with pytest.raises(siren.errors.NotSerializableError,match='stopping condition'):
        report.weighter.save(str(tmp_path/'callback'))
