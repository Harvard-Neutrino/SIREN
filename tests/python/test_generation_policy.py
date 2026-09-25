"""Strict generation rejects numerical failures at the failing attempt."""

import math

import pytest

import siren
from siren import dataclasses as dc, distributions as d, injection, interactions
from siren.errors import GenerationFailure, InjectionShortfall, WeightCalculationError
from test_failure_report import _FixedPrimary, _FixedSecondary, _ChainDecay
from test_injector_api import _forced_failure_injector


class _Source(_FixedPrimary):
    def RequiredVariables(self):
        return set()

    def SetVariables(self):
        dv = d.DistributionVariable
        return {dv.PrimaryEnergy, dv.PrimaryDirection, dv.InitialPosition,
                dv.InteractionVertex}


class _Decay(_ChainDecay):
    def DifferentialDecayWidth(self, record):
        return 1.0

    def FinalStateProbability(self, record):
        return 1.0


def _native_sampler(depth, source=None):
    source = _Source() if source is None else source
    types = [dc.ParticleType.NuMu, dc.ParticleType.NuE, dc.ParticleType.NuTau,
             dc.ParticleType.Gamma]
    models = [_Decay(types[i], types[i + 1], source) for i in range(depth + 1)]
    target = siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 100)), 1.0, 0.0)
    mixture = injection.MultiChannelPhaseSpace([
        injection.DetectorDirected2BodyChannel(target, 0, injection.DirectedMode.Cone)])
    phase_spaces = {models[-1].signature: mixture}
    inj = injection.Injector(
        detector=siren.detector.DetectorModel(), events=2, max_attempts=3, seed=331663,
        primary_type=types[0], primary_interactions=[models[0]],
        primary_injection_distributions=[source],
        primary_phase_spaces=phase_spaces if depth == 0 else {},
        secondary_interactions={types[i]: [models[i]] for i in range(1, depth + 1)},
        secondary_injection_distributions={types[i]: [_FixedSecondary()]
                                           for i in range(1, depth + 1)},
        secondary_phase_spaces={types[depth]: phase_spaces} if depth else {},
        stopping_condition=lambda tree, parent, index: False)
    return inj, source


@pytest.mark.parametrize("depth", [0, 1, 2])
@pytest.mark.parametrize("shortfall", ["warn", "raise", "ignore"])
def test_strict_sampler_failure_stops_immediately_with_report(depth, shortfall):
    inj, source = _native_sampler(depth)
    with pytest.raises(GenerationFailure, match="SamplingFailure") as exc:
        inj.generate(on_failure="raise", on_shortfall=shortfall)
    report = exc.value.report
    assert report.attempts == inj.engine.InjectionAttempts() == 1
    assert report.successes == 0
    assert report.by_vertex[0].reason == injection.FailureReason.SamplingFailure
    assert report.by_vertex[0].depth == depth
    assert len(report.last_failed_tree.tree) == max(1, depth)
    # A new run resets the ledger and succeeds with valid on-shell inputs.
    source.momentum = [10.0, 0.0, 0.0, math.sqrt(96.0)]
    trees = inj.generate(on_failure="raise", on_shortfall="raise")
    assert len(trees) == inj.engine.InjectionAttempts() == 2
    assert inj.report().by_vertex == []


def test_default_retry_policy_keeps_existing_behavior():
    inj, _ = _native_sampler(0)
    assert inj.generate(on_shortfall="ignore") == []
    assert inj.report().attempts == 3
    assert inj.report().dominant().reason == injection.FailureReason.SamplingFailure


def test_strict_policy_rejects_closed_kinematics_conservatively():
    inj, source = _native_sampler(0)
    source.mass = 0.9
    source.momentum = [10.0, 0.0, 0.0, math.sqrt(100.0 - 0.9**2)]
    with pytest.raises(GenerationFailure, match="KinematicallyForbidden"):
        inj.generate(on_failure="raise")
    assert inj.engine.InjectionAttempts() == 1


def test_strict_geometric_misses_use_attempt_budget_and_shortfall_policy():
    inj = _forced_failure_injector(events=2, max_attempts=4)
    with pytest.raises(InjectionShortfall) as exc:
        inj.generate(on_failure="raise", on_shortfall="raise")
    assert exc.value.report.attempts == inj.engine.InjectionAttempts() == 4
    assert inj.engine.FailedEvents() == 4
    assert exc.value.report.dominant().reason == injection.FailureReason.NoTargetsOnPath


def test_invalid_policy_does_not_reset_an_existing_run():
    inj = _forced_failure_injector(events=1, max_attempts=2)
    inj.generate(on_shortfall="ignore")
    with pytest.raises(ValueError, match="on_failure"):
        inj.generate(on_failure="typo")
    assert inj.engine.InjectionAttempts() == 2


def test_generate_facade_propagates_strict_failure_before_weighting():
    class NoWeighting:
        def weight_all(self, trees):
            pytest.fail("weighting must not run after a strict failure")
    inj, _ = _native_sampler(1)
    with pytest.raises(GenerationFailure, match="SamplingFailure"):
        siren.generate(inj, NoWeighting(), events=1, on_failure="raise")


def test_generate_facade_propagates_weight_errors():
    class BrokenWeighter:
        def weight_all(self, trees):
            raise WeightCalculationError("physical density undefined")
    inj, source = _native_sampler(0)
    source.momentum = [10.0, 0.0, 0.0, math.sqrt(96.0)]
    with pytest.raises(WeightCalculationError, match="physical density undefined"):
        siren.generate(inj, BrokenWeighter(), events=1, on_failure="raise")
    assert inj.engine.InjectionAttempts() == 1


def test_strict_and_retry_runs_preserve_geometric_misses_and_event_weights():
    from test_weighter_inheritance import _injector
    strict, retry = _injector(), _injector()
    strict_trees = strict.generate(5, on_failure="raise", on_shortfall="raise")
    retry_trees = retry.generate(5, on_failure="retry", on_shortfall="raise")
    assert strict.engine.InjectionAttempts() == retry.engine.InjectionAttempts()
    assert strict.engine.FailedEvents() == retry.engine.FailedEvents()
    for a, b in zip(strict_trees, retry_trees):
        assert a.tree[0].record.primary_momentum == b.tree[0].record.primary_momentum
        assert a.tree[0].record.interaction_vertex == b.tree[0].record.interaction_vertex
    import numpy as np
    np.testing.assert_array_equal(injection.Weighter(strict).weight_all(strict_trees),
                                  injection.Weighter(retry).weight_all(retry_trees))


def test_strict_secondary_zero_length_path_uses_attempt_budget():
    inj, source = _native_sampler(1)
    source.momentum = [10.0, 0.0, 0.0, math.sqrt(96.0)]
    inj.secondary_injection_distributions = {
        dc.ParticleType.NuE: [d.SecondaryBoundedVertexDistribution(0.0)]}
    inj.secondary_phase_spaces = {}
    assert inj.generate(1, on_failure="raise", on_shortfall="ignore") == []
    assert inj.engine.InjectedEvents() == 0
    assert inj.engine.InjectionAttempts() == inj.engine.FailedEvents() == 3
    assert inj.report().dominant().reason == injection.FailureReason.NoPathThroughVolume


def test_strict_generation_propagates_unexpected_python_errors():
    class BrokenSource(_Source):
        def Sample(self, *args):
            raise RuntimeError("source data corrupt")

    inj, _ = _native_sampler(0, BrokenSource())
    with pytest.raises(RuntimeError, match="source data corrupt"):
        inj.generate(on_failure="raise")


def test_generate_facade_keeps_zero_weight_events():
    from test_weighter_inheritance import _injector
    inj = _injector(chain=False)
    inj.primary.physical = [d.NormalizationConstant(0.0)]
    result = siren.generate(inj, injection.Weighter(inj), events=2,
                            on_failure="raise", on_shortfall="raise")
    assert len(result) == 2
    assert list(result.weights) == [0.0, 0.0]
    assert inj.engine.InjectedEvents() == 2
