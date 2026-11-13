"""A sampler may attach a physical factor to the vertex it fills.

Many generators do not draw exactly from the physical density. A partially
unweighted generator draws from an envelope, accepts with probability
min(1, f/c) and returns a weight max(1, f/c) for points where f exceeds the
cap c. The weight carries what the acceptance could not.

The toy here draws one variable y in [0, 1] from f(y) = 2 y with cap c = 1.5
under a uniform envelope. Accepted points have density min(f, c) / A with
A = integral of min(f, c) = 0.9375, so the factor A max(1, f/c) has mean one
and the weighted distribution of y is f. The model's final-state density is
opaque (it returns a constant FinalStateProbability), which is the case the
factor exists for.

The test checks the mean factor, the weighted y distribution and the rate
against an exact sampler, the zero-factor and invalid-factor behaviour, and
that the factor survives pickle, HepMC3 and pooled reweighting without
calling the sampler again.
"""
import os
import pickle

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import distributions as dist
from siren import injection as inj
from siren import interactions as xs
from siren import utilities
from siren.math import Vector3D

NuMu = dc.Particle.ParticleType.NuMu
CAP = 1.5
ACCEPTED_MASS = 0.9375  # integral over [0, 1] of min(2 y, CAP)
N = 4000


def f(y):
    return 2.0 * y


def _detector_model():
    from siren import detector
    from siren import _util

    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except (ValueError, RuntimeError) as e:
        pytest.skip(f"CCM detector model unavailable: {e}")
    return dm


class Toy(xs.CrossSection):
    """DummyCrossSection rates with a toy final-state variable y.

    mode "partial": partially unweighted sampling with the factor;
    mode "exact":   y drawn exactly from f by inversion, no factor;
    mode "zero":    as partial, but factor zero below y = 0.25;
    mode "negative": reports a negative factor.
    """

    def __init__(self, mode="partial"):
        super().__init__()
        self._d = xs.DummyCrossSection()
        self.mode = mode
        self.calls = 0

    def equal(self, other):
        return other is self

    def GetPossiblePrimaries(self):
        return self._d.GetPossiblePrimaries()

    def GetPossibleTargets(self):
        return self._d.GetPossibleTargets()

    def GetPossibleTargetsFromPrimary(self, p):
        return self._d.GetPossibleTargetsFromPrimary(p)

    def GetPossibleSignatures(self):
        return self._d.GetPossibleSignatures()

    def GetPossibleSignaturesFromParents(self, p, t):
        return self._d.GetPossibleSignaturesFromParents(p, t)

    def InteractionThreshold(self, r):
        return self._d.InteractionThreshold(r)

    def DensityVariables(self):
        return ["toy_y"]

    def SampleFinalState(self, r, rand):
        self.calls += 1
        self._d.SampleFinalState(r, rand)
        if self.mode == "exact":
            y = np.sqrt(rand.Uniform(0, 1))
            factor = 1.0
        else:
            while True:
                y = rand.Uniform(0, 1)
                if rand.Uniform(0, 1) * CAP < f(y):
                    break
            factor = ACCEPTED_MASS * max(1.0, f(y) / CAP)
            if self.mode == "zero" and y < 0.25:
                factor = 0.0
            if self.mode == "negative":
                factor = -factor
        params = dict(r.interaction_parameters)
        params["toy_y"] = y
        r.interaction_parameters = params
        r.sampler_factor = factor

    def TotalCrossSection(self, r):
        return self._d.TotalCrossSection(r)

    def TotalCrossSectionAllFinalStates(self, r):
        return self._d.TotalCrossSection(r)

    def DifferentialCrossSection(self, r):
        return self._d.TotalCrossSection(r)

    def FinalStateProbability(self, r):
        return 1.0

    def FinalStateDensityEvaluable(self, signature):
        return False

    def FinalStateSamplerKey(self, signature):
        return "toy-generator"


def _processes(model):
    primary = inj.PrimaryInjectionProcess()
    primary.primary_type = NuMu
    primary.interactions = xs.InteractionCollection(NuMu, [model])
    primary.distributions = [
        dist.PrimaryMass(0), dist.PowerLaw(2.0, 0.5, 5.0), dist.IsotropicDirection(),
        dist.PointSourcePositionDistribution(Vector3D(0, 0, 0), 25.0)]
    physical = inj.PhysicalProcess()
    physical.primary_type = NuMu
    physical.interactions = xs.InteractionCollection(NuMu, [model])
    physical.distributions = [dist.PrimaryMass(0), dist.IsotropicDirection()]
    return primary, physical


def _run(dm, model, n, seed):
    primary, physical = _processes(model)
    generator = inj._Injector(n, dm, primary, utilities.SIREN_random(seed))
    events = [generator.GenerateEvent() for _ in range(n)]
    return generator, physical, [e for e in events if len(e.tree)]


def _y(event):
    return event.tree[0].record.interaction_parameters["toy_y"]


@pytest.fixture(scope="module")
def partial():
    dm = _detector_model()
    model = Toy("partial")
    generator, physical, events = _run(dm, model, N, 101)
    weighter = inj._Weighter([generator], dm, physical)
    return dict(dm=dm, model=model, generator=generator, physical=physical,
                events=events, weighter=weighter)


def test_mean_factor_is_one(partial):
    factors = np.array([e.tree[0].record.sampler_factor for e in partial["events"]])
    assert np.all(factors >= ACCEPTED_MASS) and np.any(factors > ACCEPTED_MASS)
    mean, err = factors.mean(), factors.std(ddof=1) / np.sqrt(len(factors))
    assert abs(mean - 1.0) < 4 * err


def test_weighted_distribution_is_the_physical_one(partial):
    weighter, events = partial["weighter"], partial["events"]
    w = np.array([weighter.EventWeight(e) for e in events])
    y = np.array([_y(e) for e in events])
    edges = np.linspace(0, 1, 11)
    expected = np.diff(edges ** 2)  # integral of 2 y over each bin
    observed = np.array([w[(y >= a) & (y < b)].sum() for a, b in zip(edges[:-1], edges[1:])]) / w.sum()
    sigma = np.array([np.sqrt((w[(y >= a) & (y < b)] ** 2).sum()) for a, b in zip(edges[:-1], edges[1:])]) / w.sum()
    pull = (observed - expected) / sigma
    assert np.all(np.abs(pull) < 4), pull
    # Without the factor the distribution is the accepted one, min(f, c) / A,
    # which differs visibly in the top bins.
    w0 = w / np.array([e.tree[0].record.sampler_factor for e in events])
    top = (y >= 0.9)
    assert w0[top].sum() / w0.sum() < 0.9 * (w[top].sum() / w.sum())


def test_rate_agrees_with_an_exact_sampler(partial):
    dm = partial["dm"]
    exact = Toy("exact")
    generator, physical, events = _run(dm, exact, N, 202)
    exact_weighter = inj._Weighter([generator], dm, physical)
    w_exact = np.array([exact_weighter.EventWeight(e) for e in events])
    w_partial = np.array([partial["weighter"].EventWeight(e) for e in partial["events"]])
    # Each run estimates the same rate as sum(w).
    r_exact, r_partial = w_exact.sum(), w_partial.sum()
    err = np.hypot(np.sqrt((w_exact ** 2).sum()), np.sqrt((w_partial ** 2).sum()))
    assert abs(r_exact - r_partial) < 4 * err
    assert all(e.tree[0].record.sampler_factor == 1.0 for e in events)


def test_explanation_reports_the_factor(partial):
    event = partial["events"][0]
    breakdown = partial["weighter"].EventWeightWithBreakdown(event)
    vertex = breakdown.vertices[0]
    assert vertex.sampler_factor == event.tree[0].record.sampler_factor
    assert breakdown.total == pytest.approx(partial["weighter"].EventWeight(event), rel=1e-15)


def test_zero_factor_is_a_zero_weight_event_not_a_failure():
    dm = _detector_model()
    model = Toy("zero")
    generator, physical, events = _run(dm, model, 400, 303)
    zeros = [e for e in events if e.tree[0].record.sampler_factor == 0.0]
    assert len(zeros) > 0
    assert generator.ZeroFactorEvents() == len(zeros)
    assert generator.FailedEvents() == 0
    weighter = inj._Weighter([generator], dm, physical)
    for e in zeros:
        assert weighter.EventWeight(e) == 0.0
        flags = weighter.EventWeightWithBreakdown(e).vertices[0].flags
        assert "sampler factor zero (weight 0)" in flags


def test_negative_factor_is_refused():
    dm = _detector_model()
    primary, _ = _processes(Toy("negative"))
    generator = inj._Injector(4, dm, primary, utilities.SIREN_random(9))
    with pytest.raises(ValueError, match="non-negative"):
        generator.GenerateEvent()


def test_factor_survives_round_trips_without_resampling(partial, tmp_path):
    model, weighter = partial["model"], partial["weighter"]
    events = partial["events"][:200]
    calls = model.calls
    w = [weighter.EventWeight(e) for e in events]
    restored = [pickle.loads(pickle.dumps(e)) for e in events]
    assert [weighter.EventWeight(e) for e in restored] == w
    path = str(tmp_path / "toy.siren_events")
    dc.SaveInteractionTrees(events, path)
    assert [weighter.EventWeight(e) for e in dc.LoadInteractionTrees(path)] == w
    from siren import hepmc3 as sio
    out = str(tmp_path / "toy.hepmc3")
    try:
        sio.SaveInteractionTreesAsHepMC3(events, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise
    back = sio.LoadInteractionTreesFromHepMC3(out)
    # HepMC3's text format keeps 15 significant digits of a double attribute.
    np.testing.assert_allclose([e.tree[0].record.sampler_factor for e in back],
                               [e.tree[0].record.sampler_factor for e in events], rtol=1e-12)
    np.testing.assert_allclose([weighter.EventWeight(e) for e in back], w, rtol=1e-9)
    assert model.calls == calls


def test_pooled_reweighting_uses_the_stored_factor(partial):
    dm, model = partial["dm"], partial["model"]
    primary, physical = _processes(model)
    second = inj._Injector(N, dm, primary, utilities.SIREN_random(404))
    second_events = [second.GenerateEvent() for _ in range(1000)]
    calls = model.calls
    pooled = inj._Weighter([partial["generator"], second], dm, physical)
    for e in partial["events"][:200] + second_events[:200]:
        breakdown = pooled.EventWeightWithBreakdown(e)
        assert np.isfinite(breakdown.total) and breakdown.total > 0
        assert all(v.sampler_factor == e.tree[0].record.sampler_factor for v in breakdown.vertices)
    assert model.calls == calls
