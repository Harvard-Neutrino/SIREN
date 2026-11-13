"""Opaque final states: a model may declare its final-state density not evaluable.

An event generator samples a final state it cannot evaluate a density for. Its
FinalStateProbability is then a placeholder that is only correct if it cancels
between the generation and the physical side of the weight. The weighter
therefore refuses any configuration in which it would not cancel:

* a phase-space proposal registered for the opaque channel on either side;
* another interaction producing the same channel on either side;
* a different final-state sampler on the physical side, unless both models
  name the same sampler through a non-empty FinalStateSamplerKey;
* an evaluable model on one side and the opaque one on the other.

What remains allowed is reweighting everything else: the flux, the position
distribution and the total cross section (through a second model object that
names the same sampler). Each allowed case is checked against an evaluable
twin whose FinalStateProbability is the same constant, so the weights must
agree exactly.
"""
import os

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import distributions as dist
from siren import injection as inj
from siren import interactions as xs
from siren import utilities
from siren.errors import ConfigurationError
from siren.math import Vector3D

NuMu = dc.Particle.ParticleType.NuMu
N_INJECT = 48


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


class Wrapped(xs.CrossSection):
    """DummyCrossSection with a constant final-state density placeholder.

    opaque declares the density not evaluable; key names the sampler;
    scale multiplies the total cross section.
    """

    def __init__(self, opaque=True, key="", scale=1.0):
        super().__init__()
        self._d = xs.DummyCrossSection()
        self.opaque, self.key, self.scale = opaque, key, scale

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
        return self._d.DensityVariables()

    def SampleFinalState(self, r, rand):
        self._d.SampleFinalState(r, rand)

    def TotalCrossSection(self, r):
        return self.scale * self._d.TotalCrossSection(r)

    def TotalCrossSectionAllFinalStates(self, r):
        return self.scale * self._d.TotalCrossSection(r)

    def DifferentialCrossSection(self, r):
        return self.TotalCrossSection(r)

    def FinalStateProbability(self, r):
        return 1.0

    def FinalStateDensityEvaluable(self, signature):
        return not self.opaque

    def FinalStateSamplerKey(self, signature):
        return self.key


def _primary(model, alpha=2.0, radius=25.0):
    process = inj.PrimaryInjectionProcess()
    process.primary_type = NuMu
    process.interactions = xs.InteractionCollection(NuMu, [model] if not isinstance(model, list) else model)
    process.distributions = [
        dist.PrimaryMass(0),
        dist.PowerLaw(alpha, 0.5, 5.0),
        dist.IsotropicDirection(),
        dist.PointSourcePositionDistribution(Vector3D(0, 0, 0), radius),
    ]
    return process


def _physical(model, alpha=None, radius=None):
    process = inj.PhysicalProcess()
    process.primary_type = NuMu
    process.interactions = xs.InteractionCollection(NuMu, [model] if not isinstance(model, list) else model)
    distributions = [dist.PrimaryMass(0), dist.IsotropicDirection()]
    if alpha is not None:
        distributions.append(dist.PowerLaw(alpha, 0.5, 5.0))
    if radius is not None:
        distributions.append(dist.PointSourcePositionDistribution(Vector3D(0, 0, 0), radius))
    process.distributions = distributions
    return process


@pytest.fixture(scope="module")
def setup():
    dm = _detector_model()
    opaque = Wrapped(opaque=True, key="dummy-generator")
    twin = Wrapped(opaque=False)
    gen_opaque = inj._Injector(N_INJECT, dm, _primary(opaque), utilities.SIREN_random(77))
    gen_twin = inj._Injector(N_INJECT, dm, _primary(twin), utilities.SIREN_random(77))
    events, twin_events = [], []
    for _ in range(N_INJECT):
        events.append(gen_opaque.GenerateEvent())
        twin_events.append(gen_twin.GenerateEvent())
    events = [e for e in events if len(e.tree)]
    if len(events) < N_INJECT // 2:
        pytest.skip("too few events")
    return dict(dm=dm, opaque=opaque, twin=twin, gen_opaque=gen_opaque,
                gen_twin=gen_twin, events=events)


def _weights(weighter, events):
    w = np.array([weighter.EventWeight(e) for e in events])
    assert np.all(np.isfinite(w)) and np.all(w > 0)
    return w


@pytest.mark.parametrize("alpha,radius", [(None, None), (1.5, None), (None, 40.0), (1.5, 40.0)])
def test_flux_and_position_reweighting_match_an_evaluable_twin(setup, alpha, radius):
    dm = setup["dm"]
    w_opaque = _weights(inj._Weighter([setup["gen_opaque"]], dm,
                                      _physical(setup["opaque"], alpha, radius)), setup["events"])
    w_twin = _weights(inj._Weighter([setup["gen_twin"]], dm,
                                    _physical(setup["twin"], alpha, radius)), setup["events"])
    np.testing.assert_array_equal(w_opaque, w_twin)


def test_total_cross_section_reweighting_through_a_shared_sampler_key(setup):
    dm = setup["dm"]
    doubled = Wrapped(opaque=True, key="dummy-generator", scale=2.0)
    doubled_twin = Wrapped(opaque=False, scale=2.0)
    w_opaque = _weights(inj._Weighter([setup["gen_opaque"]], dm, _physical(doubled)), setup["events"])
    w_twin = _weights(inj._Weighter([setup["gen_twin"]], dm, _physical(doubled_twin)), setup["events"])
    np.testing.assert_array_equal(w_opaque, w_twin)
    base = _weights(inj._Weighter([setup["gen_opaque"]], dm, _physical(setup["opaque"])), setup["events"])
    # A doubled cross section in a thin target nearly doubles every weight.
    np.testing.assert_allclose(w_opaque / base, 2.0, rtol=1e-3)


def _refused(setup, physical, generator=None, match=None):
    with pytest.raises(ConfigurationError, match=match) as info:
        inj._Weighter([generator or setup["gen_opaque"]], setup["dm"], physical)
    message = str(info.value)
    # The error names the vertex: its primary and its channel.
    assert "NuMu" in message or "14" in message
    return message


def test_a_different_sampler_without_a_shared_key_is_refused(setup):
    _refused(setup, _physical(Wrapped(opaque=True, key="")), match="sampler")
    _refused(setup, _physical(Wrapped(opaque=True, key="another-generator")), match="sampler")


def test_an_evaluable_model_on_one_side_is_refused(setup):
    _refused(setup, _physical(setup["twin"]), match="not evaluable")
    with pytest.raises(ConfigurationError, match="not evaluable"):
        inj._Weighter([setup["gen_twin"]], setup["dm"], _physical(setup["opaque"]))


def test_a_mixture_with_another_producer_is_refused(setup):
    _refused(setup, _physical([setup["opaque"], Wrapped(opaque=False)]), match="produce")


def test_a_proposal_for_the_opaque_channel_is_refused(setup):
    model = Wrapped(opaque=True, key="dummy-generator")
    process = _primary(model)
    for signature in model.GetPossibleSignaturesFromParents(NuMu, dc.Particle.ParticleType.Ar40Nucleus):
        process.SetPhaseSpace(signature, inj.MultiChannelPhaseSpace(
            [inj.PhysicalCrossSectionChannel(model, signature)], [1.0]))
    generator = inj._Injector(N_INJECT, setup["dm"], process, utilities.SIREN_random(5))
    _refused(setup, _physical(model), generator=generator, match="proposal")


def test_the_native_defaults_are_evaluable_with_no_key():
    model = xs.DummyCrossSection()
    for signature in model.GetPossibleSignatures()[:3]:
        assert model.FinalStateDensityEvaluable(signature)
        assert model.FinalStateSamplerKey(signature) == ""
