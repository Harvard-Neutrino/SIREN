"""External generators through ExternalGeneratorCrossSection.

Two interfaces are tested.

* ExternalGenerator: the generator fills the vertex itself. A Python
  generator produces a small cascade with statuses, parents and typed
  attributes; the test checks that SIREN keeps all of it and that the
  final-state density is declared opaque under the generator's key.

* TrialGenerator + BatchResampler: an unweighting generator makes trials
  against an envelope. The toy draws a variable x with physical density
  p(x) = 2 x from a uniform proposal, so the raw trial weight is r = 2 x;
  with cap c = 1.5 a trial emits with probability min(1, r / c) and weight
  W = max(1, r / c). Then E[W 1{emit}] = 1 / c, the envelope is c * sigma,
  and the emit probability is 0.625. BatchResampler with M trials keeps one
  emitted trial with probability W_j / sum W and attaches the factor
  (sum W / M) * envelope / sigma.

The rate is checked against an independent analytic answer: a monoenergetic
beam along z, injected uniformly in a cylinder inside a uniform argon
sphere, has the thin-target expectation sum(w) = sigma * n_Ar * V, the number
of interactions per unit fluence, with V the cylinder volume. SIREN weights
already carry the 1/N of the attempt count, and lengths are in metres, so the
sum is compared in m^2.
"""
import math

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import distributions as dist
from siren import injection as inj
from siren import interactions as xs
from siren import utilities
from siren import geometry as geo
from siren.math import Vector3D

P = dc.Particle.ParticleType
NuMu, Ar40 = P.NuMu, P.Ar40Nucleus
SIGMA = 1e-36            # cm^2, constant
ENERGY = 2.0             # GeV
CAP = 1.5
EMIT = 0.625             # integral of min(1, 2 x / CAP)
RADIUS, LENGTH = 1.0, 2.0  # m, injection cylinder
DENSITY = 1.4            # g / cm^3 of argon-40
AVOGADRO = 6.02214076e23


def channel():
    s = dc.InteractionSignature()
    s.primary_type = NuMu
    s.target_type = Ar40
    s.secondary_types = [NuMu, P.Hadrons]
    return s


def argon_world():
    dm = siren.detector.DetectorModel()
    materials = dm.Materials
    materials.AddMaterial("ARGON", {1000180400: 1.0})
    dm.Materials = materials
    sector = siren.detector.DetectorSector()
    sector.name = "argon"
    sector.level = 0
    sector.material_id = materials.GetMaterialId("ARGON")
    sector.geo = geo.Sphere(geo.Placement(Vector3D(0, 0, 0)), 100.0, 0.0)
    sector.density = siren.detector.ConstantDensityDistribution(DENSITY)
    dm.AddSector(sector)
    return dm


def fill_two_body(record, x):
    """Outgoing neutrino with energy (1 - x) E along z; the rest is Hadrons."""
    e = record.primary_momentum[0]
    lep = [(1 - x) * e, 0.0, 0.0, (1 - x) * e]
    target = [record.target_mass, 0.0, 0.0, 0.0]
    had = [a + b - c for a, b, c in zip(record.primary_momentum, target, lep)]
    for i, p in enumerate((lep, had)):
        s = record.get_secondary_particle_record(i)
        s.four_momentum = p
        s.mass = 0.0 if i == 0 else math.sqrt(max(p[0] ** 2 - sum(v * v for v in p[1:]), 0.0))
        s.helicity = 0.0


class ToyTrials(xs.TrialGenerator):
    def __init__(self, key="toy-trials"):
        super().__init__()
        self.key = key
        self.last_x = None
        self.kept_x = None
        self.trials = 0

    def GetChannels(self):
        return [channel()]

    def TotalCrossSection(self, c, energy):
        return SIGMA

    def Envelope(self, c, energy):
        return CAP * SIGMA

    def Trial(self, record, rand):
        self.trials += 1
        x = rand.Uniform(0, 1)
        r = 2.0 * x
        if rand.Uniform(0, 1) >= min(1.0, r / CAP):
            return 0.0
        self.last_x = x
        return max(1.0, r / CAP)

    def Keep(self):
        self.kept_x = self.last_x

    def Fill(self, record):
        params = dict(record.interaction_parameters)
        params["x"] = self.kept_x
        record.interaction_parameters = params
        fill_two_body(record, self.kept_x)

    def SamplerKey(self):
        return self.key


class Cascade(xs.ExternalGenerator):
    """Fills mu-like neutrino + struck nucleon -> Delta -> p pi+ with attributes."""

    def GetChannels(self):
        return [channel()]

    def TotalCrossSection(self, c, energy):
        return SIGMA

    def GenerateFinalState(self, record, rand):
        record.secondary_types = [NuMu, P.Neutron, P.Hadrons, P.PPlus, P.PiPlus]
        e = record.primary_momentum[0]
        for i in range(5):
            s = record.get_secondary_particle_record(i)
            s.four_momentum = [0.2 * e, 0.0, 0.0, 0.2 * e]
            s.mass = 0.0
            s.helicity = 0.0
        record.set_secondary_status(1, 21)
        record.set_secondary_status(2, 2)
        record.set_secondary_parent(2, 1)
        record.set_secondary_parent(3, 2)
        record.set_secondary_parent(4, 2)
        record.interaction_int_parameters = {"process_id": 1200}
        record.interaction_string_parameters = {"process_name": "toy RES"}

    def SamplerKey(self):
        return "toy-cascade"


def processes(model, dm, n, seed):
    primary = inj.PrimaryInjectionProcess()
    primary.primary_type = NuMu
    primary.interactions = xs.InteractionCollection(NuMu, [model])
    primary.distributions = [
        dist.PrimaryMass(0), dist.Monoenergetic(ENERGY),
        dist.FixedDirection(Vector3D(0, 0, 1)),
        dist.CylinderVolumePositionDistribution(
            geo.Cylinder(geo.Placement(Vector3D(0, 0, 0)), RADIUS, 0.0, LENGTH)),
    ]
    physical = inj.PhysicalProcess()
    physical.primary_type = NuMu
    physical.interactions = xs.InteractionCollection(NuMu, [model])
    physical.distributions = [dist.PrimaryMass(0), dist.Monoenergetic(ENERGY),
                              dist.FixedDirection(Vector3D(0, 0, 1))]
    generator = inj._Injector(n, dm, primary, utilities.SIREN_random(seed))
    return generator, physical


@pytest.fixture(scope="module")
def batch():
    dm = argon_world()
    trials = ToyTrials()
    resampler = xs.BatchResampler(trials, 8)
    model = xs.ExternalGeneratorCrossSection(resampler)
    n = 3000
    generator, physical = processes(model, dm, n, 41)
    events = [generator.GenerateEvent() for _ in range(n)]
    assert all(len(e.tree) for e in events)
    weighter = inj._Weighter([generator], dm, physical)
    return dict(dm=dm, trials=trials, resampler=resampler, model=model,
                generator=generator, physical=physical, events=events, weighter=weighter)


def test_opaque_under_the_generator_key(batch):
    model = batch["model"]
    assert not model.FinalStateDensityEvaluable(channel())
    assert model.FinalStateSamplerKey(channel()) == "toy-trials|batch=8"
    assert list(model.GetPossibleSignatures()) == [channel()]


def test_batch_attributes_and_mean_factor(batch):
    factors = []
    for e in batch["events"]:
        r = e.tree[0].record
        assert r.interaction_int_parameters["batch_trials"] == 8
        assert 0 <= r.interaction_int_parameters["batch_emitted"] <= 8
        factors.append(r.sampler_factor)
    factors = np.array(factors)
    assert abs(factors.mean() - 1.0) < 4 * factors.std(ddof=1) / np.sqrt(len(factors))


def test_weighted_distribution_is_the_physical_one(batch):
    w = np.array([batch["weighter"].EventWeight(e) for e in batch["events"]])
    x = np.array([e.tree[0].record.interaction_parameters.get("x", -1.0) for e in batch["events"]])
    keep = w > 0
    w, x = w[keep], x[keep]
    edges = np.linspace(0, 1, 11)
    expected = np.diff(edges ** 2)
    observed, sigma = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x < b)
        observed.append(w[m].sum() / w.sum())
        sigma.append(np.sqrt((w[m] ** 2).sum()) / w.sum())
    pull = (np.array(observed) - expected) / np.array(sigma)
    assert np.all(np.abs(pull) < 4), pull


def test_rate_matches_flux_times_cross_section_times_targets(batch):
    w = np.array([batch["weighter"].EventWeight(e) for e in batch["events"]])
    n_ar = DENSITY / 39.9623831 * AVOGADRO                       # cm^-3
    volume = math.pi * (100 * RADIUS) ** 2 * (100 * LENGTH)       # cm^3
    expected = SIGMA * n_ar * volume * 1e-4                      # m^2 per unit fluence
    estimate = w.sum()
    error = np.sqrt((w ** 2).sum())
    assert abs(estimate - expected) < 4 * error, (estimate, expected, error)


def test_empty_batches_are_zero_weight_events():
    dm = argon_world()
    trials = ToyTrials()
    model = xs.ExternalGeneratorCrossSection(xs.BatchResampler(trials, 2))
    n = 2000
    generator, physical = processes(model, dm, n, 43)
    events = [generator.GenerateEvent() for _ in range(n)]
    zero = [e for e in events if e.tree[0].record.sampler_factor == 0.0]
    p_zero = (1 - EMIT) ** 2
    assert abs(len(zero) / n - p_zero) < 4 * math.sqrt(p_zero * (1 - p_zero) / n)
    assert generator.ZeroFactorEvents() == len(zero)
    assert generator.FailedEvents() == 0
    weighter = inj._Weighter([generator], dm, physical)
    assert all(weighter.EventWeight(e) == 0.0 for e in zero)
    assert all(e.tree[0].record.interaction_int_parameters["batch_emitted"] == 0 for e in zero)


def test_a_generator_fills_history_and_attributes():
    dm = argon_world()
    gen = Cascade()
    model = xs.ExternalGeneratorCrossSection(gen)
    generator, physical = processes(model, dm, 20, 47)
    events = [generator.GenerateEvent() for _ in range(20)]
    weighter = inj._Weighter([generator], dm, physical)
    for e in events:
        r = e.tree[0].record
        assert list(r.signature.secondary_types) == [NuMu, P.Hadrons]
        assert list(r.secondary_types) == [NuMu, P.Neutron, P.Hadrons, P.PPlus, P.PiPlus]
        assert list(r.secondary_statuses) == [1, 21, 2, 1, 1]
        assert list(r.secondary_parents) == [-1, -1, 1, 2, 2]
        assert r.interaction_int_parameters["process_id"] == 1200
        assert r.interaction_string_parameters["process_name"] == "toy RES"
        assert np.isfinite(weighter.EventWeight(e)) and weighter.EventWeight(e) > 0


def test_a_different_generator_is_refused_and_archiving_is_explicit(batch, tmp_path):
    other = xs.ExternalGeneratorCrossSection(xs.BatchResampler(ToyTrials(key="other"), 8))
    _, physical = processes(other, batch["dm"], 1, 1)
    from siren.errors import ConfigurationError
    with pytest.raises(ConfigurationError, match="sampler"):
        inj._Weighter([batch["generator"]], batch["dm"], physical)
    # The same trials with another batch size is another sampler, too.
    resized = xs.ExternalGeneratorCrossSection(xs.BatchResampler(batch["trials"], 4))
    _, physical = processes(resized, batch["dm"], 1, 1)
    with pytest.raises(ConfigurationError, match="sampler"):
        inj._Weighter([batch["generator"]], batch["dm"], physical)
    with pytest.raises(RuntimeError, match="cannot be archived"):
        batch["generator"].SaveInjector(str(tmp_path / "external.injector"))


def test_invalid_construction():
    with pytest.raises(ValueError):
        xs.BatchResampler(ToyTrials(), 0)

    class NoKey(Cascade):
        def SamplerKey(self):
            return ""

    with pytest.raises(ValueError, match="SamplerKey"):
        xs.ExternalGeneratorCrossSection(NoKey())
