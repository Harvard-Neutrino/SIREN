"""Sampling support and normalization of native power-law sources."""
import math

import pytest
import siren
from siren import dataclasses as dc, distributions as dist, injection as inj
from siren import interactions as xs, utilities


@pytest.mark.parametrize("index", [0.0, 1.0, 2.0])
def test_density_has_the_sampled_support(index):
    spectrum = dist.PowerLaw(index, 1.0, 2.0)
    record = dc.InteractionRecord()
    for energy in (0.0, 0.5, math.nextafter(1.0, 0.0),
                   math.nextafter(2.0, math.inf), 4.0, math.inf, math.nan):
        record.primary_momentum = [energy, 0, 0, energy]
        assert spectrum.pdf(energy) == 0.0
        assert spectrum.GenerationProbability(None, None, record) == 0.0
    for energy in (1.0, 1.5, 2.0):
        assert spectrum.pdf(energy) > 0.0


def test_degenerate_interval_only_accepts_its_energy():
    spectrum = dist.PowerLaw(2.0, 3.0, 3.0)
    assert spectrum.pdf(3.0) == 1.0
    assert spectrum.pdf(math.nextafter(3.0, 0.0)) == 0.0
    assert spectrum.pdf(math.nextafter(3.0, math.inf)) == 0.0


@pytest.mark.parametrize("index", [0.0, 1.0, 2.0])
def test_density_normalizes_and_sampling_matches_cdf(index):
    spectrum = dist.PowerLaw(index, 1.0, 2.0)
    steps = 10000
    integral = sum(spectrum.pdf(1.0 + (i + 0.5) / steps)
                   for i in range(steps)) / steps
    assert integral == pytest.approx(1.0, rel=3e-9)
    random = utilities.SIREN_random(4921)
    reference = utilities.SIREN_random(4921)
    record = dc.PrimaryDistributionRecord(dc.Particle.ParticleType.NuMu)
    for _ in range(100):
        u = reference.Uniform(0.0, 1.0)
        energy = spectrum.SampleEnergy(random, None, None, record)
        cdf = math.log(energy) / math.log(2.0) if index == 1.0 else (
            (energy ** (1 - index) - 1) / (2 ** (1 - index) - 1))
        assert cdf == pytest.approx(u, rel=2e-13, abs=2e-15)


@pytest.mark.parametrize("reference_energy", [0.5, 1.5, 4.0])
def test_reference_normalization_can_be_outside_sampling_bounds(reference_energy):
    spectrum = dist.PowerLaw(2.0, 1.0, 2.0)
    spectrum.SetNormalizationAtEnergy(7.0, reference_energy)
    physical_density = spectrum.normalization * spectrum.pdf(1.5)
    assert physical_density == pytest.approx(7.0 * (reference_energy / 1.5) ** 2)


def test_disjoint_energy_injectors_do_not_dilute_each_others_weights():
    from test_trivial_cross_section import detector_model, _constant_cross_section, _make_injector
    detector = detector_model.__wrapped__()
    model = _constant_cross_section()
    generators = []
    for seed, lo, hi in ((13, 1.0, 2.0), (14, 3.0, 4.0)):
        generator, distributions = _make_injector(detector, model, seed, 32)
        process = generator.GetPrimaryProcess()
        process.distributions = [dist.PowerLaw(2.0, lo, hi) if isinstance(d, dist.PowerLaw)
                                 else d for d in distributions]
        generators.append(generator)
    physical = inj.PhysicalProcess()
    physical.primary_type = dc.Particle.ParticleType.NuMu
    physical.interactions = xs.InteractionCollection(physical.primary_type, [model])
    physical.distributions = [dist.PowerLaw(2.0, 1.0, 4.0),
                              dist.PrimaryNeutrinoHelicityDistribution(),
                              dist.IsotropicDirection()]
    pooled = inj._Weighter(generators, detector, physical)
    batches = [[g.GenerateEvent() for _ in range(32)] for g in generators]
    for generator, events in zip(generators, batches):
        single = inj._Weighter([generator], detector, physical)
        for event in events:
            assert len(event.tree) == 1
            assert pooled.EventWeight(event) == pytest.approx(single.EventWeight(event), rel=1e-12, abs=0)
