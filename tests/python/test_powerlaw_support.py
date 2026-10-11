"""Sampling support and normalization of native power-law sources."""
import math

import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc, distributions as dist, injection as inj
from siren import interactions as xs, math as smath, utilities

NuMu = dc.Particle.ParticleType.NuMu


@pytest.fixture(scope="module")
def detector():
    # Uniform argon sphere, so no detector dataset is needed.
    dm = siren.detector.DetectorModel()
    materials = dm.Materials
    materials.AddMaterial("ARGON", {1000180400: 1.0})
    dm.Materials = materials
    sector = siren.detector.DetectorSector()
    sector.name = "argon"
    sector.level = 0
    sector.material_id = materials.GetMaterialId("ARGON")
    sector.geo = siren.geometry.Sphere(
        siren.geometry.Placement(smath.Vector3D(0, 0, 0)), 100.0, 0.0)
    sector.density = siren.detector.ConstantDensityDistribution(1.4)
    dm.AddSector(sector)
    return dm


@pytest.fixture(scope="module")
def model():
    return xs.TrivialCrossSection(1.0e-38, [NuMu], [dc.Particle.ParticleType.Nucleon])


def _injector(detector, model, seed, energy_min, energy_max, n_inject=32):
    process = inj.PrimaryInjectionProcess()
    process.primary_type = NuMu
    process.interactions = xs.InteractionCollection(NuMu, [model])
    process.distributions = [
        dist.PrimaryMass(0),
        dist.PowerLaw(2.0, energy_min, energy_max),
        dist.PrimaryNeutrinoHelicityDistribution(),
        dist.IsotropicDirection(),
        dist.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    return inj._Injector(n_inject, detector, process, utilities.SIREN_random(seed))


def _physical(model, energy_min, energy_max):
    physical = inj.PhysicalProcess()
    physical.primary_type = NuMu
    physical.interactions = xs.InteractionCollection(NuMu, [model])
    physical.distributions = [dist.PowerLaw(2.0, energy_min, energy_max),
                              dist.PrimaryNeutrinoHelicityDistribution(),
                              dist.IsotropicDirection()]
    return physical


def _disjoint_injectors(detector, model):
    return [_injector(detector, model, seed, lo, hi)
            for seed, lo, hi in ((13, 1.0, 2.0), (14, 3.0, 4.0))]


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


def test_single_energy_source_matches_monoenergetic():
    spectrum = dist.PowerLaw(2.0, 3.0, 3.0)
    reference = dist.Monoenergetic(3.0)
    for offset in (0.0, 1e-9, -1e-9, 1e-5, -1e-5):
        energy = 3.0 * (1 + offset)
        expected = 1.0 if abs(offset) < 1e-6 else 0.0
        assert spectrum.pdf(energy) == reference.pdf(energy) == expected


@pytest.mark.parametrize("index", [0.0, 1.0, 2.0])
def test_density_normalizes_and_sampling_matches_cdf(index):
    spectrum = dist.PowerLaw(index, 1.0, 2.0)
    steps = 10000
    integral = sum(spectrum.pdf(1.0 + (i + 0.5) / steps)
                   for i in range(steps)) / steps
    assert integral == pytest.approx(1.0, rel=3e-9)
    random = utilities.SIREN_random(4921)
    reference = utilities.SIREN_random(4921)
    record = dc.PrimaryDistributionRecord(NuMu)
    for _ in range(100):
        u = reference.Uniform(0.0, 1.0)
        energy = spectrum.SampleEnergy(random, None, None, record)
        cdf = math.log(energy) / math.log(2.0) if index == 1.0 else (
            (energy ** (1 - index) - 1) / (2 ** (1 - index) - 1))
        assert cdf == pytest.approx(u, rel=2e-13, abs=2e-15)


@pytest.mark.parametrize("index", [1.0, 1.01, 2.58])
def test_endpoint_rounding_stays_in_support(index):
    # Without the clamp, inverse-CDF rounding puts many samples outside this
    # four-ulp interval (all of them for index 1.01).
    energy_min = 30.0
    energy_max = energy_min + 4 * math.ulp(energy_min)
    spectrum = dist.PowerLaw(index, energy_min, energy_max)
    random = utilities.SIREN_random(4921)
    record = dc.PrimaryDistributionRecord(NuMu)
    for _ in range(1000):
        energy = spectrum.SampleEnergy(random, None, None, record)
        assert energy_min <= energy <= energy_max
        assert spectrum.pdf(energy) > 0.0


@pytest.mark.parametrize("index, energy_min, energy_max", [
    (2.0, 4.0, 1.0),          # reversed
    (2.0, -1.0, 1.0),         # negative energy
    (2.0, math.nan, 1.0),
    (2.0, 1.0, math.nan),
    (math.nan, 1.0, 2.0),
    (1.0, 0.0, 2.0),          # not normalizable at zero energy
    (2.0, 0.0, 2.0),
    (1.0, 1.0, math.inf),     # not normalizable at infinite energy
    (0.5, 1.0, math.inf),
])
def test_invalid_bounds_are_rejected(index, energy_min, energy_max):
    with pytest.raises(ValueError):
        dist.PowerLaw(index, energy_min, energy_max)


def test_normalizable_open_bounds_are_accepted():
    assert dist.PowerLaw(0.5, 0.0, 2.0).pdf(1.0) > 0.0
    assert dist.PowerLaw(2.0, 1.0, math.inf).pdf(1.0) == pytest.approx(1.0)


@pytest.mark.parametrize("reference_energy", [0.5, 1.5, 4.0])
def test_reference_normalization_can_be_outside_sampling_bounds(reference_energy):
    spectrum = dist.PowerLaw(2.0, 1.0, 2.0)
    spectrum.SetNormalizationAtEnergy(7.0, reference_energy)
    physical_density = spectrum.normalization * spectrum.pdf(1.5)
    assert physical_density == pytest.approx(7.0 * (reference_energy / 1.5) ** 2)


def test_disjoint_energy_injectors_do_not_dilute_each_others_weights(detector, model):
    injectors = _disjoint_injectors(detector, model)
    physical = _physical(model, 1.0, 4.0)
    pooled = inj._Weighter(injectors, detector, physical)
    for injector in injectors:
        single = inj._Weighter([injector], detector, physical)
        for _ in range(32):
            event = injector.GenerateEvent()
            assert pooled.EventWeight(event) == pytest.approx(single.EventWeight(event), rel=1e-12, abs=0)


def test_physical_flux_outside_its_range_gives_zero_not_nan(detector, model):
    # Above 3.5, both the first injector and the physical flux have zero density.
    injectors = _disjoint_injectors(detector, model)
    physical = _physical(model, 1.0, 3.5)
    pooled = inj._Weighter(injectors, detector, physical)
    outside = 0
    for injector in injectors:
        single = inj._Weighter([injector], detector, physical)
        for _ in range(32):
            event = injector.GenerateEvent()
            energy = event.tree[0].record.primary_momentum[0]
            if energy > 3.5:
                outside += 1
                assert single.EventWeight(event) == 0.0
                assert pooled.EventWeight(event) == 0.0
            else:
                assert single.EventWeight(event) > 0.0
                assert pooled.EventWeight(event) == pytest.approx(single.EventWeight(event), rel=1e-12, abs=0)
    assert outside > 0
