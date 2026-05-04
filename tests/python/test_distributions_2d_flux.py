"""Tests for Tabulated2DFluxDistribution and related 2D flux bindings."""
import math

import pytest

siren = pytest.importorskip("siren")

# -- Module-scoped imports ---------------------------------------------------

@pytest.fixture(scope="module")
def distributions():
    from siren import distributions
    return distributions

@pytest.fixture(scope="module")
def dataclasses():
    from siren import dataclasses
    return dataclasses

@pytest.fixture(scope="module")
def utilities():
    from siren import utilities
    return utilities

# -- CSV fixture: 3x3 constant flux table -----------------------------------

ENERGIES = [1.0, 2.0, 3.0]
COS_ZENITHS = [-1.0, 0.0, 1.0]

@pytest.fixture
def csv_2d_flux(tmp_path):
    """Create a simple 3x3 flux table with constant flux=1.0."""
    p = tmp_path / "flux2d.dat"
    lines = ["# energy coszenith flux"]
    for e in ENERGIES:
        for cz in COS_ZENITHS:
            lines.append(f"{e} {cz} 1.0")
    p.write_text("\n".join(lines) + "\n")
    return str(p)

# -- Construction ------------------------------------------------------------

class TestTabulated2DFluxDistributionConstruction:
    def test_construct_from_file(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        assert dist is not None

    def test_construct_with_energy_bounds(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(
            energyMin=1.0, energyMax=3.0, fluxTableFilename=csv_2d_flux)
        assert dist is not None

    def test_construct_with_all_bounds(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(
            energyMin=1.0, energyMax=3.0,
            cosZenithMin=-1.0, cosZenithMax=1.0,
            fluxTableFilename=csv_2d_flux)
        assert dist is not None

    def test_construct_from_vectors(self, distributions):
        energies = [e for e in ENERGIES for _ in COS_ZENITHS]
        cos_zeniths = COS_ZENITHS * len(ENERGIES)
        flux = [1.0] * 9
        dist = distributions.Tabulated2DFluxDistribution(
            energies=energies, cosZeniths=cos_zeniths, flux=flux)
        assert dist is not None

# -- Metadata ----------------------------------------------------------------

class TestTabulated2DFluxDistributionMetadata:
    def test_name(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        assert dist.Name() == "Tabulated2DFluxDistribution"

    def test_get_energy_nodes(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        expected = [e for e in ENERGIES for _ in COS_ZENITHS]
        assert dist.GetEnergyNodes() == pytest.approx(expected)

    def test_get_coszenith_nodes(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        expected = COS_ZENITHS * len(ENERGIES)
        assert dist.GetCosZenithNodes() == pytest.approx(expected)

# -- Integral ----------------------------------------------------------------

class TestTabulated2DFluxDistributionIntegral:
    def test_integral_constant_flux(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        # constant flux=1 over [1,3] x [-1,1] => integral = 2 * 2 = 4
        assert dist.GetIntegral() == pytest.approx(4.0, rel=0.1)

    def test_sample_pdf_matches_integral(self, distributions, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        integral = dist.GetIntegral()
        assert dist.SamplePDF(2.0, 0.0) == pytest.approx(1.0 / integral, rel=0.1)

# -- Sampling ----------------------------------------------------------------

def _new_record(dataclasses):
    """Helper: create a fresh PrimaryDistributionRecord for sampling."""
    rec = dataclasses.PrimaryDistributionRecord(dataclasses.ParticleType.NuMu)
    rec.energy = 1.0
    rec.mass = 0.0
    rec.direction = [0.0, 0.0, 1.0]
    return rec

class TestTabulated2DFluxDistributionSampling:
    def test_sample_energy_in_bounds(self, distributions, dataclasses,
                                     utilities, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        rand = utilities.SIREN_random()
        for _ in range(500):
            rec = _new_record(dataclasses)
            dist.Sample(rand, None, None, rec)
            assert 1.0 <= rec.energy <= 3.0

    def test_sample_direction_is_unit(self, distributions, dataclasses,
                                      utilities, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        rand = utilities.SIREN_random()
        for _ in range(500):
            rec = _new_record(dataclasses)
            dist.Sample(rand, None, None, rec)
            d = rec.direction
            mag = math.sqrt(d[0]**2 + d[1]**2 + d[2]**2)
            assert mag == pytest.approx(1.0, abs=1e-6)

    def test_sample_coszenith_in_bounds(self, distributions, dataclasses,
                                        utilities, csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        rand = utilities.SIREN_random()
        for _ in range(500):
            rec = _new_record(dataclasses)
            dist.Sample(rand, None, None, rec)
            assert -1.0 <= rec.direction[2] <= 1.0

# -- GenerationProbability ---------------------------------------------------

class TestTabulated2DFluxDistributionProbability:
    def test_in_bounds_returns_nonzero(self, distributions, dataclasses,
                                       csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        record = dataclasses.InteractionRecord()
        record.primary_momentum = [2.0, 0.0, 0.0, 1.0]
        assert dist.GenerationProbability(None, None, record) > 0.0

    def test_below_energy_returns_zero(self, distributions, dataclasses,
                                       csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        record = dataclasses.InteractionRecord()
        record.primary_momentum = [0.5, 0.0, 0.0, 1.0]
        assert dist.GenerationProbability(None, None, record) == 0.0

    def test_above_energy_returns_zero(self, distributions, dataclasses,
                                       csv_2d_flux):
        dist = distributions.Tabulated2DFluxDistribution(csv_2d_flux)
        record = dataclasses.InteractionRecord()
        record.primary_momentum = [5.0, 0.0, 0.0, 1.0]
        assert dist.GenerationProbability(None, None, record) == 0.0

# -- TabulatedFluxDistribution backward compatibility ------------------------

class TestTabulatedFluxDistributionBackwardCompat:
    def test_construct_without_romberg_arg(self, distributions):
        dist = distributions.TabulatedFluxDistribution(
            [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], False)
        assert dist is not None

    def test_construct_with_romberg_arg(self, distributions):
        dist = distributions.TabulatedFluxDistribution(
            [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], False, False)
        assert dist is not None

# -- Base class existence ----------------------------------------------------

class TestPrimaryEnergyDirectionDistribution:
    def test_base_class_exists(self, distributions):
        assert hasattr(distributions, "PrimaryEnergyDirectionDistribution")
        assert isinstance(
            distributions.PrimaryEnergyDirectionDistribution, type)
