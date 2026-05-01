"""Tests for PrimaryExternalDistribution and vertex distribution bindings."""
import os
import textwrap

import pytest

siren = pytest.importorskip("siren")


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


# ---------------------------------------------------------------------------
# CSV fixture helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def csv_basic(tmp_path):
    """CSV with energy and mass columns."""
    p = tmp_path / "basic.csv"
    p.write_text("E,m\n10.0,1.0\n20.0,2.0\n30.0,3.0\n")
    return str(p)


@pytest.fixture
def csv_full(tmp_path):
    """CSV with position, momentum, energy, and mass."""
    p = tmp_path / "full.csv"
    p.write_text("E,m,x0,y0,z0,px,py,pz\n"
                 "10.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0\n")
    return str(p)


@pytest.fixture
def csv_custom_params(tmp_path):
    """CSV with a custom interaction parameter column."""
    p = tmp_path / "custom.csv"
    p.write_text("E,Q2\n10.0,1.5\n20.0,2.5\n")
    return str(p)


@pytest.fixture
def csv_no_energy(tmp_path):
    """CSV without an E column."""
    p = tmp_path / "no_energy.csv"
    p.write_text("x0,y0,z0\n1.0,2.0,3.0\n4.0,5.0,6.0\n")
    return str(p)


@pytest.fixture
def csv_with_blanks_and_comments(tmp_path):
    """CSV with blank lines and comment lines."""
    p = tmp_path / "messy.csv"
    p.write_text("E,m\n# comment line\n10.0,1.0\n\n20.0,2.0\n")
    return str(p)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestPrimaryExternalDistributionConstruction:
    def test_construct_from_file(self, distributions, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic)
        assert dist is not None

    def test_construct_with_emin(self, distributions, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic, 15.0)
        # Only E=20 and E=30 pass the E >= 15 filter
        assert dist.GetPhysicalNumEvents() == 2

    def test_construct_missing_file_raises(self, distributions):
        with pytest.raises(RuntimeError, match="file open failed"):
            distributions.PrimaryExternalDistribution("/nonexistent/path.csv")

    def test_blank_and_comment_lines_skipped(self, distributions,
                                              csv_with_blanks_and_comments):
        dist = distributions.PrimaryExternalDistribution(
            csv_with_blanks_and_comments)
        assert dist.GetPhysicalNumEvents() == 2


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

class TestPrimaryExternalDistributionMetadata:
    def test_name(self, distributions, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic)
        assert dist.Name() == "PrimaryExternalDistribution"

    def test_density_variables(self, distributions, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic)
        assert dist.DensityVariables() == ["External"]

    def test_physical_num_events(self, distributions, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic)
        assert dist.GetPhysicalNumEvents() == 3

    def test_physical_num_events_with_emin(self, distributions, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic, 25.0)
        # Only E=30 passes
        assert dist.GetPhysicalNumEvents() == 1


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

class TestPrimaryExternalDistributionSampling:
    def test_sample_sets_energy(self, distributions, dataclasses, utilities,
                                 csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic)
        rand = utilities.SIREN_random()
        record = dataclasses.PrimaryDistributionRecord(
            dataclasses.ParticleType.NuMu)
        dist.Sample(rand, None, None, record)
        assert record.energy in (10.0, 20.0, 30.0)

    def test_sample_respects_emin(self, distributions, dataclasses, utilities,
                                   csv_basic):
        emin = 15.0
        dist = distributions.PrimaryExternalDistribution(csv_basic, emin)
        rand = utilities.SIREN_random()
        for _ in range(500):
            record = dataclasses.PrimaryDistributionRecord(
                dataclasses.ParticleType.NuMu)
            dist.Sample(rand, None, None, record)
            assert record.energy >= emin


# ---------------------------------------------------------------------------
# GenerationProbability
# ---------------------------------------------------------------------------

class TestPrimaryExternalDistributionProbability:
    def test_above_emin_returns_one(self, distributions, dataclasses,
                                     csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic, 5.0)
        record = dataclasses.InteractionRecord()
        record.primary_momentum = [15.0, 0.0, 0.0, 0.0]
        assert dist.GenerationProbability(None, None, record) == 1.0

    def test_below_emin_returns_zero(self, distributions, dataclasses,
                                      csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic, 5.0)
        record = dataclasses.InteractionRecord()
        record.primary_momentum = [3.0, 0.0, 0.0, 0.0]
        assert dist.GenerationProbability(None, None, record) == 0.0

    def test_at_emin_returns_one(self, distributions, dataclasses, csv_basic):
        dist = distributions.PrimaryExternalDistribution(csv_basic, 5.0)
        record = dataclasses.InteractionRecord()
        record.primary_momentum = [5.0, 0.0, 0.0, 0.0]
        assert dist.GenerationProbability(None, None, record) == 1.0


# ---------------------------------------------------------------------------
# Vertex distributions: construction and metadata
# ---------------------------------------------------------------------------

class TestPrimaryPhysicalVertexDistribution:
    def test_default_construction(self, distributions):
        dist = distributions.PrimaryPhysicalVertexDistribution()
        assert dist is not None

    def test_name(self, distributions):
        dist = distributions.PrimaryPhysicalVertexDistribution()
        assert dist.Name() == "PrimaryPhysicalVertexDistribution"


class TestPrimaryBoundedVertexDistribution:
    def test_default_construction(self, distributions):
        dist = distributions.PrimaryBoundedVertexDistribution()
        assert dist is not None

    def test_construct_with_max_length(self, distributions):
        dist = distributions.PrimaryBoundedVertexDistribution(100.0)
        assert dist.Name() == "PrimaryBoundedVertexDistribution"

    def test_name(self, distributions):
        dist = distributions.PrimaryBoundedVertexDistribution()
        assert dist.Name() == "PrimaryBoundedVertexDistribution"
