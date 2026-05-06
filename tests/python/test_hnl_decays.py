"""Tests for HNL decay classes."""
import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def interactions():
    from siren import interactions
    return interactions


@pytest.fixture(scope="module")
def dataclasses():
    from siren import dataclasses
    return dataclasses


class TestHNLDecay:
    def test_constructor_vector(self, interactions, dataclasses):
        decay = interactions.HNLDecay(
            0.3, [0.0, 1e-3, 0.0],
            interactions.HNLDecay.ChiralNature.Majorana)
        assert decay is not None

    def test_constructor_scalar(self, interactions):
        decay = interactions.HNLDecay(
            0.3, 1e-3, interactions.HNLDecay.ChiralNature.Dirac)
        assert decay is not None

    def test_get_hnl_mass(self, interactions):
        decay = interactions.HNLDecay(
            0.5, 1e-3, interactions.HNLDecay.ChiralNature.Majorana)
        assert decay.GetHNLMass() == pytest.approx(0.5)

    def test_signatures_non_empty(self, interactions):
        decay = interactions.HNLDecay(
            0.3, 1e-3, interactions.HNLDecay.ChiralNature.Majorana)
        assert len(decay.GetPossibleSignatures()) > 0

    def test_total_decay_width_positive(self, interactions, dataclasses):
        decay = interactions.HNLDecay(
            0.3, 1e-3, interactions.HNLDecay.ChiralNature.Majorana)
        width = decay.TotalDecayWidth(dataclasses.Particle.ParticleType.N4)
        assert width > 0

    def test_heavier_decays_faster(self, interactions, dataclasses):
        light = interactions.HNLDecay(
            0.2, 1e-3, interactions.HNLDecay.ChiralNature.Majorana)
        heavy = interactions.HNLDecay(
            0.5, 1e-3, interactions.HNLDecay.ChiralNature.Majorana)
        N4 = dataclasses.Particle.ParticleType.N4
        assert heavy.TotalDecayWidth(N4) > light.TotalDecayWidth(N4)


class TestHNLDipoleDecay:
    def test_constructor(self, interactions):
        decay = interactions.HNLDipoleDecay(
            0.3, [0.0, 1e-6, 0.0],
            interactions.HNLDipoleDecay.ChiralNature.Dirac)
        assert decay is not None

    def test_signatures_non_empty(self, interactions):
        decay = interactions.HNLDipoleDecay(
            0.3, 1e-6, interactions.HNLDipoleDecay.ChiralNature.Dirac)
        assert len(decay.GetPossibleSignatures()) > 0

    def test_total_decay_width_positive(self, interactions, dataclasses):
        decay = interactions.HNLDipoleDecay(
            0.3, 1e-6, interactions.HNLDipoleDecay.ChiralNature.Dirac)
        width = decay.TotalDecayWidth(dataclasses.Particle.ParticleType.N4)
        assert width > 0


class TestElectroweakDecay:
    def test_constructor(self, interactions, dataclasses):
        primaries = {dataclasses.Particle.ParticleType.WPlus,
                     dataclasses.Particle.ParticleType.WMinus,
                     dataclasses.Particle.ParticleType.Z0}
        decay = interactions.ElectroweakDecay(primaries)
        assert decay is not None

    def test_signatures_non_empty(self, interactions, dataclasses):
        primaries = {dataclasses.Particle.ParticleType.WPlus}
        decay = interactions.ElectroweakDecay(primaries)
        assert len(decay.GetPossibleSignatures()) > 0

    def test_w_decay_width_positive(self, interactions, dataclasses):
        primaries = {dataclasses.Particle.ParticleType.WPlus}
        decay = interactions.ElectroweakDecay(primaries)
        width = decay.TotalDecayWidth(dataclasses.Particle.ParticleType.WPlus)
        assert width > 0
