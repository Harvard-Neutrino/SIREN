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

    def test_z_decay_width_positive(self, interactions, dataclasses):
        primaries = {dataclasses.Particle.ParticleType.Z0}
        decay = interactions.ElectroweakDecay(primaries)
        width = decay.TotalDecayWidth(dataclasses.Particle.ParticleType.Z0)
        assert width > 0

    def test_w_hadronic_to_leptonic_ratio(self, interactions, dataclasses):
        # SM: Gamma(W -> hadrons) / Gamma(W -> leptons) ~ 2 (within ~10%)
        primaries = {dataclasses.Particle.ParticleType.WPlus}
        decay = interactions.ElectroweakDecay(primaries)
        sigs = decay.GetPossibleSignaturesFromParent(dataclasses.Particle.ParticleType.WPlus)
        record = dataclasses.InteractionRecord()
        # wMass ~ 80.4 GeV
        record.primary_mass = 80.4
        lep_types = {dataclasses.Particle.ParticleType.EPlus,
                     dataclasses.Particle.ParticleType.MuPlus,
                     dataclasses.Particle.ParticleType.TauPlus}
        w_lep = 0.0
        w_had = 0.0
        for sig in sigs:
            record.signature = sig
            w = decay.TotalDecayWidthForFinalState(record)
            if sig.secondary_types[0] in lep_types:
                w_lep += w
            else:
                w_had += w
        assert w_lep > 0 and w_had > 0
        ratio = w_had / w_lep
        assert 1.7 < ratio < 2.3
