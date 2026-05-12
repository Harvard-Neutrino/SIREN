"""Tests for HNL DIS cross section classes and ParticleMasses."""
import math

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


@pytest.fixture(scope="module")
def utilities():
    from siren import utilities
    return utilities


# -- ParticleMasses ----------------------------------------------------------

class TestParticleTypes:
    def test_neutrino_types_exist(self, dataclasses):
        from siren.dataclasses import Particle
        assert hasattr(Particle.ParticleType, "NuE")
        assert hasattr(Particle.ParticleType, "NuMu")
        assert hasattr(Particle.ParticleType, "NuTau")

    def test_hnl_types_exist(self, dataclasses):
        from siren.dataclasses import Particle
        assert hasattr(Particle.ParticleType, "N4")

    def test_meson_types_exist(self, dataclasses):
        from siren.dataclasses import Particle
        for name in ["Pi0", "PiPlus", "KPlus", "K0_Long", "Eta"]:
            assert hasattr(Particle.ParticleType, name)


class TestParticleMasses:
    def test_proton_mass_positive(self, dataclasses):
        from siren.dataclasses import Particle
        assert dataclasses.GetParticleMass(Particle.ParticleType.PPlus) > 0

    def test_neutrinos_massless(self, dataclasses):
        from siren.dataclasses import Particle
        for nu in [Particle.ParticleType.NuE, Particle.ParticleType.NuMu,
                   Particle.ParticleType.NuTau]:
            assert dataclasses.GetParticleMass(nu) == 0.0

    def test_antiparticle_symmetry(self, dataclasses):
        from siren.dataclasses import Particle
        assert (dataclasses.GetParticleMass(Particle.ParticleType.MuMinus) ==
                dataclasses.GetParticleMass(Particle.ParticleType.MuPlus))

    def test_electron_lighter_than_muon(self, dataclasses):
        from siren.dataclasses import Particle
        e_mass = dataclasses.GetParticleMass(Particle.ParticleType.EMinus)
        mu_mass = dataclasses.GetParticleMass(Particle.ParticleType.MuMinus)
        assert 0 < e_mass < mu_mass


# -- HNLDISFromSpline -------------------------------------------------------

class TestHNLDISFromSpline:
    def test_default_constructor(self, interactions):
        xs = interactions.HNLDISFromSpline()
        assert xs is not None

    def test_possible_primaries_empty_on_default(self, interactions):
        xs = interactions.HNLDISFromSpline()
        assert len(xs.GetPossiblePrimaries()) == 0

    def test_possible_targets_empty_on_default(self, interactions):
        xs = interactions.HNLDISFromSpline()
        assert len(xs.GetPossibleTargets()) == 0


# -- HNLDipoleDISFromSpline -------------------------------------------------

class TestHNLDipoleDISFromSpline:
    def test_default_constructor(self, interactions):
        xs = interactions.HNLDipoleDISFromSpline()
        assert xs is not None

    def test_possible_primaries_empty_on_default(self, interactions):
        xs = interactions.HNLDipoleDISFromSpline()
        assert len(xs.GetPossiblePrimaries()) == 0

    def test_possible_targets_empty_on_default(self, interactions):
        xs = interactions.HNLDipoleDISFromSpline()
        assert len(xs.GetPossibleTargets()) == 0


# -- HNLDipoleFromTable ------------------------------------------------------

class TestHNLDipoleFromTable:
    def test_class_exists(self, interactions):
        assert hasattr(interactions, "HNLDipoleFromTable")

    def test_helicity_channels_exist(self, interactions):
        assert hasattr(interactions.HNLDipoleFromTable, "HelicityChannel")
        assert hasattr(interactions.HNLDipoleFromTable.HelicityChannel,
                       "Flipping")
        assert hasattr(interactions.HNLDipoleFromTable.HelicityChannel,
                       "Conserving")
