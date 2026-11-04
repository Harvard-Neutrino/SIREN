"""Channel algebra: normalization, name resolution, PairMass, tiling,
and Mixture.validate() config-time checks.

All tests are data-free: no detector model, no data files. Anything that
genuinely needs one is skipped with pytest.skip.
"""

import math

import pytest

siren = pytest.importorskip("siren")
channels = pytest.importorskip("siren.channels")


# ------------------------------------------------------------------ #
#  Helpers                                                             #
# ------------------------------------------------------------------ #


def _box_at(z, half=1.0):
    return siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, z)),
        2.0 * half, 2.0 * half, 2.0 * half)


def _decay2body_signature():
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.N4
    sig.target_type = siren.dataclasses.ParticleType.Decay
    sig.secondary_types = [siren.particles.NuLight, siren.particles.Gamma]
    return sig


def _decay3body_signature():
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.N4
    sig.target_type = siren.dataclasses.ParticleType.Decay
    sig.secondary_types = [
        siren.particles.NuLight, siren.particles.EMinus, siren.particles.EPlus]
    return sig


def _scatter_signature():
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.NuMu
    sig.target_type = siren.dataclasses.ParticleType.O16Nucleus
    sig.secondary_types = [siren.particles.NuMu, siren.dataclasses.ParticleType.O16Nucleus]
    return sig


def _scatter_record(bjorken_y=0.3, E=1.0, m_target=14.9):
    rec = siren.dataclasses.InteractionRecord()
    rec.signature = _scatter_signature()
    rec.primary_mass = 0.0
    rec.primary_momentum = [E, 0.0, 0.0, E]
    rec.target_mass = m_target
    rec.secondary_masses = [0.0, m_target]
    rec.secondary_momenta = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    rec.secondary_helicities = [0, 0]
    rec.interaction_vertex = [0.0, 0.0, 0.0]
    rec.primary_initial_position = [0.0, 0.0, 0.0]
    rec.interaction_parameters = {"bjorken_y": bjorken_y}
    return rec


# ------------------------------------------------------------------ #
#  Normalization by construction                                       #
# ------------------------------------------------------------------ #

class TestNormalization:

    def test_add_normalizes(self):
        mix = 0.3 * channels.isotropic() + 0.7 * channels.isotropic("Gamma")
        weights = [w for w, _ in mix._entries]
        assert math.isclose(sum(weights), 1.0)
        assert math.isclose(weights[0], 0.3)
        assert math.isclose(weights[1], 0.7)

    def test_unnormalized_add_still_normalizes(self):
        mix = channels.isotropic() * 3.0 + channels.isotropic("Gamma") * 1.0
        weights = [w for w, _ in mix._entries]
        assert math.isclose(sum(weights), 1.0)
        assert math.isclose(weights[0], 0.75)
        assert math.isclose(weights[1], 0.25)


    def test_mixture_plus_channel_renormalizes(self):
        mix = 0.5 * channels.isotropic() + 0.5 * channels.isotropic("Gamma")
        mix2 = mix + channels.isotropic()
        weights = [w for w, _ in mix2._entries]
        assert len(weights) == 3
        assert math.isclose(sum(weights), 1.0)

    def test_channel_plus_mixture_renormalizes(self):
        mix = 0.5 * channels.isotropic() + 0.5 * channels.isotropic("Gamma")
        mix2 = channels.isotropic() + mix
        weights = [w for w, _ in mix2._entries]
        assert len(weights) == 3
        assert math.isclose(sum(weights), 1.0)


    def test_empty_mixture_rejected(self):
        with pytest.raises(siren.utilities.ConfigurationError):
            channels.Mixture([])

    def test_negative_weight_rejected(self):
        with pytest.raises(siren.utilities.ConfigurationError):
            channels.Mixture([(-0.5, channels.isotropic()),
                               (1.5, channels.isotropic("Gamma"))])

    def test_nonfinite_weight_rejected(self):
        with pytest.raises(siren.utilities.ConfigurationError):
            channels.Mixture([(float("nan"), channels.isotropic())])
        with pytest.raises(siren.utilities.ConfigurationError):
            channels.Mixture([(float("inf"), channels.isotropic())])




# ------------------------------------------------------------------ #
#  Per-signature name resolution                                       #
# ------------------------------------------------------------------ #

class TestNameResolution:

    def test_absent_particle_raises_configuration_error(self):
        sig = _decay2body_signature()
        ch = channels.isotropic("NotAKnownParticle")
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig)

    def test_ambiguous_particle_raises_configuration_error(self):
        sig = siren.dataclasses.InteractionSignature()
        sig.primary_type = siren.particles.N4
        sig.target_type = siren.dataclasses.ParticleType.Decay
        sig.secondary_types = [siren.particles.Gamma, siren.particles.Gamma]
        ch = channels.isotropic("Gamma")
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig)


    def test_integer_daughter_bypasses_resolution(self):
        sig = _decay2body_signature()
        ch = channels.isotropic(1)
        built = ch._build(sig)
        assert isinstance(built, siren.injection.Isotropic2BodyChannel)






# ------------------------------------------------------------------ #
#  PairMass factories map to the right ctor kwargs                     #
# ------------------------------------------------------------------ #



# ------------------------------------------------------------------ #
#  tile() disjointness                                                 #
# ------------------------------------------------------------------ #



# ------------------------------------------------------------------ #
#  Mixture.validate()                                                  #
# ------------------------------------------------------------------ #









# ------------------------------------------------------------------ #
#  describe()                                                          #
# ------------------------------------------------------------------ #

class TestDescribe:

    def test_channel_describe(self):
        assert channels.isotropic().describe() == "isotropic"

    def test_mixture_describe_contains_weights_and_labels(self):
        mix = 0.3 * channels.isotropic() + 0.7 * channels.isotropic("Gamma")
        text = mix.describe()
        assert "0.30" in text
        assert "0.70" in text
        assert "isotropic" in text
        assert "+" in text


# ------------------------------------------------------------------ #
#  physical() late binding                                             #
# ------------------------------------------------------------------ #

class TestPhysical:

    def test_physical_without_models_raises(self):
        sig = _scatter_signature()
        ch = channels.physical()
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig)

    def test_physical_with_multiple_models_raises(self):
        sig = _scatter_signature()
        xs = siren.interactions.DummyCrossSection()
        ch = channels.physical()
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig, models=[xs, xs])

    def test_physical_routes_cross_section(self):
        sig = _scatter_signature()
        xs = siren.interactions.DummyCrossSection()
        ch = channels.physical()
        built = ch._build(sig, models=[xs])
        assert isinstance(built, siren.injection.PhysicalCrossSectionChannel)

    def test_physical_compiles_via_mixture(self):
        sig = _scatter_signature()
        xs = siren.interactions.DummyCrossSection()
        mix = channels.Mixture([(1.0, channels.physical())])
        mcps = mix.compile(sig, models=[xs])
        assert isinstance(mcps.channels[0], siren.injection.PhysicalCrossSectionChannel)

    def test_physical_validates_via_mixture(self):
        sig = _scatter_signature()
        xs = siren.interactions.DummyCrossSection()
        mix = channels.Mixture([(1.0, channels.physical())])
        mcps = mix.validate(sig, models=[xs])
        assert isinstance(mcps.channels[0], siren.injection.PhysicalCrossSectionChannel)


class TestCompositionScale:

    def test_scaled_mixture_composition_splits_shares(self):
        """p * mix_a + q * channel gives each sub-channel its share of the total."""
        composed = (0.3 * (0.5 * channels.isotropic(0) + 0.5 * channels.isotropic(1))
                    + 0.7 * channels.isotropic(0))
        weights = [round(w, 6) for w, _ in composed._entries]
        assert weights == [0.15, 0.15, 0.7]

    def test_standalone_mixture_stays_normalized(self):
        """A scaled mixture used on its own is still normalized to sum 1."""
        mix = 0.3 * (0.5 * channels.isotropic(0) + 0.5 * channels.isotropic(1))
        assert math.isclose(sum(w for w, _ in mix._entries), 1.0)
