"""Vertex.compile() parity with hand-built processes, and Directed sugar.

All tests are data-free: no detector model, no data files. Anything that
genuinely needs one is skipped with pytest.skip.
"""


import gc


import math


import pickle


import weakref


import pytest


siren = pytest.importorskip("siren")


vertex = pytest.importorskip("siren.vertex")


directed_mod = pytest.importorskip("siren._directed")


channels = pytest.importorskip("siren.channels")


Vertex = vertex.Vertex


Directed = directed_mod.Directed


def _box_at(z, half=1.0):
    return siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, z)),
        2.0 * half, 2.0 * half, 2.0 * half)


def _dummy_xs():
    return siren.interactions.DummyCrossSection()


def _decay2body_signature():
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.N4
    sig.target_type = siren.dataclasses.ParticleType.Decay
    sig.secondary_types = [siren.particles.NuLight, siren.particles.Gamma]
    return sig


class _KeepAliveCrossSection(siren.interactions.CrossSection):
    """Data-free Python-subclassed CrossSection, instrumented for gc checks."""

    alive_count = 0

    def __init__(self):
        siren.interactions.CrossSection.__init__(self)
        _KeepAliveCrossSection.alive_count += 1

    def __del__(self):
        _KeepAliveCrossSection.alive_count -= 1

    def __reduce__(self):
        # This test model has no parameters; reconstructing its type suffices.
        return type(self), ()

    def TotalCrossSection(self, *args, **kwargs):
        return 1.0

    def DifferentialCrossSection(self, *args):
        return 1.0

    def InteractionThreshold(self, *args):
        return 0.0

    def SampleFinalState(self, *args):
        pass

    def GetPossibleTargets(self):
        return [siren.particles.Nucleon]

    def GetPossibleTargetsFromPrimary(self, primary_type):
        return [siren.particles.Nucleon]

    def GetPossiblePrimaries(self):
        return [siren.particles.NuMu]

    def GetPossibleSignatures(self):
        sig = siren.dataclasses.InteractionSignature()
        sig.primary_type = siren.particles.NuMu
        sig.target_type = siren.particles.Nucleon
        sig.secondary_types = [siren.particles.NuMu, siren.particles.Nucleon]
        return [sig]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        return self.GetPossibleSignatures()

    def FinalStateProbability(self, *args):
        return 1.0

    def DensityVariables(self):
        return []

    def equal(self, other):
        return isinstance(other, _KeepAliveCrossSection)


class TestDirected:

    def test_to_mixture_produces_directed_and_fallback_split(self):
        """Directed(...).to_mixture() mixes the directed term and physical() fallback."""
        box = _box_at(50.0)
        d = Directed("Gamma", box, fraction=0.8)
        mix = d.to_mixture()

        assert isinstance(mix, channels.Mixture)
        weights = sorted(w for w, _ in mix._entries)
        assert math.isclose(weights[0], 0.2, abs_tol=1e-9)
        assert math.isclose(weights[1], 0.8, abs_tol=1e-9)

        desc = mix.describe()
        assert "toward" in desc
        assert "physical" in desc

    def test_default_fraction_is_098(self):
        """Directed's own default fraction is 0.98, distinct from Simulation's 0.99 default."""
        box = _box_at(50.0)
        d = Directed("Gamma", box)
        assert d.fraction == 0.98

    def test_fraction_one_has_no_fallback_term(self):
        """fraction=1.0 drops the fallback term from the mixture entirely."""
        box = _box_at(50.0)
        d = Directed("Gamma", box, fraction=1.0)
        mix = d.to_mixture()
        assert len(mix._entries) == 1

    def test_out_of_range_fraction_raises_configuration_error(self):
        """A fraction outside [0, 1] is rejected at construction time."""
        box = _box_at(50.0)
        with pytest.raises(siren.utilities.ConfigurationError):
            Directed("Gamma", box, fraction=1.5)

    def test_toward_channel_builds_against_decay_signature(self):
        """The directed term resolves to a DetectorDirected2BodyChannel for a decay signature."""
        sig = _decay2body_signature()
        box = _box_at(50.0)
        # fraction=1.0 drops the physical() fallback, which needs a model to
        # late-bind against and is not the point of this test.
        d = Directed("Gamma", box, fraction=1.0)
        mix = d.to_mixture(sig)
        engine_channels, weights = mix._flatten(sig)
        kinds = [type(ch).__name__ for ch in engine_channels]
        assert "DetectorDirected2BodyChannel" in kinds

    def test_by_signature_maps_to_different_mixtures(self):
        """Directed.by_signature() dispatches a distinct mixture per signature."""
        sig_a = _decay2body_signature()
        sig_b = siren.dataclasses.InteractionSignature()
        sig_b.primary_type = siren.particles.N4
        sig_b.target_type = siren.dataclasses.ParticleType.Decay
        sig_b.secondary_types = [siren.particles.NuLight, siren.particles.EMinus,
                                  siren.particles.EPlus]

        box = _box_at(50.0)
        iso_mixture = channels.Mixture([(1.0, channels.isotropic(0))])
        by_sig = Directed.by_signature({
            sig_a: Directed("Gamma", box, fraction=0.9),
            sig_b: iso_mixture,
        })

        mix_a = by_sig.to_mixture(sig_a)
        mix_b = by_sig.to_mixture(sig_b)

        assert isinstance(mix_a, channels.Mixture)
        assert isinstance(mix_b, channels.Mixture)
        weights_a = sorted(w for w, _ in mix_a._entries)
        assert math.isclose(weights_a[1], 0.9, abs_tol=1e-9)
        assert mix_b is iso_mixture

    def test_by_signature_unknown_signature_raises_configuration_error(self):
        """Looking up a signature absent from the mapping is a ConfigurationError."""
        sig_a = _decay2body_signature()
        sig_b = siren.dataclasses.InteractionSignature()
        sig_b.primary_type = siren.particles.N4
        sig_b.target_type = siren.dataclasses.ParticleType.Decay
        sig_b.secondary_types = [siren.particles.NuLight, siren.particles.Gamma,
                                  siren.particles.Gamma]

        box = _box_at(50.0)
        by_sig = Directed.by_signature({sig_a: Directed("Gamma", box)})
        with pytest.raises(siren.utilities.ConfigurationError):
            by_sig.to_mixture(sig_b)

    def test_vertex_kinematics_accepts_directed(self):
        """A Vertex's kinematics= accepts a Directed interchangeably with a Mixture."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        # DummyCrossSection's signature is a Scatter2to2 topology, so the
        # directed term must be scatter_toward() (variable= selects it),
        # not the default 2-body decay toward(); fraction=1.0 drops the
        # physical() fallback, which needs a model to late-bind against.
        d = Directed(sig.secondary_types[0], _box_at(2.0), fraction=1.0,
                     variable=siren.injection.ScatteringVariable.Q2)

        v = Vertex(sig.primary_type, xs, distributions=[], kinematics=d)
        proc = v.compile(is_primary=True)
        assert proc.HasPhaseSpace(sig)

