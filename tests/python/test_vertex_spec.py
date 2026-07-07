"""Vertex.compile() parity with hand-built processes, and Directed sugar.

All tests are data-free: no detector model, no data files. Anything that
genuinely needs one is skipped with pytest.skip.
"""

import gc
import math

import pytest

siren = pytest.importorskip("siren")
vertex = pytest.importorskip("siren.vertex")
directed_mod = pytest.importorskip("siren._directed")
channels = pytest.importorskip("siren.channels")

Vertex = vertex.Vertex
Directed = directed_mod.Directed


# ------------------------------------------------------------------ #
#  Helpers                                                             #
# ------------------------------------------------------------------ #

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


# ------------------------------------------------------------------ #
#  Vertex.compile() parity with a hand-built process                   #
# ------------------------------------------------------------------ #

class TestVertexCompileParity:

    def test_primary_process_fields_match_hand_built(self):
        """Vertex.compile(is_primary=True) matches a hand-built PrimaryInjectionProcess."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        dists = [siren.distributions.PrimaryMass(0)]

        v = Vertex(sig.primary_type, xs, distributions=dists)
        proc = v.compile(is_primary=True)

        assert isinstance(proc, siren.injection.PrimaryInjectionProcess)
        assert proc.primary_type == sig.primary_type
        assert list(proc.interactions.GetCrossSections()) == [xs]
        assert list(proc.distributions) == dists

    def test_secondary_process_fields_match_hand_built(self):
        """Vertex.compile(is_primary=False) matches a hand-built SecondaryInjectionProcess."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        dists = [siren.distributions.SecondaryPhysicalVertexDistribution()]

        v = Vertex(sig.primary_type, xs, distributions=dists)
        proc = v.compile(is_primary=False)

        assert isinstance(proc, siren.injection.SecondaryInjectionProcess)
        assert proc.secondary_type == sig.primary_type
        assert list(proc.interactions.GetCrossSections()) == [xs]
        assert list(proc.distributions) == dists

    def test_kinematics_registers_phase_space_for_each_signature(self):
        """A phase space is registered for every signature the model can produce."""
        xs = _dummy_xs()
        sigs = xs.GetPossibleSignatures()
        v = Vertex(sigs[0].primary_type, xs, distributions=[],
                   kinematics=channels.isotropic(0))

        # A phase space is registered for each signature sharing this Vertex's
        # primary_type (InteractionCollection resolves by primary_type).
        matching = [s for s in sigs if s.primary_type == sigs[0].primary_type]
        proc = v.compile(is_primary=True)
        for sig in matching:
            assert proc.HasPhaseSpace(sig)
            ps = proc.GetPhaseSpace(sig)
            assert isinstance(ps, siren.injection.MultiChannelPhaseSpace)

    def test_no_kinematics_registers_no_phase_space(self):
        """kinematics=None leaves the process with no registered phase space."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        v = Vertex(sig.primary_type, xs, distributions=[])
        proc = v.compile(is_primary=True)
        assert not proc.HasAnyPhaseSpace()

    def test_colliding_signature_models_raise(self):
        """Two models producing the same signature is loud, not a silent drop.

        With kinematics set, each model registers its own phase space under the
        shared signature; the second would silently replace the first even
        though InteractionCollection still selects both.
        """
        xs_a = _dummy_xs()
        xs_b = _dummy_xs()
        # Both DummyCrossSection instances produce the identical signature.
        assert xs_a.GetPossibleSignatures()[0] == xs_b.GetPossibleSignatures()[0]
        v = Vertex(xs_a.GetPossibleSignatures()[0].primary_type,
                   [xs_a, xs_b], distributions=[],
                   kinematics=channels.isotropic(0))
        with pytest.raises(siren.utilities.ConfigurationError,
                           match="produce signature"):
            v.compile(is_primary=True)

    def test_position_appended_to_distributions(self):
        """position= is appended to the distributions list, not silently dropped."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        mass = siren.distributions.PrimaryMass(0)
        pos = siren.distributions.PointSourcePositionDistribution(
            siren.math.Vector3D(0, 0, 0), 25.0)

        v = Vertex(sig.primary_type, xs, distributions=[mass], position=pos)
        assert v.distributions == [mass, pos]
        proc = v.compile(is_primary=True)
        assert list(proc.distributions) == [mass, pos]

    def test_unknown_kwarg_raises_type_error(self):
        """A typo'd keyword argument is rejected, not silently ignored."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        with pytest.raises(TypeError):
            Vertex(sig.primary_type, xs, distirbutions=[])

    def test_unknown_attribute_assignment_raises_attribute_error(self):
        """__slots__ makes an unrecognized field name unrepresentable."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        v = Vertex(sig.primary_type, xs, distributions=[])
        with pytest.raises(AttributeError):
            v.distirbutions = []


# ------------------------------------------------------------------ #
#  Fixed weighting mode reaches the engine                             #
# ------------------------------------------------------------------ #

class TestWeightingMode:

    def test_fixed_weighting_mode_reaches_process(self):
        """weighting=Fixed() compiles to compute_interaction_probability=False."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        fixed = siren.injection.VertexWeightingMode.Fixed()

        v = Vertex(sig.primary_type, xs, distributions=[], weighting=fixed)
        proc = v.compile(is_primary=True)

        assert proc.weighting_mode.compute_interaction_probability == False
        assert proc.weighting_mode.compute_position_probability == False

    def test_default_weighting_mode_is_propagated(self):
        """No weighting= given compiles to the Propagated() preset."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        v = Vertex(sig.primary_type, xs, distributions=[])
        proc = v.compile(is_primary=True)
        propagated = siren.injection.VertexWeightingMode.Propagated()
        assert proc.weighting_mode == propagated


# ------------------------------------------------------------------ #
#  Strong refs survive gc.collect()                                    #
# ------------------------------------------------------------------ #

class TestStrongRefs:

    def test_python_model_survives_gc_after_compile(self):
        """A Python-subclassed model outlives gc.collect() while its Vertex is alive."""
        def build():
            xs = _KeepAliveCrossSection()
            v = Vertex(siren.particles.NuMu, xs, distributions=[])
            proc = v.compile(is_primary=True)
            return v, proc

        before = _KeepAliveCrossSection.alive_count
        v, proc = build()
        gc.collect()

        assert _KeepAliveCrossSection.alive_count == before + 1
        # The compiled process must still be able to reach the model: this
        # is the load-bearing claim (the local `xs` in build() is long out
        # of scope, so only the Vertex's own strong ref keeps it callable).
        recovered = proc.interactions.GetCrossSections()[0]
        assert recovered.GetPossibleSignatures()[0].primary_type == siren.particles.NuMu
        assert v.interactions[0] is recovered

    def test_vertex_retains_models_and_distributions_after_compile(self):
        """The Vertex keeps its own strong refs to models/distributions/kinematics."""
        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        dists = [siren.distributions.PrimaryMass(0)]
        mix = channels.isotropic(0)

        v = Vertex(sig.primary_type, xs, distributions=dists, kinematics=mix)
        v.compile(is_primary=True)
        gc.collect()

        assert v.interactions == [xs]
        assert v.distributions == dists
        assert v.kinematics is mix


# ------------------------------------------------------------------ #
#  as_vertex_spec()                                                    #
# ------------------------------------------------------------------ #

class TestVertexSpec:

    def test_as_vertex_spec_carries_particle_expand_continue_if(self):
        """as_vertex_spec() exposes particle/expand/continue_if for compile_expansion."""
        from siren import expand as expand_mod

        xs = _dummy_xs()
        sig = xs.GetPossibleSignatures()[0]
        rule = expand_mod.child("Nucleon")
        cont = lambda tree, parent, i: True

        v = Vertex(sig.primary_type, xs, distributions=[],
                   expand=(rule,), continue_if=cont)
        spec = v.as_vertex_spec()

        assert spec.particle == sig.primary_type
        assert spec.expand == (rule,)
        assert spec.continue_if is cont
        assert sig.primary_type in spec.secondary_types
        assert siren.particles.Nucleon in spec.secondary_types


# ------------------------------------------------------------------ #
#  Directed sugar over channels.py                                     #
# ------------------------------------------------------------------ #

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
        """Directed's own default fraction is 0.98, not the legacy 0.99 stream default."""
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
