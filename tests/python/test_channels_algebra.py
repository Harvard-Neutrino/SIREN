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

    def test_three_way_mix_normalizes(self):
        box = _box_at(50.0)
        mix = (channels.isotropic()
               + channels.isotropic("Gamma")
               + channels.toward("Gamma", box))
        weights = [w for w, _ in mix._entries]
        assert len(weights) == 3
        assert math.isclose(sum(weights), 1.0)

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

    def test_nested_rescale_composes_scale(self):
        """`p * (q * mixture)` composes to p*q, not a flat replacement by p.

        The carried composition scale must multiply through nested `*`, matching
        _MixtureNode's Channel/Tiling behavior; otherwise a double-rescaled
        mixture silently drops its inner factor.
        """
        inner = 0.5 * (channels.isotropic() + channels.isotropic("Gamma"))
        assert math.isclose(inner._scale, 0.5)
        outer = 0.6 * inner
        # The inner 0.5 must survive: 0.6 * 0.5 = 0.3, not 0.6.
        assert math.isclose(outer._scale, 0.3)
        assert math.isclose(sum(w for w, _ in outer._scaled_entries()), 0.3)

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

    def test_toward_fraction_mixes_physical_fallback(self):
        box = _box_at(50.0)
        mix = channels.toward("Gamma", box, fraction=0.3)
        assert isinstance(mix, channels.Mixture)
        weights = sorted(w for w, _ in mix._entries)
        assert math.isclose(sum(weights), 1.0)
        assert math.isclose(weights[0], 0.3, abs_tol=1e-9)
        assert math.isclose(weights[1], 0.7, abs_tol=1e-9)

    def test_toward_fraction_one_is_bare_channel(self):
        box = _box_at(50.0)
        result = channels.toward("Gamma", box, fraction=1.0)
        assert isinstance(result, channels.Channel)


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

    def test_named_particle_resolves_to_correct_index(self):
        sig = _decay2body_signature()  # [NuLight, Gamma]
        box = _box_at(50.0)
        ch = channels.toward("Gamma", box)
        built = ch._build(sig)
        assert isinstance(built, siren.injection.DetectorDirected2BodyChannel)

    def test_integer_daughter_bypasses_resolution(self):
        sig = _decay2body_signature()
        ch = channels.isotropic(1)
        built = ch._build(sig)
        assert isinstance(built, siren.injection.Isotropic2BodyChannel)

    def test_toward_3body_bad_spectator_name_raises(self):
        sig = _decay3body_signature()
        box = _box_at(50.0)
        ch = channels.toward_3body("EMinus", box, spectator="NotAParticle",
                                   strategy="recursive")
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig)

    def test_toward_3body_missing_spectator_raises(self):
        sig = _decay3body_signature()
        box = _box_at(50.0)
        ch = channels.toward_3body("EMinus", box, strategy="recursive")
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig)

    def test_scatter_toward_named_daughter_resolves_index(self, monkeypatch):
        """scatter_toward(daughter=...) resolves the named daughter to its
        secondary index and threads it to the engine channel; the default
        selects the index-0 secondary."""
        sig = _scatter_signature()  # secondary_types = [NuMu, O16Nucleus]
        box = _box_at(50.0)
        captured = {}

        def _spy(target, directed_index, *args, **kwargs):
            captured["idx"] = directed_index
            return "sentinel"

        monkeypatch.setattr(
            siren.injection, "DetectorDirectedScatteringChannel", _spy)

        channels.scatter_toward(box, daughter="O16Nucleus")._build(sig)
        assert captured["idx"] == 1

        captured.clear()
        channels.scatter_toward(box)._build(sig)
        assert captured["idx"] == 0

    def test_scatter_toward_bad_daughter_name_raises(self):
        """An unknown daughter name is a loud ConfigurationError, so the
        parameter genuinely routes through _resolve_index."""
        sig = _scatter_signature()
        box = _box_at(50.0)
        ch = channels.scatter_toward(box, daughter="NotAParticle")
        with pytest.raises(siren.utilities.ConfigurationError):
            ch._build(sig)


# ------------------------------------------------------------------ #
#  PairMass factories map to the right ctor kwargs                     #
# ------------------------------------------------------------------ #

class TestPairMass:

    def test_uniform(self):
        pm = channels.PairMass.uniform()
        assert pm.mass_mode == siren.injection.InvariantMassMode.Uniform

    def test_breit_wigner_fields(self):
        pm = channels.PairMass.breit_wigner(0.017, 0.001)
        assert pm.mass_mode == siren.injection.InvariantMassMode.BreitWigner
        assert pm.resonance_mass == 0.017
        assert pm.resonance_width == 0.001

    def test_power_law_fields(self):
        pm = channels.PairMass.power_law(0.6, offset=0.1)
        assert pm.mass_mode == siren.injection.InvariantMassMode.PowerLaw
        assert pm.power_law_nu == 0.6
        assert pm.power_law_offset == 0.1

    def test_tabulated_fields(self):
        nodes = [0.1, 0.2, 0.3]
        values = [0.0, 0.5, 1.0]
        pm = channels.PairMass.tabulated(nodes, values)
        assert pm.mass_mode == siren.injection.InvariantMassMode.Tabulated
        assert pm.mass_cdf_nodes == nodes
        assert pm.mass_cdf_values == values

    def test_breit_wigner_reaches_3body_ctor_direct(self):
        sig = _decay3body_signature()
        box = _box_at(50.0)
        pm = channels.PairMass.breit_wigner(0.017, 0.001)
        ch = channels.toward_3body("EMinus", box, strategy="direct", pair_mass=pm)
        built = ch._build(sig)
        assert isinstance(built, siren.injection.DetectorDirected3BodyChannel)
        assert built.Topology() == siren.injection.PhaseSpaceTopology.Decay3Body

    def test_power_law_reaches_3body_ctor_recursive(self):
        sig = _decay3body_signature()
        box = _box_at(50.0)
        pm = channels.PairMass.power_law(0.6, offset=0.05)
        ch = channels.toward_3body("EMinus", box, spectator="NuLight",
                                   strategy="recursive", pair_mass=pm)
        built = ch._build(sig)
        assert isinstance(built, siren.injection.DetectorDirected3BodyChannel)
        measure = built.Measure()
        assert measure.type == siren.injection.PhaseSpaceMeasureType.Recursive2Body


# ------------------------------------------------------------------ #
#  tile() disjointness                                                 #
# ------------------------------------------------------------------ #

class TestTiling:

    def test_tile_returns_tiling(self):
        box = _box_at(50.0, half=1.0)
        t = channels.tile(box, 2)
        assert isinstance(t, channels.Tiling)

    def test_tile_compiles_to_nested_mixture(self):
        sig = _decay2body_signature()
        box = _box_at(50.0, half=1.0)
        t = channels.tile(box, 2)
        built = t._build(sig)
        assert isinstance(built, siren.injection.NestedMixtureChannel)
        # 2x2x2 grid cells, all disjointified (disjointify keeps every tile)
        assert len(built.mixture.channels) == 8

    def test_tile_cells_are_disjoint(self):
        # grid_cells() partitions the box into non-overlapping sub-boxes;
        # sample points and require each falls in exactly one cell.
        box = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 0.0)),
            2.0, 2.0, 2.0)
        from siren import directed_tiling
        cells = directed_tiling.grid_cells(box, 2)
        assert len(cells) == 8

        rng = siren.utilities.SIREN_random(1234)
        n_checked = 0
        for _ in range(500):
            p = siren.math.Vector3D(
                rng.Uniform(-1.0, 1.0), rng.Uniform(-1.0, 1.0), rng.Uniform(-1.0, 1.0))
            hits = [c for c in cells if c.IsInside(p)]
            assert len(hits) <= 1, "point {} landed in {} disjoint cells".format(p, len(hits))
            if hits:
                n_checked += 1
        assert n_checked > 0

    def test_tile_participates_in_mixture_algebra(self):
        box = _box_at(50.0, half=1.0)
        t = channels.tile(box, 2)
        mix = 0.4 * channels.isotropic() + 0.6 * t
        weights = [w for w, _ in mix._entries]
        assert math.isclose(sum(weights), 1.0)
        assert any(isinstance(ch, channels.Tiling) for _, ch in mix._entries)

    def test_tile_rejects_unknown_method(self):
        box = _box_at(50.0, half=1.0)
        with pytest.raises(siren.utilities.ConfigurationError):
            channels.tile(box, 2, method="not-a-method")

    def test_grid_cells_rejects_rotated_box(self):
        """A rotated Box cannot be tiled with axis-aligned sub-cells: loud.

        The sub-cells are emitted axis-aligned in world coordinates, so a
        non-identity placement rotation would neither partition nor cover the
        true (rotated) volume.
        """
        from siren import directed_tiling
        # 90 degrees about z: quaternion (0, 0, sin(pi/4), cos(pi/4)).
        q = siren.math.Quaternion(0.0, 0.0, math.sin(math.pi / 4),
                                  math.cos(math.pi / 4))
        rotated = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 50.0), q),
            2.0, 2.0, 2.0)
        with pytest.raises(siren.utilities.ConfigurationError,
                           match="axis-aligned"):
            directed_tiling.grid_cells(rotated, 2)

    def test_grid_cells_rejects_non_box(self):
        """A non-Box region raises the same ConfigurationError fault class as
        the rotated-Box rejection, not a bare TypeError."""
        from siren import directed_tiling
        sphere = siren.geometry.Sphere(
            siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 50.0)),
            1.0, 0.0)
        with pytest.raises(siren.utilities.ConfigurationError, match="Box"):
            directed_tiling.grid_cells(sphere, 2)


# ------------------------------------------------------------------ #
#  Mixture.validate()                                                  #
# ------------------------------------------------------------------ #

class TestValidate:

    def test_topology_mismatch_is_fatal(self):
        """Isotropic2Body (Decay2Body) mixed with a scattering channel
        (Scatter2to2) is a structural topology mismatch: fatal at
        construction time, before any ConvertDensity probe runs.
        """
        box = _box_at(50.0)
        iso = channels.isotropic(0)
        sc = channels.scatter_toward(box)

        # Build engine channels directly to bypass per-signature resolution
        # differences between the two factories' expected signatures --
        # the topology mismatch is a property of the channel types
        # themselves, independent of signature.
        engine_iso = siren.injection.Isotropic2BodyChannel(0)
        engine_sc = siren.injection.DetectorDirectedScatteringChannel(box, 0)
        with pytest.raises(siren.utilities.MeasureCompatibilityError):
            siren.injection.MultiChannelPhaseSpace(
                [engine_iso, engine_sc], [0.5, 0.5])

    def test_validate_catches_topology_mismatch_mixture(self):
        """A Mixture whose Channel factories build engine channels of two
        different topologies fails Mixture.validate() with
        MeasureCompatibilityError, surfaced at configuration time.
        """
        box = _box_at(50.0)

        def iso_factory(signature, *, detector=None, models=None):
            return siren.injection.Isotropic2BodyChannel(0)

        def scatter_factory(signature, *, detector=None, models=None):
            return siren.injection.DetectorDirectedScatteringChannel(box, 0)

        bad_iso = channels.Channel(iso_factory, label="iso")
        bad_sc = channels.Channel(scatter_factory, label="scatter")
        mix = 0.5 * bad_iso + 0.5 * bad_sc

        with pytest.raises(siren.utilities.MeasureCompatibilityError):
            mix.validate(_decay2body_signature())

    def test_validate_passes_for_convertible_mixed_measure_mixture(self):
        """isotropic (SolidAngleRest) + toward (SolidAngleRest): the same
        measure, so validate() must not raise.
        """
        box = _box_at(50.0)
        mix = 0.3 * channels.isotropic("Gamma") + 0.7 * channels.toward("Gamma", box)
        mcps = mix.validate(_decay2body_signature())
        assert len(mcps.channels) == 2

    def test_validate_passes_for_genuinely_different_convertible_measures(self):
        """MandelstamQ2Phi (default scatter_toward) mixed with FixedMassYPhi
        (variable=BjorkenY): genuinely different measures within
        Scatter2to2, related by an implemented ConvertDensity conversion.
        validate() must not raise for this combination.
        """
        box = _box_at(50.0)
        q2_chan = channels.scatter_toward(
            box, variable=siren.injection.ScatteringVariable.Q2)
        y_chan = channels.scatter_toward(
            box, variable=siren.injection.ScatteringVariable.BjorkenY)
        mix = 0.5 * q2_chan + 0.5 * y_chan

        sig = _scatter_signature()
        mcps = mix.compile(sig)
        topology = mcps.CommonTopology()
        common_measure = mcps.CommonMeasure()
        y = 0.3
        m_target = 14.9
        E = 1.0
        record = _scatter_record(bjorken_y=y, E=E, m_target=m_target)
        MType = siren.injection.PhaseSpaceMeasureType

        # Q2 = 2*M*E*y + constant at fixed final-state masses, so
        # |dQ2/dy| = 2*M*E, independent of y. A unit density in y therefore
        # converts to 1/(2*M*E) in Q2, and a unit Q2 density converts to
        # 2*M*E in y.
        abs_jac = 2.0 * m_target * E  # = 29.8 for these kinematics
        checked = 0
        for ch in mcps.channels:
            measure = ch.Measure()
            if measure == common_measure:
                continue
            # This is exactly the conversion validate()'s probe performs.
            result = siren.injection.ConvertDensity(
                1.0, measure, common_measure, topology, record)
            assert math.isfinite(result) and result > 0.0
            if (measure.type == MType.FixedMassYPhi
                    and common_measure.type == MType.MandelstamQ2Phi):
                expected = 1.0 / abs_jac
            elif (measure.type == MType.MandelstamQ2Phi
                    and common_measure.type == MType.FixedMassYPhi):
                expected = abs_jac
            else:
                raise AssertionError(
                    "unexpected measure pair {} -> {}".format(
                        measure.type, common_measure.type))
            assert result == pytest.approx(expected, rel=1e-12)
            checked += 1
        assert checked >= 1

    def test_validate_returns_compiled_mixture(self):
        box = _box_at(50.0)
        mix = 0.3 * channels.isotropic() + 0.7 * channels.toward(0, box)
        mcps = mix.validate(_decay2body_signature())
        assert isinstance(mcps, siren.injection.MultiChannelPhaseSpace)
        assert math.isclose(sum(mcps.weights), 1.0)

    def test_validate_default_template_signature(self):
        """validate() with no signature argument uses the built-in
        decay-shaped template signature."""
        box = _box_at(50.0)
        mix = 0.3 * channels.isotropic("Gamma") + 0.7 * channels.toward("Gamma", box)
        mcps = mix.validate()
        assert len(mcps.channels) == 2


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
