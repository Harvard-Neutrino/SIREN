"""Channel algebra: normalization, name resolution, PairMass, tiling,
and Mixture.validate() config-time checks.

All tests are data-free: no detector model, no data files. Anything that
genuinely needs one is skipped with pytest.skip.
"""


import math


import pytest


siren = pytest.importorskip("siren")


channels = pytest.importorskip("siren.channels")


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
