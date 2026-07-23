import pytest


class TestParticlesNamespace:
    """siren.particles should expose curated particle type aliases."""

    def test_neutrinos_exist(self):
        from siren.particles import NuE, NuMu, NuTau, NuEBar, NuMuBar, NuTauBar
        assert NuMu is not None

    def test_charged_leptons_exist(self):
        from siren.particles import Electron, Muon, EMinus, MuMinus
        assert Electron == EMinus

    def test_nucleons_exist(self):
        from siren.particles import Nucleon, Neutron, PPlus
        assert Nucleon is not None

    def test_bsm_types_available(self):
        from siren import particles
        assert hasattr(particles, "N4")

    def test_all_enum_members_auto_exported(self):
        """Every member of the C++ ParticleType enum should be a
        module-level attribute of siren.particles."""
        import siren.dataclasses as dc
        import siren.particles as particles
        PT = dc.ParticleType
        sentinel = PT.NuMu
        for name in dir(PT):
            if name.startswith("_"):
                continue
            val = getattr(PT, name)
            if isinstance(val, type(sentinel)):
                assert hasattr(particles, name), (
                    f"{name} is in ParticleType enum but not in siren.particles"
                )

    def test_resolve_string(self):
        from siren.particles import resolve, NuMu
        assert resolve("NuMu") == NuMu

    def test_resolve_enum(self):
        from siren.particles import resolve, NuMu
        assert resolve(NuMu) == NuMu

    def test_resolve_alias(self):
        from siren.particles import resolve, Electron
        assert resolve("Electron") == Electron

    def test_resolve_unknown_raises(self):
        from siren.particles import resolve
        with pytest.raises(ValueError, match="Unknown particle type"):
            resolve("FooBarQuux")

    def test_resolve_wrong_type_raises(self):
        from siren.particles import resolve
        with pytest.raises(TypeError):
            resolve(42)


class TestDistNamespace:
    """siren.dist should expose distribution shorthands."""

    def test_energy_types(self):
        from siren.dist import PowerLaw, Monoenergetic
        assert PowerLaw is not None
        assert Monoenergetic is not None

    def test_direction_types(self):
        from siren.dist import IsotropicDirection, FixedDirection, Cone
        assert IsotropicDirection is not None

    def test_position_types(self):
        from siren.dist import ColumnDepth, CylinderVolume, PointSource
        assert ColumnDepth is not None

    def test_fixed_direction_accepts_list(self):
        from siren.dist import FixedDirection
        d = FixedDirection([0, 0, 1])
        assert d is not None

    def test_fixed_direction_accepts_tuple(self):
        from siren.dist import FixedDirection
        d = FixedDirection((1, 0, 0))
        assert d is not None

    def test_all_concrete_distributions_auto_exported(self):
        """Every concrete class in siren.distributions should appear in siren.dist."""
        import siren.distributions as distributions
        import siren.dist as dist
        for name in dir(distributions):
            obj = getattr(distributions, name)
            if isinstance(obj, type) and not name.startswith("_"):
                assert hasattr(dist, name), (
                    f"{name} is in siren.distributions but not in siren.dist"
                )


class TestDistributionVariableEnum:
    """DistributionVariable enum and SetVariables/RequiredVariables."""

    def test_enum_exists(self):
        import siren.distributions as d
        DV = d.DistributionVariable
        assert hasattr(DV, "PrimaryEnergy")
        assert hasattr(DV, "PrimaryDirection")
        assert hasattr(DV, "InteractionVertex")

    def test_set_variables_energy(self):
        import siren.distributions as d
        DV = d.DistributionVariable
        pl = d.PowerLaw(2, 1e3, 1e6)
        assert DV.PrimaryEnergy in pl.SetVariables()

    def test_set_variables_direction(self):
        import siren.distributions as d
        DV = d.DistributionVariable
        iso = d.IsotropicDirection()
        assert DV.PrimaryDirection in iso.SetVariables()

    def test_set_variables_position(self):
        import siren.distributions as d
        DV = d.DistributionVariable
        col = d.ColumnDepthPositionDistribution(
            600, 600.0, d.LeptonDepthFunction()
        )
        sv = col.SetVariables()
        assert DV.InteractionVertex in sv

    def test_required_variables_column_depth(self):
        """ColumnDepthPositionDistribution requires direction, energy, mass."""
        import siren.distributions as d
        DV = d.DistributionVariable
        col = d.ColumnDepthPositionDistribution(
            600, 600.0, d.LeptonDepthFunction()
        )
        rv = col.RequiredVariables()
        assert DV.PrimaryDirection in rv
        assert DV.PrimaryEnergy in rv
        assert DV.PrimaryMass in rv

    def test_required_variables_energy_empty(self):
        """Energy distributions have no required variables."""
        import siren.distributions as d
        pl = d.PowerLaw(2, 1e3, 1e6)
        assert len(pl.RequiredVariables()) == 0

    def test_delta_function_detection(self):
        """Delta functions: SetVariables has the var, DensityVariables does not."""
        import siren.distributions as d
        DV = d.DistributionVariable

        fd = d.FixedDirection([0, 0, 1])
        assert DV.PrimaryDirection in fd.SetVariables()
        assert len(fd.DensityVariables()) == 0

        mono = d.Monoenergetic(1000)
        assert DV.PrimaryEnergy in mono.SetVariables()
        assert len(mono.DensityVariables()) == 0

    def test_non_delta_has_density(self):
        """Non-delta distributions appear in both SetVariables and DensityVariables."""
        import siren.distributions as d
        DV = d.DistributionVariable

        pl = d.PowerLaw(2, 1e3, 1e6)
        assert DV.PrimaryEnergy in pl.SetVariables()
        assert "PrimaryEnergy" in pl.DensityVariables()

    def test_pidar_inherits_set_variables(self):
        """PiDARNuEDistribution inherits SetVariables from base without
        manual registration."""
        import siren.distributions as d
        DV = d.DistributionVariable
        pidar = d.PiDARNuEDistribution()
        assert DV.PrimaryEnergy in pidar.SetVariables()


class TestOrderingValidation:
    """validate_ordering should catch dependency violations."""

    def test_correct_ordering_passes(self):
        import siren.distributions as d
        from siren._validation import validate_ordering
        dists = [
            d.PrimaryMass(0),
            d.PowerLaw(2, 1e3, 1e6),
            d.IsotropicDirection(),
            d.ColumnDepthPositionDistribution(
                600, 600.0, d.LeptonDepthFunction()
            ),
        ]
        validate_ordering(dists)

    def test_wrong_ordering_raises(self):
        import siren.distributions as d
        from siren._validation import validate_ordering
        dists = [
            d.ColumnDepthPositionDistribution(
                600, 600.0, d.LeptonDepthFunction()
            ),
            d.PrimaryMass(0),
            d.PowerLaw(2, 1e3, 1e6),
            d.IsotropicDirection(),
        ]
        import pytest
        with pytest.raises(ValueError, match="requires variables"):
            validate_ordering(dists)

    def test_duplicate_set_variables_raises(self):
        import siren.distributions as d
        from siren._validation import validate_ordering
        dists = [
            d.PrimaryMass(0),
            d.PowerLaw(2, 1e3, 1e6),
            d.IsotropicDirection(),
            d.FixedDirection([0, 0, 1]),  # duplicate PrimaryDirection
            d.ColumnDepthPositionDistribution(
                600, 600.0, d.LeptonDepthFunction()
            ),
        ]
        import pytest
        with pytest.raises(ValueError, match="already been set by"):
            validate_ordering(dists)

    def test_bounded_vertex_requires_initial_position(self):
        """PrimaryBoundedVertexDistribution consumes the primary's initial
        position (record.GetInitialPosition() in SamplePosition), so it must
        declare InitialPosition required and must not declare it set; a chain
        that never sets an initial position upstream fails ordering validation.

        Regression guard: the base VertexPositionDistribution sets only
        InteractionVertex, so RequiredVariables = {InitialPosition,
        PrimaryDirection} does not overlap SetVariables and ordering is
        satisfiable given an upstream position source. Dropping InitialPosition
        from RequiredVariables to "resolve" a non-existent overlap silently
        removes a real dependency, which this test forbids.
        """
        import siren
        import siren.distributions as d
        from siren._validation import validate_ordering
        DV = d.DistributionVariable
        bounded = d.PrimaryBoundedVertexDistribution(
            siren.geometry.Box(widths=(2.0, 2.0, 2.0)))
        assert DV.InitialPosition in bounded.RequiredVariables()
        assert DV.PrimaryDirection in bounded.RequiredVariables()
        assert DV.InteractionVertex in bounded.SetVariables()
        assert DV.InitialPosition not in bounded.SetVariables()
        # Mass, energy, and direction are set, but nothing sets InitialPosition,
        # so the bounded vertex's requirement is unmet and ordering must fail.
        dists = [
            d.PrimaryMass(0),
            d.PowerLaw(2, 1e3, 1e6),
            d.IsotropicDirection(),
            bounded,
        ]
        import pytest
        with pytest.raises(ValueError, match="InitialPosition"):
            validate_ordering(dists)


class TestReweightingValidation:
    """Physical and injection measures must use compatible variables."""

    @staticmethod
    def _fixed_target_distributions():
        import siren.distributions as d
        import siren.geometry as g

        target = g.Cylinder(1.0, 0.0, 2.0)
        return (
            d.FixedTargetPositionDistribution(target, 10.0),
            d.FixedTargetAreaDistribution(target),
        )

    def test_propagated_vertex_leaves_external_area_factor(self):
        """The 1D normalized-position density is covered by a 3D vertex."""
        from siren._validation import validate_reweighting_compatibility

        vertex, _ = self._fixed_target_distributions()
        validate_reweighting_compatibility([vertex], [])

    def test_fixed_target_area_and_position_density_cancel_vertex(self):
        """Two area coordinates plus the 1D position factor cover the vertex."""
        from siren._validation import validate_reweighting_compatibility

        vertex, area = self._fixed_target_distributions()
        validate_reweighting_compatibility([vertex], [area])

    def test_downstream_mode_controls_normalized_position_measure(self):
        from siren._validation import validate_reweighting_compatibility

        _, area = self._fixed_target_distributions()
        with pytest.raises(ValueError, match="PrimaryPositionLongitudinal"):
            validate_reweighting_compatibility([area], [])

        validate_reweighting_compatibility(
            [area], [], compute_position_probability=False)

    def test_interaction_probability_is_dimensionless(self):
        from siren._validation import validate_reweighting_compatibility

        vertex, area = self._fixed_target_distributions()
        validate_reweighting_compatibility(
            [vertex], [area], compute_interaction_probability=False)

    def test_simulation_builds_weighter_only_after_compatibility_check(
        self, monkeypatch,
    ):
        import importlib

        simulation_module = importlib.import_module("siren.Simulation")
        injection_distributions = [object()]
        physical_distributions = [object()]

        def reject(injection, physical, **_mode):
            assert injection == injection_distributions
            assert physical == physical_distributions
            raise ValueError("compatibility sentinel")

        monkeypatch.setattr(
            simulation_module, "validate_physical_distributions",
            lambda _distributions: None)
        monkeypatch.setattr(
            simulation_module, "validate_reweighting_compatibility", reject)

        simulation = simulation_module.Simulation.__new__(
            simulation_module.Simulation)
        simulation._injection_distributions = injection_distributions
        simulation._physical_distributions = physical_distributions

        with pytest.raises(ValueError, match="compatibility sentinel"):
            simulation._build_weighter(object())


class TestTopLevelExports:
    """Top-level siren.* should expose key functions and classes."""

    def test_simulation_class(self):
        import siren
        assert hasattr(siren, "Simulation")

    def test_results_class(self):
        import siren
        assert hasattr(siren, "Results")

    def test_generate_events(self):
        import siren
        assert callable(siren.GenerateEvents)

    def test_save_events(self):
        import siren
        assert callable(siren.SaveEvents)

    def test_load_events(self):
        import siren
        assert callable(siren.LoadEvents)

    def test_load_detector(self):
        import siren
        assert callable(siren.load_detector)

    def test_load_flux(self):
        import siren
        assert callable(siren.load_flux)

    def test_load_processes(self):
        import siren
        assert callable(siren.load_processes)


class TestSimulationValidation:
    """Simulation should validate inputs at construction time."""

    def test_missing_energy_raises(self):
        import siren
        with pytest.raises(ValueError, match="energy"):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                targets="Nucleon",
                process="CC",
                direction=siren.dist.IsotropicDirection(),
                position=siren.dist.ColumnDepth(
                    600, 600.0, siren.distributions.LeptonDepthFunction()
                ),
            )

    def test_missing_direction_raises(self):
        import siren
        with pytest.raises(ValueError, match="direction"):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                targets="Nucleon",
                process="CC",
                energy=siren.dist.PowerLaw(2, 1e3, 1e6),
                position=siren.dist.ColumnDepth(
                    600, 600.0, siren.distributions.LeptonDepthFunction()
                ),
            )

    def test_missing_position_raises(self):
        import siren
        with pytest.raises(ValueError, match="position"):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                targets="Nucleon",
                process="CC",
                energy=siren.dist.PowerLaw(2, 1e3, 1e6),
                direction=siren.dist.IsotropicDirection(),
            )

    def test_conflicting_energy_raises(self):
        import siren
        with pytest.raises(ValueError, match="Cannot specify both"):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                targets="Nucleon",
                process="CC",
                energy=siren.dist.PowerLaw(2, 1e3, 1e6),
                injection_energy=siren.dist.Monoenergetic(1000),
                direction=siren.dist.IsotropicDirection(),
                position=siren.dist.ColumnDepth(
                    600, 600.0, siren.distributions.LeptonDepthFunction()
                ),
            )

    def test_conflicting_flux_and_physical_energy_raises(self):
        import siren
        with pytest.raises(ValueError, match="Cannot specify both.*flux"):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                targets="Nucleon",
                process="CC",
                flux=siren.dist.PowerLaw(2, 1e3, 1e6),
                physical_energy=siren.dist.Monoenergetic(1000),
                direction=siren.dist.IsotropicDirection(),
                position=siren.dist.ColumnDepth(
                    600, 600.0, siren.distributions.LeptonDepthFunction()
                ),
            )

    def test_no_interactions_raises(self):
        import siren
        with pytest.raises(TypeError):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="NuMu",
                energy=siren.dist.PowerLaw(2, 1e3, 1e6),
                direction=siren.dist.IsotropicDirection(),
                position=siren.dist.ColumnDepth(
                    600, 600.0, siren.distributions.LeptonDepthFunction()
                ),
            )

    def test_bad_particle_name_raises(self):
        import siren
        with pytest.raises(ValueError, match="Unknown particle type"):
            siren.Simulation(
                events=1,
                detector="IceCube",
                primary="InvalidParticle",
                interactions="CSMSDISSplines",
                energy=siren.dist.PowerLaw(2, 1e3, 1e6),
                direction=siren.dist.IsotropicDirection(),
                position=siren.dist.ColumnDepth(
                    600, 600.0, siren.distributions.LeptonDepthFunction()
                ),
            )


class TestSimulationConstruction:
    """Simulation should resolve all inputs correctly."""

    def test_string_detector(self):
        import siren
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        assert sim.detector_model is not None

    def test_string_primary(self):
        import siren
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        assert sim.primary_type == siren.particles.NuMu

    def test_enum_primary(self):
        import siren
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary=siren.particles.NuMu,
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        assert sim.primary_type == siren.particles.NuMu

    def test_shared_energy_goes_to_both(self):
        import siren
        e = siren.dist.PowerLaw(2, 1e3, 1e6)
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=e,
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        # Energy should appear in both injection and physical
        assert e in sim.injection_distributions
        assert e in sim.physical_distributions

    def test_split_direction(self):
        import siren
        inj_d = siren.dist.FixedDirection([0, 0, 1])
        phys_d = siren.dist.IsotropicDirection()
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            injection_direction=inj_d,
            physical_direction=phys_d,
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        assert inj_d in sim.injection_distributions
        assert phys_d in sim.physical_distributions

    def test_mass_auto_inferred(self):
        import siren
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        # Should have PrimaryMass in injection but NOT physical
        mass_dists = [
            d for d in sim.injection_distributions
            if isinstance(d, siren.distributions.PrimaryMass)
        ]
        assert len(mass_dists) == 1

        phys_mass = [
            d for d in sim.physical_distributions
            if isinstance(d, siren.distributions.PrimaryMass)
        ]
        assert len(phys_mass) == 0

    def test_repr(self):
        import siren
        sim = siren.Simulation(
            events=42,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        r = repr(sim)
        assert "42" in r
        assert "NuMu" in r


class TestSimulationRun:
    """Simulation.run() should produce valid Results."""

    @pytest.fixture
    def sim(self):
        import siren
        return siren.Simulation(
            events=5,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )

    def test_run_returns_results(self, sim):
        from siren.Results import Results
        results = sim.run()
        assert isinstance(results, Results)

    def test_correct_event_count(self, sim):
        results = sim.run()
        assert len(results) == 5

    def test_weights_are_positive(self, sim):
        results = sim.run()
        assert all(w > 0 for w in results.weights)

    def test_iteration(self, sim):
        results = sim.run()
        pairs = list(results)
        assert len(pairs) == 5
        event, weight = pairs[0]
        assert weight > 0

    def test_slicing(self, sim):
        results = sim.run()
        subset = results[:2]
        assert len(subset) == 2

    def test_indexing(self, sim):
        results = sim.run()
        event, weight = results[0]
        assert weight > 0


class TestInjectorIterator:
    """Injector should support iteration."""

    def test_iter(self):
        import siren
        injector = siren.injection.Injector()
        injector.number_of_events = 3
        injector.detector_model = siren.load_detector("IceCube")
        injector.primary_type = siren.particles.NuMu
        xs, _ = siren.load_processes(
            "CSMSDISSplines",
            primary_types=[siren.particles.NuMu],
            target_types=[siren.particles.Nucleon],
            isoscalar=True,
            process_types=["CC"],
        )
        injector.primary_interactions = xs[siren.particles.NuMu]
        injector.primary_injection_distributions = [
            siren.distributions.PrimaryMass(0),
            siren.distributions.PowerLaw(2, 1e3, 1e6),
            siren.distributions.IsotropicDirection(),
            siren.distributions.ColumnDepthPositionDistribution(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        ]
        events = list(injector)
        assert len(events) == 3

    def test_len(self):
        import siren
        injector = siren.injection.Injector()
        injector.number_of_events = 7
        assert len(injector) == 7


class TestSecondaryBiasing:
    """Phase space channel configuration for secondary biasing."""

    def test_isotropic_channel_construction(self):
        import siren
        ch = siren.injection.Isotropic2BodyChannel(0)
        assert ch.Name() == "Isotropic2Body"
        assert ch.Measure().type == siren.injection.PhaseSpaceMeasureType.SolidAngleRest

    def test_isotropic_channel_topology_measure(self):
        import siren
        ch = siren.injection.Isotropic2BodyChannel(0)
        assert ch.Topology() == siren.injection.PhaseSpaceTopology.Decay2Body
        assert ch.Measure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()

    def test_detector_directed_channel_construction(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        ch = siren.injection.DetectorDirected2BodyChannel(fid, 0)
        assert ch.Name() == "DetectorDirected2Body"
        assert ch.Measure().type == siren.injection.PhaseSpaceMeasureType.SolidAngleRest

    def test_detector_directed_channel_topology_measure(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        ch = siren.injection.DetectorDirected2BodyChannel(fid, 0)
        assert ch.Topology() == siren.injection.PhaseSpaceTopology.Decay2Body
        assert ch.Measure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()

    def test_detector_directed_3body_channel_construction(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        ch = siren.injection.DetectorDirected3BodyChannel(
            factorization=siren.injection.ThreeBodyMode.Recursive,
            target=fid,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_index=1,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
        )
        assert ch.Name() == "DetectorDirected3Body"
        assert ch.Measure().type == siren.injection.PhaseSpaceMeasureType.Recursive2Body

    def test_3body_channel_topology_measure(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        ch = siren.injection.DetectorDirected3BodyChannel(
            factorization=siren.injection.ThreeBodyMode.Direct,
            target=fid, directed_index=2)
        assert ch.Topology() == siren.injection.PhaseSpaceTopology.Decay3Body
        assert ch.Measure() == siren.injection.PhaseSpaceMeasure.Recursive2Body(2, 0, 1)

    def test_detector_directed_scattering_channel_construction(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        ch = siren.injection.DetectorDirectedScatteringChannel(
            fid,
            directed_index=0,
            variable=siren.injection.ScatteringVariable.Q2,
        )
        assert ch.Name() == "DetectorDirectedScattering"
        assert ch.Measure().type == siren.injection.PhaseSpaceMeasureType.MandelstamQ2Phi

    def test_scattering_channel_topology_measure(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        ch_q2 = siren.injection.DetectorDirectedScatteringChannel(
            fid, variable=siren.injection.ScatteringVariable.Q2)
        ch_by = siren.injection.DetectorDirectedScatteringChannel(
            fid, variable=siren.injection.ScatteringVariable.BjorkenY)
        assert ch_q2.Topology() == siren.injection.PhaseSpaceTopology.Scatter2to2
        assert ch_q2.Measure() == siren.injection.PhaseSpaceMeasure.MandelstamQ2Phi()
        assert ch_by.Topology() == siren.injection.PhaseSpaceTopology.Scatter2to2
        assert ch_by.Measure() == siren.injection.PhaseSpaceMeasure.FixedMassYPhi()

    def test_multi_channel_construction(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [
            siren.injection.Isotropic2BodyChannel(0),
            siren.injection.DetectorDirected2BodyChannel(fid, 0),
        ]
        mc.weights = [0.01, 0.99]
        assert len(mc.channels) == 2
        assert mc.CommonTopology() == siren.injection.PhaseSpaceTopology.Decay2Body
        assert mc.CommonMeasure() == siren.injection.PhaseSpaceMeasure.SolidAngleRest()
        assert mc.CommonMeasure().type == siren.injection.PhaseSpaceMeasureType.SolidAngleRest

    def test_topology_mismatch_raises(self):
        """Mixing Decay2Body with Scatter2to2 should throw."""
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [
            siren.injection.Isotropic2BodyChannel(0),
            siren.injection.DetectorDirectedScatteringChannel(fid),
        ]
        mc.weights = [0.5, 0.5]
        with pytest.raises(RuntimeError, match="[Tt]opology"):
            mc.CommonTopology()

    def test_secondary_process_phase_space(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        sip = siren.injection.SecondaryInjectionProcess()
        assert not sip.HasAnyPhaseSpace()

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [
            siren.injection.Isotropic2BodyChannel(0),
            siren.injection.DetectorDirected2BodyChannel(fid, 0),
        ]
        mc.weights = [0.01, 0.99]

        # SetPhaseSpace requires a signature
        sig = siren.dataclasses.InteractionSignature()
        sig.primary_type = siren.particles.NuMu
        sig.target_type = siren.particles.Nucleon
        sig.secondary_types = [siren.particles.MuMinus, siren.particles.PPlus]
        sip.SetPhaseSpace(sig, mc)
        assert sip.HasPhaseSpace(sig)
        assert sip.HasAnyPhaseSpace()

    def test_injector_primary_phase_space(self):
        import siren
        fid = siren.get_fiducial_volume("IceCube")
        injector = siren.injection.Injector()

        sig = siren.dataclasses.InteractionSignature()
        sig.primary_type = siren.particles.NuMu
        sig.target_type = siren.particles.Nucleon
        sig.secondary_types = [siren.particles.MuMinus, siren.particles.PPlus]

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [siren.injection.DetectorDirected2BodyChannel(fid, 0)]
        mc.weights = [1.0]

        injector.primary_phase_spaces = {sig: mc}
        assert injector.primary_phase_spaces[sig] is mc

    def test_kinematics_utilities(self):
        import siren
        assert abs(siren.injection.Kallen(1.0, 0.0, 0.0) - 1.0) < 1e-14
        p = siren.injection.TwoBodyRestMomentum(0.3, 0.106, 0.140)
        assert p > 0


class TestWeighterBatch:
    """Weighter.weight_all should batch-weight events."""

    def test_weight_all(self):
        import siren
        sim = siren.Simulation(
            events=3,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        results = sim.run()
        # weight_all should give same results as individual weighting
        batch_weights = results._weighter.weight_all(results.events)
        assert len(batch_weights) == 3
        for bw, rw in zip(batch_weights, results.weights):
            assert abs(bw - rw) < 1e-12


class TestDistributionsList:
    """Simulation.__init__ should accept a list of physical and injection distributions."""

    def test_list_of_distributions(self):
        import siren
        injection_dists = [
            siren.dist.PrimaryMass(0),
            siren.dist.PowerLaw(2, 1e3, 1e6),
            siren.dist.IsotropicDirection(),
            siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        ]
        physical_dists = [
            siren.dist.PrimaryMass(0),
            siren.dist.PowerLaw(2, 1e3, 1e6),
            siren.dist.IsotropicDirection(),
        ]
        sim = siren.Simulation(
            events=5,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            injection_distributions=injection_dists,
            physical_distributions=physical_dists,
        )
        assert len(sim.injection_distributions) == len(injection_dists)
        assert len(sim.physical_distributions) == len(physical_dists)


class TestReweight:
    """Simulation.reweight() should produce new weights without
    re-generating events."""

    @pytest.fixture
    def sim_and_results(self):
        import siren
        sim = siren.Simulation(
            events=5,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        results = sim.run()
        return sim, results

    def test_reweight_before_run_raises(self):
        import siren
        sim = siren.Simulation(
            events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )
        with pytest.raises(RuntimeError, match="Must call run"):
            sim.reweight(physical_energy=siren.dist.PowerLaw(1, 1e3, 1e6))

    def test_reweight_returns_results(self, sim_and_results):
        import siren
        from siren.Results import Results
        sim, original = sim_and_results
        reweighted = sim.reweight(
            physical_energy=siren.dist.PowerLaw(1, 1e3, 1e6),
        )
        assert isinstance(reweighted, Results)

    def test_reweight_same_event_count(self, sim_and_results):
        import siren
        sim, original = sim_and_results
        reweighted = sim.reweight(
            physical_energy=siren.dist.PowerLaw(1, 1e3, 1e6),
        )
        assert len(reweighted) == len(original)

    def test_reweight_same_events(self, sim_and_results):
        """Events should be identical objects (not copies)."""
        import siren
        sim, original = sim_and_results
        reweighted = sim.reweight(
            physical_energy=siren.dist.PowerLaw(1, 1e3, 1e6),
        )
        for orig_ev, new_ev in zip(original.events, reweighted.events):
            assert orig_ev is new_ev

    def test_reweight_different_weights(self, sim_and_results):
        """Different physical energy should produce different weights."""
        import siren
        sim, original = sim_and_results
        reweighted = sim.reweight(
            physical_energy=siren.dist.PowerLaw(1, 1e3, 1e6),
        )
        # With a different power law index, weights should differ
        any_different = any(
            abs(a - b) > 1e-15
            for a, b in zip(original.weights, reweighted.weights)
        )
        assert any_different, "All weights identical despite different energy spectrum"

    def test_reweight_no_changes_same_weights(self, sim_and_results):
        """Reweighting with no changes should reproduce original weights."""
        sim, original = sim_and_results
        reweighted = sim.reweight()
        for orig_w, new_w in zip(original.weights, reweighted.weights):
            assert abs(orig_w - new_w) < 1e-12

    def test_reweight_direction(self, sim_and_results):
        """Reweight to a cone -- must have overlapping support with
        the injection direction (isotropic) for nonzero weights."""
        import siren
        sim, original = sim_and_results
        reweighted = sim.reweight(
            physical_direction=siren.dist.Cone([0, 0, 1], 3.14159),
        )
        assert len(reweighted) == len(original)
        assert all(w > 0 for w in reweighted.weights)

    def test_reweight_multiple_calls(self, sim_and_results):
        """Multiple reweight calls should be independent."""
        import siren
        sim, original = sim_and_results
        r1 = sim.reweight(physical_energy=siren.dist.PowerLaw(1, 1e3, 1e6))
        r2 = sim.reweight(physical_energy=siren.dist.PowerLaw(3, 1e3, 1e6))
        # r1 and r2 should have different weights from each other
        any_different = any(
            abs(a - b) > 1e-15
            for a, b in zip(r1.weights, r2.weights)
        )
        assert any_different


class TestEventsAlias:
    """`n_events` is a deprecated alias for `events`."""

    def _kwargs(self, siren):
        return dict(
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            targets="Nucleon",
            process="CC",
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()
            ),
        )

    def test_n_events_warns_and_sets_count(self):
        import siren
        with pytest.warns(DeprecationWarning, match="n_events"):
            sim = siren.Simulation(n_events=7, **self._kwargs(siren))
        assert "events=7" in repr(sim)

    def test_events_does_not_warn(self):
        import siren
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            sim = siren.Simulation(events=7, **self._kwargs(siren))
        assert "events=7" in repr(sim)

    def test_both_events_and_n_events_raises(self):
        import siren
        with pytest.raises(TypeError, match="both 'events' and 'n_events'"):
            siren.Simulation(events=7, n_events=7, **self._kwargs(siren))


# ==================================================================== #
#  Measure consistency: structural and closure tests                    #
# ==================================================================== #

class TestMeasureConsistency:
    """Verify that multi-channel phase space densities are consistent
    across channel types."""

    def _make_2body_record(self):
        """Create a template InteractionRecord for a 2-body decay."""
        import siren
        record = siren.dataclasses.InteractionRecord()
        record.signature.primary_type = siren.particles.NuMu
        record.signature.target_type = siren.dataclasses.Particle.ParticleType.Decay
        record.signature.secondary_types = [
            siren.particles.EMinus, siren.particles.PiPlus
        ]
        # Parent: 300 MeV mass, 1 GeV energy, moving along z
        M = 0.300
        E = 1.0
        p = (E * E - M * M) ** 0.5
        record.primary_mass = M
        record.primary_momentum = [E, 0, 0, p]
        record.secondary_masses = [0.000511, 0.13957]
        record.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
        # Decay position far outside the target volume so the
        # directed channel has meaningful solid angle coverage
        record.interaction_vertex = [0, 0, -3000]
        return record

    def test_isotropic_self_consistency(self):
        """Isotropic channel: sample then evaluate own density.
        Should always return 1/(4*pi)."""
        import copy
        import siren
        import math

        record = self._make_2body_record()
        iso = siren.injection.Isotropic2BodyChannel(0)
        rng = siren.utilities.SIREN_random(42)

        for _ in range(20):
            r = copy.deepcopy(record)
            iso.Sample(rng, None, r)
            d = iso.Density(None, r)
            assert abs(d - 1.0 / (4 * math.pi)) < 1e-10

    def test_directed_density_positive_for_isotropic_samples(self):
        """Events from isotropic channel should get non-negative
        density from the directed channel (may be 0 if they miss
        the target, but never negative or NaN)."""
        import copy
        import siren
        import math

        record = self._make_2body_record()
        fid = siren.get_fiducial_volume("IceCube")
        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(fid, 0)
        rng = siren.utilities.SIREN_random(42)

        for _ in range(50):
            r = copy.deepcopy(record)
            iso.Sample(rng, None, r)
            d = directed.Density(None, r)
            assert d >= 0
            assert math.isfinite(d)

    def test_closure_isotropic_vs_multichannel(self):
        """Closure test: generate from multi-channel (isotropic + directed),
        weight by isotropic_density / multi_channel_density, and check
        the mean is consistent with 1.0.

        This is the gold standard test for measure consistency."""
        import copy
        import siren
        import math

        record = self._make_2body_record()
        fid = siren.get_fiducial_volume("IceCube")

        iso = siren.injection.Isotropic2BodyChannel(0)
        directed = siren.injection.DetectorDirected2BodyChannel(fid, 0)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, directed]
        mc.weights = [0.5, 0.5]

        rng = siren.utilities.SIREN_random(123)
        N = 500
        weights = []

        for _ in range(N):
            r = copy.deepcopy(record)
            mc.Sample(rng, None, r)
            d_iso = iso.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_iso / d_mc)

        assert len(weights) > N * 0.5, (
            f"Too few valid weights: {len(weights)}/{N}"
        )

        mean_w = sum(weights) / len(weights)
        # Should be 1.0 within statistical uncertainty.
        # With N=500, typical std(w)/sqrt(N) ~ 0.05
        assert 0.8 < mean_w < 1.2, (
            f"Closure test failed: mean weight = {mean_w:.4f} "
            f"(expected ~1.0)"
        )

    def test_validate_channels_no_diagnostics(self):
        """ValidateChannels should return no diagnostics for a
        well-formed multi-channel with compatible measures."""
        import siren

        record = self._make_2body_record()
        fid = siren.get_fiducial_volume("IceCube")

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [
            siren.injection.Isotropic2BodyChannel(0),
            siren.injection.DetectorDirected2BodyChannel(fid, 0),
        ]
        mc.weights = [0.5, 0.5]

        rng = siren.utilities.SIREN_random(99)
        diagnostics = mc.ValidateChannelDensities(rng, None, record, 50)

        nan_or_neg = [d for d in diagnostics
                      if "NaN" in d or "negative" in d]
        assert len(nan_or_neg) == 0, (
            f"Unexpected diagnostics: {nan_or_neg}"
        )

    # ---- 3-body channel tests ----

    @staticmethod
    def _make_3body_record():
        import siren
        import math

        record = siren.dataclasses.InteractionRecord()
        sig = siren.dataclasses.InteractionSignature()
        sig.primary_type = siren.dataclasses.ParticleType(211)
        sig.secondary_types = [
            siren.dataclasses.ParticleType(-13),
            siren.dataclasses.ParticleType(14),
            siren.dataclasses.ParticleType(5922),
        ]
        record.signature = sig
        M = 1.0
        E = 2.0
        record.primary_mass = M
        record.primary_momentum = [E, 0, 0, math.sqrt(E * E - M * M)]
        record.secondary_masses = [0.1, 0.3, 0.2]
        record.secondary_momenta = [[0, 0, 0, 0]] * 3
        record.secondary_helicities = [0, 0, 0]
        record.interaction_vertex = [0, 0, 0]
        record.primary_initial_position = [0, 0, 0]
        return record

    def test_3body_closure(self):
        """Closure test for DetectorDirected3BodyChannel.

        Mix a flat (isotropic) 3-body channel with the directed one
        via MultiChannelPhaseSpace. The importance weights
        flat_density / mixture_density must average to 1.0."""
        import copy
        import math
        import siren

        Cone = siren.injection.DirectedMode.Cone
        huge = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
            1e6, 1e6, 1e6)
        placement = siren.geometry.Placement(
            siren.math.Vector3D(0, 0, 10))
        target = siren.geometry.Box(placement, 2.0, 2.0, 2.0)

        kw = dict(
            factorization=siren.injection.ThreeBodyMode.Recursive,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_index=2,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
            mode=Cone,
        )
        flat = siren.injection.DetectorDirected3BodyChannel(target=huge, **kw)
        directed = siren.injection.DetectorDirected3BodyChannel(
            target=target, **kw)

        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [flat, directed]
        mc.weights = [0.5, 0.5]

        rng = siren.utilities.SIREN_random(999)
        N = 5000
        weights = []

        for _ in range(N):
            r = copy.deepcopy(self._make_3body_record())
            mc.Sample(rng, None, r)
            d_flat = flat.Density(None, r)
            d_mc = mc.Density(None, r)
            if d_mc > 0:
                weights.append(d_flat / d_mc)

        assert len(weights) > N * 0.5, (
            f"Too few valid weights: {len(weights)}/{N}"
        )

        mean_w = sum(weights) / len(weights)
        assert 0.85 < mean_w < 1.15, (
            f"3-body closure test failed: mean weight = {mean_w:.4f} "
            f"(expected ~1.0)"
        )

    def test_3body_density_normalization(self):
        """Check that the directed 3-body density integrates to 1.

        Sample from a flat (isotropic) reference and compute the mean
        of directed_density / flat_density. This estimates the integral
        of the directed density over the phase space."""
        import copy
        import math
        import siren

        Cone = siren.injection.DirectedMode.Cone
        huge = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
            1e6, 1e6, 1e6)
        placement = siren.geometry.Placement(
            siren.math.Vector3D(0, 0, 10))
        target = siren.geometry.Box(placement, 2.0, 2.0, 2.0)

        kw = dict(
            factorization=siren.injection.ThreeBodyMode.Recursive,
            spectator_index=0,
            pair_first_index=1,
            pair_second_index=2,
            directed_index=2,
            mass_mode=siren.injection.InvariantMassMode.Uniform,
            mode=Cone,
        )
        flat = siren.injection.DetectorDirected3BodyChannel(target=huge, **kw)
        directed = siren.injection.DetectorDirected3BodyChannel(
            target=target, **kw)

        rng = siren.utilities.SIREN_random(777)
        N = 50000
        weights = []

        for _ in range(N):
            r = copy.deepcopy(self._make_3body_record())
            flat.Sample(rng, None, r)
            d_flat = flat.Density(None, r)
            d_dir = directed.Density(None, r)
            if d_flat > 0:
                weights.append(d_dir / d_flat)

        assert len(weights) > N * 0.5
        mean_w = sum(weights) / len(weights)
        se = (sum((w - mean_w) ** 2 for w in weights)
              / len(weights)) ** 0.5 / len(weights) ** 0.5
        assert abs(mean_w - 1.0) < max(5 * se, 0.15), (
            f"Density normalization failed: mean = {mean_w:.4f} "
            f"+/- {se:.4f} (expected 1.0)"
        )
