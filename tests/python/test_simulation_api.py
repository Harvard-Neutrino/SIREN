"""
Tests for the new high-level Simulation API (Phases 1-3).
"""

import pytest


# ==================================================================== #
#  Phase 1: Sub-namespaces and top-level convenience                    #
# ==================================================================== #

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
        assert DV.InitialPosition in sv
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

    def test_top_level_count_reasonable(self):
        """Top-level namespace should not be bloated."""
        import siren
        public = [x for x in dir(siren) if not x.startswith("_")]
        # Should be around 20, definitely under 30
        assert len(public) < 30, f"Too many top-level names: {len(public)}"


# ==================================================================== #
#  Phase 2: Simulation class                                            #
# ==================================================================== #

class TestSimulationValidation:
    """Simulation should validate inputs at construction time."""

    def test_missing_energy_raises(self):
        import siren
        with pytest.raises(ValueError, match="energy"):
            siren.Simulation(
                n_events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                target="Nucleon",
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
                n_events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                target="Nucleon",
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
                n_events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                target="Nucleon",
                process="CC",
                energy=siren.dist.PowerLaw(2, 1e3, 1e6),
                direction=siren.dist.IsotropicDirection(),
            )

    def test_conflicting_energy_raises(self):
        import siren
        with pytest.raises(ValueError, match="Cannot specify both"):
            siren.Simulation(
                n_events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                target="Nucleon",
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
                n_events=1,
                detector="IceCube",
                primary="NuMu",
                interactions="CSMSDISSplines",
                target="Nucleon",
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
                n_events=1,
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
                n_events=1,
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
            n_events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=1,
            detector="IceCube",
            primary=siren.particles.NuMu,
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=1,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=42,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
            n_events=5,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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


# ==================================================================== #
#  Phase 3: Ergonomic improvements                                      #
# ==================================================================== #

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


class TestWeighterBatch:
    """Weighter.weight_all should batch-weight events."""

    def test_weight_all(self):
        import siren
        sim = siren.Simulation(
            n_events=3,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
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
