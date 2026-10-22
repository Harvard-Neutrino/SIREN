from simulation_fixtures import offline_detector, offline_fiducial
import pytest

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

class TestWeighterBatch:
    """Weighter.weight_all should batch-weight events."""

    def test_weight_all(self):
        import siren
        sim = siren.Simulation(
            events=3,
            detector=offline_detector(),
            primary="NuMu",
            interactions=[siren.interactions.DummyCrossSection()],
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

class TestEventsAlias:
    """`n_events` is a deprecated alias for `events`."""

    def _kwargs(self, siren):
        return dict(
            detector=offline_detector(),
            primary="NuMu",
            interactions=[siren.interactions.DummyCrossSection()],
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

    def test_legacy_slots_warn_even_with_events(self):
        import siren
        import warnings
        with pytest.warns(DeprecationWarning, match="convenience slots"):
            sim = siren.Simulation(events=7, **self._kwargs(siren))
        assert "events=7" in repr(sim)

    def test_both_events_and_n_events_raises(self):
        import siren
        with pytest.raises(TypeError, match="both events and n_events"):
            siren.Simulation(events=7, n_events=7, **self._kwargs(siren))
