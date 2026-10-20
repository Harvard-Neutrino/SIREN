"""The public compiler preserves native fixed-seed generation and weights."""

import pytest
import siren
from test_trivial_cross_section import (
    detector_model, _make_injector, _constant_cross_section, _physical_process,
)
from test_native_phase_space import fixed_decay_setup, record_state
from test_phase_space_decay import decay


def test_vertex_matches_existing_scattering_events_and_weights(detector_model):
    model = _constant_cross_section()
    raw, dists = _make_injector(detector_model, model, seed=7821, n_inject=100)
    physical = _physical_process(model, dists)
    vertex = siren.Vertex(
        raw.GetPrimaryProcess().primary_type, model,
        distributions=dists, physical=physical.distributions,
    )
    compiled = siren.injection.Injector(
        detector=detector_model, primary=vertex, seed=7821, events=100,
    )
    reference = [raw.GenerateEvent() for _ in range(100)]
    candidate = compiled.generate(on_failure="raise", on_shortfall="raise")
    raw_weights = siren.injection._Weighter([raw], detector_model, physical)
    weights = siren.injection.Weighter(compiled)
    assert [record_state(e) for e in candidate] == [record_state(e) for e in reference]
    assert [weights(e) for e in candidate] == [raw_weights.EventWeight(e) for e in reference]
    assert compiled.injection_attempts == raw.InjectionAttempts() == 100


def test_stationary_decay_keeps_native_direction_marginalization():
    model = decay(partial=1e-9, total=1e-9)
    mass = 0.1349768
    raw, raw_weights = fixed_decay_setup(
        model, siren.particles.Pi0, mass, 8, energy=mass,
        direction=siren.dist.IsotropicDirection(),
    )
    process = raw.GetPrimaryProcess()
    vertex = siren.Vertex(
        siren.particles.Pi0, model,
        distributions=list(process.distributions),
        physical=list(process.distributions),
        weighting=siren.Fixed(),
    )
    generator = siren.injection.Injector(
        detector=raw.GetDetectorModel(), primary=vertex, seed=6217, events=8,
    )
    expected = [raw.GenerateEvent() for _ in range(8)]
    actual = generator.generate(on_failure="raise", on_shortfall="raise")
    weights = siren.injection.Weighter(generator)
    assert [record_state(e) for e in actual] == [record_state(e) for e in expected]
    assert [weights(e) for e in actual] == pytest.approx(
        [raw_weights.EventWeight(e) for e in expected], rel=0, abs=0
    )
    assert sum(weights(e) for e in actual) == pytest.approx(1.0)
