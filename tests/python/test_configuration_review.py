"""Public configuration and generation regressions from the U04 review."""
import pickle
from types import SimpleNamespace

import numpy as np
import pytest
import siren

from siren import distributions as d, injection as inj
from siren._validation import validate_ordering
from test_injector_api import (
    _working_injector, _forced_failure_injector, _chain_injector,
    _process_with_phase_space,
)
from test_weighter_inheritance import _injector


def test_batches_share_attempt_denominator():
    split, single = _working_injector(events=5), _working_injector(events=5)
    a = split.generate(5, on_shortfall="raise")
    first_attempts = split.injection_attempts
    a += split.generate(5, on_shortfall="raise")
    b = single.generate(10, on_shortfall="raise")
    assert split.injection_attempts >= first_attempts + 5
    assert split.injection_attempts == single.injection_attempts
    np.testing.assert_array_equal(
        inj.Weighter(split).weight_all(a), inj.Weighter(single).weight_all(b))


@pytest.mark.parametrize("round_trip", ["live", "pickle", "archive"])
def test_configured_count_survives_retry_budget(round_trip, tmp_path):
    generator = _working_injector(events=5)
    generator.generate(2, on_shortfall="raise")
    if round_trip == "pickle":
        generator = pickle.loads(pickle.dumps(generator))
    elif round_trip == "archive":
        path = str(tmp_path / "count")
        generator.save(path)
        generator = inj.Injector.load(path)
    assert generator.number_of_events == 5
    assert len(generator) == 5
    assert len(generator.generate(on_shortfall="raise")) == 5
    assert generator.injection_attempts >= 7


def test_default_retry_budget_cannot_overflow_native_counter():
    generator = _forced_failure_injector(events=1)
    with pytest.raises(siren.errors.InjectionShortfall, match="min_efficiency") as exc:
        generator.generate(4_294_968, on_shortfall="raise", min_efficiency=0.5)
    assert exc.value.report.attempts == 50


def test_zero_efficiency_stops_at_warmup():
    generator = _forced_failure_injector(events=100, max_attempts=1000)
    with pytest.raises(siren.errors.InjectionShortfall, match="min_efficiency") as exc:
        generator.generate(on_shortfall="raise", min_efficiency=0.5)
    assert exc.value.report.attempts == 50
    assert exc.value.report.successes == 0
    assert generator.number_of_events == 100


@pytest.mark.parametrize("position", [
    d.FixedTargetPositionDistribution,
    d.PrimaryPhysicalVertexDistribution,
    d.PrimaryBoundedVertexDistribution,
])
def test_transport_position_requires_mass_and_energy_before_sampling(position):
    dv = d.DistributionVariable
    source = SimpleNamespace(
        RequiredVariables=lambda: set(),
        ProvidedVariables=lambda: {dv.InitialPosition, dv.PrimaryDirection},
    )
    with pytest.raises(ValueError, match="PrimaryEnergy|PrimaryMass"):
        validate_ordering([source, position(), d.PrimaryMass(1), d.Monoenergetic(2)])


@pytest.mark.parametrize("secondary", [False, True])
def test_built_injector_refuses_distribution_replacement(secondary):
    generator = _chain_injector() if secondary else _working_injector()
    generator.engine
    with pytest.raises(siren.errors.ConfigurationError, match="built|initializ"):
        if secondary:
            generator.secondary_injection_distributions = {
                siren.particles.NuMu: [d.SecondaryPhysicalVertexDistribution()]}
        else:
            generator.primary_injection_distributions = (
                generator.primary_injection_distributions)


@pytest.mark.parametrize("keyword", [
    "primary_weighting_mode", "primary_phase_spaces", "primary_decay_channels",
])
def test_vertex_injector_rejects_legacy_process_keywords(keyword):
    reference = _injector(chain=False)
    values = {"primary_weighting_mode": siren.PropagatedFromCreation(),
              "primary_phase_spaces": {}, "primary_decay_channels": []}
    with pytest.raises((TypeError, siren.errors.ConfigurationError), match="legacy|mix"):
        inj.Injector(events=2, detector=reference.detector_model,
                     primary=reference.primary, **{keyword: values[keyword]})


@pytest.mark.parametrize("surface", ["spec", "legacy", "unknown"])
def test_weighter_rejects_ignored_arguments(surface):
    generator = _working_injector()
    with pytest.raises((TypeError, siren.errors.ConfigurationError), match="legacy|mix|unknown"):
        if surface == "spec":
            inj.Weighter(generator, primary_physical_distributions=[])
        elif surface == "legacy":
            inj.Weighter(injectors=[generator], primary_physical=[])
        else:
            inj.Weighter(generator, overrides={"typo_physical": []})


@pytest.mark.parametrize("secondary", [False, True])
def test_pool_rejects_full_weighting_mode_mismatch(secondary):
    first, second = _chain_injector(), _chain_injector()
    process = (second.engine.GetSecondaryProcessMap()[siren.particles.NuMu]
               if secondary else second.engine.GetPrimaryProcess())
    process.weighting_mode = siren.Fixed()
    with pytest.raises(siren.errors.ConfigurationError, match="mode"):
        inj.Weighter(first, second).engine


def test_three_channel_sum_preserves_written_coefficients():
    a, b, c = [siren.channels.isotropic(0) for _ in range(3)]
    mixture = (0.2 * a + 0.3 * b + 0.5 * c).compile(
        siren.channels._template_signature())
    assert list(mixture.weights) == pytest.approx([0.2, 0.3, 0.5])


def test_density_probe_diagnostics_are_visible():
    generator = _forced_failure_injector()
    validation = SimpleNamespace(
        probe_channel_densities=lambda *args: ["probe template outside support"])
    with pytest.warns(RuntimeWarning, match="outside support"):
        generator._probe_phase_spaces(validation, _process_with_phase_space(), [],
                                      siren.utilities.SIREN_random(1))


def test_weighter_pickle_retains_declared_physical_target():
    generator = _injector(chain=False)
    trees = generator.generate(3, on_shortfall="raise")
    weighter = inj.Weighter(generator)
    restored = pickle.loads(pickle.dumps(weighter))
    np.testing.assert_array_equal(restored.weight_all(trees), weighter.weight_all(trees))


def test_instance_load_replaces_configuration(tmp_path):
    source = _working_injector(events=7)
    source.generate(2, on_shortfall="raise")
    filename = str(tmp_path / "instance-load")
    source.save(filename)
    target = _working_injector(events=3)
    target.engine
    assert target.load(filename) is target
    assert target.number_of_events == 7
    assert target.injection_attempts == source.injection_attempts


def test_count_is_stable_inside_progress_callback():
    generator = _working_injector(events=5)
    observed = []
    generator.generate(2, progress=lambda *args: observed.append(generator.number_of_events))
    assert observed == [5, 5]


@pytest.mark.parametrize("positional_list", [False, True])
def test_weighter_rejects_two_injector_arguments(positional_list):
    first, second = _working_injector(), _working_injector()
    with pytest.raises(TypeError, match="positional.*injectors"):
        inj.Weighter([first] if positional_list else first, injectors=[second])


def test_injector_pickle_retains_retry_policy():
    generator = _forced_failure_injector(events=1, max_attempts=3)
    restored = pickle.loads(pickle.dumps(generator))
    with pytest.raises(siren.errors.InjectionShortfall) as exc:
        restored.generate(on_shortfall="raise")
    assert exc.value.report.attempts == 3


def test_injector_accepts_legacy_three_field_pickle_state():
    generator = _working_injector(events=3)
    state = generator.__getstate__()[:3]
    restored = inj.Injector.__new__(inj.Injector)
    restored.__setstate__(state)
    assert restored.number_of_events == 3
    assert len(restored.generate(on_shortfall="raise")) == 3


def test_old_weighter_pickle_accepts_native_user_attribute():
    from pathlib import Path
    import base64
    encoded = (Path(__file__).parent / "fixtures" /
               "weighter-native-attribute-v0.b64").read_text()
    restored = pickle.loads(base64.b64decode(encoded))
    assert restored.native == "user metadata"
    assert restored.injectors is None
