"""Facade persistence coverage for per-process phase-space maps."""

import math
import os
import pickle
from pathlib import Path

import pytest

import siren
from siren import _util
from siren import dataclasses as dc
from siren import detector
from siren import distributions
from siren import injection
from siren import interactions
from siren import math as smath
from siren import utilities


_NuMu = dc.Particle.ParticleType.NuMu
_NuE = dc.Particle.ParticleType.NuE
_WEIGHTS = [0.25, 0.75]
_KP_ACCUMULATOR = [1.5, 2.5]
_KP_COUNT = 7
_KP_SUCC_SELECT = [3.0, 4.0]
_KP_FAIL_SELECT = [0.5, 1.5]


def _injection_distributions(max_distance):
    return [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(
            smath.Vector3D(0.0, 0.0, 0.0), max_distance),
    ]


def _physical_distributions():
    return [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]


def _signature_for(xs, primary_type=_NuMu):
    return next(
        sig for sig in xs.GetPossibleSignatures()
        if sig.primary_type == primary_type)


def _stateful_mixture(xs, signature):
    mixture = injection.MultiChannelPhaseSpace(
        [
            injection.PhysicalCrossSectionChannel(xs, signature),
            injection.PhysicalCrossSectionChannel(xs, signature),
        ],
        _WEIGHTS,
        True,
    )
    mixture.kp_accumulator = list(_KP_ACCUMULATOR)
    mixture.kp_count = _KP_COUNT
    mixture.kp_succ_select = list(_KP_SUCC_SELECT)
    mixture.kp_fail_select = list(_KP_FAIL_SELECT)
    return mixture


def _assert_mixture_state(mixture):
    assert list(mixture.weights) == _WEIGHTS
    assert list(mixture.kp_accumulator) == _KP_ACCUMULATOR
    assert mixture.kp_count == _KP_COUNT
    assert list(mixture.kp_succ_select) == _KP_SUCC_SELECT
    assert list(mixture.kp_fail_select) == _KP_FAIL_SELECT
    assert [channel.Name() for channel in mixture.channels] == [
        "PhysicalCrossSection",
        "PhysicalCrossSection",
    ]


def _primary_facade(events=10, seed=13, detector_model=None,
                    max_distance=0.0, phase_space_primary=_NuMu,
                    attach_phase_space=True):
    dm = detector_model or detector.DetectorModel()
    xs = interactions.DummyCrossSection()
    signature = _signature_for(xs, phase_space_primary)
    injector = injection.Injector(
        detector=dm,
        number_of_events=events,
        seed=seed,
        primary_type=_NuMu,
        primary_interactions=[xs],
        primary_injection_distributions=(
            _injection_distributions(max_distance)),
    )
    injector._build()
    mixture = _stateful_mixture(xs, signature)
    if attach_phase_space:
        injector.primary_phase_spaces = {signature: mixture}
    return injector, signature, mixture


def _raw_injector_with_shared_maps(events=10, seed=13):
    dm = detector.DetectorModel()
    xs = interactions.DummyCrossSection()
    collection = interactions.InteractionCollection(_NuMu, [xs])
    signature = _signature_for(xs)
    mixture = _stateful_mixture(xs, signature)

    primary = injection.PrimaryInjectionProcess()
    primary.primary_type = _NuMu
    primary.interactions = collection
    primary.distributions = _injection_distributions(0.0)
    primary.SetPhaseSpace(signature, mixture)

    secondary = injection.SecondaryInjectionProcess()
    secondary.primary_type = _NuMu
    secondary.interactions = collection
    secondary.distributions = [
        distributions.SecondaryPhysicalVertexDistribution()]
    secondary.SetPhaseSpace(signature, mixture)

    raw = injection._Injector(
        events, dm, primary, [secondary], utilities.SIREN_random(seed))
    return raw, signature, mixture


def _secondary_facade(tmp_path):
    raw, signature, mixture = _raw_injector_with_shared_maps()
    seed_archive = str(tmp_path / "secondary_seed.siren_injector")
    raw.SaveInjector(seed_archive)
    return injection.Injector.load(seed_archive), signature, mixture


def _assert_facade_map_identity(injector, signature, secondary=False):
    primary_map = injector.engine.GetPrimaryProcess().GetPhaseSpaceMap()
    primary_mixture = primary_map[signature]
    assert injector.primary_phase_spaces[signature] is primary_mixture

    if secondary:
        secondary_process = injector.engine.GetSecondaryProcessMap()[_NuMu]
        secondary_mixture = secondary_process.GetPhaseSpaceMap()[signature]
        assert injector.secondary_phase_spaces[_NuMu][signature] is (
            secondary_mixture)
        assert secondary_mixture is primary_mixture

    _assert_mixture_state(primary_mixture)
    return primary_mixture


def test_facade_save_load_roundtrips_shared_phase_space_maps(tmp_path):
    source, signature, _mixture = _secondary_facade(tmp_path)
    archive = str(tmp_path / "facade_roundtrip.siren_injector")

    source.save(archive)
    loaded = injection.Injector.load(archive)

    _assert_facade_map_identity(loaded, signature, secondary=True)


def test_facade_pickle_roundtrips_shared_phase_space_maps(tmp_path):
    source, signature, _mixture = _secondary_facade(tmp_path)

    restored = pickle.loads(
        pickle.dumps(source, protocol=pickle.HIGHEST_PROTOCOL))

    _assert_facade_map_identity(restored, signature, secondary=True)


def test_rebuild_after_load_reuses_loaded_phase_space_objects(tmp_path):
    source, signature, _mixture = _primary_facade()
    archive = str(tmp_path / "rebuild.siren_injector")
    source.save(archive)
    loaded = injection.Injector.load(archive)
    loaded_mixture = loaded.primary_phase_spaces[signature]

    # Lazy reassembly reuses the mixture stored in the facade configuration.
    loaded._Injector__injector = None
    rebuilt_engine = loaded.engine
    assert rebuilt_engine is not None

    rebuilt_mixture = _assert_facade_map_identity(loaded, signature)
    assert rebuilt_mixture is loaded_mixture


def test_facade_weighter_save_load_accepts_channel_carrying_injector(tmp_path):
    source, signature, _mixture = _primary_facade()
    weighter = injection.Weighter(
        source, primary_physical=_physical_distributions())
    archive = str(tmp_path / "facade_weighter")

    weighter.save(archive)
    loaded = injection.Weighter()
    loaded.load(archive)

    loaded_injectors = loaded.engine.GetInjectors()
    assert len(loaded_injectors) == 1
    loaded_mixture = (
        loaded_injectors[0].GetPrimaryProcess().GetPhaseSpaceMap()[signature])
    _assert_mixture_state(loaded_mixture)


def test_forced_failure_resume_preserves_rng_and_phase_space_map(tmp_path):
    source, signature, _mixture = _primary_facade(events=10, seed=13)
    for _ in range(3):
        assert len(source.generate_event().tree) == 0
    assert source.injection_attempts == 3

    engine = source.engine.GetRandom()
    for _ in range(20):
        engine.Uniform(0.0, 1.0)
    archive = str(tmp_path / "forced_resume.siren_injector")
    source.save(archive)
    continuation = [engine.Uniform(0.0, 1.0) for _ in range(8)]

    loaded = injection.Injector.load(archive)
    resumed = [
        loaded.engine.GetRandom().Uniform(0.0, 1.0) for _ in range(8)]
    assert resumed == continuation
    assert loaded.seed == 13
    assert loaded.injection_attempts == 3
    assert loaded.injected_events == 0
    assert loaded.engine.FailedEvents() == 3
    assert loaded.injection_attempts == (
        loaded.injected_events + loaded.engine.FailedEvents())
    _assert_facade_map_identity(loaded, signature)


def test_old_physical_process_v1_fixture_loads_without_phase_spaces():
    fixture = (
        Path(__file__).resolve().parents[1]
        / "fixtures" / "injector_physical_process_v1.siren_injector")
    assert fixture.exists()

    loaded = injection.Injector.load(str(fixture))

    assert loaded.number_of_events == 3
    assert loaded.injection_attempts == 2
    assert loaded.injected_events == 0
    assert loaded.engine.FailedEvents() == 2
    assert loaded.injection_attempts == (
        loaded.injected_events + loaded.engine.FailedEvents())
    assert not loaded.engine.GetPrimaryProcess().HasAnyPhaseSpace()
    assert all(
        not process.HasAnyPhaseSpace()
        for process in loaded.engine.GetSecondaryProcessMap().values())
    assert loaded.primary_phase_spaces == {}
    assert loaded.secondary_phase_spaces == {}


def _load_ccm_detector():
    """Load CCM or skip only when its detector files are unavailable."""
    try:
        detector_dir = _util.get_detector_model_path("CCM")
    except ValueError as exc:
        pytest.skip("CCM detector model path unavailable: {}".format(exc))
    materials = os.path.join(detector_dir, "materials.dat")
    densities = os.path.join(detector_dir, "densities.dat")
    if not (os.path.exists(materials) and os.path.exists(densities)):
        pytest.skip("CCM detector data files not present")
    try:
        dm = detector.DetectorModel()
        dm.LoadMaterialModel(materials)
        dm.LoadDetectorModel(densities)
        return dm
    except RuntimeError as exc:
        pytest.skip("CCM detector model unavailable: {}".format(exc))


def _next_nonempty_event(injector, attempts=300):
    for _ in range(attempts):
        event = injector.generate_event()
        if event.tree:
            return event
    raise AssertionError("no event generated within {} attempts".format(attempts))


def _physics_bytes(event):
    """Serialize event physics after removing process-global ID allocation."""
    normalized = pickle.loads(
        pickle.dumps(event, protocol=pickle.HIGHEST_PROTOCOL))
    for node_index, datum in enumerate(normalized.tree):
        record = datum.record
        id_base = 1000 * node_index
        record.primary_id = dc.ParticleID(0, id_base)
        record.target_id = dc.ParticleID(0, id_base + 1)
        record.secondary_ids = [
            dc.ParticleID(0, id_base + 2 + secondary_index)
            for secondary_index in range(len(record.secondary_ids))]
    return pickle.dumps(normalized, protocol=pickle.HIGHEST_PROTOCOL)


def test_ccm_weighter_event_weight_matches_after_save_load(tmp_path):
    dm = _load_ccm_detector()
    source, signature, _mixture = _primary_facade(
        events=50, seed=1234, detector_model=dm, max_distance=25.0,
        attach_phase_space=False)
    event = _next_nonempty_event(source)
    source.primary_phase_spaces = {signature: _mixture}
    weighter = injection.Weighter(
        source, primary_physical=_physical_distributions())
    reference = weighter(event)
    assert math.isfinite(reference) and reference > 0.0

    archive = str(tmp_path / "ccm_weighter")
    weighter.save(archive)
    loaded = injection.Weighter()
    loaded.load(archive)

    assert loaded(event) == pytest.approx(reference, rel=1e-12)
    loaded_mixture = (
        loaded.engine.GetInjectors()[0]
        .GetPrimaryProcess().GetPhaseSpaceMap()[signature])
    _assert_mixture_state(loaded_mixture)


def test_ccm_generation_continues_exactly_with_phase_space_map(tmp_path):
    dm = _load_ccm_detector()
    source, signature, _mixture = _primary_facade(
        events=100, seed=4321, detector_model=dm, max_distance=25.0,
        phase_space_primary=_NuE)
    for _ in range(3):
        source.generate_event()

    archive = str(tmp_path / "ccm_continuation.siren_injector")
    source.save(archive)

    expected = []
    for _ in range(100):
        event = source.generate_event()
        expected.append(_physics_bytes(event))
        if event.tree:
            break
    assert expected and pickle.loads(expected[-1]).tree

    loaded = injection.Injector.load(archive)
    actual = [
        _physics_bytes(loaded.generate_event())
        for _ in range(len(expected))]
    assert actual == expected
    _assert_facade_map_identity(loaded, signature)
