"""Injector serialization and Fixed-vertex channel-factor pins.

Data-free and fixed-seed throughout:

  (a) The raw C++ injector (siren.injection._Injector) round-trips
      InjectionAttempts, InjectedEvents, FailedEvents, and each process's
      weighting mode across SaveInjector/LoadInjector.

  (b) FixedVertexChannelSelectionProbability returns exactly 1.0 for a
      single-channel process (an exact float no-op) and
      selected_rate/total_rate when multiple channels compete, so a Fixed
      vertex charges the same channel-selection factor the injector's rate
      selection applied.
"""

import math

import pytest

import siren
from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import math as smath
from siren import utilities

_NuMu = dc.Particle.ParticleType.NuMu


def _distributions(max_distance):
    return [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(
            smath.Vector3D(0, 0, 0), max_distance),
    ]


def _raw_forced_failure_injector(events, seed, weighting_mode):
    """A raw _Injector whose every attempt fails (bare detector, zero source
    distance means every path misses its target), with a Fixed-mode primary
    process. No detector data files are needed."""
    dm = siren.detector.DetectorModel()
    collection = interactions.InteractionCollection(
        _NuMu, [interactions.DummyCrossSection()])
    primary = injection.PrimaryInjectionProcess(_NuMu, collection)
    primary.distributions = _distributions(0.0)
    primary.weighting_mode = weighting_mode
    random = utilities.SIREN_random(seed)
    return injection._Injector(events, dm, primary, random)


# ------------------------------------------------------------------ #
#  (a) SaveInjector/LoadInjector round-trip                           #
# ------------------------------------------------------------------ #

def test_save_load_preserves_counters_and_weighting_mode(tmp_path):
    """Counters and the primary process's weighting mode survive a raw
    SaveInjector/LoadInjector round-trip."""
    events = 6
    fixed = injection.VertexWeightingMode.Fixed()
    inj = _raw_forced_failure_injector(events=events, seed=13, weighting_mode=fixed)

    # Every attempt fails; GenerateEvent catches the failure internally and
    # returns an empty tree, so drive it up to the configured attempt cap.
    for _ in range(events):
        tree = inj.GenerateEvent()
        assert len(tree.tree) == 0

    assert inj.InjectionAttempts() == events
    assert inj.InjectedEvents() == 0
    assert inj.FailedEvents() == events
    # The accounting invariant that motivated serializing FailedEvents.
    assert inj.InjectionAttempts() == inj.InjectedEvents() + inj.FailedEvents()

    path = str(tmp_path / "injector_roundtrip")
    inj.SaveInjector(path)

    reloaded = injection._Injector(1, path, utilities.SIREN_random(0))

    assert reloaded.InjectionAttempts() == inj.InjectionAttempts()
    assert reloaded.InjectedEvents() == inj.InjectedEvents()
    assert reloaded.FailedEvents() == inj.FailedEvents()
    assert reloaded.InjectionAttempts() == (
        reloaded.InjectedEvents() + reloaded.FailedEvents())

    # The archive round-trip preserves the Fixed weighting mode.
    assert reloaded.GetPrimaryProcess().GetWeightingMode() == fixed
    assert reloaded.GetPrimaryProcess().GetWeightingMode() != (
        injection.VertexWeightingMode.Propagated())


def test_failed_events_default_when_not_yet_generated(tmp_path):
    """A freshly-built injector round-trips with all counters at zero."""
    inj = _raw_forced_failure_injector(
        events=3, seed=1, weighting_mode=injection.VertexWeightingMode.Propagated())
    path = str(tmp_path / "injector_fresh")
    inj.SaveInjector(path)
    reloaded = injection._Injector(1, path, utilities.SIREN_random(0))
    assert reloaded.InjectionAttempts() == 0
    assert reloaded.InjectedEvents() == 0
    assert reloaded.FailedEvents() == 0
    assert reloaded.GetPrimaryProcess().GetWeightingMode() == (
        injection.VertexWeightingMode.Propagated())


def test_save_load_resumes_the_rng_stream(tmp_path):
    """A version-2 archive restores the RNG engine state, so a reloaded injector
    RESUMES its generation stream where the saved one left off -- it does not
    restart from the seed, and the fallback engine handed to the loading
    constructor is discarded. The saved seed survives only as a label."""
    inj = _raw_forced_failure_injector(
        events=10, seed=13, weighting_mode=injection.VertexWeightingMode.Fixed())
    engine = inj.GetRandom()
    # Advance the engine, then snapshot the injector mid-stream.
    for _ in range(20):
        engine.Uniform(0.0, 1.0)
    path = str(tmp_path / "injector_rng")
    inj.SaveInjector(path)
    # What the original engine produces next is the continuation to reproduce.
    continuation = [engine.Uniform(0.0, 1.0) for _ in range(8)]

    # Reload with a DIFFERENT fallback seed; a version-2 archive ignores it.
    reloaded = injection._Injector(1, path, utilities.SIREN_random(999))
    resumed = [reloaded.GetRandom().Uniform(0.0, 1.0) for _ in range(8)]
    assert resumed == continuation

    # It is a resume, not a restart-from-seed: seed 13 replayed from the start
    # gives the pre-snapshot draws, which must not match the continuation.
    restart = [utilities.SIREN_random(13).Uniform(0.0, 1.0) for _ in range(8)]
    assert resumed != restart
    # The originating seed is preserved as a label (not the fallback 999).
    assert reloaded.GetRandom().get_seed() == 13


# ------------------------------------------------------------------ #
#  FixedVertexChannelSelectionProbability                         #
# ------------------------------------------------------------------ #

class _ConstWidthDecay(siren.interactions.Decay):
    """Minimal data-free decay double with a single fixed signature and a
    constant total width, so its selection rate is proportional to `width`.
    The two secondary types choose the signature, so two instances with
    different secondaries present as two competing channels."""

    def __init__(self, secondary_types, width):
        siren.interactions.Decay.__init__(self)
        self._width = float(width)
        self._sig = siren.dataclasses.InteractionSignature()
        self._sig.primary_type = _NuMu
        self._sig.target_type = dc.Particle.ParticleType.Decay
        self._sig.secondary_types = list(secondary_types)

    def equal(self, other):
        return self is other

    def TotalDecayWidthAllFinalStates(self, record):
        return self._width

    def TotalDecayWidth(self, arg):
        # Both the (ParticleType) and (record) overloads dispatch here.
        return self._width

    def DifferentialDecayWidth(self, record):
        return self._width

    def SampleFinalState(self, csdr, random):
        pass

    def GetPossibleSignatures(self):
        return [self._sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        return [self._sig]

    def FinalStateProbability(self, record):
        return 1.0

    def DensityVariables(self):
        return []


def _moving_decay_record(signature, mass=0.1, energy=1.0):
    record = siren.dataclasses.InteractionRecord()
    record.signature = signature
    record.primary_mass = mass
    pz = math.sqrt(max(energy**2 - mass**2, 0.0))
    record.primary_momentum = [energy, 0.0, 0.0, pz]
    record.primary_initial_position = [0.0, 0.0, 0.0]
    record.interaction_vertex = [0.0, 0.0, 0.0]
    record.secondary_masses = []
    record.secondary_momenta = []
    return record


def test_single_channel_fixed_factor_is_exactly_one():
    """One candidate signature -> exact 1.0 (single-channel Fixed configs are
    unaffected)."""
    dm = siren.detector.DetectorModel()
    decay = _ConstWidthDecay(
        [dc.Particle.ParticleType.NuE, dc.Particle.ParticleType.Gamma], width=1.0)
    collection = interactions.InteractionCollection(_NuMu, [decay])
    record = _moving_decay_record(decay.GetPossibleSignatures()[0])

    factor = injection.FixedVertexChannelSelectionProbability(dm, collection, record)
    assert factor == 1.0


def test_two_channel_fixed_factor_is_branching_fraction():
    """Two equal-width channels -> the observed channel carries exactly 1/2,
    matching the injector's rate selection; unequal widths give the width
    fraction. The local decay references keep the python doubles alive for the
    duration of the C++ calls."""
    dm = siren.detector.DetectorModel()

    decay_a = _ConstWidthDecay(
        [dc.Particle.ParticleType.NuE, dc.Particle.ParticleType.Gamma], width=1.0)
    decay_b = _ConstWidthDecay(
        [dc.Particle.ParticleType.NuMu, dc.Particle.ParticleType.Gamma], width=1.0)
    collection = interactions.InteractionCollection(_NuMu, [decay_a, decay_b])

    # A record on channel A: equal widths -> selected/total = 1/2.
    record_a = _moving_decay_record(decay_a.GetPossibleSignatures()[0])
    factor_a = injection.FixedVertexChannelSelectionProbability(dm, collection, record_a)
    assert factor_a == pytest.approx(0.5)

    # The same record through the single-channel collection has factor 1.0, so
    # the two-channel factor is a genuine 1/2 relative to the single-channel
    # case.
    single = interactions.InteractionCollection(_NuMu, [decay_a])
    assert injection.FixedVertexChannelSelectionProbability(
        dm, single, record_a) == 1.0

    # Unequal widths (3 vs 1) -> width fraction 3/4 on the wider channel.
    decay_wide = _ConstWidthDecay(
        [dc.Particle.ParticleType.NuE, dc.Particle.ParticleType.Gamma], width=3.0)
    decay_narrow = _ConstWidthDecay(
        [dc.Particle.ParticleType.NuMu, dc.Particle.ParticleType.Gamma], width=1.0)
    collection2 = interactions.InteractionCollection(_NuMu, [decay_wide, decay_narrow])
    record_wide = _moving_decay_record(decay_wide.GetPossibleSignatures()[0])
    factor_wide = injection.FixedVertexChannelSelectionProbability(
        dm, collection2, record_wide)
    assert factor_wide == pytest.approx(0.75)
