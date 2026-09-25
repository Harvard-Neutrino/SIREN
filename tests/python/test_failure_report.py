"""FailureReason/FailureLedger cross the pybind boundary and are filterable.

FailureLedger.entries() keys on (depth, parent_pdg, FailureReason); reasons
must be usable as dict keys/values (hashable, comparable) from Python, a
fresh ledger must start empty, and a forced GenerateEvent failure must land
in the ledger under its FailureReason so callers can filter by reason.
"""
from __future__ import annotations

import math

import pytest

from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import geometry
from siren import math as smath
from siren import utilities

_NuMu = dc.Particle.ParticleType.NuMu


# --------------------------------------------------------------------------- #
# (a) FailureReason enum: importable, hashable, usable as a dict key           #
# --------------------------------------------------------------------------- #

def test_failure_reason_enum_members_present():
    """All documented FailureReason members are importable."""
    expected = [
        "Unspecified", "NoPathThroughVolume", "NoTargetsOnPath",
        "NoColumnDepthSolution", "KinematicallyForbidden",
        "UnregisteredSecondaryType", "PrimaryVertexFailure", "TopLevelCatch",
        "SamplingFailure",
    ]
    for name in expected:
        assert hasattr(injection.FailureReason, name), (
            f"siren.injection.FailureReason missing member {name!r}")


def test_failure_reason_hashable_as_dict_key():
    """FailureReason members compare and hash correctly as dict keys."""
    d = {injection.FailureReason.NoTargetsOnPath: 1}
    assert injection.FailureReason.NoTargetsOnPath in d
    assert d[injection.FailureReason.NoTargetsOnPath] == 1
    assert injection.FailureReason.NoTargetsOnPath != injection.FailureReason.TopLevelCatch
    assert injection.FailureReason.TopLevelCatch == injection.FailureReason.TopLevelCatch


# --------------------------------------------------------------------------- #
# (b) A fresh injector's ledger starts empty                                   #
# --------------------------------------------------------------------------- #

def _bare_primary_process(max_distance):
    """A data-free DummyCrossSection primary process using
    PointSourcePositionDistribution with the given max_distance. Returns
    (process, keepalive).
    """
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(_NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = _NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), max_distance),
    ]
    return primary_inj, (xs, int_col)


def test_fresh_injector_ledger_is_empty():
    """A freshly constructed injector reports an empty failure ledger
    before any GenerateEvent call."""
    dm = detector.DetectorModel()  # bare: no materials/sectors loaded
    primary_inj, _keepalive = _bare_primary_process(max_distance=25.0)
    rand = utilities.SIREN_random(11)
    inj = injection._Injector(5, dm, primary_inj, [], rand)

    assert inj.GetFailureLedger().entries() == {}


# --------------------------------------------------------------------------- #
# (c) A forced, data-free GenerateEvent failure lands in the ledger under its  #
#     FailureReason, and entries are filterable by reason.                    #
#                                                                              #
# Forcing mechanism: PointSourcePositionDistribution with max_distance=0.0    #
# clips the sampling path to a single point, so                              #
# Path.GetInteractionDepthInBounds(...) == 0 unconditionally (no detector     #
# materials required) and the distribution throws InjectionFailure tagged    #
# FailureReason.NoTargetsOnPath. This needs no detector data files at all --  #
# a bare, unloaded DetectorModel() is enough.                                #
# --------------------------------------------------------------------------- #

def _make_forced_failure_injector(seed=99, n_inject=5):
    dm = detector.DetectorModel()  # bare: no materials/sectors loaded
    primary_inj, keepalive = _bare_primary_process(max_distance=0.0)
    rand = utilities.SIREN_random(seed)
    inj = injection._Injector(n_inject, dm, primary_inj, [], rand)
    return inj, (dm, primary_inj, rand) + keepalive


def test_forced_failure_populates_ledger_with_reason():
    """A forced NoTargetsOnPath failure appears in the ledger keyed by that
    FailureReason, with a positive count and a non-empty exemplar message."""
    inj, _keepalive = _make_forced_failure_injector()

    ev = inj.GenerateEvent()
    assert len(ev.tree) == 0, "forced failure must yield an empty tree"

    entries = inj.GetFailureLedger().entries()
    assert len(entries) == 1
    (depth, parent_pdg, reason), (count, exemplar) = next(iter(entries.items()))
    assert reason == injection.FailureReason.NoTargetsOnPath
    assert depth == 0
    assert parent_pdg == int(_NuMu)
    assert count == 1
    assert inj.InjectionAttempts() == inj.FailedEvents() == 1
    assert inj.InjectedEvents() == 0
    assert isinstance(exemplar, str) and len(exemplar) > 0


def test_forced_failure_repeated_aggregates_count_not_entries():
    """Repeating the same forced failure increments the existing entry's
    count rather than adding new entries (aggregated accounting)."""
    inj, _keepalive = _make_forced_failure_injector(n_inject=4)

    for _ in range(4):
        ev = inj.GenerateEvent()
        assert len(ev.tree) == 0

    entries = inj.GetFailureLedger().entries()
    assert len(entries) == 1
    (_, _, reason), (count, _) = next(iter(entries.items()))
    assert reason == injection.FailureReason.NoTargetsOnPath
    assert count == 4

    # Aggregating over the ledger's entries reproduces the total.
    assert sum(count for count, _ in entries.values()) == 4


def test_ledger_entries_filterable_by_reason():
    """Entries can be filtered down to a single FailureReason, as a caller
    building a per-reason report would do."""
    inj, _keepalive = _make_forced_failure_injector(n_inject=3)

    for _ in range(3):
        inj.GenerateEvent()

    entries = inj.GetFailureLedger().entries()
    matching = {k: v for k, v in entries.items()
                if k[2] == injection.FailureReason.NoTargetsOnPath}
    not_matching = {k: v for k, v in entries.items()
                    if k[2] == injection.FailureReason.TopLevelCatch}

    assert len(matching) == 1
    assert len(not_matching) == 0


def test_ledger_clear_resets_to_empty():
    """FailureLedger.Clear() empties the live ledger held by the injector."""
    inj, _keepalive = _make_forced_failure_injector(n_inject=2)

    inj.GenerateEvent()
    assert inj.GetFailureLedger().entries() != {}

    inj.GetFailureLedger().Clear()
    assert inj.GetFailureLedger().entries() == {}


class _FixedPrimary(distributions.VertexPositionDistribution):
    """Controlled inputs for native sampler failure tests, without detector data."""

    def __init__(self):
        super().__init__()
        self.mass = 2.0
        self.momentum = [10.0, 0.0, 0.0, 10.0]

    def Sample(self, random, detector_model, interactions, record):
        record.mass = self.mass
        record.four_momentum = self.momentum
        record.initial_position = record.interaction_vertex = [0.0, 0.0, 0.0]
        record.helicity = record.initial_time = record.interaction_time = 0.0


class _FixedSecondary(distributions.SecondaryVertexPositionDistribution):
    def SampleVertex(self, random, detector_model, interactions, record):
        record.length = 0.0


class _ChainDecay(interactions.Decay):
    """Test double delivering controlled inputs to a downstream native sampler.

    Intermediate final states deliberately copy the primary fixture; they are
    plumbing fixtures, not a physical decay model. Only the last vertex uses
    the native two-body sampler under test (with daughter masses 0.5, 0.5).
    """

    def __init__(self, primary, secondary, source):
        super().__init__()
        self.source = source
        self.signature = dc.InteractionSignature()
        self.signature.primary_type = primary
        self.signature.target_type = dc.ParticleType.Decay
        self.signature.secondary_types = [secondary, dc.ParticleType.Gamma]

    def equal(self, other):
        return self is other

    def TotalDecayLengthAllFinalStates(self, record):
        return 1.0

    def TotalDecayLength(self, record):
        return 1.0

    def TotalDecayWidthAllFinalStates(self, record):
        return 1.0

    def TotalDecayWidth(self, record):
        return 1.0

    def SampleFinalState(self, record, random):
        for daughter in record.secondary_particle_records:
            daughter.mass = self.source.mass
            daughter.four_momentum = self.source.momentum
            daughter.helicity = 0.0

    def SecondaryMasses(self, types):
        return [0.5] * len(types)

    def GetPossibleSignatures(self):
        return [self.signature]

    def GetPossibleSignaturesFromParent(self, primary):
        return [self.signature] if primary == self.signature.primary_type else []

    def DensityVariables(self):
        return ["cos_theta"]

    def Topology(self):
        return injection.PhaseSpaceTopology.Decay2Body

    def Measure(self):
        return injection.PhaseSpaceMeasure.SolidAngleRest()


@pytest.mark.parametrize("depth", [0, 1, 2], ids=["primary", "secondary", "grandchild"])
@pytest.mark.parametrize("angular_sector", [False, True], ids=["cone", "sector"])
def test_native_sampling_failures_preserve_reason_and_attempt_accounting(depth, angular_sector):
    """Native failures reach the ledger unchanged through every generation path.

    As in the C++ DirectedRejectionExhaustion tests, beta=1 with a finite supplied
    mass deterministically exhausts the inverse solver. A subthreshold mass
    instead has physical zero support; on-shell inputs succeed in the same run.
    """
    source = _FixedPrimary()
    secondary_vertex = _FixedSecondary()
    types = [dc.ParticleType.NuMu, dc.ParticleType.NuE, dc.ParticleType.NuTau,
             dc.ParticleType.Gamma]
    models = [_ChainDecay(types[i], types[i + 1], source) for i in range(depth + 1)]
    processes = []
    for i, model in enumerate(models):
        process = (injection.PrimaryInjectionProcess() if i == 0
                   else injection.SecondaryInjectionProcess())
        process.primary_type = types[i]
        process.interactions = interactions.InteractionCollection(types[i], [model])
        process.distributions = [source if i == 0 else secondary_vertex]
        processes.append(process)

    target = geometry.Sphere(geometry.Placement(smath.Vector3D(0, 0, 100)), 1.0, 0.0)
    channel = (injection.DetectorDirectedAngularSectorChannel(
        target, 0.0, 1.0, 0.0, 2 * math.pi, 0) if angular_sector
        else injection.DetectorDirected2BodyChannel(target, 0, injection.DirectedMode.Cone))
    mixture = injection.MultiChannelPhaseSpace([channel])
    processes[-1].SetPhaseSpace(models[-1].signature, mixture)
    inj = injection._Injector(5, detector.DetectorModel(), processes[0],
                              processes[1:], utilities.SIREN_random(331663))
    inj.SetStoppingCondition(lambda tree, parent, secondary_index: False)
    reason = injection.FailureReason
    # At secondary depths the key identifies the producing parent, which has
    # a different PDG from the particle whose sampler failed.
    parent_pdg = int(types[max(0, depth - 1)])
    sampling_key = (depth, parent_pdg, reason.SamplingFailure)
    forbidden_key = (depth, parent_pdg, reason.KinematicallyForbidden)

    for attempt in (1, 2):
        assert len(inj.GenerateEvent().tree) == 0
        entries = inj.GetFailureLedger().entries()
        assert set(entries) == {sampling_key}
        assert entries[sampling_key][0] == attempt
        assert inj.InjectionAttempts() == inj.FailedEvents() == attempt
        assert inj.InjectedEvents() == 0
        assert len(inj.GetLastFailedTree().tree) == max(1, depth)
        if attempt == 1:
            exemplar = entries[sampling_key][1]
            assert exemplar
        else:
            assert entries[sampling_key][1] == exemplar
        if depth:
            assert f"secondary pdg {int(types[depth])}" in exemplar

    source.mass = 0.9
    source.momentum = [10.0, 0.0, 0.0, math.sqrt(100.0 - source.mass**2)]
    assert len(inj.GenerateEvent().tree) == 0
    entries = inj.GetFailureLedger().entries()
    assert set(entries) == {sampling_key, forbidden_key}
    assert entries[sampling_key] == (2, exemplar)
    assert entries[forbidden_key][0] == 1
    assert inj.GetLastFailureReason()

    source.mass = 2.0
    source.momentum = [10.0, 0.0, 0.0, math.sqrt(96.0)]
    for event_number in (0, 1):
        event = inj.GenerateEvent()
        assert len(event.tree) == depth + 1
        assert event.header.event_number == event_number
    assert inj.GetFailureLedger().entries() == entries
    assert inj.InjectionAttempts() == 5
    assert inj.InjectedEvents() == 2
    assert inj.FailedEvents() == sum(count for count, _ in entries.values()) == 3
    assert inj.InjectionAttempts() == inj.InjectedEvents() + inj.FailedEvents()
    assert inj.UnregisteredSecondaryCount() > 0
    with pytest.raises(RuntimeError, match="maximum number of injection attempts"):
        inj.GenerateEvent()
    assert inj.InjectionAttempts() == 5
    assert inj.GetFailureLedger().entries() == entries

    inj.ResetInjectedEvents()
    assert inj.EventsToInject() == 5
    assert inj.InjectionAttempts() == inj.InjectedEvents() == inj.FailedEvents() == 0
    assert inj.UnregisteredSecondaryCount() == 0
    assert inj.GetFailureLedger().entries() == {}
    assert inj.GetLastFailureReason() == ""
    assert len(inj.GetLastFailedTree().tree) == 0
    assert inj.GenerateEvent().header.event_number == 0

    inj.ResetInjectedEvents(1)
    assert inj.EventsToInject() == 1
    assert inj.InjectionAttempts() == inj.InjectedEvents() == inj.FailedEvents() == 0


# --------------------------------------------------------------------------- #
# (d) The rendered InjectionReport over the ledger.                           #
# --------------------------------------------------------------------------- #

def test_injection_report_renders_table_and_dominant():
    """InjectionReport.from_ledger renders an attrition table and a dominant()."""
    from siren.report import InjectionReport

    inj, _keepalive = _make_forced_failure_injector(n_inject=4)
    for _ in range(4):
        inj.GenerateEvent()

    report = InjectionReport.from_ledger(
        inj.GetFailureLedger(),
        attempts=inj.InjectionAttempts(),
        successes=inj.InjectedEvents(),
        last_failed_tree=inj.GetLastFailedTree())

    text = str(report)
    assert "InjectionReport" in text
    assert "NoTargetsOnPath" in text
    assert "count" in text

    dominant = report.dominant()
    assert dominant is not None
    assert dominant.reason_name == "NoTargetsOnPath"
    assert dominant.count == 4
    assert dominant.hint  # a non-empty one-line remedy


def test_injection_report_efficiency_and_particle_name():
    """The report exposes efficiency and a resolved particle name per bucket."""
    from siren.report import InjectionReport

    inj, _keepalive = _make_forced_failure_injector(n_inject=3)
    for _ in range(3):
        inj.GenerateEvent()

    report = InjectionReport.from_ledger(
        inj.GetFailureLedger(),
        attempts=inj.InjectionAttempts(),
        successes=inj.InjectedEvents())

    assert report.successes == 0
    assert report.efficiency == 0.0
    bucket = report.dominant()
    # The parent PDG (14, NuMu) resolves to a readable particle name.
    assert bucket.particle != str(bucket.pdg) or bucket.pdg == 14
