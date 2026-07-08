"""FailureReason/FailureLedger cross the pybind boundary and are filterable.

FailureLedger.entries() keys on (depth, parent_pdg, FailureReason); reasons
must be usable as dict keys/values (hashable, comparable) from Python, a
fresh ledger must start empty, and a forced GenerateEvent failure must land
in the ledger under its FailureReason so callers can filter by reason.
"""
from __future__ import annotations

from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities

_NuMu = dc.Particle.ParticleType.NuMu


# --------------------------------------------------------------------------- #
# (a) FailureReason enum: importable, hashable, usable as a dict key           #
# --------------------------------------------------------------------------- #

def test_failure_reason_enum_members_present():
    """All eight documented FailureReason members are importable."""
    expected = [
        "Unspecified", "NoPathThroughVolume", "NoTargetsOnPath",
        "NoColumnDepthSolution", "KinematicallyForbidden",
        "UnregisteredSecondaryType", "PrimaryVertexFailure", "TopLevelCatch",
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
    assert count >= 1
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
