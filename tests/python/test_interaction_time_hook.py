"""Trampoline test for the post-selection vertex-time hook.

A Python CrossSection subclass overrides SampleInteractionTime to shift the
vertex interaction time. Running an injection through the Injector, the override
must land on the root record's interaction_time and propagate to the daughter
production times (secondary_times) via the CrossSectionDistributionRecord
back-sync. A default (unoverridden) cross section leaves the flight-time value
untouched, so the shifted run is observably different from the default run.
"""
import os

import pytest

siren = pytest.importorskip("siren")

SHIFTED_TIME = 987.5


def _detector_model():
    from siren import detector
    from siren import _util

    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except Exception as e:  # detector model files not available in this env
        pytest.skip(f"CCM detector model unavailable: {e}")
    return dm


def _make_cross_section(interactions, override_time):
    """A CrossSection wrapping DummyCrossSection that optionally overrides the
    vertex time via the SampleInteractionTime hook."""

    class TimeShiftCrossSection(interactions.CrossSection):
        def __init__(self):
            super().__init__()
            self._d = interactions.DummyCrossSection()

        def equal(self, other):
            return False

        def GetPossiblePrimaries(self):
            return self._d.GetPossiblePrimaries()

        def GetPossibleTargets(self):
            return self._d.GetPossibleTargets()

        def GetPossibleTargetsFromPrimary(self, p):
            return self._d.GetPossibleTargetsFromPrimary(p)

        def GetPossibleSignatures(self):
            return self._d.GetPossibleSignatures()

        def GetPossibleSignaturesFromParents(self, p, t):
            return self._d.GetPossibleSignaturesFromParents(p, t)

        def InteractionThreshold(self, r):
            return self._d.InteractionThreshold(r)

        def DensityVariables(self):
            return self._d.DensityVariables()

        def SampleFinalState(self, r, rand):
            self._d.SampleFinalState(r, rand)

        def TotalCrossSection(self, r):
            return self._d.TotalCrossSection(r)

        def TotalCrossSectionAllFinalStates(self, r):
            return self._d.TotalCrossSection(r)

        def DifferentialCrossSection(self, r):
            return self._d.TotalCrossSection(r)

        def FinalStateProbability(self, r):
            return 1.0

        def SampleInteractionTime(self, r, rand):
            if override_time:
                return SHIFTED_TIME
            # Defer to the base identity default.
            return super().SampleInteractionTime(r, rand)

    return TimeShiftCrossSection()


def _generate(override_time):
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions
    from siren import math as smath
    from siren import utilities

    NuMu = dc.Particle.ParticleType.NuMu
    dm = _detector_model()

    xs = _make_cross_section(interactions, override_time)
    int_col = interactions.InteractionCollection(NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    rand = utilities.SIREN_random(42)
    inj = injection._Injector(40, dm, primary_inj, rand)

    events = []
    for _ in range(40):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            break
        if len(ev.tree) > 0:
            events.append(ev)
    if not events:
        pytest.skip("no non-empty events generated")
    # Keep the Python cross section alive alongside the events.
    return events, xs


def test_cross_section_time_hook_shifts_vertex_and_daughter_times():
    events, _keepalive = _generate(override_time=True)
    for ev in events:
        root = ev.tree[0].record
        assert root.interaction_time == pytest.approx(SHIFTED_TIME)
        # The daughters inherit the overridden vertex time as their production
        # time through the CrossSectionDistributionRecord back-sync.
        assert len(root.secondary_times) > 0
        for t in root.secondary_times:
            assert t == pytest.approx(SHIFTED_TIME)


def test_default_time_hook_leaves_flight_time():
    events, _keepalive = _generate(override_time=False)
    # The default hook is the identity, so no vertex is pinned to SHIFTED_TIME.
    for ev in events:
        root = ev.tree[0].record
        assert root.interaction_time != pytest.approx(SHIFTED_TIME)
