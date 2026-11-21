"""Numeric reweighting closure for the HepMC3 round trip (C5).

Generate events with a real injector + weighter (a data-free DummyCrossSection
setup), compute each event weight, write the events to a HepMC3 file, read them
back, recompute the weights from the reconstructed trees, and assert they match.
This is the end-to-end guarantee that the HepMC3 round trip preserves everything
the Weighter reads.

The primary flux uses a power-law energy spectrum and a neutrino helicity
distribution so the recomputed weight actually depends on the round-tripped
energy and helicity (a monoenergetic, helicity-free setup would pass even if
those fields were corrupted). A second test injects a secondary interaction so
the events are multi-vertex cascades, whose reweighting reads a non-root
interaction's primary_initial_position through the secondary vertex-position
distribution.

Skips cleanly when SIREN was built without HepMC3, or when the CCM detector
material files are unavailable in the environment.
"""
import pytest
import numpy as np

siren = pytest.importorskip("siren")


def _build(inject_secondaries):
    """Build (weighter, events) for a data-free DummyCrossSection setup.

    When inject_secondaries is True the injector also simulates one generation of
    secondary interactions, producing multi-vertex cascade events.
    """
    import os
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions
    from siren import detector
    from siren import math as smath
    from siren import utilities
    from siren import _util

    NuMu = dc.Particle.ParticleType.NuMu

    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except (ValueError, RuntimeError) as e:
        # ValueError: get_detector_model_path found no CCM resource folder.
        # RuntimeError: LoadMaterialModel/LoadDetectorModel could not open a
        # materials.dat/densities.dat file. Both mean the CCM data files are
        # genuinely missing from this environment; anything else is a real bug.
        pytest.skip(f"CCM detector model unavailable: {e}")

    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),          # spectrum -> weight depends on energy
        distributions.PrimaryNeutrinoHelicityDistribution(),  # -> weight depends on helicity
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = NuMu
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    n_inject = 40
    rand = utilities.SIREN_random(1234)

    secondary_inj_processes = []
    secondary_phys_processes = []
    if inject_secondaries:
        sec_xs = interactions.DummyCrossSection()
        sec_col = interactions.InteractionCollection(NuMu, [sec_xs])

        sec_inj = injection.SecondaryInjectionProcess()
        sec_inj.primary_type = NuMu
        sec_inj.interactions = sec_col
        sec_inj.distributions = [distributions.SecondaryPhysicalVertexDistribution()]
        secondary_inj_processes = [sec_inj]

        sec_phys = injection.PhysicalProcess()
        sec_phys.primary_type = NuMu
        sec_phys.interactions = sec_col
        sec_phys.distributions = [distributions.SecondaryPhysicalVertexDistribution()]
        secondary_phys_processes = [sec_phys]

        inj = injection._Injector(n_inject, dm, primary_inj, secondary_inj_processes, rand)
        # Simulate exactly one generation of secondaries (stop at depth >= 1).
        inj.SetStoppingCondition(lambda tree, datum, i: datum.depth(tree) >= 1)
        weighter = injection._Weighter([inj], dm, primary_phys, secondary_phys_processes)
    else:
        inj = injection._Injector(n_inject, dm, primary_inj, rand)
        weighter = injection._Weighter([inj], dm, primary_phys)

    events = []
    for _ in range(n_inject):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            # The only RuntimeError GenerateEvent() raises is the max-attempts
            # sentinel ("Injector has already made the maximum number of
            # injection attempts!"); InjectionFailure is caught internally and
            # surfaces as an empty tree, so this cannot mask another failure.
            break
        if len(ev.tree) > 0:
            events.append(ev)
    if not events:
        pytest.skip("no non-empty events generated")

    return weighter, events


def _closure(weighter, events, tmp_path, rel=1e-9):
    from siren import hepmc3 as sio

    weights_before = [weighter.EventWeight(ev) for ev in events]
    for w in weights_before:
        assert np.isfinite(w) and w > 0

    out = str(tmp_path / "closure.hepmc3")
    try:
        sio.SaveInteractionTreesAsHepMC3(events, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise

    loaded = sio.LoadInteractionTreesFromHepMC3(out)
    assert len(loaded) == len(events)

    weights_after = [weighter.EventWeight(ev) for ev in loaded]
    for before, after in zip(weights_before, weights_after):
        assert after == pytest.approx(before, rel=rel), (before, after)
    return loaded


def test_reweighting_numeric_closure(tmp_path):
    """Single-vertex events; weight depends on energy, helicity, geometry."""
    weighter, events = _build(inject_secondaries=False)
    _closure(weighter, events, tmp_path)


def test_reweighting_numeric_closure_cascade(tmp_path):
    """Multi-vertex cascades; reweighting reads a non-root primary_initial_position."""
    weighter, events = _build(inject_secondaries=True)
    loaded = _closure(weighter, events, tmp_path)
    # Ensure the test actually exercised cascades (not just single-vertex events).
    assert any(len(ev.tree) > 1 for ev in loaded), "expected at least one multi-vertex event"


def _build_params():
    """Setup whose physical weight depends on a sampled interaction parameter.

    A small Python cross section makes its final-state probability a function of
    bjorken_x (which DummyCrossSection stores in interaction_parameters). It is
    used as the physical process; generation uses the flat DummyCrossSection so
    the parameter dependence does not cancel in the weight ratio. Reweighting an
    event therefore reads interaction_parameters back off the record.
    """
    import os
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions
    from siren import detector
    from siren import math as smath
    from siren import utilities
    from siren import _util

    NuMu = dc.Particle.ParticleType.NuMu

    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except Exception as e:
        pytest.skip(f"CCM detector model unavailable: {e}")

    class ParamCrossSection(interactions.CrossSection):
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
            params = r.interaction_parameters
            x = params["bjorken_x"] if "bjorken_x" in params else 0.5
            return self._d.TotalCrossSection(r) * (0.5 + x)

        def FinalStateProbability(self, r):
            t = self._d.TotalCrossSection(r)
            return (self.DifferentialCrossSection(r) / t) if t > 0 else 0.0

    gen_xs = interactions.DummyCrossSection()
    phys_xs = ParamCrossSection()
    gen_col = interactions.InteractionCollection(NuMu, [gen_xs])
    phys_col = interactions.InteractionCollection(NuMu, [phys_xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = gen_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = NuMu
    primary_phys.interactions = phys_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.IsotropicDirection(),
    ]

    rand = utilities.SIREN_random(7)
    inj = injection._Injector(40, dm, primary_inj, rand)
    weighter = injection._Weighter([inj], dm, primary_phys)

    events = []
    for _ in range(40):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            # Max-attempts sentinel only; see the identical guard in _build().
            break
        if len(ev.tree) > 0:
            events.append(ev)
    if not events:
        pytest.skip("no non-empty events generated")

    # Return the Python cross sections to keep them alive alongside the weighter.
    return weighter, events, (gen_xs, phys_xs)


def test_reweighting_numeric_closure_interaction_parameters(tmp_path):
    """The physical weight depends on the sampled interaction_parameters, so the
    closure exercises their round trip numerically (not just field-by-field)."""
    weighter, events, _keepalive = _build_params()
    _closure(weighter, events, tmp_path)
