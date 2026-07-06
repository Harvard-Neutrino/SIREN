"""EventWeightWithBreakdown().total equals EventWeight() over real generated events.

Pins that the per-vertex breakdown decomposition is a lossless refactor of the
scalar weight: for a normal (non-degenerate) tree, summing/multiplying the
recorded per-vertex factors must reproduce the same double EventWeight()
returns, to floating-point exactness (both compute from the same factors, so
the tolerance is 1e-12, not a physics-level rtol). Uses the same data-free
DummyCrossSection + CCM-detector multi-vertex assembly as the golden
regression pin, but generates a fresh (unarchived) sample -- this is a
same-run invariant check, not a frozen numeric anchor.
"""
from __future__ import annotations

import math
import os

import pytest

from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util

_NuMu = dc.Particle.ParticleType.NuMu

CHAIN_SEED = 5678
N_EVENTS = 50
MAX_ATTEMPTS = 2000


# --------------------------------------------------------------------------- #
# Assembly (copied from test_golden_regression.py -- importing a test module   #
# is fragile, so this file carries its own minimal copy).                      #
# --------------------------------------------------------------------------- #

def _skip_unless_ccm_data():
    """Skip only when the CCM detector data files are absent.

    Missing files are an environment condition; any error raised while
    PARSING present files is a real regression and must propagate.
    """
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except ValueError as e:
        pytest.skip(f"CCM detector model path unavailable: {e}")
    missing = [p for p in (os.path.join(det_dir, "materials.dat"),
                           os.path.join(det_dir, "densities.dat"))
               if not os.path.exists(p)]
    if missing:
        pytest.skip("CCM detector data missing: " + ", ".join(missing))


def _load_ccm_detector():
    """Load the CCM detector model; errors propagate to the caller."""
    dm = detector.DetectorModel()
    det_dir = _util.get_detector_model_path("CCM")
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    return dm


def _build_chain(detector_model, n_inject, seed):
    """Assemble a real multi-vertex Injector + Weighter over a data-free
    DummyCrossSection. One generation of secondaries is simulated so every
    event is a multi-vertex cascade. Returns (injector, weighter, keepalive)
    where keepalive holds Python objects that must outlive the C++ engines.
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
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = _NuMu
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    sec_xs = interactions.DummyCrossSection()
    sec_col = interactions.InteractionCollection(_NuMu, [sec_xs])

    sec_inj = injection.SecondaryInjectionProcess()
    sec_inj.primary_type = _NuMu
    sec_inj.interactions = sec_col
    sec_inj.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    sec_phys = injection.PhysicalProcess()
    sec_phys.primary_type = _NuMu
    sec_phys.interactions = sec_col
    sec_phys.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    rand = utilities.SIREN_random(seed)
    inj = injection._Injector(n_inject, detector_model, primary_inj, [sec_inj], rand)
    # Simulate exactly one generation of secondaries (stop at depth >= 1).
    inj.SetStoppingCondition(lambda tree, datum, i: datum.depth(tree) >= 1)
    weighter = injection._Weighter([inj], detector_model, primary_phys, [sec_phys])

    keepalive = (xs, int_col, sec_xs, sec_col, primary_inj, primary_phys,
                 sec_inj, sec_phys, rand)
    return inj, weighter, keepalive


def _generate_events(inj, n_events, max_attempts):
    """Generate n_events accepted (non-empty) events."""
    events = []
    for _ in range(max_attempts):
        if len(events) >= n_events:
            break
        try:
            ev = inj.GenerateEvent()
        except RuntimeError as err:
            if "maximum number of injection attempts" not in str(err):
                raise
            break
        if len(ev.tree) == 0:
            continue
        events.append(ev)
    if len(events) < n_events:
        raise RuntimeError(
            f"chain produced only {len(events)} accepted events "
            f"(< n_events={n_events}) within {max_attempts} attempts")
    return events


# --------------------------------------------------------------------------- #
# Tests                                                                       #
# --------------------------------------------------------------------------- #

def test_breakdown_total_matches_event_weight():
    """EventWeightWithBreakdown(ev).total == EventWeight(ev) to 1e-12 over a
    real multi-vertex generated sample."""
    _skip_unless_ccm_data()
    detector_model = _load_ccm_detector()
    inj, weighter, _keepalive = _build_chain(detector_model, MAX_ATTEMPTS, CHAIN_SEED)

    # Weight AFTER generation completes: EventWeight normalizes by the
    # realized injected count, so all events must share the same (final)
    # normalization -- weighting mid-generation would bake the running count
    # into each weight.
    events = _generate_events(inj, N_EVENTS, MAX_ATTEMPTS)

    for ev in events:
        bd = weighter.EventWeightWithBreakdown(ev)
        w = weighter.EventWeight(ev)
        assert math.isclose(bd.total, w, rel_tol=1e-12, abs_tol=0.0), (
            f"breakdown total {bd.total!r} != EventWeight {w!r}")


def test_breakdown_vertices_cover_the_tree_with_no_flags():
    """Each vertex in the breakdown carries a finite generation/physical
    factor and no error flag for a normal (non-degenerate) multi-vertex tree."""
    _skip_unless_ccm_data()
    detector_model = _load_ccm_detector()
    inj, weighter, _keepalive = _build_chain(detector_model, MAX_ATTEMPTS, CHAIN_SEED)

    events = _generate_events(inj, N_EVENTS, MAX_ATTEMPTS)

    for ev in events:
        bd = weighter.EventWeightWithBreakdown(ev)
        assert len(bd.vertices) == len(ev.tree)
        assert len(bd.vertices) >= 2, "one secondary generation must yield >=2 vertices"
        for v in bd.vertices:
            assert math.isfinite(v.generation)
            assert math.isfinite(v.physical)
            assert v.generation > 0.0
            assert "generation density zero" not in v.flags
            assert "outside physical support (weight 0)" not in v.flags
