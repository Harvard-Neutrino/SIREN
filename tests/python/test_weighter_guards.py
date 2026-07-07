"""Weighter fail-loud + realized-count normalization behavior.

Built as Python tests because Weighter_TEST.cxx has no injector fixture and the
C++ CCM_HNL fixture needs detector data files that are absent on some machines.
This uses the same data-free chain the golden-regression and HepMC3 closure
tests use: a DummyCrossSection over the (checked-in) CCM detector model, a
PowerLaw energy spectrum, and a point-source vertex -- so a genuine
injection + weighting round-trip runs without external tables.

Pinned:
  * the event weight normalizes by the REALIZED injected count
    (InjectedEvents()), not the requested count. Weighting the same event at two
    different realized counts M and N scales as w_M * M == w_N * N.

Skipped (not cheaply constructible from Python; see notes in each test). No
C++ test exercises these two guard behaviors either -- Weighter_TEST.cxx covers
only the OneMinusExp/LogOneMinusExp numeric helpers -- so the guard behaviors
are currently unpinned; only the WeightCalculationError *type* is pinned (in
test_errors.py):
  * generation_probability <= 0 -> WeightCalculationError
  * physical_probability == 0 -> weight exactly 0.0 (no raise)

ASCII only. No network. Each test stays well under ~30s.
"""
import math
import os

import pytest

siren = pytest.importorskip("siren")

from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util


_NuMu = siren.dataclasses.Particle.ParticleType.NuMu


# ------------------------------------------------------------------ #
#  Data-free single-vertex chain over the CCM detector                #
# ------------------------------------------------------------------ #

def _load_ccm_detector():
    """Load the CCM detector model or skip if its data files are missing."""
    try:
        det_dir = _util.get_detector_model_path("CCM")
        materials = os.path.join(det_dir, "materials.dat")
        densities = os.path.join(det_dir, "densities.dat")
        if not (os.path.exists(materials) and os.path.exists(densities)):
            pytest.skip("CCM detector data files not present")
        dm = detector.DetectorModel()
        dm.LoadMaterialModel(materials)
        dm.LoadDetectorModel(densities)
        return dm
    except (ValueError, RuntimeError) as exc:
        pytest.skip("CCM detector model unavailable: {}".format(exc))


def _build_chain(detector_model, n_inject, seed):
    """A single-vertex Injector + Weighter over a DummyCrossSection.

    The generation flux is a genuine PowerLaw energy spectrum (so the weight
    depends on energy), plus helicity/direction/vertex distributions. Returns
    (injector, weighter, keepalive) -- keepalive holds Python objects that must
    outlive the C++ engines.
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

    rand = utilities.SIREN_random(seed)
    inj = injection._Injector(n_inject, detector_model, primary_inj, [], rand)
    weighter = injection._Weighter([inj], detector_model, primary_phys, [])

    keepalive = (xs, int_col, primary_inj, primary_phys, rand)
    return inj, weighter, keepalive


def _generate_until(inj, target):
    """Generate until InjectedEvents() reaches `target`; return the last
    non-empty event (a valid tree to weight)."""
    last = None
    guard = 0
    while inj.InjectedEvents() < target:
        ev = inj.GenerateEvent()
        guard += 1
        if len(ev.tree) > 0:
            last = ev
        if guard > 20 * max(target, 1):
            break
    return last


# ------------------------------------------------------------------ #
#  weight normalizes by realized injected count                    #
# ------------------------------------------------------------------ #

def test_weight_scales_with_realized_injected_count():
    """The same event, weighted when the injector has realized M injected events
    vs N, must satisfy w_M * M == w_N * N (the generation seed is the realized
    count InjectedEvents(), so w ~ 1/InjectedEvents and the count-independent
    physical/generation-density product is w * InjectedEvents).
    """
    dm = _load_ccm_detector()
    inj, weighter, _keep = _build_chain(dm, n_inject=200, seed=1234)

    M = 10
    event = _generate_until(inj, M)
    assert event is not None and len(event.tree) > 0, "no event generated"
    iM = inj.InjectedEvents()
    assert iM >= M
    w_M = weighter.EventWeight(event)
    assert math.isfinite(w_M) and w_M > 0.0

    # Generate more, then re-weight the SAME stored event at a larger realized
    # count. EventWeight reads the injector's CURRENT InjectedEvents().
    N = 80
    _generate_until(inj, N)
    iN = inj.InjectedEvents()
    assert iN >= N and iN > iM
    w_N = weighter.EventWeight(event)
    assert math.isfinite(w_N) and w_N > 0.0

    # Contract: w_M * M == w_N * N (both equal the realized-count-independent
    # physical/generation-density product).
    lhs = w_M * iM
    rhs = w_N * iN
    assert math.isclose(lhs, rhs, rel_tol=1e-9), (
        "weight did not scale with realized injected count: "
        "w_M={} (M={}), w_N={} (N={}); w_M*M={} vs w_N*N={} "
        "(if these are equal only because w_M==w_N, the realized-count seed "
        "was not applied)".format(w_M, iM, w_N, iN, lhs, rhs))

    # Guard against count-independent normalization (that would give
    # w_M == w_N, hence w_M*M != w_N*N for M != N).
    assert not math.isclose(w_M, w_N, rel_tol=1e-9), (
        "weight is independent of realized count (w_M == w_N); the "
        "realized-count normalization (InjectedEvents seeding) is not active")


def test_weight_falls_back_to_events_to_inject_before_generation():
    """Before any generation InjectedEvents()==0, so the weight seed falls back
    to EventsToInject(); after generation it uses the realized count. The two
    weights of the same event differ by exactly the count ratio.
    """
    dm = _load_ccm_detector()
    inj, weighter, _keep = _build_chain(dm, n_inject=100, seed=99)

    # Generate a handful, keep one event.
    M = 7
    event = _generate_until(inj, M)
    assert event is not None and len(event.tree) > 0
    realized = inj.InjectedEvents()
    w_realized = weighter.EventWeight(event)

    # Reset the run counters (zero-arg reset preserves the quota) so
    # InjectedEvents()==0 and the seed falls back to EventsToInject().
    inj.ResetInjectedEvents()
    assert inj.InjectedEvents() == 0
    quota = inj.EventsToInject()
    assert quota == 100
    w_fallback = weighter.EventWeight(event)

    assert math.isfinite(w_realized) and w_realized > 0.0
    assert math.isfinite(w_fallback) and w_fallback > 0.0
    # w ~ 1/seed, so w_realized * realized == w_fallback * quota.
    assert math.isclose(w_realized * realized, w_fallback * quota, rel_tol=1e-9), (
        "fallback seed mismatch: w_realized*realized={} vs w_fallback*quota={}"
        .format(w_realized * realized, w_fallback * quota))


# ------------------------------------------------------------------ #
#  Not-cheaply-constructible guard paths (documented skips)            #
# ------------------------------------------------------------------ #

@pytest.mark.skip(
    reason="generation_probability<=0 is not cheaply constructible from Python: "
           "the injector's distributions return positive densities for every "
           "in-support generated event, and doctoring a tree datum's record "
           "does not propagate (InteractionTreeDatum.record returns a copy), so "
           "there is no reliable Python-only way to drive a distribution's "
           "GenerationProbability to 0 for a tree the weighter's bounds "
           "machinery still accepts. This guard behavior is NOT covered by any "
           "C++ test: Weighter_TEST.cxx tests only the OneMinusExp/"
           "LogOneMinusExp numeric helpers. Only the WeightCalculationError "
           "type (existence + RuntimeError subclassing) is pinned, in "
           "test_errors.py; the raise-on-nonpositive-generation behavior is "
           "currently unpinned.")
def test_generation_density_nonpositive_raises_weight_calc_error():
    # Intended contract:
    #   with a tree whose generation density is 0/negative,
    #   weighter.EventWeight(tree) raises siren.utilities.WeightCalculationError
    #   (a RuntimeError subclass), message carrying '[siren-docs: errors#weight-calc]'.
    raise NotImplementedError


@pytest.mark.skip(
    reason="physical_probability==0 (out-of-physical-support -> weight exactly "
           "0.0 WITHOUT raising) is not cheaply constructible from Python: "
           "giving the physical process an energy distribution disjoint from "
           "the injection range did not zero the physical probability in the "
           "single-vertex chain (the physical energy term is not evaluated the "
           "way a disjoint-support trick assumes), and record doctoring does "
           "not propagate through the tree. The guard polarity (physical==0 -> "
           "0.0, no throw; generation<=0 -> throw) is NOT covered by any C++ "
           "test either: Weighter_TEST.cxx tests only the OneMinusExp/"
           "LogOneMinusExp numeric helpers, so this behavior is currently "
           "unpinned.")
def test_physical_density_zero_gives_zero_weight_without_raising():
    # Intended contract:
    #   with a tree outside the physical support (physical_probability == 0),
    #   weighter.EventWeight(tree) == 0.0 exactly and does NOT raise.
    raise NotImplementedError
