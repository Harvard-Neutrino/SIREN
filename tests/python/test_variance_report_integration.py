"""End-to-end variance_report()/ESS through a real biased-decay Injector chain.

The Results contract tests (test_results_contract.py) drive variance_report()
only through stub engines with pre-seeded accumulators; no test drives it
through a genuinely assembled Injector + Weighter with a MultiChannelPhaseSpace
mixing a physical channel and a directed biased channel. This test closes that
gap: it builds a data-free 2-body decay (a minimal siren.DecayModel, no
DarkNews and no external tables), registers a [physical decay @ 0.1, directed
@ 0.9] mixture on the primary injection process, injects 200 events through the
real siren.Injector wrapper, weights them through the real siren.Weighter
wrapper, wraps them in a Results, and exercises the report surface end to end.

The chain is fixed-seed and deterministic (identical seed -> identical weights),
uses the checked-in CCM detector model (skipped when its data files are absent),
and runs in well under a second.

Asserted invariants:
  (a) variance_report() runs and every reported W and share is finite;
  (b) 0 < ess() <= N (a real effective sample size);
  (c) each vertex mixture's per-channel shares reconcile to the report's own
      contract: shares == W_i / sum_j W_j and sum(shares) == 1 per vertex
      (the within-vertex variance-fraction decomposition VarianceReport
      promises), and the mixture is non-degenerate (a finite-weight event
      reached it);
  (d) every event weight is finite and non-negative.

ASCII only. No network. No external tables.
"""
import math
import os

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import injection
from siren import distributions
from siren import detector
from siren import math as smath
from siren import _util
from siren.Injector import Injector
from siren.Weighter import Weighter
from siren.Results import Results


P = siren.dataclasses.ParticleType

# Fixed configuration -- the seed and sizes are part of the determinism check.
SEED = 1234
N_EVENTS = 200
M_N4 = 0.02                 # parent (N4) mass [GeV]
DIRECTED_FRACTION = 0.9     # weight on the directed biased channel


class _IsoDecay(siren.DecayModel):
    """Minimal data-free 2-body decay N4 -> NuLight Gamma on the rest-frame
    solid-angle measure (the closure-by-construction default sampler)."""

    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)

    def density_variables(self):
        return "cost"


def _load_ccm_detector():
    """Load the CCM detector model or skip if its data files are missing."""
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except (ValueError, RuntimeError) as exc:
        pytest.skip("CCM detector model path unavailable: {}".format(exc))
    materials = os.path.join(det_dir, "materials.dat")
    densities = os.path.join(det_dir, "densities.dat")
    if not (os.path.exists(materials) and os.path.exists(densities)):
        pytest.skip("CCM detector data files not present")
    dm = detector.DetectorModel()
    dm.LoadMaterialModel(materials)
    dm.LoadDetectorModel(densities)
    return dm


def _build_results(detector_model, seed=SEED, n_events=N_EVENTS):
    """Assemble the real biased-decay chain and return its weighted Results.

    A [physical decay, directed 2-body] mixture is registered on the primary
    injection process so the injector samples through the biased multi-channel
    phase space and the weighter reweights back to the physical density.
    """
    decay = _IsoDecay()
    sig = decay.GetPossibleSignatures()[0]

    box = siren.geometry.Box(
        siren.geometry.Placement(smath.Vector3D(0.0, 0.0, 8.0)),
        6.0, 6.0, 6.0)
    mixture = injection.MultiChannelPhaseSpace()
    mixture.channels = [
        injection.PhysicalDecayChannel(decay, sig),
        injection.DetectorDirected2BodyChannel(box, 0),
    ]
    mixture.weights = [1.0 - DIRECTED_FRACTION, DIRECTED_FRACTION]

    inj = Injector(
        number_of_events=n_events,
        detector_model=detector_model,
        seed=seed,
        primary_type=P.N4,
        primary_interactions=[decay],
        primary_injection_distributions=[
            distributions.PrimaryMass(M_N4),
            distributions.PowerLaw(2.0, M_N4 * 1.05, 0.05),
            distributions.IsotropicDirection(),
            distributions.PointSourcePositionDistribution(
                smath.Vector3D(0.0, 0.0, 0.0), 25.0),
        ],
        primary_phase_spaces={sig: mixture},
    )
    trees = inj.generate(n_events, on_shortfall="raise")

    weighter = Weighter(
        inj,
        primary_physical=[
            distributions.PrimaryMass(M_N4),
            distributions.IsotropicDirection(),
        ],
    )
    weights = np.array([weighter(tree) for tree in trees], dtype=float)
    gen_times = [0.0] * len(trees)
    return Results(trees, weights, gen_times, weighter, inj)


def test_chain_registers_a_biased_multichannel_mixture():
    """The assembled injector exposes exactly one multi-channel decay mixture
    (physical + directed) to the variance machinery."""
    dm = _load_ccm_detector()
    results = _build_results(dm)

    mixtures = results._injector.engine.GetPhaseSpaces()
    assert len(mixtures) == 1
    labels = [ch.Name() for ch in mixtures[0].channels]
    assert any(name.startswith("Physical") for name in labels)
    assert any(name.startswith("DetectorDirected") for name in labels)


def test_variance_report_runs_and_reconciles():
    """variance_report() runs end to end and its per-vertex shares reconcile
    to W_i / sum_j W_j summing to 1 (the report's own decomposition contract)."""
    dm = _load_ccm_detector()
    results = _build_results(dm)

    report = results.variance_report()

    # (a) The report runs and yields the one decay mixture with finite numbers.
    assert len(report.vertices) == 1
    vertex = report.vertices[0]
    assert len(vertex.W) == 2 and len(vertex.shares) == 2
    assert all(math.isfinite(w) for w in vertex.W)
    assert all(math.isfinite(s) for s in vertex.shares)

    # (c) A finite-weight event reached the mixture, so it is non-degenerate...
    assert vertex.degenerate is False
    assert vertex.count == N_EVENTS
    # ...and the shares reconcile to the report's stated contract:
    # shares_i == W_i / sum_j W_j, summing to 1 within the vertex.
    total_W = sum(vertex.W)
    assert total_W > 0.0
    expected_shares = [w / total_W for w in vertex.W]
    assert vertex.shares == pytest.approx(expected_shares, abs=1e-12)
    assert sum(vertex.shares) == pytest.approx(1.0, abs=1e-12)
    assert all(s >= 0.0 for s in vertex.shares)


def test_ess_is_a_real_effective_sample_size():
    """(b) ess() returns a finite effective sample size in (0, N] and equals
    the closed-form (sum w)^2 / sum(w^2) over the run's own finite weights."""
    dm = _load_ccm_detector()
    results = _build_results(dm)

    ess = results.ess()
    assert math.isfinite(ess)
    assert 0.0 < ess <= len(results)

    # ess() is defined over the finite (non-inf, non-nan) entries of the run's
    # own per-event weights -- exactly the population results.weights exposes.
    w = np.asarray(results.weights)
    finite = w[np.isfinite(w)]
    ess_expected = float(np.sum(finite)) ** 2 / float(np.sum(finite * finite))
    assert ess == pytest.approx(ess_expected, rel=1e-9)


def test_weights_are_finite_and_nonnegative():
    """(d) Every reweighted event weight is finite and non-negative."""
    dm = _load_ccm_detector()
    results = _build_results(dm)

    weights = np.asarray(results.weights)
    assert weights.shape == (N_EVENTS,)
    assert np.all(np.isfinite(weights))
    assert np.all(weights >= 0.0)


def test_chain_is_deterministic():
    """The fixed seed reproduces the same weights and ESS exactly."""
    dm = _load_ccm_detector()
    r1 = _build_results(dm)
    r2 = _build_results(dm)

    np.testing.assert_array_equal(np.asarray(r1.weights), np.asarray(r2.weights))
    assert r1.ess() == r2.ess()
