"""Slow-rescaling (xi, y) sampling tests for QuarkDISFromSpline.

These were previously two orphan smoke scripts under tests/slow_rescaling/
(smoke_quarkdis_100.py, smoke_quarkdis_10k.py) that were never collected by
pytest (testpaths = tests/python) and aborted the session with module-level
sys.exit() calls. They are now real, collectable pytest tests.

The charm slow-rescaling FITS splines are LHAPDF-derived and not committed to
the repository, so every spline-dependent test is gated behind the
SIREN_CHARM_SPLINE_DIR environment variable. Point it at a directory containing

    dsdxidy_nu-N-cc-charm-CT14nlo_central.fits
    sigma_nu-N-cc-charm-CT14nlo_central.fits

(e.g. the FASRC Maboi_M_Muon_SR spline set) to run the physics asserts. When the
variable is unset, or the FITS files are absent, the whole module skips cleanly.

Number of events for the differential-positivity test is configurable via
SIREN_CHARM_NEVENTS (default 2000, kept modest so the suite stays fast; set it
to 10000 to reproduce the original 10k smoke run).
"""
import math
import os

import pytest

siren = pytest.importorskip("siren")

# ---------------------------------------------------------------------------
# Constants (mirror C++ siren::utilities::Constants values exactly)
# ---------------------------------------------------------------------------
M_C = 1.27                          # Constants::charmMass
M_D0 = 1.86484                      # Constants::D0Mass
M_N = (0.938272 + 0.939565) / 2     # isoscalar nucleon mass
M_MU = 0.105658                     # muon mass
Q2MIN = 1.0                         # GeV^2
E_NU = 100.0                        # neutrino energy in GeV

# Number of kinematic-bound events (the cheap test) and re-evaluation events
# (the heavier differential-positivity test).
N_BOUNDS = 100
N_DIFF = int(os.environ.get("SIREN_CHARM_NEVENTS", "2000"))

# Per-event resample budget when the sampler raises a recoverable
# InjectionFailure (surfaces in Python as RuntimeError). The production injector
# framework retries automatically; a direct SampleFinalState caller must do so.
MAX_RETRIES = 100

# ---------------------------------------------------------------------------
# Spline-file gating (resolves the old hardcoded /n/holylfs05 cluster path)
# ---------------------------------------------------------------------------
_SPLINE_DIR = os.environ.get("SIREN_CHARM_SPLINE_DIR")
_DIFF_FILE = (
    os.path.join(_SPLINE_DIR, "dsdxidy_nu-N-cc-charm-CT14nlo_central.fits")
    if _SPLINE_DIR
    else None
)
_TOTAL_FILE = (
    os.path.join(_SPLINE_DIR, "sigma_nu-N-cc-charm-CT14nlo_central.fits")
    if _SPLINE_DIR
    else None
)
_have_splines = bool(
    _SPLINE_DIR
    and _DIFF_FILE
    and _TOTAL_FILE
    and os.path.exists(_DIFF_FILE)
    and os.path.exists(_TOTAL_FILE)
)

pytestmark = pytest.mark.skipif(
    not _have_splines,
    reason=(
        "set SIREN_CHARM_SPLINE_DIR to a directory containing "
        "dsdxidy_nu-N-cc-charm-*.fits and sigma_nu-N-cc-charm-*.fits "
        "to run the charm slow-rescaling tests"
    ),
)


# ---------------------------------------------------------------------------
# Expected kinematic bounds (mirror C++ SampleFinalState sampling bounds)
# ---------------------------------------------------------------------------
Y_MAX = 1.0 - M_MU / E_NU
W2_THR = (M_N + M_D0) ** 2
Y_MIN = (W2_THR - M_N ** 2 + Q2MIN) / (2.0 * M_N * E_NU)
XI_MIN = (M_C ** 2 + Q2MIN) / (2.0 * M_N * E_NU * Y_MAX)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def charm_xs():
    """Build the QuarkDISFromSpline cross section once for the module."""
    import siren.interactions
    import siren.dataclasses

    PT = siren.dataclasses.Particle.ParticleType
    xs = siren.interactions.QuarkDISFromSpline(
        _DIFF_FILE,
        _TOTAL_FILE,
        int(1),                 # interaction_type: CC
        M_N,                    # isoscalar target mass
        int(1),                 # minimum Q2 (GeV^2)
        [PT.NuMu],              # primary types
        [PT.O16Nucleus],        # target types
        "m",                    # units
    )
    return xs


@pytest.fixture(scope="module")
def signature(charm_xs):
    sigs = list(charm_xs.GetPossibleSignatures())
    assert len(sigs) > 0, "QuarkDISFromSpline returned no signatures"
    return sigs[0]


@pytest.fixture
def rng():
    import siren.utilities
    return siren.utilities.SIREN_random(1234)


def _make_neutrino_record(sig):
    """A fresh input InteractionRecord for a 100 GeV nu_mu along +z."""
    import siren.dataclasses
    ir = siren.dataclasses.InteractionRecord()
    ir.signature = sig
    ir.primary_momentum = [E_NU, 0.0, 0.0, E_NU]   # massless nu along +z
    ir.primary_mass = 0.0
    ir.target_mass = M_N
    return ir


def _sample_with_retries(charm_xs, sig, rng):
    """Sample one event, retrying on recoverable InjectionFailure.

    Returns the populated CrossSectionDistributionRecord, or None if the
    per-event resample budget was exhausted (an expected, rare NaN-guard
    rejection, not a test failure). InjectionFailure derives from
    std::runtime_error, so it surfaces in Python as RuntimeError.
    """
    import siren.dataclasses
    ir = _make_neutrino_record(sig)
    for retry in range(MAX_RETRIES + 1):
        cdr = siren.dataclasses.CrossSectionDistributionRecord(ir)
        try:
            charm_xs.SampleFinalState(cdr, rng)
            return cdr
        except RuntimeError:
            # Recoverable per-event sampling failure; resample with fresh
            # RNG draws. Hard configuration errors (bad signature, units, or
            # spline) would already have fired in the fixtures, not here.
            if retry < MAX_RETRIES:
                continue
            return None
    return None


# ---------------------------------------------------------------------------
# Test 1: kinematic bounds on the sampled (xi, y) over 100 events
# ---------------------------------------------------------------------------
def test_quarkdis_kinematic_bounds(charm_xs, signature, rng):
    """Every sampled event must satisfy the slow-rescaling kinematic bounds."""
    assert XI_MIN < 1.0, "degenerate bounds: xi_min >= 1"
    assert Y_MIN < Y_MAX, "degenerate bounds: y_min >= y_max"

    failures = []
    rejected = 0
    sampled = 0

    for event_idx in range(N_BOUNDS):
        cdr = _sample_with_retries(charm_xs, signature, rng)
        if cdr is None:
            rejected += 1
            continue
        sampled += 1

        params = dict(cdr.interaction_parameters)
        xi = params["bjorken_xi"]
        y = params["bjorken_y"]
        x = params["bjorken_x"]

        # Derived quantities (mirror C++ slow-rescaling relations).
        Q2 = 2.0 * M_N * E_NU * y * xi - M_C ** 2
        W2 = M_N ** 2 + 2.0 * M_N * E_NU * y * (1.0 - xi) + M_C ** 2

        event_failures = []
        if not (XI_MIN <= xi <= 1.0):
            event_failures.append(f"xi={xi:.6g} not in [{XI_MIN:.6g}, 1.0]")
        if not (Y_MIN <= y <= Y_MAX):
            event_failures.append(
                f"y={y:.6g} not in [{Y_MIN:.6g}, {Y_MAX:.6g}]"
            )
        if not (x > 0.0):
            event_failures.append(f"x={x:.6g} not > 0")
        if not (Q2 >= Q2MIN):
            event_failures.append(f"Q2={Q2:.6g} < Q2MIN={Q2MIN:.6g}")
        if not (W2 > W2_THR):
            event_failures.append(f"W2={W2:.6g} not > W2_thr={W2_THR:.6g}")

        if event_failures:
            failures.append(
                f"event {event_idx}: xi={xi:.6g} y={y:.6g} x={x:.6g} "
                f"Q2={Q2:.6g} W2={W2:.6g} | " + "; ".join(event_failures)
            )

    # The vast majority of slots must produce a sample; a few NaN-guard
    # rejections are tolerated but a wholesale failure indicates a real bug.
    assert sampled > 0, "no events were sampled at all"
    assert rejected <= N_BOUNDS // 10, (
        f"{rejected}/{N_BOUNDS} events were unrecoverable "
        f"(expected at most {N_BOUNDS // 10})"
    )
    assert not failures, (
        f"{len(failures)} of {sampled} sampled events violated kinematic "
        "bounds:\n" + "\n".join(failures[:10])
    )


# ---------------------------------------------------------------------------
# Test 2: differential cross section is finite-positive on the production path
# ---------------------------------------------------------------------------
def test_quarkdis_differential_positive(charm_xs, signature, rng):
    """Re-evaluate the spline on finalized records via the production path.

    This drives the REAL weighting density: SampleFinalState populates the
    CrossSectionDistributionRecord, cdr.finalize(ir_out) materializes the
    secondary momenta, and DifferentialCrossSection(ir_out) is evaluated on the
    finalized record. That exercises the primary-momentum Q2 branch of
    DifferentialCrossSection -- exactly the path the Weighter runs -- instead
    of the deliberately-broken zero-momenta fallback used by the old smoke
    script. We assert a high finite-positive fraction.
    """
    import siren.dataclasses

    positive = []
    nonpositive = 0
    eval_errors = 0
    rejected = 0
    sampled = 0

    n_secondaries = len(signature.secondary_types)

    for event_idx in range(N_DIFF):
        cdr = _sample_with_retries(charm_xs, signature, rng)
        if cdr is None:
            rejected += 1
            continue
        sampled += 1

        # Production path: materialize the sampled state into an output record.
        ir_out = siren.dataclasses.InteractionRecord()
        ir_out.signature = signature
        ir_out.primary_momentum = [E_NU, 0.0, 0.0, E_NU]
        ir_out.primary_mass = 0.0
        cdr.finalize(ir_out)

        # finalize must populate the secondary momenta from the sampled state.
        assert len(ir_out.secondary_momenta) == n_secondaries, (
            f"event {event_idx}: finalize wrote "
            f"{len(ir_out.secondary_momenta)} secondary momenta, "
            f"expected {n_secondaries}"
        )

        try:
            val = charm_xs.DifferentialCrossSection(ir_out)
        except Exception:
            eval_errors += 1
            continue

        if math.isfinite(val) and val > 0.0:
            positive.append(val)
        else:
            nonpositive += 1

    assert sampled > 0, "no events were sampled at all"
    assert eval_errors == 0, (
        f"{eval_errors} of {sampled} DifferentialCrossSection evaluations "
        "raised on the production (finalized) path"
    )

    positive_fraction = len(positive) / sampled
    assert positive_fraction > 0.95, (
        f"only {positive_fraction:.4f} of {sampled} finalized records had a "
        f"finite-positive DifferentialCrossSection ({nonpositive} non-positive)"
    )

    mean_log_xs = sum(math.log10(v) for v in positive) / len(positive)
    assert math.isfinite(mean_log_xs), (
        f"mean log10 differential cross section is not finite: {mean_log_xs}"
    )


# ---------------------------------------------------------------------------
# Test 3: Sample == Density closure on the finalized record
# ---------------------------------------------------------------------------
def test_quarkdis_sample_density_closure(charm_xs, signature, rng):
    """The two Q2 derivations in DifferentialCrossSection must agree.

    DifferentialCrossSection(InteractionRecord) takes the primary-momentum Q2
    branch when the finalized record carries real secondary momenta, and the
    stored-(xi, y) fallback branch when momenta are absent. Both must yield the
    same differential value for a self-consistently sampled event; FinalState-
    Probability (= dxs/txs, the physical density the Weighter uses) must be
    finite and non-negative.
    """
    import siren.dataclasses

    checked = 0
    for event_idx in range(N_BOUNDS):
        cdr = _sample_with_retries(charm_xs, signature, rng)
        if cdr is None:
            continue

        params = dict(cdr.interaction_parameters)

        # Finalized record -> primary-momentum Q2 branch.
        ir_out = siren.dataclasses.InteractionRecord()
        ir_out.signature = signature
        ir_out.primary_momentum = [E_NU, 0.0, 0.0, E_NU]
        ir_out.primary_mass = 0.0
        cdr.finalize(ir_out)
        dxs_primary = charm_xs.DifferentialCrossSection(ir_out)

        if not (math.isfinite(dxs_primary) and dxs_primary > 0.0):
            # Skip rare edge points where the spline returns 0; closure is only
            # meaningful where the differential value is positive on both paths.
            continue

        # Same record stripped of secondary momenta -> stored-(xi, y) fallback.
        ir_fb = siren.dataclasses.InteractionRecord()
        ir_fb.signature = signature
        ir_fb.primary_momentum = [E_NU, 0.0, 0.0, E_NU]
        ir_fb.primary_mass = 0.0
        ir_fb.target_mass = M_N
        ir_fb.secondary_momenta = [[0.0, 0.0, 0.0, 0.0]] * len(signature.secondary_types)
        ir_fb.secondary_masses = [0.0] * len(signature.secondary_types)
        ir_fb.interaction_parameters = {
            "energy": E_NU,
            "bjorken_xi": params["bjorken_xi"],
            "bjorken_y": params["bjorken_y"],
        }
        dxs_fallback = charm_xs.DifferentialCrossSection(ir_fb)

        # Both branches read the same spline at the same (E, xi, y); they should
        # agree to within float roundoff in the Q2 derivation.
        rel = abs(dxs_primary - dxs_fallback) / max(dxs_primary, dxs_fallback)
        assert rel < 1e-3, (
            f"event {event_idx}: primary-path dxs={dxs_primary:.6g} disagrees "
            f"with fallback dxs={dxs_fallback:.6g} (rel={rel:.3g}); the two Q2 "
            "derivations in DifferentialCrossSection are inconsistent"
        )

        # Physical density used by the Weighter.
        fsp = charm_xs.FinalStateProbability(ir_out)
        assert math.isfinite(fsp), f"event {event_idx}: FinalStateProbability not finite: {fsp}"
        assert fsp >= 0.0, f"event {event_idx}: FinalStateProbability negative: {fsp}"

        checked += 1

    assert checked > 0, "no positive-dxs events were available for closure check"
