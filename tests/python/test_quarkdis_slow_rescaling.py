"""Slow-rescaling (xi, y) sampling tests for QuarkDISFromSpline.

The charm FITS splines are LHAPDF-derived and not committed, so the whole module
is gated behind SIREN_CHARM_SPLINE_DIR (a directory containing
dsdxidy_nu-N-cc-charm-CT14nlo_central.fits and sigma_nu-N-cc-charm-CT14nlo_central.fits);
it skips cleanly when unset or the files are absent. SIREN_CHARM_NEVENTS (default
2000) sizes the differential-positivity test.
"""
import math
import os

import pytest

siren = pytest.importorskip("siren")

# Constants (mirror C++ siren::utilities::Constants values exactly).
M_C = 1.27                          # Constants::charmMass
M_D0 = 1.86484                      # Constants::D0Mass
M_N = (0.938272 + 0.939565) / 2     # isoscalar nucleon mass
M_MU = 0.105658                     # muon mass
Q2MIN = 1.0                         # GeV^2
E_NU = 100.0                        # neutrino energy in GeV

N_BOUNDS = 100   # cheap kinematic-bounds test
N_DIFF = int(os.environ.get("SIREN_CHARM_NEVENTS", "2000"))   # differential test

# Per-event resample budget on recoverable InjectionFailure (RuntimeError in
# Python); a direct SampleFinalState caller must retry as the injector does.
MAX_RETRIES = 100

# Spline-file gating: resolve spline paths from SIREN_CHARM_SPLINE_DIR.
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


# Expected kinematic bounds (mirror C++ SampleFinalState sampling bounds).
Y_MAX = 1.0 - M_MU / E_NU
W2_THR = (M_N + M_D0) ** 2
Y_MIN = (W2_THR - M_N ** 2 + Q2MIN) / (2.0 * M_N * E_NU)
XI_MIN = (M_C ** 2 + Q2MIN) / (2.0 * M_N * E_NU * Y_MAX)


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
    """Sample one event, retrying on recoverable InjectionFailure (RuntimeError).

    Returns the populated record, or None if the resample budget was exhausted
    (a rare, expected NaN-guard rejection, not a test failure).
    """
    import siren.dataclasses
    ir = _make_neutrino_record(sig)
    for retry in range(MAX_RETRIES + 1):
        cdr = siren.dataclasses.CrossSectionDistributionRecord(ir)
        try:
            charm_xs.SampleFinalState(cdr, rng)
            return cdr
        except RuntimeError:
            # Recoverable per-event failure; resample. Config errors fire in fixtures.
            if retry < MAX_RETRIES:
                continue
            return None
    return None


# Test 1: kinematic bounds on the sampled (xi, y) over 100 events.
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

        # Mirror C++ slow-rescaling relations.
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

    # Most slots must sample; a few NaN-guard rejections are tolerated.
    assert sampled > 0, "no events were sampled at all"
    assert rejected <= N_BOUNDS // 10, (
        f"{rejected}/{N_BOUNDS} events were unrecoverable "
        f"(expected at most {N_BOUNDS // 10})"
    )
    assert not failures, (
        f"{len(failures)} of {sampled} sampled events violated kinematic "
        "bounds:\n" + "\n".join(failures[:10])
    )


# Test 2: differential cross section is finite-positive on the production path.
def test_quarkdis_differential_positive(charm_xs, signature, rng):
    """Re-evaluate the spline on finalized records via the production path.

    Drives the weighting density exactly as the Weighter does: SampleFinalState
    -> cdr.finalize(ir_out) -> DifferentialCrossSection on the finalized record
    (primary-momentum Q2 branch). Asserts a high finite-positive fraction.
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


# Test 3: Sample == Density closure on the finalized record.
def test_quarkdis_sample_density_closure(charm_xs, signature, rng):
    """The two Q2 derivations in DifferentialCrossSection must agree.

    Primary-momentum branch (real secondary momenta) vs stored-(xi, y) fallback
    branch (momenta absent) must yield the same value for a self-consistent event;
    FinalStateProbability (= dxs/txs, the Weighter's density) must be finite >= 0.
    """
    import siren.dataclasses

    checked = 0
    for event_idx in range(N_BOUNDS):
        cdr = _sample_with_retries(charm_xs, signature, rng)
        if cdr is None:
            continue

        params = dict(cdr.interaction_parameters)

        # Finalized record: primary-momentum Q2 branch.
        ir_out = siren.dataclasses.InteractionRecord()
        ir_out.signature = signature
        ir_out.primary_momentum = [E_NU, 0.0, 0.0, E_NU]
        ir_out.primary_mass = 0.0
        cdr.finalize(ir_out)
        dxs_primary = charm_xs.DifferentialCrossSection(ir_out)

        if not (math.isfinite(dxs_primary) and dxs_primary > 0.0):
            # Closure is only meaningful where dxs is positive on both paths.
            continue

        # Same record without secondary momenta: stored-(xi, y) fallback branch.
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

        # Same spline at the same (E, xi, y): agree to within Q2-derivation roundoff.
        rel = abs(dxs_primary - dxs_fallback) / max(dxs_primary, dxs_fallback)
        assert rel < 1e-3, (
            f"event {event_idx}: primary-path dxs={dxs_primary:.6g} disagrees "
            f"with fallback dxs={dxs_fallback:.6g} (rel={rel:.3g}); the two Q2 "
            "derivations in DifferentialCrossSection are inconsistent"
        )

        # Weighter's physical density.
        fsp = charm_xs.FinalStateProbability(ir_out)
        assert math.isfinite(fsp), f"event {event_idx}: FinalStateProbability not finite: {fsp}"
        assert fsp >= 0.0, f"event {event_idx}: FinalStateProbability negative: {fsp}"

        checked += 1

    assert checked > 0, "no positive-dxs events were available for closure check"


# Test 4: absolute charm-DIS normalization (charm fraction vs literature/Pythia).
def test_quarkdis_charm_fraction_normalization():
    """The inclusive slow-rescaling charm-CC cross section must have the right
    absolute magnitude: the charm fraction sigma_charm / sigma_CC must land in the
    literature band (~4-7% over 100 GeV - 1 TeV), cross-validating the independent
    PythiaDISCrossSection charm fraction (~6.5% at 100 GeV). The spline is read in
    cm so TotalCrossSection returns cm^2 (the per-nucleon reference below is cm^2).
    """
    import siren.interactions
    import siren.dataclasses

    PT = siren.dataclasses.Particle.ParticleType
    xs = siren.interactions.QuarkDISFromSpline(
        _DIFF_FILE, _TOTAL_FILE, int(1), M_N, int(1),
        [PT.NuMu], [PT.O16Nucleus], "cm")

    sigma_cc_per_gev = 0.677e-38       # textbook nu-N CC sigma/E [cm^2/GeV]
    prev = None
    n_checked = 0
    for E in (100.0, 200.0, 300.0):
        try:
            s = xs.TotalCrossSection(PT.NuMu, E)
        except RuntimeError:
            continue                   # energy outside the provided spline range
        assert s > 0.0
        frac = s / (sigma_cc_per_gev * E)
        assert 0.02 < frac < 0.10, (
            f"charm fraction {frac:.3f} at {E:.0f} GeV outside the literature band "
            "[0.02, 0.10]")
        if prev is not None:
            assert s > prev, "charm cross section must increase with energy"
        prev = s
        n_checked += 1
    assert n_checked >= 1, "no in-range energy point was available to normalize"
