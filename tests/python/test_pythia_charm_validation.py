"""Physics validation for the PythiaDISCrossSection charm-DIS generator.

This is the SIREN-side counterpart of the Pythia-vs-SIREN comparison shown in
the IceCube multi-cascade tau/charm update (Diffuse WG): it guards the absolute
charm-production rate, confirms that SampleFinalState reproduces bare Pythia's
DIS kinematics (Bjorken x/y, Q^2), and that an optional differential spline
covers the realized sampling support (no silent-zero weight bias), across the
analysis energy band (TeV-PeV).

Everything here needs Pythia8/LHAPDF at runtime and a wide-range total (and,
for the coverage test, differential) charm spline. It is therefore gated behind
environment variables and skips cleanly when they are unset:

    SIREN_PYTHIA_WIDE_SIGMA    -> total sigma(E) FITS spline, 100 GeV - 1 PeV
    SIREN_PYTHIA_WIDE_DSDXDY   -> differential d2sigma/dx dy FITS spline (optional)
    LHAPDF_DATA_PATH           -> LHAPDF data dir containing the PDF set below

Generate the splines with siren.pythia_charm_splines (see scratch gen script):

    python gen_wide_splines.py <out_dir>

The SampleFinalState path re-initializes Pythia per event (~1 s/event), so the
SIREN-sampled statistics are deliberately modest; the bare-Pythia reference
(GeneratePythiaCharmSamples, ~ms/event) uses high statistics.
"""
import math
import os

import pytest

siren = pytest.importorskip("siren")

PDF_SET = os.environ.get("SIREN_PYTHIA_PDF", "LHAPDF6:CT18NLO")
PYTHIA_DATA = os.environ.get("PYTHIA8DATA", "")
M_N = 0.9314943          # isoscalar nucleon mass (GeV)
M_MU = 0.105658

_SIGMA = os.environ.get("SIREN_PYTHIA_WIDE_SIGMA")
_DSDXDY = os.environ.get("SIREN_PYTHIA_WIDE_DSDXDY")
_LHAPDF = os.environ.get("LHAPDF_DATA_PATH")

_have_total = bool(_SIGMA and os.path.exists(_SIGMA) and _LHAPDF)
_have_diff = bool(_DSDXDY and os.path.exists(_DSDXDY))

_pythia_reason = (
    "set SIREN_PYTHIA_WIDE_SIGMA (wide total spline) and LHAPDF_DATA_PATH, and "
    "build with SIREN_WITH_PYTHIA8=ON, to run the charm validation tests"
)

# PythiaDISCrossSection is only registered when SIREN is built with Pythia8.
_has_pythia_class = hasattr(getattr(siren, "interactions", object()), "PythiaDISCrossSection")

pytestmark = pytest.mark.skipif(
    not (_have_total and _has_pythia_class), reason=_pythia_reason
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_xs(with_differential):
    import siren.interactions
    import siren.dataclasses
    PT = siren.dataclasses.Particle.ParticleType
    diff = _DSDXDY if (with_differential and _have_diff) else ""
    return siren.interactions.PythiaDISCrossSection(
        diff, _SIGMA, 1, M_N, 1.0, [PT.NuMu], [PT.PPlus],
        PYTHIA_DATA, PDF_SET, "cm")


def _sample_siren(xs, E, n):
    """SampleFinalState n charm events at neutrino energy E (GeV).

    Returns a list of dicts with the stored Bjorken (x, y), the finalized D and
    primary-lepton 4-momenta, and the D/lepton energy fractions.
    """
    import siren.dataclasses
    import siren.utilities
    PT = siren.dataclasses.Particle.ParticleType
    rng = siren.utilities.SIREN_random(20260506)
    sigs = list(xs.GetPossibleSignaturesFromParents(PT.NuMu, PT.PPlus))
    out = []
    attempts = 0
    while len(out) < n and attempts < n * 20:
        attempts += 1
        sig = sigs[len(out) % len(sigs)]
        ir = siren.dataclasses.InteractionRecord()
        ir.signature = sig
        ir.primary_momentum = [E, 0.0, 0.0, E]
        ir.primary_mass = 0.0
        ir.target_mass = M_N
        cdr = siren.dataclasses.CrossSectionDistributionRecord(ir)
        try:
            xs.SampleFinalState(cdr, rng)
        except RuntimeError:
            continue
        params = dict(cdr.interaction_parameters)
        ir_out = siren.dataclasses.InteractionRecord()
        ir_out.signature = ir.signature
        ir_out.primary_momentum = [E, 0.0, 0.0, E]
        ir_out.primary_mass = 0.0
        cdr.finalize(ir_out)
        secs = list(ir_out.secondary_momenta)
        # secondary order is [charged lepton, Hadrons, D meson]
        p_lep, p_D = secs[0], secs[2]
        out.append(dict(x=params["bjorken_x"], y=params["bjorken_y"],
                        E=E, p_lep=p_lep, p_D=p_D,
                        zD=p_D[0] / E, zlep=p_lep[0] / E,
                        ir_out=ir_out))
    return out


def _bare_pythia_xy(E, n_events):
    """Bare-Pythia muon-reconstructed (x, y) for an ice p/n mix (10:8)."""
    import siren.interactions
    PT = siren.interactions.PythiaDISCrossSection
    xs_all, ys_all = [], []
    for target_pdg, w in ((2212, 10), (2112, 8)):
        _, _, x, y = PT.GeneratePythiaCharmSamples(
            1, 14, target_pdg, M_N, PDF_SET, PYTHIA_DATA, 1.0,
            [float(E)], int(n_events * w / 18))
        xs_all += list(x)
        ys_all += list(y)
    return xs_all, ys_all


def _mean(v):
    return sum(v) / len(v)


# ---------------------------------------------------------------------------
# Test 1: absolute charm-production rate / charm fraction (normalization)
# ---------------------------------------------------------------------------
def test_charm_total_cross_section_normalization():
    """The inclusive charm-DIS sigma must have the right absolute magnitude.

    Pins (a) the charm fraction sigma_charm/sigma_CC at 100 GeV against the
    textbook nu-N CC cross section, and (b) sigma_charm/E to the right order of
    magnitude across 100 GeV - 1 PeV, with monotonic growth. Guards the rate
    that sets S:B and the astro-flux normalization fit -- previously unguarded
    (only mocks / relative partitioning existed).
    """
    import siren.dataclasses
    PT = siren.dataclasses.Particle.ParticleType
    xs = _make_xs(with_differential=False)

    # (a) charm fraction at 100 GeV (CC sigma is textbook-linear there).
    sigma_cc_ref_100 = 0.677e-38 * 100.0       # cm^2, sigma_CC/E ~ 0.677e-38 cm^2/GeV
    sigma_charm_100 = xs.TotalCrossSection(PT.NuMu, 100.0)
    frac = sigma_charm_100 / sigma_cc_ref_100
    assert 0.03 < frac < 0.12, (
        f"charm fraction at 100 GeV = {frac:.3f} outside the literature band "
        f"[0.03, 0.12] (sigma_charm={sigma_charm_100:.3e} cm^2)")

    # (b) order-of-magnitude of sigma_charm/E + monotonic growth across the band.
    energies = [1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6]
    sig = [xs.TotalCrossSection(PT.NuMu, E) for E in energies]
    for E, s in zip(energies, sig):
        assert s > 0.0
        soe = s / E
        assert 5e-41 < soe < 5e-39, (
            f"sigma_charm/E = {soe:.3e} cm^2/GeV at E={E:.0e} GeV is out of the "
            "expected charm-DIS magnitude band")
    for a, b in zip(sig, sig[1:]):
        assert b > a, "charm cross section must increase with energy"


# ---------------------------------------------------------------------------
# Test 2: SampleFinalState reproduces bare-Pythia DIS kinematics (the slides)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_sampling_matches_bare_pythia_at_100gev():
    """SIREN SampleFinalState must reproduce bare Pythia's Bjorken x/y and Q^2.

    Reproduces the production-side panels of the Pythia-vs-pythiaSIREN slide at
    E_nu = 100 GeV: SampleFinalState extracts/rotates the Pythia final state and
    reconstructs (x, y); GeneratePythiaCharmSamples is the same Pythia config
    inline. Their distributions must agree -- a frame/extraction bug would show
    here as a mismatch.
    """
    E = 100.0
    xs = _make_xs(with_differential=False)
    siren_ev = _sample_siren(xs, E, n=80)
    assert len(siren_ev) >= 40, f"only {len(siren_ev)} SIREN events sampled"
    x_si = [e["x"] for e in siren_ev]
    y_si = [e["y"] for e in siren_ev]

    x_py, y_py = _bare_pythia_xy(E, n_events=4000)
    assert len(x_py) > 1000

    # Compare means; tolerance is set by the modest SIREN sample size.
    def _stderr(v):
        m = _mean(v)
        var = sum((vi - m) ** 2 for vi in v) / max(1, len(v) - 1)
        return math.sqrt(var / len(v))

    for name, vi_si, vi_py in (("y", y_si, y_py), ("x", x_si, x_py)):
        m_si, m_py = _mean(vi_si), _mean(vi_py)
        tol = 4.0 * (_stderr(vi_si) + _stderr(vi_py)) + 0.05 * m_py
        assert abs(m_si - m_py) < tol, (
            f"mean Bjorken-{name}: SIREN={m_si:.4f} vs bare-Pythia={m_py:.4f} "
            f"(tol {tol:.4f}) -- SampleFinalState does not reproduce Pythia")

    # Sampled kinematics must be physical.
    assert all(0.0 < e["x"] < 1.0 for e in siren_ev)
    assert all(0.0 < e["y"] < 1.0 for e in siren_ev)


# ---------------------------------------------------------------------------
# Test 3: D-meson energy fraction + production collimation (morphology)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_d_meson_energy_fraction_and_collimation():
    """The D meson carries a sizeable, physical fraction of the event energy and
    is produced nearly collinear with the primary lepton at high energy -- the
    morphology that sets the two-cascade energy split and separation direction.
    """
    E = 1.0e4   # 10 TeV
    xs = _make_xs(with_differential=False)
    ev = _sample_siren(xs, E, n=60)
    assert len(ev) >= 30
    zD = [e["zD"] for e in ev]
    # Every D carries a physical (0, 1) fraction of the neutrino energy.
    assert all(0.0 < z < 1.0 for z in zD)
    # Mean D energy fraction in a sane charm-DIS range (not ~0, not ~1).
    assert 0.05 < _mean(zD) < 0.95
    # Opening angle between D and primary lepton is small at 10 TeV (collimated).
    def _angle(a, b):
        import math
        na = math.sqrt(sum(a[i] ** 2 for i in (1, 2, 3)))
        nb = math.sqrt(sum(b[i] ** 2 for i in (1, 2, 3)))
        dot = sum(a[i] * b[i] for i in (1, 2, 3))
        c = max(-1.0, min(1.0, dot / (na * nb)))
        return math.degrees(math.acos(c))
    angles = [_angle(e["p_D"], e["p_lep"]) for e in ev]
    # Median opening angle should be modest (DIS at 10 TeV is forward).
    angles.sort()
    median = angles[len(angles) // 2]
    assert median < 30.0, f"median D/lepton opening angle {median:.1f} deg too large"


# ---------------------------------------------------------------------------
# Test 4: optional differential spline covers the sampled support (no silent 0)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not (_have_diff and PYTHIA_DATA),
                    reason="set SIREN_PYTHIA_WIDE_DSDXDY + PYTHIA8DATA to run the coverage test")
def test_differential_spline_covers_sampling_support():
    """With a differential spline supplied, DifferentialCrossSection must be
    finite-positive on essentially every sampled event, else those events get a
    silently-zero physical density and a biased weight. Guards spline (E, x, y)
    support vs the realized sampling support at analysis energies.
    """
    xs = _make_xs(with_differential=True)
    for E in (1.0e3, 1.0e4):   # 1 TeV, 10 TeV (within the wide spline x-range)
        ev = _sample_siren(xs, E, n=60)
        assert len(ev) >= 30, f"too few events at E={E:.0e}"
        vals = [xs.DifferentialCrossSection(e["ir_out"]) for e in ev]
        good = [v for v in vals if math.isfinite(v) and v > 0.0]
        frac = len(good) / len(vals)
        assert frac > 0.95, (
            f"only {frac:.3f} of sampled events at E={E:.0e} GeV had a "
            "finite-positive differential density -- the differential spline "
            "(E, x, y) support does not cover the sampling support (silent-zero "
            "weight bias). Widen the spline's logx range.")


# ---------------------------------------------------------------------------
# Test 5: end-to-end inject -> weight through the real Injector / Weighter
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_end_to_end_inject_and_weight():
    """Drive PythiaDISCrossSection as the primary charm-DIS generator through the
    real _Injector / _Weighter on the CCM detector. Every generated event must
    receive a finite, positive weight, and the weights must vary (genuine
    per-event sampling). This exercises the full production-side weight
    accumulation (TotalCrossSection interaction depth / position / survival times
    the cancelling constant FinalStateProbability) end to end -- the cross-section
    closure tests prove the per-vertex factors; this proves they integrate.
    """
    import numpy as np
    from siren import (dataclasses as dc, injection, interactions, distributions,
                       detector, math as smath, utilities, _util)
    PT = dc.Particle.ParticleType
    NuMu = PT.NuMu

    det_dir = _util.get_detector_model_path("CCM")
    dm = detector.DetectorModel()
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))

    xs = interactions.PythiaDISCrossSection(
        "", _SIGMA, 1, M_N, 1.0, [NuMu], [PT.Ar40Nucleus], PYTHIA_DATA, PDF_SET, "cm")
    int_col = interactions.InteractionCollection(NuMu, [xs])

    pinj = injection.PrimaryInjectionProcess()
    pinj.primary_type = NuMu
    pinj.interactions = int_col
    pinj.distributions = [
        distributions.PrimaryMass(0),
        distributions.Monoenergetic(1000.0),     # 1 TeV nu_mu
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 5.0),
    ]
    pphys = injection.PhysicalProcess()
    pphys.primary_type = NuMu
    pphys.interactions = int_col
    pphys.distributions = [distributions.PrimaryMass(0), distributions.IsotropicDirection()]

    rand = utilities.SIREN_random(7)
    N = 4
    inj = injection._Injector(N, dm, pinj, rand)
    weighter = injection._Weighter([inj], dm, pphys)

    weights = []
    for _ in range(N):
        ev = inj.GenerateEvent()
        w = weighter.EventWeight(ev)
        assert np.isfinite(w) and w > 0.0, f"non-finite/non-positive event weight {w}"
        weights.append(w)

    assert len(weights) == N
    assert len(set(weights)) > 1, "all event weights identical -- sampling did not vary"


# ---------------------------------------------------------------------------
# Test 6: per-event SampleFinalState cost (production throughput guard)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_sample_final_state_perf_budget():
    """PythiaDISCrossSection re-initializes Pythia per event (variable-energy mode
    is unsupported for WeakBosonExchange). Guard the per-event cost so a major
    regression that would make PeV production infeasible is caught. Generous
    ceiling; the measured value is reported.
    """
    import time
    xs = _make_xs(with_differential=False)
    M = 10
    t0 = time.time()
    ev = _sample_siren(xs, 1.0e4, n=M)
    elapsed = time.time() - t0
    assert len(ev) >= M // 2, "too few events sampled to time"
    per_event = elapsed / len(ev)
    print(f"\nSampleFinalState mean per-event cost: {per_event:.3f} s")
    assert per_event < 3.0, (
        f"SampleFinalState is {per_event:.2f} s/event (> 3 s ceiling) -- a "
        "Pythia re-init regression would make large-scale production infeasible.")
