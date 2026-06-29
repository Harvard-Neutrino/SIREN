"""Physics validation for the PythiaDISCrossSection charm-DIS generator.

Guards the absolute charm-production rate, that SampleFinalState reproduces bare
Pythia's DIS kinematics (Bjorken x/y, Q^2), and that an optional differential
spline covers the realized sampling support (no silent-zero weight bias) across
the TeV-PeV band. Needs Pythia8/LHAPDF at runtime plus charm splines; gated behind
env vars, skips cleanly when unset:

    SIREN_PYTHIA_WIDE_SIGMA    -> total sigma(E) FITS spline, 100 GeV - 1 PeV
    SIREN_PYTHIA_WIDE_DSDXDY   -> differential d2sigma/dx dy FITS spline (optional)
    LHAPDF_DATA_PATH           -> LHAPDF data dir containing the PDF set below

SampleFinalState re-inits Pythia per event (~1 s/event) so SIREN statistics are
modest; the bare-Pythia reference (GeneratePythiaCharmSamples, ~ms/event) is high-stat.
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
        # secondary order: [charged lepton, Hadrons, D meson]
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


# Test 1: absolute charm-production rate / charm fraction (normalization).
def test_charm_total_cross_section_normalization():
    """The inclusive charm-DIS sigma must have the right absolute magnitude.

    Pins (a) charm fraction sigma_charm/sigma_CC at 100 GeV vs the textbook nu-N CC
    cross section, and (b) sigma_charm/E to the right order of magnitude across
    100 GeV - 1 PeV with monotonic growth.
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


# Test 2: SampleFinalState reproduces bare-Pythia DIS kinematics.
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_sampling_matches_bare_pythia_at_100gev():
    """SIREN SampleFinalState must reproduce bare Pythia's Bjorken x/y at 100 GeV.

    SampleFinalState extracts/rotates the Pythia final state and reconstructs
    (x, y); GeneratePythiaCharmSamples is the same Pythia config inline. A
    frame/extraction bug would show here as a distribution mismatch.
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

    # Physical kinematics.
    assert all(0.0 < e["x"] < 1.0 for e in siren_ev)
    assert all(0.0 < e["y"] < 1.0 for e in siren_ev)


# Test 3: D-meson energy fraction + production collimation (morphology).
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_d_meson_energy_fraction_and_collimation():
    """The D meson carries a sizeable, physical energy fraction and is nearly
    collinear with the primary lepton at high energy -- the morphology setting the
    two-cascade energy split and separation direction.
    """
    E = 1.0e4   # 10 TeV
    xs = _make_xs(with_differential=False)
    ev = _sample_siren(xs, E, n=60)
    assert len(ev) >= 30
    zD = [e["zD"] for e in ev]
    # Every D carries a physical (0, 1) energy fraction, mean in a sane range.
    assert all(0.0 < z < 1.0 for z in zD)
    assert 0.05 < _mean(zD) < 0.95
    # D/lepton opening angle is small at 10 TeV (collimated).
    def _angle(a, b):
        import math
        na = math.sqrt(sum(a[i] ** 2 for i in (1, 2, 3)))
        nb = math.sqrt(sum(b[i] ** 2 for i in (1, 2, 3)))
        dot = sum(a[i] * b[i] for i in (1, 2, 3))
        c = max(-1.0, min(1.0, dot / (na * nb)))
        return math.degrees(math.acos(c))
    angles = [_angle(e["p_D"], e["p_lep"]) for e in ev]
    # Modest median (DIS at 10 TeV is forward).
    angles.sort()
    median = angles[len(angles) // 2]
    assert median < 30.0, f"median D/lepton opening angle {median:.1f} deg too large"


# Test 4: optional differential spline covers the sampled support (no silent 0).
@pytest.mark.skipif(not (_have_diff and PYTHIA_DATA),
                    reason="set SIREN_PYTHIA_WIDE_DSDXDY + PYTHIA8DATA to run the coverage test")
def test_differential_spline_covers_sampling_support():
    """With a differential spline, DifferentialCrossSection must be finite-positive
    on essentially every sampled event, else those events get a silently-zero
    density and a biased weight. Guards spline (E, x, y) support vs realized support.
    """
    xs = _make_xs(with_differential=True)
    for E in (1.0e3, 1.0e4):   # 1 TeV, 10 TeV (within the wide spline x-range)
        ev = _sample_siren(xs, E, n=60)
        assert len(ev) >= 30, f"too few events at E={E:.0e}"
        in_range = 0
        out_of_range = 0
        silent_zero = 0
        for e in ev:
            try:
                v = xs.DifferentialCrossSection(e["ir_out"])
            except RuntimeError:
                out_of_range += 1   # out-of-support correctly RAISES
                continue
            if math.isfinite(v) and v > 0.0:
                in_range += 1
            else:
                silent_zero += 1
        # Contract: out-of-range raises; nothing returns a silent zero.
        assert silent_zero == 0, (
            f"{silent_zero} sampled events at E={E:.0e} GeV returned a silent-zero "
            "differential density instead of raising")
        frac = in_range / len(ev)
        assert frac > 0.95, (
            f"only {frac:.3f} of sampled events at E={E:.0e} GeV fell inside the "
            "differential spline (E, x, y) support; widen the spline's logx/logy range.")


# Test 5: end-to-end inject -> weight through the real Injector / Weighter.
@pytest.mark.skipif(not PYTHIA_DATA, reason="set PYTHIA8DATA to run the Pythia-sampling tests")
def test_end_to_end_rate_closure():
    """Quantitative inject->weight rate closure for charm DIS through the real
    _Injector / _Weighter on CCM.

    With injection==physical the shared distributions, cross-section probability,
    and PointSource position propagator all cancel, leaving EventWeight_i =
    InteractionProbability_i / N. Therefore:
      - per event, EventWeight == GetInteractionProbabilities/N to machine
        precision; they are independent code paths, so agreement IS the
        unbiasedness proof -- the PR#74 SampleFinalState/FinalStateProbability/
        TotalCrossSection closure-break class would break it, and
      - sum(EventWeight) == mean(P_int), bounded by a thin-target sigma*n_Ar*L_eff.
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
    E0 = 1000.0   # 1 TeV monoenergetic
    pinj.distributions = [
        distributions.PrimaryMass(0),
        distributions.Monoenergetic(E0),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 5.0),
    ]
    N = 10
    inj = injection._Injector(N, dm, pinj, rand)
    weighter = injection._Weighter([inj], dm, pphys)

    weights, int_probs = [], []
    for _ in range(N):
        ev = inj.GenerateEvent()
        w = weighter.EventWeight(ev)
        ip = weighter.GetInteractionProbabilities(ev, 0)[0]
        assert np.isfinite(w) and w > 0.0, f"non-finite/non-positive event weight {w}"
        weights.append(w)
        int_probs.append(ip)

    # (1) Per-event unbiasedness: EventWeight == single-event interaction probability / N.
    for w, ip in zip(weights, int_probs):
        assert abs(w - ip / N) <= 1e-9 * (ip / N), f"EventWeight {w} != P_int/N {ip / N}"

    # (2) Sum of weights is the unbiased physical-rate estimator == mean(P_int).
    sum_w = sum(weights)
    mean_ip = sum(int_probs) / N
    assert abs(sum_w - mean_ip) <= 1e-9 * mean_ip

    # (3) finite, positive, genuinely varying.
    assert len(set(weights)) > 1, "all event weights identical -- sampling did not vary"

    # (4) Analytic thin-target anchor: sum_w ~ sigma_charm * n_Ar * L_eff, with the
    # effective LAr path L_eff between ~0 and the 5 m injection cap (isotropic rays
    # exit the volume). sigma is the per-nucleon charm spline value at E0.
    sigma_cm2 = xs.TotalCrossSection(NuMu, E0)
    n_Ar = 1.396 * 6.022e23 / 40.0          # CCM liquid-argon number density [cm^-3]
    L_max_cm = 500.0
    upper = sigma_cm2 * n_Ar * L_max_cm
    assert 0.0 < sum_w < upper, f"sum_w {sum_w:.3e} not in (0, sigma*n*L_max={upper:.3e})"
    assert sum_w > 0.02 * upper, f"sum_w {sum_w:.3e} implausibly small vs thin-target {upper:.3e}"
