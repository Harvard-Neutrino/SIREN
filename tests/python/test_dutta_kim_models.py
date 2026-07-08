"""Regression tests for the self-contained Dutta-Kim vector portal models."""

import math
import os

import pytest

siren = pytest.importorskip("siren")
from siren.Injector import Injector as _Injector
from siren.Weighter import Weighter as _Weighter


@pytest.fixture(scope="module")
def vector_portal(processes_dir):
    from siren import _util

    return _util.load_module(
        "test_dutta_kim_VectorPortal",
        str(processes_dir / "DarkNewsTables" / "VectorPortal.py"),
    )


@pytest.fixture(scope="module")
def meson_production_module(processes_dir):
    from siren import _util

    return _util.load_module(
        "test_dutta_kim_MesonProduction",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )


@pytest.fixture(scope="module")
def beam_decays(processes_dir):
    from siren import _util

    return _util.load_module(
        "test_dutta_kim_BeamDecays",
        str(processes_dir / "BeamDecays" / "processes.py"),
    )


@pytest.fixture(scope="module")
def darknews_processes(processes_dir):
    from siren import _util

    return _util.load_module(
        "test_dutta_kim_processes",
        str(processes_dir / "DarkNewsTables" / "processes.py"),
    )


def _decay_record(signature, parent_mass, parent_momentum, secondary_masses):
    record = siren.dataclasses.InteractionRecord()
    record.signature = signature
    record.primary_mass = parent_mass
    record.primary_momentum = parent_momentum
    record.primary_initial_position = [0.0, 0.0, 0.0]
    record.interaction_vertex = [0.0, 0.0, 0.0]
    record.secondary_masses = secondary_masses
    record.secondary_momenta = [[0.0, 0.0, 0.0, 0.0]
                                for _ in secondary_masses]
    return record


def _record(signature, primary_mass, primary_energy, vertex=(0.0, 0.0, -100.0)):
    record = siren.dataclasses.InteractionRecord()
    record.signature = signature
    record.primary_mass = primary_mass
    pz = math.sqrt(max(primary_energy**2 - primary_mass**2, 0.0))
    record.primary_momentum = [primary_energy, 0.0, 0.0, pz]
    record.primary_initial_position = [0.0, 0.0, 0.0]
    record.interaction_vertex = list(vertex)
    record.secondary_masses = []
    record.secondary_momenta = []
    record.secondary_helicities = []
    return record


def _zero_secondaries(record, masses):
    record.secondary_masses = list(masses)
    record.secondary_momenta = [[0.0, 0.0, 0.0, 0.0] for _ in masses]
    record.secondary_helicities = [0 for _ in masses]
    return record


def test_vector_portal_offshell_loader_is_self_contained(darknews_processes):
    primary, secondary = darknews_processes.load_vector_portal_offshell(
        m_chi=0.008,
        m_chi_prime=0.035,
        m_V1=0.017,
        m_V2=0.05,
        g_D=1.0,
        epsilon_1=1.0e-4,
        epsilon_2=1.0e-4,
        detector_model=None,
    )

    chi_type = siren.dataclasses.ParticleType(5917)
    v1_type = siren.dataclasses.ParticleType(5922)

    assert chi_type in primary
    assert v1_type in secondary
    assert primary[chi_type][0].DensityVariables() == ["s_pair", "cos_theta_sub"]
    assert secondary[v1_type][0].DensityVariables() == ["cos_theta"]


def test_vector_portal_onshell_loader_preserves_topology(darknews_processes):
    primary, secondary = darknews_processes.load_vector_portal_onshell(
        m_chi=0.008,
        m_chi_prime=0.050,
        m_V1=0.017,
        m_V2=0.200,
        g_D=1.0,
        epsilon_1=7.0e-5,
        epsilon_2=1.0e-4,
        detector_model=None,
    )

    chi_type = siren.dataclasses.ParticleType(5917)
    chi_prime_type = siren.dataclasses.ParticleType(5918)
    v1_prod_type = siren.dataclasses.ParticleType(5922)
    v1_signal_type = siren.dataclasses.ParticleType(5923)

    assert primary[chi_type][0].DensityVariables() == ["Q2"]
    assert secondary[v1_prod_type][0].DensityVariables() == ["cos_theta"]
    assert secondary[chi_prime_type][0].DensityVariables() == ["cos_theta"]
    assert secondary[v1_signal_type][0].DensityVariables() == ["cos_theta"]


def test_dark_photon_to_chi_decay_channel_samples_rest_frame_density(
        vector_portal):
    decay = vector_portal.DarkPhotonToChiDecay(0.017, 0.008, 1.0)
    signature = decay.GetPossibleSignatures()[0]
    record = _decay_record(
        signature,
        parent_mass=0.017,
        parent_momentum=[0.017, 0.0, 0.0, 0.0],
        secondary_masses=[0.008, 0.008],
    )

    channel = siren.injection.PhysicalDecayChannel(decay, signature)
    assert channel.Measure().type == (
        siren.injection.PhaseSpaceMeasureType.SolidAngleRest
    )

    channel.Sample(siren.utilities.SIREN_random(1), None, record)

    assert channel.Density(None, record) == pytest.approx(1.0 / (4.0 * math.pi))
    assert record.secondary_masses == pytest.approx([0.008, 0.008])
    assert record.secondary_momenta[0][0] == pytest.approx(0.017 / 2.0)
    assert record.secondary_momenta[1][0] == pytest.approx(0.017 / 2.0)


def test_cone_biased_dark_photon_to_chi_decay_infers_custom(vector_portal):
    decay = vector_portal.BiasedDarkPhotonToChiDecay(
        0.017,
        0.008,
        1.0,
        detector_position=(0.0, 0.0, 10.0),
        detector_radius=1.0,
    )
    signature = decay.GetPossibleSignatures()[0]
    channel = siren.injection.PhysicalDecayChannel(decay, signature)

    assert decay.DensityVariables() == ["lab_cone_chi"]
    assert channel.Measure().type == siren.injection.PhaseSpaceMeasureType.Unspecified


def test_vector_portal_offshell_cross_section_channel_samples(vector_portal):
    cross_section = vector_portal.VectorPortalOffShellXS(
        0.008,
        0.035,
        0.017,
        0.05,
        1.0,
        1.0e-4,
    )
    signature = cross_section.GetPossibleSignatures()[0]

    record = siren.dataclasses.InteractionRecord()
    record.signature = signature
    record.primary_mass = 0.008
    energy = 1.0
    pz = math.sqrt(energy * energy - record.primary_mass * record.primary_mass)
    record.primary_momentum = [energy, 0.0, 0.0, pz]
    record.primary_initial_position = [0.0, 0.0, 0.0]
    record.interaction_vertex = [0.0, 0.0, 0.0]
    record.target_mass = 37.215
    record.secondary_masses = [0.008, 0.017, 37.215]
    record.secondary_momenta = [[0.0, 0.0, 0.0, 0.0],
                                [0.0, 0.0, 0.0, 0.0],
                                [0.0, 0.0, 0.0, 0.0]]

    channel = siren.injection.PhysicalCrossSectionChannel(
        cross_section,
        signature,
    )
    assert channel.Measure().type == siren.injection.PhaseSpaceMeasureType.Recursive2Body

    channel.Sample(siren.utilities.SIREN_random(2), None, record)
    assert record.secondary_masses == pytest.approx([0.008, 0.017, 37.215])
    assert all(momentum[0] > 0.0 for momentum in record.secondary_momenta)


def test_offshell_s_pair_sample_matches_density(vector_portal):
    """s_pair is sampled from EXACTLY the shared Breit-Wigner map whose
    Density FinalStateProbability reports (Contract C1).
    """
    import numpy as np

    cross_section = vector_portal.VectorPortalOffShellXS(
        0.008, 0.035, 0.017, 0.05, 1.0, 1.0e-4)
    E_chi = 1.0

    bw_map = cross_section._s_pair_mapping(E_chi)
    assert bw_map is not None

    s_min = bw_map.Forward(0.0)
    s_max = bw_map.Forward(1.0)
    assert s_max > s_min

    # Reported density is normalized over [s_min, s_max]: its CDF spans 0..1.
    # (Computed via the analytic CDF -- the narrow resonance defeats a uniform
    # grid quadrature.)
    assert (bw_map.Inverse(s_max) - bw_map.Inverse(s_min)) == \
        pytest.approx(1.0, abs=1e-9)

    # The sampler draws from this SAME map, so pushing samples back through the
    # map's CDF (Inverse) yields a uniform distribution on [0, 1].  This is the
    # Sample == Density check: if _sample_s_pair drifted from the reported
    # density, the transformed samples would not be uniform.
    rng = siren.utilities.SIREN_random(7)
    u = np.array([bw_map.Inverse(cross_section._sample_s_pair(E_chi, rng))
                  for _ in range(400)])
    assert np.all((u >= -1e-9) & (u <= 1.0 + 1e-9))
    assert u.mean() == pytest.approx(0.5, abs=0.08)
    assert u.min() < 0.2 and u.max() > 0.8


def test_vector_portal_upscattering_cross_section(vector_portal):
    cross_section = vector_portal.VectorPortalUpscatteringXS(
        m_chi=0.008,
        m_chi_prime=0.050,
        m_V2=0.200,
        g_D=1.0,
        epsilon=1.0e-4,
    )
    signature = cross_section.GetPossibleSignatures()[0]
    record = _record(signature, 0.008, 1.0)
    record.target_mass = 37.215
    _zero_secondaries(record, cross_section.SecondaryMasses(signature.secondary_types))

    q2min = cross_section.Q2Min(record)
    q2max = cross_section.Q2Max(record)
    q2 = 0.5 * (q2min + q2max)

    assert q2max > q2min
    assert cross_section.TotalCrossSection(record) > 0.0
    assert cross_section.DifferentialCrossSection(
        signature.primary_type,
        energy=1.0,
        Q2=q2,
    ) >= 0.0
    assert cross_section.DensityVariables() == ["Q2"]

    channel = siren.injection.PhysicalCrossSectionChannel(
        cross_section,
        signature,
    )
    channel.Sample(siren.utilities.SIREN_random(4), None, record)
    assert channel.Measure().type == siren.injection.PhaseSpaceMeasureType.MandelstamQ2
    assert channel.Density(None, record) > 0.0


def test_dutta_kim_directed_phase_space_topology(vector_portal, processes_dir):
    from siren import _util

    meson = _util.load_module(
        "test_dutta_kim_MesonProduction_for_topology",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )

    fiducial = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
        0.1, 0.1, 0.1)
    random = siren.utilities.SIREN_random(1337)

    pion_decay = meson.MesonThreeBodySIRENDecay(m_mediator=0.017, g_mu=1.0e-3)
    pion_record = _record(
        pion_decay.GetPossibleSignatures()[0],
        0.13957039,
        0.13957039,
    )
    _zero_secondaries(pion_record, [0.10565837, 0.0, 0.017])
    pion_channel = siren.injection.DetectorDirected3BodyChannel(
        factorization=siren.injection.ThreeBodyMode.Recursive,
        target=fiducial,
        spectator_index=0,
        pair_first_index=1,
        pair_second_index=2,
        directed_index=2,
        mass_mode=siren.injection.InvariantMassMode.Uniform,
    )
    pion_channel.Sample(random, None, pion_record)
    assert pion_channel.Density(None, pion_record) > 0.0
    assert pion_decay.FinalStateProbability(pion_record) >= 0.0

    v1_decay = vector_portal.DarkPhotonToChiDecay(0.017, 0.008, 1.0)
    v1_record = _record(v1_decay.GetPossibleSignatures()[0], 0.017, 0.08)
    _zero_secondaries(v1_record, [0.008, 0.008])
    v1_channel = siren.injection.DetectorDirected2BodyChannel(fiducial, 0)
    v1_channel.Sample(random, None, v1_record)
    assert v1_channel.Density(None, v1_record) > 0.0
    assert v1_decay.FinalStateProbability(v1_record) > 0.0

    scatter = vector_portal.VectorPortalUpscatteringXS(
        0.008,
        0.050,
        0.200,
        1.0,
        1.0e-4,
    )
    scatter_sig = scatter.GetPossibleSignatures()[0]
    scatter_record = _record(scatter_sig, 0.008, 1.0)
    scatter_record.target_mass = 37.215
    _zero_secondaries(
        scatter_record,
        scatter.SecondaryMasses(scatter_sig.secondary_types),
    )
    scatter_channel = siren.injection.DetectorDirectedScatteringChannel(
        fiducial,
        directed_index=0,
        variable=siren.injection.ScatteringVariable.Q2,
    )
    scatter_channel.Sample(random, None, scatter_record)
    assert scatter_channel.Density(None, scatter_record) > 0.0
    assert scatter.DifferentialCrossSection(scatter_record) >= 0.0

    chi_prime_decay = vector_portal.ChiPrimeDecay(0.008, 0.050, 0.017, 1.0)
    chi_prime_record = _record(
        chi_prime_decay.GetPossibleSignatures()[0],
        0.050,
        0.12,
    )
    _zero_secondaries(chi_prime_record, [0.008, 0.017])
    chi_prime_channel = siren.injection.DetectorDirected2BodyChannel(fiducial, 1)
    chi_prime_channel.Sample(random, None, chi_prime_record)
    assert chi_prime_channel.Density(None, chi_prime_record) > 0.0
    assert chi_prime_decay.FinalStateProbability(chi_prime_record) > 0.0


def test_meson_three_body_recursive2body_jacobian_matches_sampler(processes_dir):
    """MesonThreeBodySIRENDecay: FinalStateProbability equals the density its own
    SampleFinalState draws, verifying the Recursive2Body-to-Dalitz Jacobian.

    The rejection sampler draws (E_nu, E_phi) proportional to |M|^2 in the parent
    rest frame, then places the overall orientation isotropically. So its
    generated density in the (s_pair, cos_theta_sub) marginal is

        g(s_pair, cos_sub) = |M|^2 / Z * |d(E_nu,E_phi)/d(s_pair,cos_sub)|

    with Z = integral |M|^2 dE_nu dE_phi = 64 pi^3 m_M Gamma (from total_width),
    and the change-of-variables Jacobian |d(E_nu,E_phi)/d(s_pair,cos_sub)| in
    closed form is (1/(2 m_M)) (p_pair/sqrt_s) p_nu_star -- an INDEPENDENT
    kinematic identity, NOT the Recursive2Body Jacobian under test.

    FinalStateProbability returns the full R2B density; its marginal over the two
    flat isotropic angles (dOmega_pair = 4 pi, azimuth = 2 pi) multiplies by
    8 pi^2. The per-sample importance weight w = f_marg / g must therefore be a
    CONSTANT (== 1) for every sample when the R2B Jacobian matches the sampler.
    p_spectator = sqrt(lambda_parent)/(2 sqrt(s_pair)). Substituting /(2 m_M)
    makes w vary sample to sample (CV ~ 0.12) and mis-scales it, so this test
    fails on an incorrect Jacobian.
    """
    import numpy as np
    from siren import _util

    meson = _util.load_module(
        "test_dutta_kim_MesonProduction_jacobian",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )

    decay = meson.MesonThreeBodySIRENDecay(
        m_meson=0.13957039, m_lepton=0.10565837, m_mediator=0.017, g_mu=1.0e-3)
    core = decay._decay
    m_M = core.m_M
    m_l = core.m_l
    m_phi = core.m_phi
    # Z = integral |M|^2 dE_nu dE_phi, from Gamma = Z / (64 pi^3 m_M).
    Z = 64.0 * math.pi ** 3 * m_M * decay._total_width

    signature = decay.GetPossibleSignatures()[0]

    def _rest_energies(rec):
        """(E_nu, E_phi) in the parent rest frame from the record secondaries."""
        idx = {int(st): i for i, st in enumerate(signature.secondary_types)}
        i_nu = idx[decay.pdgid_neutrino]
        i_phi = idx[decay.pdgid_mediator]
        p_parent = np.array(rec.primary_momentum, dtype=float)
        e0 = p_parent[0]
        mom0 = np.linalg.norm(p_parent[1:])
        mass0 = math.sqrt(max(e0 ** 2 - mom0 ** 2, 0.0))
        p_nu = np.array(rec.secondary_momenta[i_nu], dtype=float)
        p_phi = np.array(rec.secondary_momenta[i_phi], dtype=float)
        if mom0 < 1e-12:
            return p_nu[0], p_phi[0]
        beta = mom0 / e0
        gamma = e0 / mass0
        bhat = p_parent[1:] / mom0
        e_nu = gamma * (p_nu[0] - beta * float(np.dot(p_nu[1:], bhat)))
        e_phi = gamma * (p_phi[0] - beta * float(np.dot(p_phi[1:], bhat)))
        return e_nu, e_phi

    def _sampler_density_marginal(e_nu, e_phi):
        """g(s_pair, cos_sub): |M|^2/Z times an independent change-of-vars Jacobian."""
        matel = core._matel_sq(e_nu, e_phi)
        if matel <= 0.0:
            return 0.0
        e_l = m_M - e_nu - e_phi
        s_pair = m_M ** 2 + m_l ** 2 - 2.0 * m_M * e_l
        if s_pair <= 0.0:
            return 0.0
        sqrt_s = math.sqrt(s_pair)
        lam_parent = (m_M ** 2 - (m_l + sqrt_s) ** 2) * (m_M ** 2 - (m_l - sqrt_s) ** 2)
        if lam_parent <= 0.0:
            return 0.0
        p_pair = math.sqrt(lam_parent) / (2.0 * m_M)
        p_nu_star = (s_pair - m_phi ** 2) / (2.0 * sqrt_s)
        d_e_d_sc = (1.0 / (2.0 * m_M)) * (p_pair / sqrt_s) * p_nu_star
        return matel / Z * d_e_d_sc

    # Test both a rest-frame pion and a strongly boosted one; the R2B density is
    # frame-independent, so the weight must be flat in either frame.
    for e_pi in (m_M, 2.0):
        pz = math.sqrt(max(e_pi ** 2 - m_M ** 2, 0.0))
        template = siren.dataclasses.InteractionRecord()
        template.signature = signature
        template.primary_mass = m_M
        template.primary_momentum = [e_pi, 0.0, 0.0, pz]
        template.secondary_masses = list(decay.SecondaryMasses(signature.secondary_types))
        template.secondary_momenta = [[0.0, 0.0, 0.0, 0.0]
                                      for _ in template.secondary_masses]
        template.secondary_helicities = [0, 0, 0]

        rng = siren.utilities.SIREN_random(2024)
        weights = []
        for _ in range(20000):
            csdr = siren.dataclasses.CrossSectionDistributionRecord(template)
            decay.SampleFinalState(csdr, rng)
            out = siren.dataclasses.InteractionRecord()
            out.signature = signature
            out.primary_mass = template.primary_mass
            out.primary_momentum = template.primary_momentum
            out.secondary_masses = list(template.secondary_masses)
            out.secondary_momenta = [[0.0, 0.0, 0.0, 0.0]
                                     for _ in template.secondary_masses]
            out.secondary_helicities = [0, 0, 0]
            csdr.finalize(out)
            e_nu, e_phi = _rest_energies(out)
            g = _sampler_density_marginal(e_nu, e_phi)
            if g <= 0.0:
                continue
            f_marg = decay.FinalStateProbability(out) * 8.0 * math.pi ** 2
            weights.append(f_marg / g)

        weights = np.array(weights)
        assert weights.size > 15000
        mean = weights.mean()
        cv = weights.std() / mean
        # Sample == Density: the ratio is a constant equal to 1. The wrong
        # /(2 m_M) Jacobian gives mean ~ 0.2 and CV ~ 0.12, far outside these.
        assert mean == pytest.approx(1.0, abs=1e-3)
        assert cv < 1.0e-2


def test_meson_three_body_sampler_matches_dalitz_marginal(processes_dir):
    """_sample_energies draws |M|^2 on the Dalitz plane: the E_nu marginal
    matches the band-integrated physical marginal and rejects the
    band-averaged alternative that a uniform-in-band proposal with
    unweighted acceptance would produce."""
    import numpy as np
    from siren import _util

    meson = _util.load_module(
        "test_dutta_kim_MesonProduction_marginal",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )

    d = meson.MesonThreeBodyDecay(0.13957039, 0.10565837, 0.017, 1.0e-3, "scalar")

    # Deterministic bound on the acceptance weight |M|^2 * band(E_nu).
    max_weight = 0.0
    for e_nu in np.linspace(0.0, d.E_nu_max, 200):
        lims = d._E_phi_limits(e_nu)
        if lims[0] is None:
            continue
        band = lims[1] - lims[0]
        for e_phi in np.linspace(lims[0], lims[1], 40):
            max_weight = max(max_weight, d._matel_sq(e_nu, e_phi) * band)
    max_weight *= 1.2

    rng = siren.utilities.SIREN_random(123)
    n = 20000
    samples = np.array(
        [meson._sample_energies(d, max_weight, rng)[0] for _ in range(n)])

    grid = np.linspace(1e-6, d.E_nu_max * (1 - 1e-6), 100)
    p_phys, p_avg = [], []
    for e_nu in grid:
        lims = d._E_phi_limits(e_nu)
        if lims[0] is None or lims[1] <= lims[0]:
            p_phys.append(0.0)
            p_avg.append(0.0)
            continue
        xs = np.linspace(lims[0], lims[1], 300)
        integ = np.trapezoid([d._matel_sq(e_nu, x) for x in xs], xs)
        p_phys.append(integ)                        # int |M|^2 dE_phi
        p_avg.append(integ / (lims[1] - lims[0]))   # band-averaged |M|^2
    p_phys = np.array(p_phys) / np.trapezoid(p_phys, grid)
    p_avg = np.array(p_avg) / np.trapezoid(p_avg, grid)

    hist, edges = np.histogram(samples, bins=24, range=(0.0, d.E_nu_max),
                               density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    err = np.sqrt(np.maximum(hist, 1e-12) / (n * (edges[1] - edges[0])))
    chi2_phys = np.mean(
        ((hist - np.interp(centers, grid, p_phys)) / np.maximum(err, 1e-12)) ** 2)
    chi2_avg = np.mean(
        ((hist - np.interp(centers, grid, p_avg)) / np.maximum(err, 1e-12)) ** 2)

    assert chi2_phys < 3.0, (
        f"E_nu marginal disagrees with the physical Dalitz marginal "
        f"(chi2/dof {chi2_phys:.2f})")
    assert chi2_avg > 10.0, (
        f"test lost discriminating power (chi2/dof vs band-averaged "
        f"alternative {chi2_avg:.2f})")


def test_meson_three_body_check_closure_passes(processes_dir):
    """siren.check_closure certifies MesonThreeBodySIRENDecay: its Recursive2Body
    sampler and FinalStateProbability are the same distribution.

    The gauge reads the fallback cos(theta) in the parent rest frame for a
    Recursive2Body measure; the lab coordinate is degenerate at ~1 for a
    boosted parent, collapsing the flatness test to a single bin.
    """
    from siren import _util

    meson = _util.load_module(
        "test_dutta_kim_MesonProduction_closure",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )

    decay = meson.MesonThreeBodySIRENDecay(m_mediator=0.017, g_mu=1.0e-3)
    report = siren.check_closure(decay, primary_energy=2.0, samples=6000, seed=0)
    assert report.ok, str(report)
    fmean, _ferr = report.flatness
    assert fmean == pytest.approx(1.0, abs=0.1)


def test_meson_three_body_width_absolute_anchors(meson_production_module):
    """Pin the absolute three-body widths pi/K -> l nu phi in the PDG f_M
    convention, anchored to the explicit Dirac-spinor |M|^2 and to the
    Dutta-Kim PRL 129, 111803 Table 2 branching ratios per g^2."""
    mp = meson_production_module

    m_pi = 0.13957039
    m_K = 0.49368
    m_mu = 0.10565837

    # total_width() in the PDG f_M convention (g_mu = 1.0).
    expected_width = {
        ("pi", "scalar"): 1.2247070252317424e-18,
        ("pi", "pseudoscalar"): 2.773199175624213e-21,
        ("K", "scalar"): 4.9727779507952735e-18,
        ("K", "pseudoscalar"): 5.362959599511565e-19,
    }

    m_meson = {"pi": m_pi, "K": m_K}
    widths = {}
    for meson, m_M in m_meson.items():
        for ptype in ("scalar", "pseudoscalar"):
            decay = mp.MesonThreeBodyDecay(m_M, m_mu, 0.001, 1.0, ptype)
            w = decay.total_width()
            widths[(meson, ptype)] = w
            assert w == pytest.approx(expected_width[(meson, ptype)], rel=1e-6)

    # Dutta-Kim Table 2 sanity band: abs BR per g^2 = total_width / Gamma_tot
    # / g^2 (g_mu = 1.0 above). Gamma_tot from the PDG mean lifetimes.
    hbar = 6.582119569e-25  # GeV s
    gamma_tot = {"pi": hbar / 2.6033e-8, "K": hbar / 1.2380e-8}
    table_br = {
        ("pi", "scalar"): 4.8e-2,
        ("pi", "pseudoscalar"): 1.1e-4,
        ("K", "scalar"): 9.2e-2,
        ("K", "pseudoscalar"): 1.0e-2,
    }
    for key, w in widths.items():
        meson, _ptype = key
        br_per_g2 = w / gamma_tot[meson]
        assert br_per_g2 == pytest.approx(table_br[key], rel=0.12)


def test_vector_meson_br_ratio_anchors(vector_portal):
    """Pin the M -> l nu V1 production branching ratio Gamma3/Gamma2.

    Two independent anchors fix it: the six-figure spinor-sum values (explicit
    Dirac spinor sums agreeing with the Carlson-Rislow closed form), and the
    Dutta-Kim PRL 129, 111803 Table 2 vector entries once the paper's stated
    BR(V1 -> observed) = 50% is folded out (its ee rows equal its invisible-X
    rows, fixing that reading).
    """
    vp = vector_portal
    m_pi = 0.13957039
    m_K = 0.49368
    m_mu = 0.10565837
    m_V1 = 0.017
    epsilon = 7e-5

    br3_pi = vp._meson_to_v1_branching_ratio(m_pi, m_mu, m_V1, epsilon)
    br3_K = vp._meson_to_v1_branching_ratio(m_K, m_mu, m_V1, epsilon)
    assert br3_pi == pytest.approx(3.80806e-13, rel=1e-4)
    assert br3_K == pytest.approx(2.21940e-09, rel=1e-4)

    # Absolute BR against 2x the Dutta-Kim Table 2 vector entries (the table
    # folds BR(V1 -> 2chi) = 0.5).  BR2 are the PDG M -> mu nu branchings.
    BR2_pi = 0.999877
    BR2_K = 0.6356
    assert br3_pi * BR2_pi == pytest.approx(2.0 * 0.17e-12, rel=0.20)
    assert br3_K * BR2_K == pytest.approx(2.0 * 680e-12, rel=0.08)

    # C_V^2 scales as epsilon^2, so the ratio scales as epsilon^2 exactly.
    br3_pi_6 = vp._meson_to_v1_branching_ratio(m_pi, m_mu, m_V1, 6e-5)
    assert br3_pi_6 == pytest.approx(br3_pi * (6.0 / 7.0)**2, rel=1e-9)


def test_dark_sector_width_anchors(vector_portal):
    """Pin the dark-sector partial widths to explicit spinor-sum
    values, the HNL limit, and the shared V1 -> ee Lorentz structure.
    """
    vp = vector_portal

    # V1 -> chi chi vector-pair width (explicit Dirac spinor sum).
    assert vp.DarkPhotonToChiDecay(0.017, 0.008, 1.0).total_width() == \
        pytest.approx(2.19869e-04, rel=1e-4)

    # chi' -> chi V1 fermion -> fermion + vector width.  Signature is
    # ChiPrimeDecay(m_chi, m_chi_prime, m_V1, g_D).
    assert vp.ChiPrimeDecay(0.008, 0.040, 0.017, 1.0).total_width() == \
        pytest.approx(2.77278e-03, rel=1e-4)
    assert vp.ChiPrimeDecay(0.008, 0.050, 0.017, 1.0).total_width() == \
        pytest.approx(6.81945e-03, rel=1e-4)

    # HNL N -> nu V limit: as m_chi -> 0 the width collapses to the closed
    # form g^2 m'^3/(16 pi m_V^2) (1 - m_V^2/m'^2)^2 (1 + 2 m_V^2/m'^2).
    mp, mv, g = 0.040, 0.017, 1.0
    hnl_expected = (g**2 * mp**3 / (16.0 * math.pi * mv**2)
                    * (1.0 - mv**2 / mp**2)**2 * (1.0 + 2.0 * mv**2 / mp**2))
    assert vp._chi_prime_to_chi_v1_width(0.040, 1e-12, 0.017, 1.0) == \
        pytest.approx(hnl_expected, rel=1e-6)

    # V1 -> ee shares the vector-pair implementation with V1 -> chi chi;
    # it reproduces the (alpha eps^2 m_V/3) beta (1 + 2 r^2) inline form
    # bit-for-bit.
    me, m_V1 = 0.000511, 0.017
    beta_e = math.sqrt(1.0 - (2.0 * me / m_V1)**2)
    ee_expected = ((1.0 / 137.036) * m_V1 / 3.0
                   * beta_e * (1.0 + 2.0 * (me / m_V1)**2))
    assert vp.DarkPhotonDecay(0.017, 1.0).total_width() == \
        pytest.approx(ee_expected, rel=1e-9)


# ------------------------------------------------------------------ #
#  End-to-end chain test                                               #
# ------------------------------------------------------------------ #

M_CHI = 8e-3
M_CHI_PRIME = 50e-3
M_V1 = 17e-3
M_V2 = 200e-3
M_PION = 0.13957039
M_MUON = 0.10565837
M_ARGON40 = 37.215

G_D = 1.0
EPSILON_1 = 7e-5
EPSILON_2 = 1e-4
G_MU = 1e-3

PDGID_PION = 211
PDGID_MUPLUS = -13
PDGID_NUMU = 14
PDGID_V1_PROD = 5922
PDGID_CHI = 5917
PDGID_CHI_PRIME = 5918
PDGID_V1_SIGNAL = 5923
PDGID_ARGON40 = 1000180400


def _pt(pdgid):
    return siren.dataclasses.Particle.ParticleType(pdgid)


def _make_mc(channels, weights):
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = channels
    mc.weights = weights
    return mc


def test_dutta_kim_end_to_end_chain(vector_portal, processes_dir):
    """Generate events through the full 5-vertex on-shell chain and weight them.

    Uses the SBND detector (LAr TPC) from the SBN loader. The pion is
    the primary (monoenergetic at 2 GeV, fixed +z direction from the
    BNB target). Secondaries are:
      V1_prod -> chi chi
      chi Ar40 -> chi' Ar40
      chi' -> chi V1_signal
      V1_signal -> e+ e-

    Each secondary vertex has detector-directed biasing toward a box
    approximating the SBND fiducial volume. The test asserts that
    generated events reach the full chain depth and that weights are
    finite and positive.
    """
    from siren import _util

    meson = _util.load_module(
        "test_e2e_MesonProduction",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )

    # -- Detector --
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    det_origin = detector_model.DetectorOrigin.get()
    det_x, det_y, det_z = det_origin.GetX(), det_origin.GetY(), det_origin.GetZ()

    placement = siren.geometry.Placement(
        siren.math.Vector3D(det_x, det_y, det_z))
    fiducial = siren.geometry.Box(placement, 2.0, 2.0, 2.5)

    # -- Physics models --
    pion_decay = meson.MesonThreeBodySIRENDecay(
        m_meson=M_PION,
        m_lepton=M_MUON,
        m_mediator=M_V1,
        g_mu=G_MU,
        pdgid_meson=PDGID_PION,
        pdgid_lepton=PDGID_MUPLUS,
        pdgid_neutrino=PDGID_NUMU,
        pdgid_mediator=PDGID_V1_PROD,
    )

    v1_to_chi = vector_portal.DarkPhotonToChiDecay(
        M_V1, M_CHI, G_D,
        pdgid_V1=PDGID_V1_PROD,
        pdgid_chi=PDGID_CHI,
    )

    upscatter = vector_portal.VectorPortalUpscatteringXS(
        M_CHI, M_CHI_PRIME, M_V2, G_D, EPSILON_2,
        pdgid_chi=PDGID_CHI,
        pdgid_chi_prime=PDGID_CHI_PRIME,
        nuclear_pdgid=PDGID_ARGON40,
        nuclear_mass=M_ARGON40,
        A=40, Z=18,
    )

    chi_prime_decay = vector_portal.ChiPrimeDecay(
        M_CHI, M_CHI_PRIME, M_V1, G_D,
        pdgid_chi_prime=PDGID_CHI_PRIME,
        pdgid_chi=PDGID_CHI,
        pdgid_V1=PDGID_V1_SIGNAL,
    )

    visible_decay = vector_portal.DarkPhotonDecay(
        M_V1, EPSILON_1,
        pdgid_V1=PDGID_V1_SIGNAL,
    )

    # -- Primary injection distributions --
    primary_mass = siren.distributions.PrimaryMass(M_PION)
    primary_energy = siren.distributions.Monoenergetic(2.0)
    primary_direction = siren.distributions.FixedDirection(
        siren.math.Vector3D(0, 0, 1)
    )
    primary_position = siren.distributions.PointSourcePositionDistribution(
        [0.0, 0.0, 0.0],
        300.0,
    )

    # -- Secondary distributions --
    secondary_vertex = siren.distributions.SecondaryPhysicalVertexDistribution()
    secondary_bounded_vertex = siren.distributions.SecondaryBoundedVertexDistribution(fiducial, 1000.0)

    # -- Biased phase spaces for each secondary --
    v1_prod_sig = v1_to_chi.GetPossibleSignatures()[0]
    chi_sig = upscatter.GetPossibleSignatures()[0]
    chi_prime_sig = chi_prime_decay.GetPossibleSignatures()[0]
    v1_signal_sig = visible_decay.GetPossibleSignatures()[0]

    secondary_phase_spaces = {
        _pt(PDGID_V1_PROD): {
            v1_prod_sig: _make_mc(
                [siren.injection.PhysicalDecayChannel(v1_to_chi, v1_prod_sig),
                 siren.injection.DetectorDirected2BodyChannel(fiducial, 0)],
                [0.01, 0.99],
            ),
        },
        _pt(PDGID_CHI): {
            chi_sig: _make_mc(
                [siren.injection.PhysicalCrossSectionChannel(upscatter, chi_sig),
                 siren.injection.DetectorDirectedScatteringChannel(
                    fiducial,
                    directed_index=0,
                    variable=siren.injection.ScatteringVariable.Q2,
                )],
                [0.01, 0.99],
            ),
        },
        _pt(PDGID_CHI_PRIME): {
            chi_prime_sig: _make_mc(
                [siren.injection.PhysicalDecayChannel(chi_prime_decay, chi_prime_sig),
                 siren.injection.DetectorDirected2BodyChannel(fiducial, 1)],
                [0.01, 0.99],
            ),
        },
        _pt(PDGID_V1_SIGNAL): {
            v1_signal_sig: _make_mc(
                [siren.injection.Isotropic2BodyChannel(0)],
                [1.0],
            ),
        },
    }

    # -- Build Injector --
    N_EVENTS = 5
    injector = _Injector(
        number_of_events=N_EVENTS,
        detector_model=detector_model,
        seed=42,
        primary_type=_pt(PDGID_PION),
        primary_interactions=[pion_decay],
        primary_injection_distributions=[
            primary_mass,
            primary_energy,
            primary_direction,
            primary_position,
        ],
        secondary_interactions={
            _pt(PDGID_V1_PROD): [v1_to_chi],
            _pt(PDGID_CHI): [upscatter],
            _pt(PDGID_CHI_PRIME): [chi_prime_decay],
            _pt(PDGID_V1_SIGNAL): [visible_decay],
        },
        secondary_injection_distributions={
            _pt(PDGID_V1_PROD): [secondary_vertex],
            _pt(PDGID_CHI): [secondary_vertex],
            _pt(PDGID_CHI_PRIME): [secondary_vertex],
            _pt(PDGID_V1_SIGNAL): [secondary_bounded_vertex],
        },
        secondary_phase_spaces=secondary_phase_spaces,
        stopping_condition=lambda tree, datum, i: not (
            int(datum.record.signature.secondary_types[i]) == PDGID_V1_PROD
            or (int(datum.record.signature.secondary_types[i]) == PDGID_CHI
                and int(datum.record.signature.primary_type) == PDGID_V1_PROD)
            or int(datum.record.signature.secondary_types[i]) == PDGID_CHI_PRIME
            or int(datum.record.signature.secondary_types[i]) == PDGID_V1_SIGNAL
        ),
    )

    # -- Generate events --
    events = []
    for event in injector:
        if event.tree:
            events.append(event)
        if len(events) >= N_EVENTS:
            break

    assert len(events) > 0, "No events were generated"

    for event in events:
        records = list(event.tree)
        assert len(records) >= 2, (
            f"Event tree has only {len(records)} records, "
            f"expected at least the primary + one secondary"
        )

    # -- Build Weighter and compute weights --
    weighter = _Weighter(
        injectors=[injector],
        detector_model=detector_model,
        primary_type=_pt(PDGID_PION),
        primary_interactions=[pion_decay],
        primary_physical_distributions=[
            primary_energy,
            primary_direction,
        ],
        secondary_interactions={
            _pt(PDGID_V1_PROD): [v1_to_chi],
            _pt(PDGID_CHI): [upscatter],
            _pt(PDGID_CHI_PRIME): [chi_prime_decay],
            _pt(PDGID_V1_SIGNAL): [visible_decay],
        },
    )

    for event in events:
        w = weighter(event)
        assert math.isfinite(w), f"Weight is not finite: {w}"
        assert w > 0.0, f"Weight is not positive: {w}"


def _model_family(vector_portal, meson_production_module, beam_decays):
    """One instance of every Dutta-Kim model class at test-scale parameters,
    plus the BeamDecays SM two-body decay the beam chains inject from."""
    vp, mp, beam = vector_portal, meson_production_module, beam_decays
    return [
        vp.DarkPhotonToChiDecay(m_V1=0.017, m_chi=0.008, g_D=1.0),
        vp.ChiPrimeDecay(m_chi_prime=0.050, m_chi=0.008, m_V1=0.017, g_D=1.0),
        vp.VectorPortalUpscatteringXS(
            m_chi=0.008, m_chi_prime=0.050, m_V2=0.200, g_D=1.0,
            epsilon=1.0e-4),
        vp.VectorPortalOffShellXS(0.008, 0.035, 0.017, 0.05, 1.0, 1.0e-4),
        vp.BiasedDarkPhotonToChiDecay(
            0.017, 0.008, 1.0,
            detector_position=(0.0, 0.0, 10.0), detector_radius=1.0),
        beam.MesonTwoBodyLeptonicDecay(211),
        mp.BiasedMesonThreeBodyDecay(
            m_mediator=0.017, g_mu=1.0e-3,
            detector_position=(0.0, 0.0, 10.0), detector_radius=1.0),
    ]


def test_audit_overrides_passes_over_model_family(
        vector_portal, meson_production_module, beam_decays):
    """Every model satisfies the trampoline override contract: no required
    virtual resolves only to the abstract root, and models relying on the
    default sampler declare a self-contained measure."""
    from siren import _validation

    _validation.audit_overrides(
        _model_family(vector_portal, meson_production_module, beam_decays))


def test_check_closure_family_coverage(vector_portal, meson_production_module,
                                       beam_decays):
    """check_closure certifies (or honestly declines) every model:

    * on-shell decays and the Q2 upscatter get a full shape verdict;
    * the cone-biased decay is gauged over its empirical angular range;
    * the off-shell 2->3 scatter is out of the single-coordinate shape
      check's jurisdiction and must say so rather than fail, since its
      mixture-level E[f/g] closure is pinned elsewhere in this suite.
    """
    import siren

    vp = vector_portal
    family = _model_family(vector_portal, meson_production_module, beam_decays)
    for model in family:
        report = siren.check_closure(
            model, primary_energy=1.0, samples=4000, seed=0)
        assert report.ok, "%s: %s" % (type(model).__name__, report)

    offshell = vp.VectorPortalOffShellXS(0.008, 0.035, 0.017, 0.05, 1.0, 1.0e-4)
    report = siren.check_closure(
        offshell, primary_energy=1.0, samples=4000, seed=0)
    assert any("not applicable" in n for n in report.notes), report.notes
    assert any("kinematically pinned" in n for n in report.notes), report.notes


def test_chi_box_edges_match_two_body_kinematics(vector_portal):
    """_chi_box_edges reproduces the analytic two-body box for equal daughters."""
    vp = vector_portal
    m_V1, m_chi, E_V = 0.017, 0.008, 5.0
    E_rf = m_V1 / 2.0
    p_rf = math.sqrt(E_rf**2 - m_chi**2)
    gamma = E_V / m_V1
    beta = math.sqrt(1.0 - 1.0 / gamma**2)
    lo, hi = vp._chi_box_edges(E_V, m_V1, m_chi, m_chi)
    assert lo == pytest.approx(gamma * (E_rf - beta * p_rf), rel=1e-12)
    assert hi == pytest.approx(gamma * (E_rf + beta * p_rf), rel=1e-12)


def test_chi_flux_closed_channel_raises(vector_portal):
    """Pins the typed closed-channel guard at the paper-forbidden V1 -> chi chi' configuration."""
    import numpy as np

    vp = vector_portal

    with pytest.raises(siren.utilities.ConfigurationError):
        vp.compute_chi_flux(
            0.13957, 0.10566, 0.017, 0.008, 0.050, 1.0, 7e-5,
            "pion_numu", 0.05, 3.0, physically_normalized=False)

    dk2nu_data = {
        "ptype": np.array([211, 211, 211]),
        "E": np.array([1.0, 1.5, 2.0]),
        "nimpwt": np.ones(3),
        "pot": 1.0,
    }
    with pytest.raises(siren.utilities.ConfigurationError):
        vp.compute_chi_flux_from_dk2nu(
            dk2nu_data, 211, 0.13957, 0.10566, 0.017, 0.008, 0.050,
            1.0, 7e-5, 0.05, 3.0)

    with pytest.raises(ValueError):
        vp._chi_box_edges(1.0, 0.017, 0.008, 0.050)


def test_chi_flux_integral_pin(vector_portal):
    """compute_chi_flux absolute normalization pinned end-to-end on the PionKaon table.

    Re-derives against an independent from-scratch continuous fold (anchor 1.3606835e+06).
    """
    import numpy as np

    vp = vector_portal
    flux = vp.compute_chi_flux(
        0.13957, 0.10566, 0.017, 0.008, 0.008, 1.0, 1.0e-4,
        "pion_numu", 0.05, 3.0, physically_normalized=False)
    nodes = np.array(flux.GetEnergyNodes())
    vals = np.array([flux.SampleUnnormedPDF(float(e)) for e in nodes])
    integral = float(np.trapezoid(vals, nodes))
    assert integral == pytest.approx(1360683.525557477, rel=1e-5)


def test_chi_flux_absolute_scale_and_multiplicity(vector_portal):
    """Locks the 2 x BR(V1 -> chi chi) multiplicity and the on-axis forward two-body parent-energy inversion in the chi-flux builder."""
    import numpy as np

    vp = vector_portal
    f1 = vp.compute_chi_flux(
        0.13957, 0.10566, 0.017, 0.008, 0.008, 1.0, 1.0e-4,
        "pion_numu", 0.05, 3.0, physically_normalized=False)

    # Both dark-matter daughters count; the width model gives BR(V1 -> chi chi) = 1 here.
    assert vp._v1_to_chi_chi_branching(0.017, 0.008, 1.0, 1.0e-4) == \
        pytest.approx(1.0, abs=1e-8)

    nodes = np.array(f1.GetEnergyNodes())
    vals = np.array([f1.SampleUnnormedPDF(float(e)) for e in nodes])
    assert float(np.trapezoid(vals, nodes)) == pytest.approx(
        1360683.525557477, rel=1e-5)

    # Forward two-body inversion E_M = 0.5 m_M (k + 1/k); k <= 1 (E_nu <= E_nu_rf)
    # has no forward solution.
    m_pi = 0.13957039
    E_nu_rf = (m_pi**2 - 0.10565837**2) / (2.0 * m_pi)
    assert vp._meson_energy_from_forward_nu(E_nu_rf - 1e-6, m_pi, E_nu_rf) is None
    k = 2.0 / E_nu_rf
    assert vp._meson_energy_from_forward_nu(2.0, m_pi, E_nu_rf) == \
        pytest.approx(0.5 * m_pi * (k + 1.0 / k))


def test_coherent_total_cross_section_magnitude_pin(vector_portal):
    """Pins the coherent argon magnitude (Z^2 = 324) of the Dutta-Kim M_0,SF
    amplitude with exact Kallen flux; anchored against the paper-formula
    integral and the forward slope 4 pi alpha alpha_D eps^2 Z^2 / m_V^4."""
    xs = vector_portal.VectorPortalUpscatteringXS(
        m_chi=0.008, m_chi_prime=0.050, m_V2=0.200, g_D=1.0, epsilon=1.0e-4)
    signature = xs.GetPossibleSignatures()[0]
    record = _record(signature, 0.008, 1.0)
    record.target_mass = xs._ups.MA
    assert xs.TotalCrossSection(record) == pytest.approx(
        4.2364767575893256e-35, rel=1e-9)


def test_biased_meson_adaptive_cone_pointwise_closure(meson_production_module):
    """Adaptive-cone biased decay: FinalStateProbability equals the physical
    density times 4pi/Omega(event), the invariant the single-coordinate gauge
    cannot certify when the cone varies with energy."""
    mp = meson_production_module
    biased = mp.BiasedMesonThreeBodyDecay(
        m_mediator=0.017, g_mu=1.0e-3,
        detector_position=(3.0, 1.0, 100.0), detector_radius=1.0)
    physical = mp.MesonThreeBodySIRENDecay(m_mediator=0.017, g_mu=1.0e-3)
    signature = biased.GetPossibleSignatures()[0]

    E_pi = 0.5
    m_pi = biased.m_meson
    random = siren.utilities.SIREN_random(7)
    channel = siren.injection.PhysicalDecayChannel(biased, signature)

    checked = 0
    for _ in range(200):
        record = _decay_record(
            signature,
            parent_mass=m_pi,
            parent_momentum=[E_pi, 0.0, 0.0, math.sqrt(E_pi**2 - m_pi**2)],
            secondary_masses=list(biased.SecondaryMasses(signature.secondary_types)),
        )
        channel.Sample(random, None, record)
        fb = biased.FinalStateProbability(record)
        fp = physical.FinalStateProbability(record)
        assert fb > 0.0
        E_phi_rf = biased._rest_frame_mediator_energy(record)
        gamma = E_pi / m_pi
        _, half_angle = biased._get_cone_params(
            record.interaction_vertex, gamma * E_phi_rf)
        omega = biased._cone_solid_angle(half_angle)
        assert fb == pytest.approx(fp * 4.0 * math.pi / omega, rel=1e-9)
        checked += 1
    assert checked == 200
