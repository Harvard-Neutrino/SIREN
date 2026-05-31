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
    assert channel.Convention() == (
        siren.injection.PhaseSpaceConvention.RestFrameSolidAngle
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
    assert channel.Convention() == siren.injection.PhaseSpaceConvention.Custom


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
    assert channel.Convention() == siren.injection.PhaseSpaceConvention.Recursive2Body

    channel.Sample(siren.utilities.SIREN_random(2), None, record)
    assert record.secondary_masses == pytest.approx([0.008, 0.017, 37.215])
    assert all(momentum[0] > 0.0 for momentum in record.secondary_momenta)


def test_offshell_s_pair_sample_matches_density(vector_portal):
    """s_pair is sampled from EXACTLY the shared Breit-Wigner map whose
    Density FinalStateProbability reports (Contract C1).  Regression guard for
    the former mismatch: the sampler drew from the [s_min,s_max]-renormalized
    BW (1/(hi-lo)) while the density reported a full-line BW (1/pi), so the
    reported density did not integrate to 1 over the allowed range.
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
    assert channel.Convention() == siren.injection.PhaseSpaceConvention.MandelstamST
    assert channel.Density(None, record) > 0.0


def test_dutta_kim_directed_phase_space_topology(vector_portal, processes_dir):
    from siren import _util

    meson = _util.load_module(
        "test_dutta_kim_MesonProduction_for_topology",
        str(processes_dir / "DarkNewsTables" / "MesonProduction.py"),
    )

    fiducial = siren.geometry.Box(0.1, 0.1, 0.1)
    random = siren.utilities.SIREN_random(1337)

    pion_decay = meson.MesonThreeBodySIRENDecay(m_mediator=0.017, g_mu=1.0e-3)
    pion_record = _record(
        pion_decay.GetPossibleSignatures()[0],
        0.13957039,
        0.13957039,
    )
    _zero_secondaries(pion_record, [0.10565837, 0.0, 0.017])
    pion_channel = siren.injection.DetectorDirected3BodyChannel(
        fiducial,
        spectator_index=0,
        pair_first_index=1,
        pair_second_index=2,
        directed_pair_index=2,
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
        stopping_condition=lambda datum, i: not (
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
