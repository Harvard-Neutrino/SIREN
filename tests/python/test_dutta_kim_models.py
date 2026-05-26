"""Regression tests for the self-contained Dutta-Kim vector portal models."""

import math

import pytest

siren = pytest.importorskip("siren")


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
    assert primary[chi_type][0].DensityVariables() == ["Q2"]
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
    assert channel.Convention() == siren.injection.PhaseSpaceConvention.MandelstamST

    channel.Sample(siren.utilities.SIREN_random(2), None, record)

    assert channel.Density(None, record) == pytest.approx(1.0)
    assert record.secondary_masses == pytest.approx([0.008, 0.017, 37.215])
    assert all(momentum[0] > 0.0 for momentum in record.secondary_momenta)


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
