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
