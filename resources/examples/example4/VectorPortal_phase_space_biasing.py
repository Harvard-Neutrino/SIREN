"""
Detector-directed phase-space setup for the Dutta-Kim vector-portal chain.

This example is independent of dk2nu ROOT files. It builds the SIREN
processes and multi-channel phase spaces for:

    pi+ -> mu+ nu V1
    V1 -> chi chi
    chi Ar -> chi' Ar
    chi' -> chi V1_signal
    V1_signal -> e+ e-

The pion production channel is sampled manually because current SIREN
phase-space biasing is attached to secondary injection processes. The
returned secondary_phase_spaces dictionary can be passed to siren.Injector
for a full meson-source injection once a primary pion distribution is chosen.
"""

import math
import os

import siren
from siren import _util as _siren_util


_DARK_NEWS_TABLES = os.path.join(
    _siren_util.resource_package_dir(),
    "processes",
    "DarkNewsTables",
)
_MESON = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.MesonProduction",
    os.path.join(_DARK_NEWS_TABLES, "MesonProduction.py"),
)
_VECTOR = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.VectorPortal",
    os.path.join(_DARK_NEWS_TABLES, "VectorPortal.py"),
)


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


def particle_type(pdgid):
    return siren.dataclasses.Particle.ParticleType(pdgid)


def build_processes():
    pion_decay = _MESON.MesonThreeBodySIRENDecay(
        m_meson=M_PION,
        m_lepton=M_MUON,
        m_mediator=M_V1,
        g_mu=G_MU,
        pdgid_meson=PDGID_PION,
        pdgid_lepton=PDGID_MUPLUS,
        pdgid_neutrino=PDGID_NUMU,
        pdgid_mediator=PDGID_V1_PROD,
    )

    v1_to_chi = _VECTOR.DarkPhotonToChiDecay(
        M_V1,
        M_CHI,
        G_D,
        pdgid_V1=PDGID_V1_PROD,
        pdgid_chi=PDGID_CHI,
    )

    upscatter = _VECTOR.VectorPortalUpscatteringXS(
        M_CHI,
        M_CHI_PRIME,
        M_V2,
        G_D,
        EPSILON_2,
        pdgid_chi=PDGID_CHI,
        pdgid_chi_prime=PDGID_CHI_PRIME,
        nuclear_pdgid=PDGID_ARGON40,
        nuclear_mass=M_ARGON40,
        nuclear_name="Ar40",
        A=40,
        Z=18,
    )

    chi_prime_decay = _VECTOR.ChiPrimeDecay(
        M_CHI,
        M_CHI_PRIME,
        M_V1,
        G_D,
        pdgid_chi_prime=PDGID_CHI_PRIME,
        pdgid_chi=PDGID_CHI,
        pdgid_V1=PDGID_V1_SIGNAL,
    )

    visible_decay = _VECTOR.DarkPhotonDecay(
        M_V1,
        EPSILON_1,
        pdgid_V1=PDGID_V1_SIGNAL,
    )

    secondary_interactions = {
        particle_type(PDGID_V1_PROD): [v1_to_chi],
        particle_type(PDGID_CHI): [upscatter],
        particle_type(PDGID_CHI_PRIME): [chi_prime_decay],
        particle_type(PDGID_V1_SIGNAL): [visible_decay],
    }
    return pion_decay, secondary_interactions


def _single_channel(channel):
    phase_space = siren.injection.MultiChannelPhaseSpace()
    phase_space.channels = [channel]
    phase_space.weights = [1.0]
    return phase_space


def build_secondary_phase_spaces(fiducial):
    return {
        particle_type(PDGID_V1_PROD): _single_channel(
            siren.injection.DetectorDirected2BodyChannel(fiducial, 0)
        ),
        particle_type(PDGID_CHI): _single_channel(
            siren.injection.DetectorDirectedScatteringChannel(
                fiducial,
                directed_index=0,
                variable=siren.injection.ScatteringVariable.Q2,
            )
        ),
        particle_type(PDGID_CHI_PRIME): _single_channel(
            siren.injection.DetectorDirected2BodyChannel(fiducial, 1)
        ),
        particle_type(PDGID_V1_SIGNAL): _single_channel(
            siren.injection.Isotropic2BodyChannel(0)
        ),
    }


def _record(signature, primary_mass, primary_energy, vertex=(0.0, 0.0, -110.0)):
    record = siren.dataclasses.InteractionRecord()
    record.signature = signature
    record.primary_mass = primary_mass
    pz = math.sqrt(max(primary_energy**2 - primary_mass**2, 0.0))
    record.primary_momentum = [primary_energy, 0.0, 0.0, pz]
    record.interaction_vertex = list(vertex)
    record.secondary_momenta = []
    record.secondary_masses = []
    record.secondary_helicities = []
    return record


def _set_secondaries(record, masses):
    record.secondary_masses = list(masses)
    record.secondary_momenta = [[0.0, 0.0, 0.0, 0.0] for _ in masses]
    record.secondary_helicities = [0 for _ in masses]
    return record


def sample_representative_chain(fiducial):
    pion_decay, secondary_interactions = build_processes()
    phase_spaces = build_secondary_phase_spaces(fiducial)
    random = siren.utilities.SIREN_random(12345)
    densities = {}

    pion_record = _record(pion_decay.GetPossibleSignatures()[0], M_PION, M_PION)
    _set_secondaries(pion_record, [M_MUON, 0.0, M_V1])
    pion_channel = siren.injection.DetectorDirected3BodyChannel(
        fiducial,
        spectator_index=0,
        pair_first_index=1,
        pair_second_index=2,
        directed_pair_index=2,
        mass_mode=siren.injection.InvariantMassMode.Uniform,
    )
    pion_channel.Sample(random, None, pion_record)
    densities["pi -> mu nu V1"] = pion_channel.Density(None, pion_record)

    v1_decay = secondary_interactions[particle_type(PDGID_V1_PROD)][0]
    v1_record = _record(v1_decay.GetPossibleSignatures()[0], M_V1, 0.08)
    _set_secondaries(v1_record, [M_CHI, M_CHI])
    phase_spaces[particle_type(PDGID_V1_PROD)].Sample(random, None, v1_record)
    densities["V1 -> chi chi"] = phase_spaces[particle_type(PDGID_V1_PROD)].Density(
        None,
        v1_record,
    )

    upscatter = secondary_interactions[particle_type(PDGID_CHI)][0]
    scatter_record = _record(upscatter.GetPossibleSignatures()[0], M_CHI, 1.0)
    scatter_record.target_mass = M_ARGON40
    _set_secondaries(
        scatter_record,
        upscatter.SecondaryMasses(scatter_record.signature.secondary_types),
    )
    phase_spaces[particle_type(PDGID_CHI)].Sample(random, None, scatter_record)
    densities["chi Ar -> chi' Ar"] = phase_spaces[particle_type(PDGID_CHI)].Density(
        None,
        scatter_record,
    )

    chi_prime_decay = secondary_interactions[particle_type(PDGID_CHI_PRIME)][0]
    chi_prime_record = _record(
        chi_prime_decay.GetPossibleSignatures()[0],
        M_CHI_PRIME,
        0.12,
    )
    _set_secondaries(chi_prime_record, [M_CHI, M_V1])
    phase_spaces[particle_type(PDGID_CHI_PRIME)].Sample(
        random,
        None,
        chi_prime_record,
    )
    densities["chi' -> chi V1_signal"] = phase_spaces[
        particle_type(PDGID_CHI_PRIME)
    ].Density(None, chi_prime_record)

    visible_decay = secondary_interactions[particle_type(PDGID_V1_SIGNAL)][0]
    visible_record = _record(visible_decay.GetPossibleSignatures()[0], M_V1, 0.06)
    _set_secondaries(visible_record, [0.000511, 0.000511])
    phase_spaces[particle_type(PDGID_V1_SIGNAL)].Sample(random, None, visible_record)
    densities["V1_signal -> e+ e-"] = phase_spaces[
        particle_type(PDGID_V1_SIGNAL)
    ].Density(None, visible_record)

    return densities


if __name__ == "__main__":
    fiducial_volume = siren.geometry.Box(widths=(4.0, 4.0, 5.0))
    for label, density in sample_representative_chain(fiducial_volume).items():
        print(f"{label}: density = {density:.6e}")
