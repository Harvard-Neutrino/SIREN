"""
Small end-to-end Dutta-Kim vector-portal injection smoke example.

This keeps the physical topology from the paper but uses a tiny in-memory
primary pion sample instead of dk2nu ROOT files:

    pi+ -> mu+ nu V1
    V1 -> chi chi
    chi Ar -> chi' Ar
    chi' -> chi V1_signal
    V1_signal -> e+ e-

The example is meant to verify SIREN plumbing: in-memory external primaries,
primary and secondary phase-space channels, on-shell upscattering, and
recursive secondary injection. It uses the bundled IceCube detector and O16
targets so it can run without downloading SBN GDML files.
"""

import math
import os

import siren
from siren import _util as _siren_util


N_EVENTS = int(os.environ.get("SIREN_DK_N_EVENTS", "3"))
SOURCE_Z = float(os.environ.get("SIREN_DK_SOURCE_Z", "-10000000.0"))

M_CHI = 8e-3
M_CHI_PRIME = 50e-3
M_V1 = 17e-3
M_V2 = 200e-3
M_PION = 0.13957039
M_MUON = 0.10565837
M_TARGET = 14.89916863650645

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
PDGID_TARGET = 1000080160
TARGET_A = 16
TARGET_Z = 8


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


def particle_type(pdgid):
    return siren.dataclasses.Particle.ParticleType(pdgid)


def _single_channel(channel):
    phase_space = siren.injection.MultiChannelPhaseSpace()
    phase_space.channels = [channel]
    phase_space.weights = [1.0]
    return phase_space


def _pion_source():
    energy = 5.0
    pz = math.sqrt(energy * energy - M_PION * M_PION)
    return siren.distributions.PrimaryExternalDistribution(
        ["E", "m", "x", "y", "z", "px", "py", "pz"],
        [
            [energy, M_PION, 0.0, 0.0, SOURCE_Z, 0.0, 0.0, pz],
            [6.0, M_PION, 0.0, 0.0, SOURCE_Z, 0.0, 0.0,
             math.sqrt(6.0 * 6.0 - M_PION * M_PION)],
        ],
    )


def build_injector():
    detector_model = siren.load_detector("IceCube")
    fiducial = siren.get_fiducial_volume("IceCube")

    pion_type = particle_type(PDGID_PION)
    v1_prod_type = particle_type(PDGID_V1_PROD)
    chi_type = particle_type(PDGID_CHI)
    chi_prime_type = particle_type(PDGID_CHI_PRIME)
    v1_signal_type = particle_type(PDGID_V1_SIGNAL)

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
        nuclear_pdgid=PDGID_TARGET,
        nuclear_mass=M_TARGET,
        nuclear_name="O16",
        A=TARGET_A,
        Z=TARGET_Z,
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

    primary_phase_spaces = {
        pion_decay.GetPossibleSignatures()[0]: _single_channel(
            siren.injection.DetectorDirected3BodyChannel(
                fiducial,
                spectator_index=0,
                pair_first_index=1,
                pair_second_index=2,
                directed_pair_index=2,
                mass_mode=siren.injection.InvariantMassMode.Uniform,
            )
        )
    }
    secondary_phase_spaces = {
        v1_prod_type: {
            v1_to_chi.GetPossibleSignatures()[0]: _single_channel(
                siren.injection.DetectorDirected2BodyChannel(fiducial, 0)
            )
        },
        chi_type: {
            upscatter.GetPossibleSignatures()[0]: _single_channel(
                siren.injection.DetectorDirectedScatteringChannel(
                    fiducial,
                    directed_index=0,
                    variable=siren.injection.ScatteringVariable.Q2,
                )
            )
        },
        chi_prime_type: {
            chi_prime_decay.GetPossibleSignatures()[0]: _single_channel(
                siren.injection.DetectorDirected2BodyChannel(fiducial, 1)
            )
        },
        v1_signal_type: {
            visible_decay.GetPossibleSignatures()[0]: _single_channel(
                siren.injection.Isotropic2BodyChannel(0)
            )
        },
    }

    secondary_interactions = {
        v1_prod_type: [v1_to_chi],
        chi_type: [upscatter],
        chi_prime_type: [chi_prime_decay],
        v1_signal_type: [visible_decay],
    }
    secondary_distributions = {
        v1_prod_type: [siren.distributions.SecondaryBoundedVertexDistribution()],
        chi_type: [siren.distributions.SecondaryBoundedVertexDistribution(fiducial)],
        chi_prime_type: [siren.distributions.SecondaryBoundedVertexDistribution()],
        v1_signal_type: [siren.distributions.SecondaryBoundedVertexDistribution()],
    }

    def stop(tree, datum, index):
        secondary_type = datum.record.signature.secondary_types[index]
        parent_type = datum.record.signature.primary_type
        if parent_type == v1_prod_type and secondary_type == chi_type:
            return index != 0
        if parent_type == chi_prime_type and secondary_type == chi_type:
            return True
        if secondary_type in {
            v1_prod_type,
            chi_type,
            chi_prime_type,
            v1_signal_type,
        }:
            return False
        return True

    injector = siren.injection.Injector(
        events=N_EVENTS,
        detector_model=detector_model,
        seed=12345,
        primary_type=pion_type,
        primary_interactions=[pion_decay],
        primary_injection_distributions=[_pion_source()],
        primary_phase_spaces=primary_phase_spaces,
        secondary_interactions=secondary_interactions,
        secondary_injection_distributions=secondary_distributions,
        secondary_phase_spaces=secondary_phase_spaces,
        stopping_condition=stop,
    )
    return injector


def summarize(events):
    counts = {}
    for event in events:
        for datum in event.tree:
            key = int(datum.record.signature.primary_type)
            counts[key] = counts.get(key, 0) + 1
    return counts


if __name__ == "__main__":
    injector = build_injector()
    events = list(injector)
    counts = summarize(events)
    nonempty = sum(1 for event in events if len(event.tree) > 0)
    print(f"Generated {nonempty}/{len(events)} non-empty events")
    for pdgid, count in sorted(counts.items()):
        print(f"primary {pdgid}: {count}")
