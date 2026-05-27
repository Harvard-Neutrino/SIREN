"""
Full Dutta-Kim vector-portal chain at SBND with detector-directed biasing.

Generates weighted events through the 5-vertex on-shell chain:

  1. pi+ -> mu+ nu V1        (3-body meson decay)
  2. V1 -> chi chi            (dark photon 2-body decay)
  3. chi Ar -> chi' Ar        (upscattering via t-channel V2)
  4. chi' -> chi V1_signal    (de-excitation 2-body decay)
  5. V1_signal -> e+ e-       (visible signal)

Each secondary vertex uses detector-directed phase-space biasing toward
the SBND LAr TPC fiducial volume.

Usage:
    python DuttaKim_SBND_full_chain.py [--n-events N] [--seed S]
"""

import argparse
import math
import os

import siren
from siren import _util, dataclasses, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter

# ------------------------------------------------------------------ #
#  Load physics models                                                 #
# ------------------------------------------------------------------ #

_PROC_DIR = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
_MESON = _util.load_module("DuttaKim_MesonProduction",
                           os.path.join(_PROC_DIR, "MesonProduction.py"))
_VP = _util.load_module("DuttaKim_VectorPortal",
                        os.path.join(_PROC_DIR, "VectorPortal.py"))

# ------------------------------------------------------------------ #
#  Constants and particle types                                        #
# ------------------------------------------------------------------ #

M_PION = 0.13957039
M_MUON = 0.10565837

M_CHI = 8e-3
M_CHI_PRIME = 50e-3
M_V1 = 17e-3
M_V2 = 200e-3
M_ARGON40 = 37.215

G_D = 1.0
EPSILON_1 = 7e-5
EPSILON_2 = 1e-4
G_MU = 1e-3

PT = lambda pdg: dataclasses.Particle.ParticleType(pdg)
PION = PT(211)
V1_PROD = PT(5922)
CHI = PT(5917)
CHI_PRIME = PT(5918)
V1_SIGNAL = PT(5923)


def build_models():
    """Construct all physics models for the 5-vertex chain."""
    pion_decay = _MESON.MesonThreeBodySIRENDecay(
        m_mediator=M_V1, g_mu=G_MU,
        pdgid_meson=211, pdgid_lepton=-13,
        pdgid_neutrino=14, pdgid_mediator=5922,
    )
    v1_to_chi = _VP.DarkPhotonToChiDecay(
        M_V1, M_CHI, G_D, pdgid_V1=5922, pdgid_chi=5917,
    )
    upscatter = _VP.VectorPortalUpscatteringXS(
        M_CHI, M_CHI_PRIME, M_V2, G_D, EPSILON_2,
        pdgid_chi=5917, pdgid_chi_prime=5918,
        nuclear_pdgid=1000180400, nuclear_mass=M_ARGON40, A=40, Z=18,
    )
    chi_prime_decay = _VP.ChiPrimeDecay(
        M_CHI, M_CHI_PRIME, M_V1, G_D,
        pdgid_chi_prime=5918, pdgid_chi=5917, pdgid_V1=5923,
    )
    visible_decay = _VP.DarkPhotonDecay(
        M_V1, EPSILON_1, pdgid_V1=5923,
    )
    return pion_decay, v1_to_chi, upscatter, chi_prime_decay, visible_decay


def build_biased_phase_spaces(fiducial, models):
    """Build MultiChannelPhaseSpace for each secondary vertex.

    Each secondary gets a single detector-directed channel that biases
    one daughter toward the fiducial volume.
    """
    _, v1_to_chi, upscatter, chi_prime_decay, visible_decay = models

    def mc(channels, weights):
        m = injection.MultiChannelPhaseSpace()
        m.channels = channels
        m.weights = weights
        return m

    return {
        V1_PROD: {
            v1_to_chi.GetPossibleSignatures()[0]: mc(
                [injection.DetectorDirected2BodyChannel(fiducial, 0)], [1.0]),
        },
        CHI: {
            upscatter.GetPossibleSignatures()[0]: mc(
                [injection.DetectorDirectedScatteringChannel(
                    fiducial, directed_index=0,
                    variable=injection.ScatteringVariable.Q2)], [1.0]),
        },
        CHI_PRIME: {
            chi_prime_decay.GetPossibleSignatures()[0]: mc(
                [injection.DetectorDirected2BodyChannel(fiducial, 1)], [1.0]),
        },
        V1_SIGNAL: {
            visible_decay.GetPossibleSignatures()[0]: mc(
                [injection.Isotropic2BodyChannel(0)], [1.0]),
        },
    }


def stopping_condition(datum, i):
    """Continue processing only particles in the chain.

    Only process chi particles from V1_prod decay (not from chi' decay,
    which would cause infinite recursion).
    """
    sec_type = int(datum.record.signature.secondary_types[i])
    parent_type = int(datum.record.signature.primary_type)
    if sec_type == 5922:       # V1_prod: always process
        return False
    if sec_type == 5917:       # chi: only from V1 -> chi chi
        return parent_type != 5922
    if sec_type == 5918:       # chi': always process
        return False
    if sec_type == 5923:       # V1_signal: always process
        return False
    return True                # everything else: stop


def run(n_events=10, seed=42):
    """Generate and weight events through the full chain."""

    # -- Detector --
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(2.0, 2.0, 2.5)

    # -- Physics --
    models = build_models()
    pion_decay, v1_to_chi, upscatter, chi_prime_decay, visible_decay = models

    # -- Primary distributions (monoenergetic pion along +z from BNB target) --
    primary_mass = distributions.PrimaryMass(M_PION)
    primary_energy = distributions.Monoenergetic(2.0)
    primary_dir = distributions.FixedDirection(siren.math.Vector3D(0, 0, 1))
    primary_pos = distributions.PointSourcePositionDistribution(
        [0.0, 0.0, 0.0], 300.0)

    # -- Secondary distributions --
    sv = distributions.SecondaryPhysicalVertexDistribution()

    # -- Biased phase spaces --
    phase_spaces = build_biased_phase_spaces(fiducial, models)

    # -- Build Injector --
    injector = Injector(
        number_of_events=n_events,
        detector_model=detector_model,
        seed=seed,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_injection_distributions=[
            primary_mass, primary_energy, primary_dir, primary_pos],
        secondary_interactions={
            V1_PROD: [v1_to_chi], CHI: [upscatter],
            CHI_PRIME: [chi_prime_decay], V1_SIGNAL: [visible_decay],
        },
        secondary_injection_distributions={
            V1_PROD: [sv], CHI: [sv], CHI_PRIME: [sv], V1_SIGNAL: [sv],
        },
        secondary_phase_spaces=phase_spaces,
        stopping_condition=stopping_condition,
    )

    # -- Generate events --
    events = []
    for event in injector:
        if event.tree:
            events.append(event)
        if len(events) >= n_events:
            break

    print(f"Generated {len(events)} / {n_events} events")

    # -- Build Weighter --
    weighter = Weighter(
        injectors=[injector],
        detector_model=detector_model,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_physical_distributions=[primary_energy, primary_dir],
        secondary_interactions={
            V1_PROD: [v1_to_chi], CHI: [upscatter],
            CHI_PRIME: [chi_prime_decay], V1_SIGNAL: [visible_decay],
        },
    )

    # -- Compute weights --
    print(f"\n{'Event':>5}  {'Records':>7}  {'Weight':>14}  Chain")
    print("-" * 65)

    for i, event in enumerate(events):
        w = weighter(event)
        records = list(event.tree)
        chain = " -> ".join(
            str(int(d.record.signature.primary_type)) for d in records)
        status = "OK" if math.isfinite(w) and w > 0 else "BAD"
        print(f"{i:5d}  {len(records):7d}  {w:14.4e}  {chain}  [{status}]")

    weights = [weighter(ev) for ev in events]
    finite_weights = [w for w in weights if math.isfinite(w) and w > 0]
    print(f"\n{len(finite_weights)} / {len(events)} events have finite positive weights")

    if finite_weights:
        import numpy as np
        w = np.array(finite_weights)
        print(f"Weight range: [{w.min():.4e}, {w.max():.4e}]")
        print(f"Weight mean:  {w.mean():.4e}")
        print(f"Weight std:   {w.std():.4e}")

    return events, weights


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dutta-Kim vector-portal chain at SBND")
    parser.add_argument("--n-events", type=int, default=10,
                        help="Number of events to generate")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random number seed")
    args = parser.parse_args()
    run(n_events=args.n_events, seed=args.seed)
