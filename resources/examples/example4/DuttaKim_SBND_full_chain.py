"""
Full Dutta-Kim vector-portal chain at SBND with multi-channel biasing.

Generates weighted events through the 5-vertex on-shell chain:

  1. pi+ -> mu+ nu V1        (3-body meson decay)
  2. V1 -> chi chi            (dark photon 2-body decay)
  3. chi Ar -> chi' Ar        (upscattering via t-channel V2)
  4. chi' -> chi V1_signal    (de-excitation 2-body decay)
  5. V1_signal -> e+ e-       (visible signal)

Each secondary vertex uses a MultiChannelPhaseSpace combining a
physical channel (correct angular distribution) with a detector-directed
biasing channel. The --optimize flag runs iterative Kleiss-Pittau
weight optimization to minimize total event weight variance.

Usage:
    python DuttaKim_SBND_full_chain.py [--n-events N] [--seed S] [--optimize]
"""

import argparse
import math
import os

import numpy as np

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


def _mc(channels, weights):
    m = injection.MultiChannelPhaseSpace()
    m.channels = channels
    m.weights = weights
    return m


def build_multichannel_phase_spaces(fiducial, models):
    """Build MultiChannelPhaseSpace for each secondary vertex.

    Each vertex gets a physical channel (correct density, isotropic
    sampling) and a detector-directed biasing channel. The initial
    weights are 10% physical / 90% biased.
    """
    _, v1_to_chi, upscatter, chi_prime_decay, visible_decay = models

    v1_sig = v1_to_chi.GetPossibleSignatures()[0]
    chi_sig = upscatter.GetPossibleSignatures()[0]
    chip_sig = chi_prime_decay.GetPossibleSignatures()[0]
    vis_sig = visible_decay.GetPossibleSignatures()[0]

    return {
        V1_PROD: {
            v1_sig: _mc([
                injection.PhysicalDecayChannel(v1_to_chi, v1_sig),
                injection.DetectorDirected2BodyChannel(fiducial, 0),
            ], [0.10, 0.90]),
        },
        CHI: {
            chi_sig: _mc([
                injection.PhysicalCrossSectionChannel(upscatter, chi_sig),
                injection.DetectorDirectedScatteringChannel(
                    fiducial, directed_index=0,
                    variable=injection.ScatteringVariable.Q2),
            ], [0.10, 0.90]),
        },
        CHI_PRIME: {
            chip_sig: _mc([
                injection.PhysicalDecayChannel(chi_prime_decay, chip_sig),
                injection.DetectorDirected2BodyChannel(fiducial, 1),
            ], [0.10, 0.90]),
        },
        V1_SIGNAL: {
            vis_sig: _mc([
                injection.PhysicalDecayChannel(visible_decay, vis_sig),
                injection.Isotropic2BodyChannel(0),
            ], [0.50, 0.50]),
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


def effective_sample_fraction(weights):
    """Compute ESS / N for an array of weights."""
    w = np.asarray(weights)
    w = w[(np.isfinite(w)) & (w > 0)]
    if len(w) == 0:
        return 0.0
    return (w.sum() ** 2) / (len(w) * (w ** 2).sum())


def run(n_events=100, seed=42, optimize=False, opt_iterations=5,
        opt_batch=200):
    """Generate and weight events through the full chain."""

    # -- Detector --
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(2.0, 2.0, 2.5)

    # -- Physics --
    models = build_models()
    pion_decay, v1_to_chi, upscatter, chi_prime_decay, visible_decay = models

    # -- Primary distributions --
    primary_dists = [
        distributions.PrimaryMass(M_PION),
        distributions.Monoenergetic(2.0),
        distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        distributions.PointSourcePositionDistribution(
            [0.0, 0.0, 0.0], 300.0),
    ]

    # -- Secondary distributions --
    sv = distributions.SecondaryPhysicalVertexDistribution()

    # -- Multi-channel phase spaces --
    phase_spaces = build_multichannel_phase_spaces(fiducial, models)

    # -- Total events budget: optimization batches + production --
    total_budget = n_events
    if optimize:
        total_budget += opt_iterations * opt_batch * 10

    # -- Build Injector --
    injector = Injector(
        number_of_events=total_budget,
        detector_model=detector_model,
        seed=seed,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_injection_distributions=primary_dists,
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

    # Force initialization
    for ev in injector:
        break
    injector._Injector__injector.ResetInjectedEvents()

    # -- Build Weighter --
    weighter = Weighter(
        injectors=[injector],
        detector_model=detector_model,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_physical_distributions=[
            distributions.Monoenergetic(2.0),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        ],
        secondary_interactions={
            V1_PROD: [v1_to_chi], CHI: [upscatter],
            CHI_PRIME: [chi_prime_decay], V1_SIGNAL: [visible_decay],
        },
    )

    # -- Optimize channel weights --
    if optimize:
        from siren.optimize import optimize_chain_weights
        print("Optimizing multi-channel weights...")
        optimize_chain_weights(
            injector, weighter,
            n_iterations=opt_iterations,
            batch_size=opt_batch,
            damping=0.5,
            verbose=True,
        )
        print()

    # -- Generate production events --
    injector._Injector__injector.ResetInjectedEvents()
    events = []
    for event in injector:
        if event.tree:
            events.append(event)
        if len(events) >= n_events:
            break

    print(f"Generated {len(events)} / {n_events} events")

    # -- Compute weights --
    weights = np.array([weighter(event) for event in events])
    valid_mask = np.isfinite(weights) & (weights > 0)
    valid_weights = weights[valid_mask]

    print(f"\n{'Event':>5}  {'Records':>7}  {'Weight':>14}  Chain")
    print("-" * 70)
    for i, (event, w) in enumerate(zip(events, weights)):
        records = list(event.tree)
        chain = " -> ".join(
            str(int(d.record.signature.primary_type)) for d in records)
        status = "OK" if valid_mask[i] else "BAD"
        print(f"{i:5d}  {len(records):7d}  {w:14.4e}  {chain}  [{status}]")

    # -- Summary statistics --
    print(f"\n{'='*50}")
    print(f"Valid events: {valid_mask.sum()} / {len(events)}")

    if len(valid_weights) > 0:
        eff = effective_sample_fraction(valid_weights) * 100
        cv = valid_weights.std() / valid_weights.mean()
        print(f"Weight range:  [{valid_weights.min():.4e}, "
              f"{valid_weights.max():.4e}]")
        print(f"Weight mean:   {valid_weights.mean():.4e}")
        print(f"Weight CV:     {cv:.2f}")
        print(f"Eff. sample:   {eff:.1f}%")

    return events, weights


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dutta-Kim vector-portal chain at SBND")
    parser.add_argument("--n-events", type=int, default=100,
                        help="Number of production events")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random number seed")
    parser.add_argument("--optimize", action="store_true",
                        help="Run iterative weight optimization before production")
    parser.add_argument("--opt-iterations", type=int, default=5,
                        help="Number of optimization iterations")
    parser.add_argument("--opt-batch", type=int, default=200,
                        help="Events per optimization iteration")
    args = parser.parse_args()
    run(n_events=args.n_events, seed=args.seed, optimize=args.optimize,
        opt_iterations=args.opt_iterations, opt_batch=args.opt_batch)
