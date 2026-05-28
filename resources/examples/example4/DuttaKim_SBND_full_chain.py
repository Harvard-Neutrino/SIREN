"""
Full Dutta-Kim vector-portal chain at SBND with multi-channel biasing.

Supports two chain topologies:

  On-shell (default, 5 vertices):
    pi+ -> mu+ nu V1 -> V1 -> chi chi -> chi Ar -> chi' Ar
    -> chi' -> chi V1_signal -> V1_signal -> e+ e-

  Off-shell (--offshell, 4 vertices, virtual chi'):
    pi+ -> mu+ nu V1 -> V1 -> chi chi -> chi Ar -> chi V1_signal Ar
    -> V1_signal -> e+ e-

Each secondary vertex uses a MultiChannelPhaseSpace combining a
physical channel with a detector-directed biasing channel. The
--optimize flag runs iterative Kleiss-Pittau weight optimization.

Usage:
    python DuttaKim_SBND_full_chain.py [OPTIONS]

Options:
    --dk2nu-dir DIR     Directory with dk2nu ROOT files
    --monoenergetic     Use fixed 2 GeV pion instead of dk2nu
    --offshell          Use off-shell chi' (single 2->3 scattering vertex)
    --n-events N        Number of production events (default: 100)
    --seed S            Random number seed (default: 42)
    --optimize          Run weight optimization before production
    --opt-iterations N  Optimization iterations (default: 5)
    --opt-batch N       Events per optimization iteration (default: 200)
"""

import argparse
import glob
import math
import os
import sys

import numpy as np

import siren
from siren import _util, dataclasses, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter

# ------------------------------------------------------------------ #
#  Load physics models and dk2nu reader                                #
# ------------------------------------------------------------------ #

_PROC_DIR = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
_MESON = _util.load_module("DuttaKim_MesonProduction",
                           os.path.join(_PROC_DIR, "MesonProduction.py"))
_VP = _util.load_module("DuttaKim_VectorPortal",
                        os.path.join(_PROC_DIR, "VectorPortal.py"))
_DK = _util.load_module("DuttaKim_Dk2nuReader",
                        os.path.join(_PROC_DIR, "Dk2nuReader.py"))

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

# chi' decay width (for BreitWigner sampling in off-shell mode)
_P_STAR = injection.TwoBodyRestMomentum(M_CHI_PRIME, M_CHI, M_V1)
CHI_PRIME_WIDTH = G_D**2 * _P_STAR**3 / (6.0 * math.pi * M_CHI_PRIME**2)


def load_dk2nu_pions(dk2nu_dir, detector_model):
    """Read pi+ kinematics from dk2nu files."""
    dk2nu_files = sorted(glob.glob(os.path.join(dk2nu_dir, "*dk2nu*.root")))
    if not dk2nu_files:
        dk2nu_files = sorted(glob.glob(os.path.join(dk2nu_dir, "nubeam*.root")))
    if not dk2nu_files:
        print(f"No dk2nu files found in {dk2nu_dir}")
        sys.exit(1)

    print(f"Reading {len(dk2nu_files)} dk2nu file(s) from {dk2nu_dir}")
    dk2nu_data = _DK.read_dk2nu(dk2nu_files, parent_pdg=[_DK.PTYPE_PIPLUS])
    _DK.print_summary(dk2nu_data)

    n_pions = len(dk2nu_data["E"])
    pot = dk2nu_data["pot"]
    print(f"  {n_pions} pi+ entries, {pot:.3e} POT")

    primary_dist = _DK.dk2nu_to_primary_distribution(
        dk2nu_data, detector_model, parent_pdg=_DK.PTYPE_PIPLUS)
    return primary_dist, pot


# ------------------------------------------------------------------ #
#  On-shell chain (5 vertices)                                         #
# ------------------------------------------------------------------ #

def build_onshell_models():
    """Physics models for the on-shell chain."""
    pion_decay = _MESON.MesonThreeBodySIRENDecay(
        m_mediator=M_V1, g_mu=G_MU,
        pdgid_meson=211, pdgid_lepton=-13,
        pdgid_neutrino=14, pdgid_mediator=5922)
    v1_to_chi = _VP.DarkPhotonToChiDecay(
        M_V1, M_CHI, G_D, pdgid_V1=5922, pdgid_chi=5917)
    upscatter = _VP.VectorPortalUpscatteringXS(
        M_CHI, M_CHI_PRIME, M_V2, G_D, EPSILON_2,
        pdgid_chi=5917, pdgid_chi_prime=5918,
        nuclear_pdgid=1000180400, nuclear_mass=M_ARGON40, A=40, Z=18)
    chi_prime_decay = _VP.ChiPrimeDecay(
        M_CHI, M_CHI_PRIME, M_V1, G_D,
        pdgid_chi_prime=5918, pdgid_chi=5917, pdgid_V1=5923)
    visible_decay = _VP.DarkPhotonDecay(M_V1, EPSILON_1, pdgid_V1=5923)
    return {
        "pion_decay": pion_decay,
        "secondary_interactions": {
            V1_PROD: [v1_to_chi],
            CHI: [upscatter],
            CHI_PRIME: [chi_prime_decay],
            V1_SIGNAL: [visible_decay],
        },
        "models": {
            "v1_to_chi": v1_to_chi,
            "upscatter": upscatter,
            "chi_prime_decay": chi_prime_decay,
            "visible_decay": visible_decay,
        },
    }


def build_onshell_phase_spaces(fiducial, models):
    """Multi-channel phase spaces for on-shell chain."""
    m = models["models"]
    v1_sig = m["v1_to_chi"].GetPossibleSignatures()[0]
    chi_sig = m["upscatter"].GetPossibleSignatures()[0]
    chip_sig = m["chi_prime_decay"].GetPossibleSignatures()[0]
    vis_sig = m["visible_decay"].GetPossibleSignatures()[0]

    return {
        V1_PROD: {v1_sig: _mc([
            injection.PhysicalDecayChannel(m["v1_to_chi"], v1_sig),
            injection.DetectorDirected2BodyChannel(fiducial, 0),
        ], [0.10, 0.90])},
        CHI: {chi_sig: _mc([
            injection.PhysicalCrossSectionChannel(m["upscatter"], chi_sig),
            injection.DetectorDirectedScatteringChannel(
                fiducial, directed_index=0,
                variable=injection.ScatteringVariable.Q2),
        ], [0.10, 0.90])},
        CHI_PRIME: {chip_sig: _mc([
            injection.PhysicalDecayChannel(m["chi_prime_decay"], chip_sig),
            injection.DetectorDirected2BodyChannel(fiducial, 1),
        ], [0.10, 0.90])},
        V1_SIGNAL: {vis_sig: _mc([
            injection.PhysicalDecayChannel(m["visible_decay"], vis_sig),
            injection.Isotropic2BodyChannel(0),
        ], [0.50, 0.50])},
    }


def onshell_stopping_condition(datum, i):
    sec = int(datum.record.signature.secondary_types[i])
    parent = int(datum.record.signature.primary_type)
    if sec == 5922:   return False          # V1_prod
    if sec == 5917:   return parent != 5922 or i != 0 # chi only from V1 decay
    if sec == 5918:   return False          # chi'
    if sec == 5923:   return False          # V1_signal
    return True


# ------------------------------------------------------------------ #
#  Off-shell chain (4 vertices, virtual chi')                          #
# ------------------------------------------------------------------ #

def build_offshell_models():
    """Physics models for the off-shell chain."""
    pion_decay = _MESON.MesonThreeBodySIRENDecay(
        m_mediator=M_V1, g_mu=G_MU,
        pdgid_meson=211, pdgid_lepton=-13,
        pdgid_neutrino=14, pdgid_mediator=5922)
    v1_to_chi = _VP.DarkPhotonToChiDecay(
        M_V1, M_CHI, G_D, pdgid_V1=5922, pdgid_chi=5917)
    offshell_xs = _VP.VectorPortalOffShellXS(
        M_CHI, M_CHI_PRIME, M_V1, M_V2, G_D, EPSILON_2,
        pdgid_V1=5923)
    visible_decay = _VP.DarkPhotonDecay(M_V1, EPSILON_1, pdgid_V1=5923)
    return {
        "pion_decay": pion_decay,
        "secondary_interactions": {
            V1_PROD: [v1_to_chi],
            CHI: [offshell_xs],
            V1_SIGNAL: [visible_decay],
        },
        "models": {
            "v1_to_chi": v1_to_chi,
            "offshell_xs": offshell_xs,
            "visible_decay": visible_decay,
        },
    }


def build_offshell_phase_spaces(fiducial, models):
    """Multi-channel phase spaces for off-shell chain."""
    m = models["models"]
    v1_sig = m["v1_to_chi"].GetPossibleSignatures()[0]
    offshell_sig = m["offshell_xs"].GetPossibleSignatures()[0]
    vis_sig = m["visible_decay"].GetPossibleSignatures()[0]

    return {
        V1_PROD: {v1_sig: _mc([
            injection.PhysicalDecayChannel(m["v1_to_chi"], v1_sig),
            injection.DetectorDirected2BodyChannel(fiducial, 0),
        ], [0.10, 0.90])},
        CHI: {offshell_sig: _mc([
            injection.PhysicalCrossSectionChannel(m["offshell_xs"], offshell_sig),
            injection.DetectorDirected3BodyChannel(
                fiducial,
                spectator_index=2, pair_first_index=0, pair_second_index=1,
                directed_pair_index=0,
                mass_mode=injection.InvariantMassMode.BreitWigner,
                resonance_mass=M_CHI_PRIME,
                resonance_width=CHI_PRIME_WIDTH,
                topology=injection.PhaseSpaceTopology.Scatter2to3),
        ], [0.10, 0.90])},
        V1_SIGNAL: {vis_sig: _mc([
            injection.PhysicalDecayChannel(m["visible_decay"], vis_sig),
            injection.Isotropic2BodyChannel(0),
        ], [0.50, 0.50])},
    }


def offshell_stopping_condition(datum, i):
    sec = int(datum.record.signature.secondary_types[i])
    parent = int(datum.record.signature.primary_type)
    if sec == 5922:   return False          # V1_prod
    if sec == 5917:   return parent != 5922 or i != 0 # chi only from V1 decay (only one of two produced chi)
    if sec == 5923:   return False          # V1_signal
    return True


# ------------------------------------------------------------------ #
#  Shared utilities                                                    #
# ------------------------------------------------------------------ #

def _mc(channels, weights):
    m = injection.MultiChannelPhaseSpace()
    m.channels = channels
    m.weights = weights
    return m


def effective_sample_fraction(weights):
    w = np.asarray(weights)
    w = w[(np.isfinite(w)) & (w > 0)]
    if len(w) == 0:
        return 0.0
    return (w.sum() ** 2) / (len(w) * (w ** 2).sum())


# ------------------------------------------------------------------ #
#  Main                                                                #
# ------------------------------------------------------------------ #

def run(dk2nu_dir, n_events=100, seed=42, optimize=False,
        opt_iterations=5, opt_batch=200, monoenergetic=False,
        offshell=False):

    # -- Detector --
    print("Loading SBND detector model ...")
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(4.0, 4.0, 5.0)

    # -- Chain topology --
    if offshell:
        print("\nUsing off-shell chi' chain (4 vertices, virtual chi')")
        chain = build_offshell_models()
        phase_spaces = build_offshell_phase_spaces(fiducial, chain)
        stop_fn = offshell_stopping_condition
    else:
        print("\nUsing on-shell chain (5 vertices)")
        chain = build_onshell_models()
        phase_spaces = build_onshell_phase_spaces(fiducial, chain)
        stop_fn = onshell_stopping_condition

    pion_decay = chain["pion_decay"]
    secondary_interactions = chain["secondary_interactions"]

    # -- Primary distributions --
    pot = 0.0
    if monoenergetic:
        print("Using monoenergetic 2 GeV pion along +z")
        primary_dists = [
            distributions.PrimaryMass(M_PION),
            distributions.Monoenergetic(2.0),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
            distributions.PointSourcePositionDistribution(
                [0.0, 0.0, 0.0], 300.0),
        ]
        physical_dists = [
            distributions.Monoenergetic(2.0),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        ]
        primary_mode = injection.VertexWeightingMode.Propagated()
    else:
        print()
        pion_dist, pot = load_dk2nu_pions(dk2nu_dir, detector_model)
        primary_dists = [pion_dist]
        physical_dists = [pion_dist]
        primary_mode = injection.VertexWeightingMode.Fixed()

    # -- Secondary distributions --
    sv = distributions.SecondaryPhysicalVertexDistribution()
    sec_dists = {pt: [sv] for pt in secondary_interactions}

    # -- Budget --
    total_budget = n_events
    if optimize:
        total_budget += opt_iterations * opt_batch * 10

    # -- Build Injector --
    print(f"\nBuilding injector ({total_budget} event budget) ...")
    injector = Injector(
        number_of_events=total_budget,
        detector_model=detector_model,
        seed=seed,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_injection_distributions=primary_dists,
        primary_weighting_mode=primary_mode,
        secondary_interactions=secondary_interactions,
        secondary_injection_distributions=sec_dists,
        secondary_phase_spaces=phase_spaces,
        stopping_condition=stop_fn,
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
        primary_physical_distributions=physical_dists,
        secondary_interactions=secondary_interactions,
    )

    # -- Optimize --
    if optimize:
        from siren.optimize import optimize_chain_weights
        print("\nOptimizing multi-channel weights ...")
        optimize_chain_weights(
            injector, weighter,
            n_iterations=opt_iterations,
            batch_size=opt_batch,
            damping=0.5,
            verbose=True,
        )

    # -- Generate production events --
    print(f"\nGenerating {n_events} production events ...")
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

    # -- Print events --
    print(f"\n{'Event':>5}  {'Records':>7}  {'Weight':>14}  Chain")
    print("-" * 70)
    for i, (event, w) in enumerate(zip(events, weights)):
        records = list(sorted(event.tree, key=lambda r: r.depth()))
        chain_str = " -> ".join(
            f"{int(d.record.signature.primary_type)}({d.depth()})" for d in records)
        status = "OK" if valid_mask[i] else "BAD"
        if i < 20 or not valid_mask[i]:
            print(f"{i:5d}  {len(records):7d}  {w:14.4e}  {chain_str}  [{status}]")
    if len(events) > 20:
        print(f"  ... ({len(events) - 20} more events)")

    # -- Summary --
    print(f"\n{'='*50}")
    print(f"Chain:         {'off-shell' if offshell else 'on-shell'}")
    print(f"Valid events:  {valid_mask.sum()} / {len(events)}")
    print(f"Simulated POT: {pot:.3e}")

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
    default_dk2nu = os.environ.get(
        "DK2NU_DIR",
        "/Users/aschneider/workspaces/g4bnb/sources/G4BNB",
    )

    parser = argparse.ArgumentParser(
        description="Dutta-Kim vector-portal chain at SBND")
    parser.add_argument("--dk2nu-dir", type=str, default=default_dk2nu,
                        help="Directory containing dk2nu ROOT files")
    parser.add_argument("--monoenergetic", action="store_true",
                        help="Use monoenergetic 2 GeV pion instead of dk2nu")
    parser.add_argument("--offshell", action="store_true",
                        help="Use off-shell chi' (single 2->3 vertex)")
    parser.add_argument("--n-events", type=int, default=100,
                        help="Number of production events")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random number seed")
    parser.add_argument("--optimize", action="store_true",
                        help="Run iterative weight optimization")
    parser.add_argument("--opt-iterations", type=int, default=5,
                        help="Optimization iterations")
    parser.add_argument("--opt-batch", type=int, default=200,
                        help="Events per optimization iteration")
    args = parser.parse_args()
    run(dk2nu_dir=args.dk2nu_dir, n_events=args.n_events, seed=args.seed,
        optimize=args.optimize, opt_iterations=args.opt_iterations,
        opt_batch=args.opt_batch, monoenergetic=args.monoenergetic,
        offshell=args.offshell)
