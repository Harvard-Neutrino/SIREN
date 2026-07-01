"""Diagnose why ~10% of chi injection attempts fail.

Captures partial trees from failed events and classifies each failure:
  - dead V1 direction: chi collimation cone doesn't reach Argon
  - dead pion: no V1 from this pion's 3-body decay can reach Argon
  - wrong channel: V1 was directed away from Argon but pion was viable
  - marginal miss: chi barely missed the Argon volume

Usage:
    python diagnose_chi_failures.py [--n-events N]
"""
import argparse
import math
import os
import sys

import numpy as np

import siren
from siren import _util, dataclasses, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter

from DuttaKim_SBND_full_chain import (
    build_offshell_models, build_offshell_phase_spaces,
    build_geometric_targets, build_primary_phase_spaces,
    offshell_stopping_condition, load_dk2nu_pions, default_pion_bias,
    _mc, M_PION, M_V1, M_CHI, PION, GAMMA_PION_SM, M_MUON,
)

ARGON_CENTER = np.array([0.0, 0.0, 0.0])
ARGON_HALF = np.array([5.2, 2.9, 4.3])
ARGON_BOUNDING_RADIUS = np.linalg.norm(ARGON_HALF)

PDG_NAMES = {
    211: "pi+", 5922: "V1_prod", 5917: "chi", 5918: "chi'",
    5923: "V1_sig", -13: "mu+", 14: "nu_mu",
}


def normalize(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def angle_between(a, b):
    """Angle in radians between two unit vectors."""
    dot = np.clip(np.dot(a, b), -1.0, 1.0)
    return np.arccos(dot)


def argon_half_angle(position):
    """Half-angle subtended by the Argon bounding sphere from position."""
    dist = np.linalg.norm(position - ARGON_CENTER)
    if dist < ARGON_BOUNDING_RADIUS:
        return math.pi
    return math.atan2(ARGON_BOUNDING_RADIUS, dist)


def max_chi_opening_angle(V1_momentum):
    """Max lab-frame opening angle of chi from V1 direction.

    For V1 -> chi chi, the rest-frame momentum is:
      p* = sqrt(M_V1^2/4 - m_chi^2)
    The max lab angle is:
      sin(theta_max) = p* / (gamma * beta * E*_chi)
    where E*_chi = M_V1/2.
    """
    E_V1 = V1_momentum[0]
    px, py, pz = V1_momentum[1], V1_momentum[2], V1_momentum[3]
    p_V1 = math.sqrt(px**2 + py**2 + pz**2)

    if E_V1 <= M_V1 or p_V1 < 1e-15:
        return math.pi

    gamma = E_V1 / M_V1
    beta = p_V1 / E_V1

    E_chi_rf = M_V1 / 2.0
    p_chi_rf = math.sqrt(max(E_chi_rf**2 - M_CHI**2, 0.0))

    if gamma * beta * E_chi_rf < 1e-15:
        return math.pi

    sin_max = p_chi_rf / (gamma * beta * E_chi_rf)
    if sin_max >= 1.0:
        return math.pi
    return math.asin(sin_max)


def max_V1_opening_angle(pion_momentum):
    """Max lab-frame opening angle of V1 from pion direction.

    In the pion rest frame, V1 can go in any direction. The maximum
    lab-frame angle is determined by the V1's rest-frame momentum and
    the pion boost.

    For pi -> mu nu V1, V1's maximum rest-frame momentum occurs when
    the (mu,nu) system has minimum invariant mass = m_mu (neutrino massless).
    Then E_V1_max = (m_pi^2 + m_V1^2 - m_mu^2) / (2*m_pi),
    p_V1_max = sqrt(E_V1_max^2 - m_V1^2).
    """
    E_pi = pion_momentum[0]
    px, py, pz = pion_momentum[1], pion_momentum[2], pion_momentum[3]
    p_pi = math.sqrt(px**2 + py**2 + pz**2)

    if E_pi <= M_PION or p_pi < 1e-15:
        return math.pi

    gamma = E_pi / M_PION
    beta = p_pi / E_pi

    E_V1_max_rf = (M_PION**2 + M_V1**2 - M_MUON**2) / (2.0 * M_PION)
    p_V1_max_rf = math.sqrt(max(E_V1_max_rf**2 - M_V1**2, 0.0))

    if gamma * beta * E_V1_max_rf < 1e-15:
        return math.pi

    sin_max = p_V1_max_rf / (gamma * beta * E_V1_max_rf)
    if sin_max >= 1.0:
        return math.pi
    return math.asin(sin_max)


def extract_kinematics(tree_data):
    """Extract pion and V1/chi kinematics from a tree's datum list."""
    result = {}

    for datum in tree_data:
        r = datum.record
        pid = int(r.signature.primary_type)

        if datum.depth() == 0:
            result["pion_pos"] = np.array(r.interaction_vertex[:3])
            result["pion_mom"] = np.array(r.primary_momentum)
            for si, stype in enumerate(r.signature.secondary_types):
                if int(stype) == 5922:
                    result["V1_mom"] = np.array(r.secondary_momenta[si])

        elif pid == 5922:
            result["V1_decay_pos"] = np.array(r.interaction_vertex[:3])
            for si, stype in enumerate(r.signature.secondary_types):
                if int(stype) == 5917:
                    result["chi_mom"] = np.array(r.secondary_momenta[si])
                    break

    return result


def classify_failure(kin):
    """Classify a failed event based on its kinematics.

    Returns one of: 'dead_pion', 'wrong_channel', 'marginal_miss', 'unknown'
    """
    if "chi_mom" not in kin or "V1_mom" not in kin:
        return "unknown"

    chi_pos = kin.get("V1_decay_pos", kin["pion_pos"])
    chi_dir = normalize(kin["chi_mom"][1:4])
    to_argon = normalize(ARGON_CENTER - chi_pos)

    alpha = angle_between(chi_dir, to_argon)
    theta_argon = argon_half_angle(chi_pos)
    theta_chi_max = max_chi_opening_angle(kin["V1_mom"])

    if alpha <= theta_argon + theta_chi_max:
        return "marginal_miss"

    # V1 direction can't produce chi hitting Argon.
    # Now check: could ANY V1 direction from this pion work?
    pion_dir = normalize(kin["pion_mom"][1:4])
    pion_to_argon = angle_between(pion_dir, to_argon)
    theta_V1_max = max_V1_opening_angle(kin["pion_mom"])

    # V1 can deviate up to theta_V1_max from pion direction.
    # Chi can deviate up to theta_chi_max from V1 direction.
    # Total max deviation from pion direction: theta_V1_max + theta_chi_max
    # The Argon subtends theta_argon from chi's position.
    # Pion needs to be within theta_V1_max + theta_chi_max + theta_argon
    # of the Argon direction for ANY viable V1+chi.
    if pion_to_argon > theta_V1_max + theta_chi_max + theta_argon:
        return "dead_pion"
    else:
        return "wrong_channel"


def run_diagnostic(n_events, dk2nu_dir):
    print("Loading detector...")
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")

    fiducial = siren.geometry.Box((4.5 + 2.022) * 2, 4.074645 * 2, 5.010 * 2)
    targets = build_geometric_targets(detector_model, fiducial)

    chain = build_offshell_models()
    pion_decay = chain["pion_decay"]
    secondary_interactions = chain["secondary_interactions"]

    pion_dist, pot = load_dk2nu_pions(dk2nu_dir, detector_model,
                                       sampling_bias=default_pion_bias)

    primary_mode = injection.VertexWeightingMode.Fixed()
    sv = distributions.SecondaryPhysicalVertexDistribution()
    sec_dists = {pt: [sv] for pt in secondary_interactions}
    phase_spaces = build_offshell_phase_spaces(targets, chain)
    primary_ps = build_primary_phase_spaces(targets, pion_decay)

    injector = Injector(
        number_of_events=n_events,
        detector_model=detector_model,
        seed=42,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_injection_distributions=[pion_dist],
        primary_weighting_mode=primary_mode,
        secondary_interactions=secondary_interactions,
        secondary_injection_distributions=sec_dists,
        secondary_phase_spaces=phase_spaces,
        primary_phase_spaces=primary_ps,
        stopping_condition=offshell_stopping_condition,
    )

    for ev in injector:
        break
    cpp_inj = injector._Injector__injector
    cpp_inj.ResetInjectedEvents(n_events)

    # Collect kinematics for success and failure
    success_kin = []
    failure_kin = []
    classifications = {"dead_pion": 0, "wrong_channel": 0,
                       "marginal_miss": 0, "unknown": 0}

    print(f"\nGenerating {n_events} events...")
    for i in range(n_events):
        event = cpp_inj.GenerateEvent()
        if event.tree:
            kin = extract_kinematics(event.tree)
            kin["success"] = True
            success_kin.append(kin)
        else:
            failed_tree = cpp_inj.GetLastFailedTree()
            if failed_tree.tree:
                kin = extract_kinematics(failed_tree.tree)
                kin["success"] = False
                cat = classify_failure(kin)
                kin["category"] = cat
                classifications[cat] += 1
                failure_kin.append(kin)

    n_success = len(success_kin)
    n_fail = len(failure_kin)
    print(f"\n{'='*60}")
    print(f"Results: {n_success} success, {n_fail} failures "
          f"({100*n_fail/(n_success+n_fail):.1f}% failure rate)")
    print(f"{'='*60}")

    print(f"\nFailure classification:")
    for cat, count in sorted(classifications.items(), key=lambda x: -x[1]):
        if count > 0:
            pct = 100 * count / n_fail if n_fail > 0 else 0
            print(f"  {cat:20s}: {count:4d} ({pct:5.1f}%)")

    # Angular analysis
    print(f"\n--- Angular analysis ---")
    if failure_kin:
        fail_alphas = []
        fail_theta_argon = []
        fail_theta_chi = []
        for kin in failure_kin:
            if "chi_mom" not in kin:
                continue
            chi_pos = kin.get("V1_decay_pos", kin["pion_pos"])
            chi_dir = normalize(kin["chi_mom"][1:4])
            to_argon = normalize(ARGON_CENTER - chi_pos)
            alpha = angle_between(chi_dir, to_argon)
            fail_alphas.append(np.degrees(alpha))
            fail_theta_argon.append(np.degrees(argon_half_angle(chi_pos)))
            fail_theta_chi.append(np.degrees(max_chi_opening_angle(kin["V1_mom"])))

        fail_alphas = np.array(fail_alphas)
        fail_theta_argon = np.array(fail_theta_argon)
        fail_theta_chi = np.array(fail_theta_chi)
        print(f"Failed events (n={len(fail_alphas)}):")
        print(f"  Chi angular offset from Argon: "
              f"mean={fail_alphas.mean():.2f} deg, "
              f"range=[{fail_alphas.min():.2f}, {fail_alphas.max():.2f}] deg")
        print(f"  Argon half-angle:              "
              f"mean={fail_theta_argon.mean():.2f} deg")
        print(f"  Max chi deflection from V1:    "
              f"mean={fail_theta_chi.mean():.3f} deg")

    if success_kin:
        succ_alphas = []
        for kin in success_kin:
            if "chi_mom" not in kin:
                continue
            chi_pos = kin.get("V1_decay_pos", kin["pion_pos"])
            chi_dir = normalize(kin["chi_mom"][1:4])
            to_argon = normalize(ARGON_CENTER - chi_pos)
            alpha = angle_between(chi_dir, to_argon)
            succ_alphas.append(np.degrees(alpha))

        succ_alphas = np.array(succ_alphas)
        print(f"Successful events (n={len(succ_alphas)}):")
        print(f"  Chi angular offset from Argon: "
              f"mean={succ_alphas.mean():.2f} deg, "
              f"range=[{succ_alphas.min():.2f}, {succ_alphas.max():.2f}] deg")

    # Pion energy comparison
    print(f"\n--- Pion energy comparison ---")
    if success_kin:
        succ_E = np.array([k["pion_mom"][0] for k in success_kin
                           if "pion_mom" in k])
        print(f"Success: E_pi mean={succ_E.mean():.3f} GeV, "
              f"range=[{succ_E.min():.3f}, {succ_E.max():.3f}] GeV")
    if failure_kin:
        fail_E = np.array([k["pion_mom"][0] for k in failure_kin
                           if "pion_mom" in k])
        print(f"Failure: E_pi mean={fail_E.mean():.3f} GeV, "
              f"range=[{fail_E.min():.3f}, {fail_E.max():.3f}] GeV")

    # Pion position comparison
    print(f"\n--- Pion decay position comparison ---")
    if success_kin:
        succ_z = np.array([k["pion_pos"][2] for k in success_kin
                           if "pion_pos" in k])
        succ_r = np.array([np.sqrt(k["pion_pos"][0]**2 + k["pion_pos"][1]**2)
                           for k in success_kin if "pion_pos" in k])
        print(f"Success: z mean={succ_z.mean():.1f} m, "
              f"r_trans mean={succ_r.mean():.2f} m")
    if failure_kin:
        fail_z = np.array([k["pion_pos"][2] for k in failure_kin
                           if "pion_pos" in k])
        fail_r = np.array([np.sqrt(k["pion_pos"][0]**2 + k["pion_pos"][1]**2)
                           for k in failure_kin if "pion_pos" in k])
        print(f"Failure: z mean={fail_z.mean():.1f} m, "
              f"r_trans mean={fail_r.mean():.2f} m")

    # V1 direction analysis for wrong_channel events
    wrong_channel = [k for k in failure_kin
                     if k.get("category") == "wrong_channel"]
    if wrong_channel:
        print(f"\n--- Wrong channel events (n={len(wrong_channel)}) ---")
        for i, kin in enumerate(wrong_channel[:5]):
            chi_pos = kin.get("V1_decay_pos", kin["pion_pos"])
            V1_dir = normalize(kin["V1_mom"][1:4])
            to_argon = normalize(ARGON_CENTER - chi_pos)
            V1_alpha = np.degrees(angle_between(V1_dir, to_argon))
            theta_argon = np.degrees(argon_half_angle(chi_pos))
            print(f"  Event {i}: V1 offset={V1_alpha:.2f} deg, "
                  f"Argon half-angle={theta_argon:.2f} deg, "
                  f"gap={V1_alpha - theta_argon:.2f} deg")


if __name__ == "__main__":
    default_dk2nu = os.environ.get(
        "DK2NU_DIR", "/Users/aschneider/workspaces/g4bnb/sources/G4BNB")

    parser = argparse.ArgumentParser(description="Diagnose chi injection failures")
    parser.add_argument("--dk2nu-dir", default=default_dk2nu)
    parser.add_argument("--n-events", type=int, default=2000)
    args = parser.parse_args()

    run_diagnostic(args.n_events, args.dk2nu_dir)
