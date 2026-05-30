"""Diagnose on-shell chain efficiency: where does the V1_sig vertex land?

The on-shell chain has 5 vertices:
    pi+ -> mu+ nu V1   (primary)
    V1  -> chi chi      (V1_prod, 5922)
    chi + Ar -> chi' Ar (upscatter, 5917 -> 5918)
    chi' -> chi V1_sig  (chi' decay, 5918)
    V1_sig -> e+ e-     (visible decay, 5923)

Unlike diagnose_chi_failures.py (which studies why chi injection *fails*),
this script studies the *completed* on-shell events and asks where the
signal vertex (V1_sig -> e+ e-, PDG 5923 interaction_vertex) lands
relative to the fiducial box.  The fiducial-volume metric only counts
events whose V1_sig vertex is inside the box, so events whose signal
vertex leaks outside are dead weight that depress effective sample size.

It also classifies injection failures by how far down the chain they got.

Usage:
    python diagnose_onshell_failures.py [--n-events N]
"""
import argparse
import os
import sys

import numpy as np

import siren
from siren import _util, dataclasses, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter

from DuttaKim_SBND_full_chain import (
    build_onshell_models, build_onshell_phase_spaces,
    build_geometric_targets, build_primary_phase_spaces,
    onshell_stopping_condition, load_dk2nu_pions, default_pion_bias,
    make_fiducial_metric, effective_sample_fraction,
    GAMMA_PION_SM, PION, CHI, CHI_PRIME, V1_PROD, V1_SIGNAL,
)

# Fiducial box half-widths (detector coordinates, metres).
FID_HALF = np.array([4.5 + 2.022, 4.074645, 5.010])

PDG_NAMES = {
    211: "pi+", 5922: "V1_prod", 5917: "chi", 5918: "chi'",
    5923: "V1_sig", -13: "mu+", 14: "nu_mu", 11: "e-", -11: "e+",
}


def box_signed_offset(pos):
    """Per-axis distance outside the fiducial box (0 if inside on that axis)."""
    return np.maximum(np.abs(pos) - FID_HALF, 0.0)


def inside_box(pos):
    return np.all(np.abs(pos) <= FID_HALF)


def vertex_of(tree, pdg):
    """interaction_vertex of the (first) datum whose primary_type == pdg."""
    for datum in tree:
        r = datum.record
        if int(r.signature.primary_type) == pdg:
            return np.array(r.interaction_vertex[:3])
    return None


def depth_reached(tree):
    """Return the set of primary PDGs present in a (partial) tree."""
    return set(int(d.record.signature.primary_type) for d in tree)


def classify_failure(tree):
    """Classify how far an injection failure got down the chain."""
    pdgs = depth_reached(tree)
    if not pdgs:
        return "empty"
    if 5923 in pdgs:
        return "reached_V1sig"      # got to signal but final decay failed
    if 5918 in pdgs:
        return "reached_chiprime"   # upscattered but chi' decay failed
    if 5917 in pdgs:
        return "reached_chi"        # V1->chi ok but upscatter failed
    if 5922 in pdgs:
        return "reached_V1prod"     # pion decayed but V1->chi failed
    if 211 in pdgs:
        return "reached_pion"
    return "other"


def run_diagnostic(n_events, dk2nu_dir):
    print("Loading detector...")
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")

    fiducial = siren.geometry.Box(FID_HALF[0] * 2, FID_HALF[1] * 2, FID_HALF[2] * 2)
    targets = build_geometric_targets(detector_model, fiducial)

    chain = build_onshell_models()
    pion_decay = chain["pion_decay"]
    secondary_interactions = chain["secondary_interactions"]

    pion_dist, pot = load_dk2nu_pions(dk2nu_dir, detector_model,
                                      sampling_bias=default_pion_bias)

    primary_mode = injection.VertexWeightingMode.Fixed()
    sv = distributions.SecondaryPhysicalVertexDistribution()
    sv_bounded = distributions.SecondaryBoundedVertexDistribution(fiducial)
    sec_dists = {pt: [sv] for pt in secondary_interactions}
    sec_dists[CHI] = [sv_bounded]   # current baseline: only chi upscatter bounded

    phase_spaces = build_onshell_phase_spaces(targets, chain)
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
        stopping_condition=onshell_stopping_condition,
    )
    for ev in injector:
        break
    cpp_inj = injector._Injector__injector
    cpp_inj.ResetInjectedEvents(n_events)

    weighter = Weighter(
        injectors=[injector],
        detector_model=detector_model,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_physical_distributions=[pion_dist,
                                        distributions.NormalizationConstant(
                                            pion_decay._total_width / GAMMA_PION_SM)],
        secondary_interactions=secondary_interactions,
    )
    fid_metric = make_fiducial_metric(fiducial)

    n_complete = 0
    n_fail = 0
    fail_classes = {}
    v1sig_inside = 0
    v1sig_outside = 0
    outside_offsets = []          # per-axis excess (m) for outside V1_sig
    upscatter_to_v1sig = []       # distance chi'-decay-point produced V1_sig -> V1_sig decay
    weights = []

    print(f"\nGenerating {n_events} events...")
    for i in range(n_events):
        event = cpp_inj.GenerateEvent()
        if event.tree:
            n_complete += 1
            v1sig_vtx = vertex_of(event.tree, 5923)
            upscatter_vtx = vertex_of(event.tree, 5917)   # chi continues? no; chi is produced
            chip_vtx = vertex_of(event.tree, 5918)
            if v1sig_vtx is not None:
                if inside_box(v1sig_vtx):
                    v1sig_inside += 1
                else:
                    v1sig_outside += 1
                    outside_offsets.append(box_signed_offset(v1sig_vtx))
                if chip_vtx is not None:
                    upscatter_to_v1sig.append(np.linalg.norm(v1sig_vtx - chip_vtx))
            w = weighter(event)
            weights.append(fid_metric(event, w))
        else:
            n_fail += 1
            failed = cpp_inj.GetLastFailedTree()
            cat = classify_failure(failed.tree) if failed.tree else "empty"
            fail_classes[cat] = fail_classes.get(cat, 0) + 1

    weights = np.array(weights)
    valid = weights[(np.isfinite(weights)) & (weights > 0)]

    print(f"\n{'='*64}")
    print(f"Completed: {n_complete}, Failures: {n_fail} "
          f"({100*n_fail/max(n_complete+n_fail,1):.1f}% injection failure)")
    print(f"{'='*64}")

    print("\nInjection-failure classification (how far the chain got):")
    for cat, cnt in sorted(fail_classes.items(), key=lambda x: -x[1]):
        print(f"  {cat:18s}: {cnt:4d} ({100*cnt/max(n_fail,1):5.1f}%)")

    print("\nV1_sig (signal) vertex vs fiducial box:")
    tot = v1sig_inside + v1sig_outside
    print(f"  inside  fiducial: {v1sig_inside:4d} ({100*v1sig_inside/max(tot,1):5.1f}%)")
    print(f"  outside fiducial: {v1sig_outside:4d} ({100*v1sig_outside/max(tot,1):5.1f}%)")
    if outside_offsets:
        off = np.array(outside_offsets)
        dist = np.linalg.norm(off, axis=1)
        print(f"  outside offset (m): mean={dist.mean():.3f} "
              f"median={np.median(dist):.3f} max={dist.max():.3f}")
        print(f"    per-axis mean excess: x={off[:,0].mean():.3f} "
              f"y={off[:,1].mean():.3f} z={off[:,2].mean():.3f}")
        # how many are "marginal" (within 0.5 m of the box on all axes)?
        marginal = np.sum(dist < 0.5)
        print(f"    marginal (<0.5 m outside): {marginal}/{len(dist)} "
              f"({100*marginal/len(dist):.1f}%)")

    if upscatter_to_v1sig:
        d = np.array(upscatter_to_v1sig)
        print(f"\nchi'-decay-point -> V1_sig-decay-point distance (m): "
              f"mean={d.mean():.4f} median={np.median(d):.4f} max={d.max():.4f}")
        print("  (confirms how far V1_sig travels before its e+e- decay)")

    if len(valid) > 0:
        eff = effective_sample_fraction(weights) * 100
        sm = valid.std() / valid.mean()
        print(f"\nMetric-passing: {len(valid)}/{n_complete} "
              f"({100*len(valid)/max(n_complete,1):.1f}%)")
        print(f"Weight stddev/mean: {sm:.2f}")
        print(f"Weight range: [{valid.min():.3e}, {valid.max():.3e}] "
              f"(dynamic range 1e{np.log10(valid.max()/valid.min()):.0f})")
        print(f"Eff. sample (unoptimized): {eff:.1f}%")


if __name__ == "__main__":
    default_dk2nu = os.environ.get(
        "DK2NU_DIR", "/Users/aschneider/workspaces/g4bnb/sources/G4BNB")
    parser = argparse.ArgumentParser(description="Diagnose on-shell efficiency")
    parser.add_argument("--dk2nu-dir", default=default_dk2nu)
    parser.add_argument("--n-events", type=int, default=2000)
    args = parser.parse_args()
    run_diagnostic(args.n_events, args.dk2nu_dir)
