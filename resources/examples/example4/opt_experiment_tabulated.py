"""A/B/C: Tabulated vs PowerLaw vs Uniform s_X mapping at the primary vertex.

Uses the canonical build_primary_phase_spaces (which builds the s_X CDF
table internally for Tabulated mode) so this exercises exactly the code
path shipped in the example.  Reports effective sample fraction, weight
std/mean, and sum(weights) (the physics estimate, which must be invariant
across mappings).

Usage:
    python opt_experiment_tabulated.py [--n N] [--opt-batch N] [--opt-iters N]
"""
import argparse
import os

import numpy as np

import siren
from siren import distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter
from siren.tune import optimize_chain_weights

from DuttaKim_SBND_full_chain import (
    build_onshell_models, build_onshell_phase_spaces, build_geometric_targets,
    build_primary_phase_spaces, onshell_stopping_condition,
    load_dk2nu_pions, default_pion_bias, make_fiducial_metric,
    effective_sample_fraction, GAMMA_PION_SM, PION, CHI,
)


def run_config(name, dk2nu_dir, n, opt_batch, opt_iters, mode, nu=-8.6, seed=42):
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(
        widths=((4.5 + 2.022) * 2, 4.074645 * 2, 5.010 * 2))
    targets = build_geometric_targets(detector_model, fiducial)

    chain = build_onshell_models()
    pion_decay = chain["pion_decay"]
    secondary_interactions = chain["secondary_interactions"]

    pion_dist, pot = load_dk2nu_pions(dk2nu_dir, detector_model,
                                      sampling_bias=default_pion_bias)
    br_dist = distributions.NormalizationConstant(
        pion_decay._total_width / GAMMA_PION_SM)

    sv = distributions.SecondaryPhysicalVertexDistribution()
    sv_bounded = distributions.SecondaryBoundedVertexDistribution(fiducial)
    sec_dists = {pt: [sv] for pt in secondary_interactions}
    sec_dists[CHI] = [sv_bounded]

    phase_spaces = build_onshell_phase_spaces(targets, chain)
    primary_ps = build_primary_phase_spaces(
        targets, pion_decay, mass_mode=mode, power_law_nu=nu)

    injector = Injector(
        events=n, detector_model=detector_model, seed=seed,
        primary_type=PION, primary_interactions=[pion_decay],
        primary_injection_distributions=[pion_dist],
        primary_weighting_mode=injection.VertexWeightingMode.Fixed(),
        secondary_interactions=secondary_interactions,
        secondary_injection_distributions=sec_dists,
        secondary_phase_spaces=phase_spaces,
        primary_phase_spaces=primary_ps,
        stopping_condition=onshell_stopping_condition,
    )
    for ev in injector:
        break
    injector._Injector__injector.ResetInjectedEvents(n)

    weighter = Weighter(
        injectors=[injector], detector_model=detector_model,
        primary_type=PION, primary_interactions=[pion_decay],
        primary_physical_distributions=[pion_dist, br_dist],
        secondary_interactions=secondary_interactions,
    )
    fid_metric = make_fiducial_metric(fiducial)

    optimize_chain_weights(
        injector, weighter, n_iterations=opt_iters, batch_size=opt_batch,
        damping=0.5, metric=fid_metric, verbose=False, min_weight=1e-4)

    events = []
    for event in injector:
        if event.tree:
            events.append(event)
        if len(events) >= n:
            break

    weights = np.array([weighter(event) for event in events])
    weights = np.array([fid_metric(event, w) for event, w in zip(events, weights)])
    valid = weights[(np.isfinite(weights)) & (weights > 0)]
    eff = effective_sample_fraction(weights) * 100
    sm = valid.std() / valid.mean() if len(valid) else float("nan")
    print(f"\n### {name} (mode={mode})")
    print(f"    completed={len(events)} pass={len(valid)} "
          f"({100*len(valid)/max(len(events),1):.1f}%)  "
          f"sd/mu={sm:.3f}  EFF={eff:.1f}%  sum(w)={valid.sum():.4e}")
    return {"name": name, "eff": eff, "sm": sm, "pass": len(valid),
            "completed": len(events), "sumw": valid.sum()}


if __name__ == "__main__":
    default_dk2nu = os.environ.get(
        "DK2NU_DIR", "/Users/aschneider/workspaces/g4bnb/sources/G4BNB")
    ap = argparse.ArgumentParser()
    ap.add_argument("--dk2nu-dir", default=default_dk2nu)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--opt-batch", type=int, default=1000)
    ap.add_argument("--opt-iters", type=int, default=5)
    args = ap.parse_args()

    IM = injection.InvariantMassMode
    results = [
        run_config("Uniform", args.dk2nu_dir, args.n, args.opt_batch,
                   args.opt_iters, IM.Uniform),
        run_config("PowerLaw_-8.6", args.dk2nu_dir, args.n, args.opt_batch,
                   args.opt_iters, IM.PowerLaw, nu=-8.6),
        run_config("Tabulated", args.dk2nu_dir, args.n, args.opt_batch,
                   args.opt_iters, IM.Tabulated),
    ]

    print("\n" + "=" * 60)
    print(f"{'config':<16} {'pass%':>6} {'sd/mu':>7} {'EFF%':>6} {'sum(w)':>12}")
    print("-" * 60)
    base = results[0]["sumw"]
    for r in results:
        pr = 100 * r["pass"] / max(r["completed"], 1)
        print(f"{r['name']:<16} {pr:6.1f} {r['sm']:7.3f} {r['eff']:6.1f} "
              f"{r['sumw']:12.4e}")
    print("\nsum(w) relative to Uniform (should be ~1, unbiased):")
    for r in results:
        print(f"  {r['name']:<16} {r['sumw']/base:.3f}")
