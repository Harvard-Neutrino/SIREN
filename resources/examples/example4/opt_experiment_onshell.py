"""Experiment: which on-shell biasing configuration minimizes weight variance?

The on-shell chain (5 vertices) has lower effective sampling than the
off-shell chain (4 vertices).  The diagnostic shows on-shell has a *higher*
fiducial pass rate but much larger weight variance.  The extra variance
comes from two directed vertices the off-shell chain does not have:

    upscatter   chi + Ar -> chi' + Ar   (directs chi' at a target)
    chi' decay  chi' -> chi V1_sig      (directs V1_sig at a target)

Both occur inside the bounded fiducial region: the upscatter vertex is
confined to the fiducial, and chi' decays sub-micron away, so V1_sig is
produced inside the fiducial and decays ~cm later -- regardless of which
direction chi' or V1_sig point.  Directing them therefore adds importance
weight variance without improving acceptance.  Using the physical channel
alone at such a vertex makes its weight ratio identically 1.

This script runs the same optimize+produce loop under several
configurations and reports metric-pass rate, weight stddev/mean, and
effective sample fraction so we can pick the best.

Usage:
    python opt_experiment_onshell.py [--n N] [--opt-batch N] [--opt-iters N]
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
    build_onshell_models, build_geometric_targets, build_primary_phase_spaces,
    onshell_stopping_condition, load_dk2nu_pions, default_pion_bias,
    make_fiducial_metric, effective_sample_fraction,
    _build_2body_channels, _build_scatter_channels, _mc,
    GAMMA_PION_SM, PION, CHI, CHI_PRIME, V1_PROD, V1_SIGNAL,
)


def build_phase_spaces(targets, models, *, bias_upscatter, bias_chiprime):
    """On-shell phase spaces with toggleable biasing per vertex."""
    m = models["models"]
    v1_sig = m["v1_to_chi"].GetPossibleSignatures()[0]
    chi_sig = m["upscatter"].GetPossibleSignatures()[0]
    chip_sig = m["chi_prime_decay"].GetPossibleSignatures()[0]
    vis_sig = m["visible_decay"].GetPossibleSignatures()[0]
    geo = list(targets.values())

    ps = {}
    # V1 -> chi chi: always directed (chi must reach the fiducial so the
    # bounded upscatter can occur).  Direct index 0 (the continuing chi).
    ps[V1_PROD] = {v1_sig: _build_2body_channels(geo, 0, m["v1_to_chi"], v1_sig)}

    # chi + Ar -> chi' + Ar
    if bias_upscatter:
        ps[CHI] = {chi_sig: _build_scatter_channels(geo, 0, m["upscatter"], chi_sig)}
    else:
        ps[CHI] = {chi_sig: _mc(
            [injection.PhysicalCrossSectionChannel(m["upscatter"], chi_sig)], [1.0])}

    # chi' -> chi V1_sig
    if bias_chiprime:
        ps[CHI_PRIME] = {chip_sig: _build_2body_channels(
            geo, 1, m["chi_prime_decay"], chip_sig)}
    else:
        ps[CHI_PRIME] = {chip_sig: _mc(
            [injection.PhysicalDecayChannel(m["chi_prime_decay"], chip_sig)], [1.0])}

    # V1_sig -> e+ e-: physical + isotropic (unchanged from baseline)
    ps[V1_SIGNAL] = {vis_sig: _mc([
        injection.PhysicalDecayChannel(m["visible_decay"], vis_sig),
        injection.Isotropic2BodyChannel(0),
    ], [0.50, 0.50])}
    return ps


def run_config(name, dk2nu_dir, n, opt_batch, opt_iters,
               *, bias_upscatter, bias_chiprime, seed=42):
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

    phase_spaces = build_phase_spaces(
        targets, chain, bias_upscatter=bias_upscatter, bias_chiprime=bias_chiprime)
    primary_ps = build_primary_phase_spaces(targets, pion_decay)

    injector = Injector(
        events=n,
        detector_model=detector_model,
        seed=seed,
        primary_type=PION,
        primary_interactions=[pion_decay],
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
        injectors=[injector],
        detector_model=detector_model,
        primary_type=PION,
        primary_interactions=[pion_decay],
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
    dyn = (np.log10(valid.max() / valid.min()) if len(valid) else float("nan"))
    print(f"\n### CONFIG: {name}")
    print(f"    bias_upscatter={bias_upscatter} bias_chiprime={bias_chiprime}")
    print(f"    completed={len(events)}  pass_metric={len(valid)} "
          f"({100*len(valid)/max(len(events),1):.1f}%)")
    print(f"    stddev/mean={sm:.2f}  dyn_range=1e{dyn:.0f}  EFF={eff:.1f}%")
    return {"name": name, "eff": eff, "sm": sm, "pass": len(valid),
            "completed": len(events)}


if __name__ == "__main__":
    default_dk2nu = os.environ.get(
        "DK2NU_DIR", "/Users/aschneider/workspaces/g4bnb/sources/G4BNB")
    parser = argparse.ArgumentParser()
    parser.add_argument("--dk2nu-dir", default=default_dk2nu)
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--opt-batch", type=int, default=1000)
    parser.add_argument("--opt-iters", type=int, default=5)
    args = parser.parse_args()

    configs = [
        ("A_baseline_both_directed",  dict(bias_upscatter=True,  bias_chiprime=True)),
        ("B_chiprime_physical",       dict(bias_upscatter=True,  bias_chiprime=False)),
        ("C_both_physical",           dict(bias_upscatter=False, bias_chiprime=False)),
        ("D_upscatter_physical",      dict(bias_upscatter=False, bias_chiprime=True)),
    ]
    results = []
    for name, flags in configs:
        results.append(run_config(name, args.dk2nu_dir, args.n,
                                  args.opt_batch, args.opt_iters, **flags))

    print("\n" + "=" * 60)
    print(f"{'config':<28} {'pass%':>6} {'sd/mu':>7} {'EFF%':>6}")
    print("-" * 60)
    for r in results:
        pr = 100 * r["pass"] / max(r["completed"], 1)
        print(f"{r['name']:<28} {pr:6.1f} {r['sm']:7.2f} {r['eff']:6.1f}")
