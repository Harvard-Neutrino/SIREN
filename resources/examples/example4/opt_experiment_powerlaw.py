"""A/B test: PowerLaw vs Uniform s_X mapping at the primary pion vertex.

The biased DetectorDirected3BodyChannel at the primary pi -> mu nu V1
vertex samples the (mu,nu) pair invariant mass s_X uniformly, while the
physical marginal dGamma/dE_V1 (lepton-propagator shaped) rises toward
high s_X.  A PowerLaw mapping g(s) ~ (s - m2)^(-nu) with negative nu
concentrates the proposal at high s_X.  A 1D fit minimizing int p^2/g
(the importance-weight second moment) gives nu ~ -8.6, m2 = 0, dropping
this factor's std/mean from 0.73 to ~0.35.

This script runs the full optimize+produce loop with the on-shell chain
under Uniform and PowerLaw primary mappings and reports effective sample
fraction, so we can confirm the per-vertex win survives end to end.

Usage:
    python opt_experiment_powerlaw.py [--n N] [--opt-batch N] [--opt-iters N]
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
    onshell_stopping_condition, load_dk2nu_pions, default_pion_bias,
    make_fiducial_metric, effective_sample_fraction, _mc,
    GAMMA_PION_SM, PION, CHI,
)


def build_primary_ps(targets, pion_decay, mode, nu=-8.6, offset=0.0):
    """Primary 3-body phase space with selectable s_X mapping."""
    sig = pion_decay.GetPossibleSignatures()[0]
    geo_list = list(targets.values())
    channels = [injection.PhysicalDecayChannel(pion_decay, sig)]
    for target in geo_list:
        channels.append(
            injection.DetectorDirected3BodyChannel(
                target,
                directed_index=2,
                mass_mode=mode,
                resonance_mass=0.0,
                resonance_width=0.0,
                power_law_nu=nu,
                power_law_offset=offset,
                topology=injection.PhaseSpaceTopology.Decay3Body))
    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return {sig: _mc(channels, weights)}


def run_config(name, dk2nu_dir, n, opt_batch, opt_iters, mode, nu, seed=42):
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
    primary_ps = build_primary_ps(targets, pion_decay, mode, nu=nu)

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
    tot = valid.sum()
    print(f"\n### CONFIG: {name}  (mode={mode}, nu={nu})")
    print(f"    completed={len(events)}  pass_metric={len(valid)} "
          f"({100*len(valid)/max(len(events),1):.1f}%)")
    print(f"    stddev/mean={sm:.3f}  EFF={eff:.1f}%  sum(w)={tot:.4e}")
    return {"name": name, "eff": eff, "sm": sm, "pass": len(valid),
            "completed": len(events), "sumw": tot}


if __name__ == "__main__":
    default_dk2nu = os.environ.get(
        "DK2NU_DIR", "/Users/aschneider/workspaces/g4bnb/sources/G4BNB")
    parser = argparse.ArgumentParser()
    parser.add_argument("--dk2nu-dir", default=default_dk2nu)
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--opt-batch", type=int, default=1000)
    parser.add_argument("--opt-iters", type=int, default=5)
    parser.add_argument("--nu", type=float, default=-8.6)
    args = parser.parse_args()

    PL = injection.InvariantMassMode.PowerLaw
    UN = injection.InvariantMassMode.Uniform
    results = [
        run_config("Uniform_baseline", args.dk2nu_dir, args.n,
                   args.opt_batch, args.opt_iters, UN, 0.8),
        run_config("PowerLaw_fit", args.dk2nu_dir, args.n,
                   args.opt_batch, args.opt_iters, PL, args.nu),
    ]

    print("\n" + "=" * 58)
    print(f"{'config':<22} {'pass%':>6} {'sd/mu':>7} {'EFF%':>6} {'sum(w)':>11}")
    print("-" * 58)
    for r in results:
        pr = 100 * r["pass"] / max(r["completed"], 1)
        print(f"{r['name']:<22} {pr:6.1f} {r['sm']:7.3f} {r['eff']:6.1f} "
              f"{r['sumw']:11.4e}")
    # physics-invariance check: sum(w) must agree within MC error
    if all(r["sumw"] > 0 for r in results):
        ratio = results[1]["sumw"] / results[0]["sumw"]
        print(f"\nsum(w) ratio PowerLaw/Uniform = {ratio:.3f} "
              f"(should be ~1 if the mapping is unbiased)")
