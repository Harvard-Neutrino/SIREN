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
import glob
import sys
import os

import numpy as np

import siren
from siren import _util, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter
from siren.tune import optimize_chain_weights

# DuttaKim_SBND_full_chain.py was rewritten to a spec-form (Vertex/channels)
# module; it still exports the on-shell model factory and the particle-type
# constants (V1_SIGNAL was renamed V1_SIG), but no longer exports the
# geometric-target/phase-space/stopping-condition helpers this experiment
# script needs to build its own toggleable per-config biasing -- those are
# reproduced locally below from the pre-rewrite implementation.
from DuttaKim_SBND_full_chain import (
    build_onshell_models, build_sX_cdf,
    PION, CHI, CHI_PRIME, V1_PROD, V1_SIG as V1_SIGNAL,
)

_PROC_DIR = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
_DK = _util.load_module("DuttaKim_Dk2nuReader",
                        os.path.join(_PROC_DIR, "Dk2nuReader.py"))

GAMMA_PION_SM = 2.5281e-17  # GeV, PDG total width of pi+


def default_pion_bias(E, px, py, pz, vx, vy, vz):
    """Bias toward pions likely to produce observable signal events.

    Combines three factors (all in geometry/beam coordinates):
      - Energy: higher-energy pions produce more collimated V1/chi,
        increasing both the upscattering cross section and the
        geometric acceptance for the downstream detector. Uses E^2
        because the chi opening angle scales as 1/gamma ~ m/E,
        and acceptance scales with the solid angle ratio.
      - Forward angle: pions decaying forward along the beam axis
        are more likely to send secondaries toward the detector.
      - Transverse position: pions closer to the beam axis have
        better geometric acceptance for a downstream detector.
    """
    p = np.sqrt(px**2 + py**2 + pz**2)
    cos_theta = np.where(p > 0, pz / p, 0.0)
    r_trans = np.sqrt(vx**2 + vy**2)
    return E**2 * np.maximum(cos_theta, 0.01) * np.exp(-r_trans / 200.0)


def load_dk2nu_pions(dk2nu_dir, detector_model, sampling_bias=None):
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
        dk2nu_data, detector_model, parent_pdg=_DK.PTYPE_PIPLUS,
        sampling_bias=sampling_bias)
    if sampling_bias is not None:
        print("  Biased pion selection enabled")
    return primary_dist, pot


def build_geometric_targets(detector_model, fiducial):
    """Build a set of biasing target geometries of increasing size.

    Returns a dict with named targets: a tight fiducial box, spheres of
    increasing radius centered on the detector origin, and cylinder
    segments along the beam axis at two radii covering four z-ranges
    (near-target, mid-range, near-detector, downstream).
    """
    from siren.math import Vector3D
    from siren.geometry import Placement, Sphere, Cylinder

    det_placement = Placement(Vector3D(0, 0, 0))

    # Cylinder segments along the beam axis (+z in detector coordinates).
    # Each segment is (center_z, full_length); Cylinder(placement, radius,
    # inner_radius, z) uses the FULL height z, so a cylinder with z=50
    # extends from center_z-25 to center_z+25.
    cyl_segments = {
        "target":   (-88.0, 50.0),   # z = -113 to -63  (near BNB target)
        "mid":      (-38.0, 50.0),   # z = -63  to -13  (mid-range)
        "near_det": ( -3.0, 40.0),   # z = -23  to +17  (near detector)
        "down":     ( 52.0, 70.0),   # z = +17  to +87  (downstream)
    }
    cyl_radii = [2.0, 5.0]

    targets = {
        "fiducial": fiducial,
        "sphere_2m": Sphere(det_placement, 2.0, 0.0).create(),
        "sphere_5m": Sphere(det_placement, 5.0, 0.0).create(),
        "sphere_10m": Sphere(det_placement, 10.0, 0.0).create(),
        "sphere_20m": Sphere(det_placement, 20.0, 0.0).create(),
    }

    for seg_name, (center_z, full_len) in cyl_segments.items():
        for radius in cyl_radii:
            name = f"cyl_r{radius:.0f}m_{seg_name}"
            placement = Placement(Vector3D(0, 0, center_z))
            targets[name] = Cylinder(placement, radius, 0.0, full_len).create()

    return targets


def _mc(channels, weights):
    m = injection.MultiChannelPhaseSpace()
    m.channels = channels
    m.weights = weights
    return m


def _build_2body_channels(targets, daughter_index, decay, sig):
    """Physical channel + one directed channel per target geometry."""
    channels = [injection.PhysicalDecayChannel(decay, sig)]
    for target in targets:
        channels.append(
            injection.DetectorDirected2BodyChannel(target, daughter_index))
    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return _mc(channels, weights)


def _build_scatter_channels(targets, directed_index, xs, sig):
    """Physical channel + one directed channel per target geometry."""
    channels = [injection.PhysicalCrossSectionChannel(xs, sig)]
    for target in targets:
        channels.append(
            injection.DetectorDirectedScatteringChannel(
                target, directed_index=directed_index,
                variable=injection.ScatteringVariable.Q2))
    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return _mc(channels, weights)


def build_primary_phase_spaces(targets, pion_decay):
    """Multi-channel phase space for the primary pion 3-body decay.

    Directs V1 (secondary index 2) toward the target geometries; the
    complementary (mu, nu) pair invariant mass s_X is drawn from the
    physical marginal's tabulated CDF (build_sX_cdf), which keeps f/g ~ 1
    at this vertex for any masses.
    """
    sig = pion_decay.GetPossibleSignatures()[0]
    geo_list = list(targets.values())
    cdf_nodes, cdf_values = build_sX_cdf(pion_decay)

    channels = [injection.PhysicalDecayChannel(pion_decay, sig)]
    for target in geo_list:
        channels.append(
            injection.DetectorDirected3BodyChannel(
                target,
                directed_index=2,
                mass_mode=injection.InvariantMassMode.Tabulated,
                topology=injection.PhaseSpaceTopology.Decay3Body,
                mass_cdf_nodes=cdf_nodes,
                mass_cdf_values=cdf_values))
    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return {sig: _mc(channels, weights)}


def onshell_stopping_condition(tree, datum, i):
    sec = int(datum.record.signature.secondary_types[i])
    parent = int(datum.record.signature.primary_type)
    if sec == 5922:   return False          # V1_prod
    if sec == 5917:   return parent != 5922 or i != 0  # chi only from V1 decay
    if sec == 5918:   return False          # chi'
    if sec == 5923:   return False          # V1_signal
    return True


def effective_sample_fraction(weights):
    w = np.asarray(weights)
    w = w[(np.isfinite(w)) & (w > 0)]
    if len(w) == 0:
        return 0.0
    return (w.sum() ** 2) / (len(w) * (w ** 2).sum()) * len(w) / len(weights)


def make_fiducial_metric(fiducial, signal_pdgids=(5923,)):
    """Metric selecting events with a signal vertex inside the fiducial volume.

    Returns a callable metric(event, weight) -> float suitable for
    optimize_chain_weights(metric=...): the weight if any interaction
    vertex for a signal particle (default: V1_signal, PDG 5923) is inside
    the fiducial geometry, else 0.
    """
    from siren.math import Vector3D

    def metric(event, w):
        for datum in event.tree:
            r = datum.record
            if int(r.signature.primary_type) in signal_pdgids:
                vtx = Vector3D(r.interaction_vertex[0],
                               r.interaction_vertex[1],
                               r.interaction_vertex[2])
                if fiducial.IsInside(vtx):
                    return w
        return 0.0

    return metric


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
