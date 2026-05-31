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
GAMMA_PION_SM = 2.5281e-17  # GeV, PDG total width of pi+

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


def build_geometric_targets(detector_model, fiducial, fiducial_only=False):
    """Build a set of biasing target geometries of increasing size.

    Returns a dict with named targets:
      - fiducial: tight box matching the SBND TPC active volume
      - sphere_5m, sphere_10m, sphere_20m: spheres centered on the
        detector origin with increasing radius
      - Cylinder segments along the beam axis at two radii (2m, 5m),
        covering four z-ranges in detector coordinates:
          near-target (z=-113 to -63m), mid-range (-63 to -13m),
          near-detector (-23 to +17m), downstream (+17 to +87m)

    If ``fiducial_only`` is True, return just the fiducial target.  The
    optimizer cannot tell the 13 near-degenerate geometric targets apart --
    it drives them to near-uniform weights (the Kleiss-Pittau degeneracy
    signature) -- so one physical + one fiducial-directed channel per biased
    vertex gives the same accuracy at far fewer density evaluations (gap R2).
    """
    if fiducial_only:
        return {"fiducial": fiducial}

    from siren.math import Vector3D
    from siren.geometry import Placement, Sphere, Cylinder

    det_placement = Placement(Vector3D(0, 0, 0))

    # Cylinder segments along the beam axis (+z in detector coordinates).
    # Positions are in detector coordinates where the detector center is
    # at the origin and the BNB target is at z ~ -113m.
    #
    # Each segment is defined by (center_z, full_length) in detector coords.
    # Note: Cylinder(placement, radius, inner_radius, z) uses FULL height z,
    # so a cylinder with z=50 extends from center_z-25 to center_z+25.
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


def _build_2body_channels(targets, daughter_index, decay, sig):
    """Build multi-channel for a 2-body decay vertex.

    Physical channel + one directed channel per target geometry.
    """
    channels = [injection.PhysicalDecayChannel(decay, sig)]
    for target in targets:
        channels.append(
            injection.DetectorDirected2BodyChannel(target, daughter_index))

    n = len(channels)
    # Physical gets a small weight, rest split equally
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return _mc(channels, weights)


def _build_scatter_channels(targets, directed_index, xs, sig):
    """Build multi-channel for a 2->2 scattering vertex."""
    channels = [injection.PhysicalCrossSectionChannel(xs, sig)]
    for target in targets:
        channels.append(
            injection.DetectorDirectedScatteringChannel(
                target, directed_index=directed_index,
                variable=injection.ScatteringVariable.Q2))

    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return _mc(channels, weights)


def _build_scatter_3body_channels(targets, xs, sig):
    """Build multi-channel for a 2->3 off-shell scattering vertex.

    Signature: chi -> [chi, V1_sig, Ar].  Direct V1_sig (index 1)
    toward the detector.  BW on the pair mass M(chi+V1_sig) targets
    the chi' resonance.  Recursive mode is fine here because the
    CM boost is small (beta ~ 0.03).
    """
    channels = [injection.PhysicalCrossSectionChannel(xs, sig)]
    for target in targets:
        channels.append(
            injection.DetectorDirected3BodyChannel(
                target,
                spectator_index=2, pair_first_index=0, pair_second_index=1,
                directed_pair_index=1,
                mass_mode=injection.InvariantMassMode.BreitWigner,
                resonance_mass=M_CHI_PRIME,
                resonance_width=CHI_PRIME_WIDTH,
                topology=injection.PhaseSpaceTopology.Scatter2to3))

    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return _mc(channels, weights)


def build_sX_cdf_table(pion_decay, n_nodes=257):
    """Tabulate the cumulative of the physical s_X = M^2(mu, nu) marginal.

    The directed primary channel samples the (mu, nu) pair invariant mass
    s_X; its physical marginal is the Carlson-Rislow dGamma differential,
    which the decay model already exposes as dGamma/dE_V1 (E_V1 = the V1
    rest-frame energy).  E_V1 maps linearly to s_X via
        E_V1 = (m_M^2 + m_V1^2 - s_X) / (2 m_M),
    a constant Jacobian, so dGamma/ds_X is proportional to dGamma/dE_V1.

    Returns (s_nodes, cdf_nodes): ascending s_X grid and the cumulative
    (trapezoid) integral of the marginal, for InvariantMassMode.Tabulated.
    The cumulative is left unnormalized (the C++ mapping renormalizes over
    the channel's allowed [s_min, s_max]); only its shape matters, so the
    result is independent of couplings and overall constants.

    This is exact for any meson / lepton / mediator masses -- no tuned
    exponent -- and drives the per-vertex importance ratio f/g to ~1.
    """
    engine = pion_decay._decay   # underlying MesonThreeBodyDecay
    m_M = engine.m_M
    m_l = engine.m_l
    m_V1 = engine.m_phi

    s_min = m_l * m_l
    s_max = (m_M - m_V1) ** 2
    # Pad the endpoints slightly inward: the marginal -> 0 at both edges
    # and the matrix element has an integrable edge, so sampling the open
    # interval avoids a zero-density node.
    eps = 1e-9 * (s_max - s_min)
    s_nodes = np.linspace(s_min + eps, s_max - eps, n_nodes)

    E_V1 = (m_M ** 2 + m_V1 ** 2 - s_nodes) / (2.0 * m_M)
    dens = np.asarray(engine.differential_decay_rate(E_V1), dtype=float)
    dens = np.clip(dens, 0.0, None)

    # Cumulative trapezoid integral -> CDF nodes (ascending, starts at 0).
    cdf = np.concatenate([[0.0], np.cumsum(
        0.5 * (dens[1:] + dens[:-1]) * np.diff(s_nodes))])
    if cdf[-1] <= 0.0:
        raise RuntimeError("s_X marginal integrated to zero; cannot tabulate")
    return s_nodes.tolist(), cdf.tolist()


def build_primary_phase_spaces(targets, pion_decay,
                               mass_mode=injection.InvariantMassMode.Tabulated,
                               power_law_nu=-8.6, power_law_offset=0.0):
    """Multi-channel phase space for the primary pion 3-body decay.

    Uses direct lab-frame biasing: V1 (secondary index 2) is
    directed toward the target geometries using the pion's boost.
    The complementary (mu, nu) system supplies the pair invariant mass
    s_X = M^2(mu, nu), which the directed channel must sample.

    s_X mapping: the Carlson-Rislow matrix element carries a charged-
    lepton propagator 1/D^2 (D = t - m_mu^2), so the physical marginal
    dGamma/dE_V1 is not flat -- it rises toward high s_X (low rest-frame
    E_V1).  Sampling s_X uniformly leaves a residual f/g spread
    (per-vertex std/mean ~0.73) that, after the upscatter Q^2 fix, is the
    largest remaining variance source in the chain.

    Default mode is Tabulated: the proposal is the physical s_X marginal
    itself (a CDF built once from the decay model, see build_sX_cdf_table),
    so f/g -> 1 at this vertex for *any* masses, with no tuned parameter
    and no coupling dependence (the cumulative shape is coupling-free).
    This lifts end-to-end effective sampling from ~40% (Uniform) to ~50%.

    The optimal-exponent PowerLaw (nu ~ -8.6 for the pi->mu benchmark)
    remains available via mass_mode=InvariantMassMode.PowerLaw, and
    InvariantMassMode.Uniform recovers the old flat proposal -- but note
    a fixed PowerLaw exponent is mass-specific (it is unsafe to reuse for
    a different meson/lepton), whereas the Tabulated map adapts itself.
    """
    sig = pion_decay.GetPossibleSignatures()[0]
    geo_list = list(targets.values())

    cdf_nodes, cdf_values = [], []
    if mass_mode == injection.InvariantMassMode.Tabulated:
        cdf_nodes, cdf_values = build_sX_cdf_table(pion_decay)

    channels = [injection.PhysicalDecayChannel(pion_decay, sig)]
    for target in geo_list:
        channels.append(
            injection.DetectorDirected3BodyChannel(
                target,
                directed_index=2,
                mass_mode=mass_mode,
                resonance_mass=0.0,
                resonance_width=0.0,
                power_law_nu=power_law_nu,
                power_law_offset=power_law_offset,
                topology=injection.PhaseSpaceTopology.Decay3Body,
                mass_cdf_nodes=cdf_nodes,
                mass_cdf_values=cdf_values))
    n = len(channels)
    weights = [0.02] + [(1.0 - 0.02) / (n - 1)] * (n - 1)
    return {sig: _mc(channels, weights)}


def build_onshell_phase_spaces(targets, models):
    """Multi-channel phase spaces for on-shell chain.

    Only two vertices are directed toward the detector:

      - V1 -> chi chi: direct the continuing chi (index 0) so it reaches
        the fiducial volume, where the bounded upscatter vertex confines
        the chi + Ar -> chi' + Ar interaction.
      - V1_sig -> e+ e-: a physical/isotropic mix (location is what the
        metric cares about, and that is fixed upstream).

    The intermediate upscatter (chi + Ar -> chi' + Ar) and the chi' decay
    (chi' -> chi V1_sig) use the *physical* channel only -- they are NOT
    directed.  This is the key difference from the naive 5-vertex setup:
    once the upscatter vertex is confined to the fiducial (see sec_dists
    in run()), chi' is produced there, decays sub-micron later (lab
    decay length ~1e-10 m for M_CHI_PRIME=50 MeV, G_D=1), and the
    resulting V1_sig travels only ~1 cm before its e+ e- decay (measured
    median 4.5 mm, max 0.43 m).  The signal vertex therefore lands inside
    the fiducial regardless of the chi' or V1_sig directions, so directing
    them cannot improve acceptance -- it only injects large importance
    weights (the directed Q2 scattering channel spans >60 orders of
    magnitude against the physical 1/(Q2+m_V2^2)^2).  Using the physical
    channel alone makes each of these vertices contribute a weight ratio
    of exactly 1.  Empirically this lifts the on-shell effective sample
    fraction from ~8% to ~40%, with the bulk of the gain coming from the
    upscatter (the dominant variance source).
    """
    m = models["models"]
    v1_sig = m["v1_to_chi"].GetPossibleSignatures()[0]
    chi_sig = m["upscatter"].GetPossibleSignatures()[0]
    chip_sig = m["chi_prime_decay"].GetPossibleSignatures()[0]
    vis_sig = m["visible_decay"].GetPossibleSignatures()[0]

    geo_list = list(targets.values())

    return {
        V1_PROD: {v1_sig: _build_2body_channels(
            geo_list, 0, m["v1_to_chi"], v1_sig)},
        CHI: {chi_sig: _mc([
            injection.PhysicalCrossSectionChannel(m["upscatter"], chi_sig),
        ], [1.0])},
        CHI_PRIME: {chip_sig: _mc([
            injection.PhysicalDecayChannel(m["chi_prime_decay"], chip_sig),
        ], [1.0])},
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


def build_offshell_phase_spaces(targets, models):
    """Multi-channel phase spaces for off-shell chain."""
    m = models["models"]
    v1_sig = m["v1_to_chi"].GetPossibleSignatures()[0]
    offshell_sig = m["offshell_xs"].GetPossibleSignatures()[0]
    vis_sig = m["visible_decay"].GetPossibleSignatures()[0]

    geo_list = list(targets.values())

    return {
        V1_PROD: {v1_sig: _build_2body_channels(
            geo_list, 0, m["v1_to_chi"], v1_sig)},
        CHI: {offshell_sig: _build_scatter_3body_channels(
            geo_list, m["offshell_xs"], offshell_sig)},
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
    return (w.sum() ** 2) / (len(w) * (w ** 2).sum()) * len(w) / len(weights)


def make_fiducial_metric(fiducial, signal_pdgids=(5923,)):
    """Build a metric that selects events with a signal vertex in the fiducial volume.

    Returns a callable ``metric(event, weight) -> float`` suitable for
    ``optimize_chain_weights(metric=...)``.  Returns the weight if any
    interaction vertex for a signal particle (default: V1_signal, PDG 5923)
    is inside the fiducial geometry; returns 0 otherwise.
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


# ------------------------------------------------------------------ #
#  Main                                                                #
# ------------------------------------------------------------------ #

def run(dk2nu_dir, n_events=100, seed=42, optimize=False,
        opt_iterations=5, opt_batch=200, monoenergetic=False,
        offshell=False, target_set="all_13", tile_n=2):

    # -- Detector --
    print("Loading SBND detector model ...")
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")

    # -- Biasing targets of increasing size --
    x_half_width = 4.5 + 2.022
    y_half_width = 4.074645
    z_half_width = 5.010
    fiducial = siren.geometry.Box(x_half_width*2, y_half_width*2, z_half_width*2)
    targets = build_geometric_targets(detector_model, fiducial)

    # Opt-in volume-aware target sets (Option B; defaults unchanged).  Because
    # the per-vertex mixture is built as "physical + one directed channel per
    # target", swapping the target set swaps the directed-channel basis:
    #   all_13        -- the 13 nested/overlapping probes (default; degenerate)
    #   fiducial_only -- one fiducial-directed channel (gap R2)
    #   tiled         -- a disjoint grid of the fiducial into tile_n^3 cells
    #   union         -- a single BooleanGeometry union of the 13 probes
    # Tiling/union are a non-degenerate basis: the optimizer tunes the cells
    # cleanly and overlapping probes stop being redundant.  See
    # docs/volume_aware_directed_channel_plan.md.
    if target_set == "fiducial_only":
        targets = {"fiducial": fiducial}
    elif target_set == "tiled":
        from siren import directed_tiling as dt
        cells = dt.grid_cells(fiducial, tile_n)
        targets = {f"cell_{i}": c for i, c in enumerate(cells)}
    elif target_set == "union":
        from siren import directed_tiling as dt
        targets = {"union": dt.make_union(list(targets.values()))}
    elif target_set != "all_13":
        raise ValueError(f"unknown target_set {target_set!r}")
    print(f"  Built {len(targets)} biasing targets ({target_set}): "
          f"{list(targets.keys())}")

    # -- Chain topology --
    if offshell:
        print("\nUsing off-shell chi' chain (4 vertices, virtual chi')")
        chain = build_offshell_models()
        phase_spaces = build_offshell_phase_spaces(targets, chain)
        stop_fn = offshell_stopping_condition
    else:
        print("\nUsing on-shell chain (5 vertices)")
        chain = build_onshell_models()
        phase_spaces = build_onshell_phase_spaces(targets, chain)
        stop_fn = onshell_stopping_condition

    pion_decay = chain["pion_decay"]
    secondary_interactions = chain["secondary_interactions"]

    # -- Primary distributions --
    pot = 0.0
    bsm_width = pion_decay._total_width
    br_bsm = bsm_width / GAMMA_PION_SM
    print(f"  BSM pion decay width: {bsm_width:.4e} GeV")
    print(f"  SM pion total width:  {GAMMA_PION_SM:.4e} GeV")
    print(f"  BSM branching ratio:  {br_bsm:.4e}")

    if monoenergetic:
        print("\nUsing monoenergetic 2 GeV pion along +z")
        # For Propagated mode, the InteractionProbability already
        # accounts for the tiny BSM decay rate (long decay length),
        # so we do NOT add the branching ratio here.
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
        pion_dist, pot = load_dk2nu_pions(dk2nu_dir, detector_model,
                                          sampling_bias=default_pion_bias)
        primary_dists = [pion_dist]
        # For Fixed mode (dk2nu), the pion already decayed via SM.
        # The BSM branching ratio must be included explicitly in the
        # physical probability since InteractionProbability is skipped.
        br_dist = distributions.NormalizationConstant(br_bsm)
        physical_dists = [pion_dist, br_dist]
        primary_mode = injection.VertexWeightingMode.Fixed()

    # -- Secondary distributions --
    # Most secondaries use physical vertex sampling (decay/scatter
    # anywhere along the path).  For chi, use bounded vertex sampling
    # to confine the scattering to the fiducial volume -- the LAr
    # cryostat extends well beyond the TPC, and without bounding,
    # ~50% of chi scatters land outside the fiducial.
    #
    # Only the chi upscatter vertex is bounded.  For the on-shell chain
    # the chi' produced there decays essentially in place (lab decay
    # length ~1e-10 m for M_CHI_PRIME=50 MeV), and the resulting V1_sig
    # travels only ~1 cm before decaying, so bounding the upscatter is
    # enough to keep the whole downstream signal inside the fiducial --
    # CHI_PRIME does not need its own bounded vertex distribution.
    sv = distributions.SecondaryPhysicalVertexDistribution()
    sv_bounded = distributions.SecondaryBoundedVertexDistribution(fiducial)
    sec_dists = {pt: [sv] for pt in secondary_interactions}
    sec_dists[CHI] = [sv_bounded]

    # -- Primary phase space biasing --
    primary_ps = build_primary_phase_spaces(targets, pion_decay)
    print(f"  Primary phase space: {len(primary_ps)} signature(s), "
          f"{len(list(primary_ps.values())[0].channels)} channels")

    # -- Budget --
    total_budget = n_events
    if optimize:
        total_budget += opt_iterations * opt_batch

    # -- Build Injector --
    print(f"\nBuilding injector ({total_budget} event budget) ...")
    injector = Injector(
        number_of_events=n_events,
        detector_model=detector_model,
        seed=seed,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_injection_distributions=primary_dists,
        primary_weighting_mode=primary_mode,
        secondary_interactions=secondary_interactions,
        secondary_injection_distributions=sec_dists,
        secondary_phase_spaces=phase_spaces,
        primary_phase_spaces=primary_ps,
        stopping_condition=stop_fn,
    )

    # Force initialization
    for ev in injector:
        break
    injector._Injector__injector.ResetInjectedEvents(n_events)

    # -- Build Weighter --
    weighter = Weighter(
        injectors=[injector],
        detector_model=detector_model,
        primary_type=PION,
        primary_interactions=[pion_decay],
        primary_physical_distributions=physical_dists,
        secondary_interactions=secondary_interactions,
    )

    # Build the fiducial volume in detector coordinates for the metric.
    # Use the fiducial box itself (NOT targets[0]): the signal-in-fiducial
    # metric is independent of the directed-channel target basis, which may be
    # a tiling/union (Option B) rather than the fiducial.
    fid_metric = make_fiducial_metric(fiducial)

    # -- Optimize --
    if optimize:
        from siren.optimize import optimize_chain_weights

        print("\nOptimizing multi-channel weights ...")
        print("  Metric: signal vertex inside fiducial volume")
        optimize_chain_weights(
            injector, weighter,
            n_iterations=opt_iterations,
            batch_size=opt_batch,
            damping=0.5,
            metric=fid_metric,
            verbose=True,
            min_weight=1e-4,
            # Canonical Kleiss-Pittau rule: turns off a directed channel's
            # isotropic fallback in favor of the physical channel at
            # non-isotropic vertices (the memoryless sqrt_W cannot).
            update_rule="alpha_sqrt_W",
        )

    # -- Generate production events --
    print(f"\nGenerating {n_events} production events ...")
    events = []
    for event in injector:
        if event.tree:
            events.append(event)
        if len(events) >= n_events:
            break

    print(f"Generated {len(events)} / {n_events} events")

    # -- Compute weights --
    weights = np.array([weighter(event) for event in events])
    weights = np.array([fid_metric(event, w) for event, w in zip(events, weights)])
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
        eff = effective_sample_fraction(weights) * 100
        sigma_over_mu = valid_weights.std() / valid_weights.mean()
        print(f"Weight range:  [{valid_weights.min():.4e}, "
              f"{valid_weights.max():.4e}]")
        print(f"Weight mean:   {valid_weights.mean():.4e}")
        print(f"Weight stddev/mean:     {sigma_over_mu:.2f}")
        print(f"Eff. sample:   {eff:.1f}%")
        print(f"Total weight:   {valid_weights.sum():.4e} per POT ({valid_weights.sum() * 1e21:.3e} per 1e21 POT)")

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
    parser.add_argument("--target-set", type=str, default="all_13",
                        choices=["all_13", "fiducial_only", "tiled", "union"],
                        help="Directed-channel target basis (Option B)")
    parser.add_argument("--tile-n", type=int, default=2,
                        help="Grid cells per axis for --target-set tiled")
    args = parser.parse_args()
    run(dk2nu_dir=args.dk2nu_dir, n_events=args.n_events, seed=args.seed,
        optimize=args.optimize, opt_iterations=args.opt_iterations,
        opt_batch=args.opt_batch, monoenergetic=args.monoenergetic,
        offshell=args.offshell, target_set=args.target_set, tile_n=args.tile_n)
