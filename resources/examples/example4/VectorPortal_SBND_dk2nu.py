"""
Vector-portal dark matter at SBND, driven by BNB dk2nu pion files.

Off-shell chain (4 vertices, virtual chi'):

    pi+ -> mu+ nu V1;   V1 -> chi chi;   chi Ar -> chi V1_sig Ar;
    V1_sig -> e+ e-

The pion vertex takes its kinematics from a PrimaryExternalDistribution
built by the SIREN dk2nu reader (decay point, momentum, decay time, and
per-POT importance weight in a "weight" column) and is weighted in Fixed
mode, since the file already accounts for the decay having happened
there. Every vertex is a siren.Vertex carrying its interaction model, its
vertex distributions, a channels-algebra phase space (detector-directed
mixed with the physical channel), and declarative expand rules replacing
a hand-written stopping condition. The off-shell scattering aims the
signal V1 at the SBND active volume with a Breit-Wigner pair mass on the
virtual chi' resonance, and its vertex is bounded to the active volume;
the Weighter removes all of the biasing. The on-shell (5-vertex) variant
with an in-tree monoenergetic pion source is DuttaKim_SBND_full_chain.py.

Reference: Dutta et al., PRL 129, 111803 (2022) [arXiv:2110.11944]

Requires dk2nu beam-simulation ROOT files (not shipped with SIREN); set
DK2NU_DIR to their directory.

Usage:
    python VectorPortal_SBND_dk2nu.py [--events N] [--seed S]
        [--dk2nu-dir DIR]
"""

import argparse
import glob
import os
import sys

import numpy as np

import siren
from siren import _util, channels, dataclasses, distributions, expand
from siren import dk2nu as DK2NU
from siren.Injector import Injector
from siren.Weighter import Weighter

_PROC = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
_MESON = _util.load_module("DuttaKim_MesonProduction",
                           os.path.join(_PROC, "MesonProduction.py"))
_VP = _util.load_module("DuttaKim_VectorPortal",
                        os.path.join(_PROC, "VectorPortal.py"))

# Model parameters (Dutta et al. Table I, double-mediator).
M_PION = 0.13957039
M_CHI, M_CHI_PRIME = 8e-3, 50e-3
M_V1, M_V2 = 17e-3, 200e-3
G_D, EPSILON_1, EPSILON_2, G_MU = 1.0, 7e-5, 1e-4, 1e-3

PION = dataclasses.Particle.ParticleType(211)
V1_PROD = siren.particles.define("V1_prod", 5922, M_V1)
CHI = siren.particles.define("chi", 5917, M_CHI)
V1_SIG = siren.particles.define("V1_sig", 5923, M_V1)

# SBND active argon volume in detector coordinates, from the volTPCActive
# sectors of sbnd_v02_06.gdml as placed by load_detector("SBN", "SBND").
# The BNB proton target sits at (-0.74, 0.59, -112.92) m in this frame.
ACTIVE_CENTER = (0.0, 0.59, -0.415)
ACTIVE_WIDTHS = (4.026, 4.074645, 5.01)


def build_models():
    """Physics models for the off-shell chain.

    All four models are the unbiased physical ones; the detector-directed
    biasing lives in the channels algebra on each Vertex, where the
    Weighter can see and remove it.
    """
    pion_decay = _MESON.MesonThreeBodySIRENDecay(
        m_meson=M_PION, m_mediator=M_V1, g_mu=G_MU,
        pdgid_meson=211, pdgid_lepton=-13,
        pdgid_neutrino=14, pdgid_mediator=5922)
    v1_to_chi = _VP.DarkPhotonToChiDecay(
        M_V1, M_CHI, G_D, pdgid_V1=5922, pdgid_chi=5917)
    offshell_xs = _VP.VectorPortalOffShellXS(
        M_CHI, M_CHI_PRIME, M_V1, M_V2, G_D, EPSILON_2, pdgid_V1=5923)
    visible_decay = _VP.DarkPhotonDecay(M_V1, EPSILON_1, pdgid_V1=5923)
    return pion_decay, v1_to_chi, offshell_xs, visible_decay


def build_sX_cdf(pion_decay, n_nodes=257):
    """CDF of the physical s_X = M^2(mu, nu) marginal for the directed
    primary channel's pair-mass proposal (exact for any masses)."""
    d = pion_decay._decay
    s_min, s_max = d.m_l ** 2, (d.m_M - d.m_phi) ** 2
    eps = 1e-9 * (s_max - s_min)
    s = np.linspace(s_min + eps, s_max - eps, n_nodes)
    E_V1 = (d.m_M ** 2 + d.m_phi ** 2 - s) / (2.0 * d.m_M)
    dens = np.clip(np.asarray(d.differential_decay_rate(E_V1), float), 0.0, None)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(s))])
    return s.tolist(), cdf.tolist()


def build_vertices(models, external, fiducial, max_length=150.0):
    """The chain as spec-form Vertex objects with channels algebra."""
    pion_decay, v1_to_chi, offshell_xs, visible_decay = models
    sx = channels.PairMass.tabulated(*build_sX_cdf(pion_decay))
    # The biasing Breit-Wigner mirrors the physical chi' resonance width.
    chip_width = offshell_xs._chi_prime_width()
    sv = distributions.SecondaryPhysicalVertexDistribution

    primary = siren.Vertex(
        PION, pion_decay,
        distributions=[external],
        physical=[external],
        # The dk2nu row already accounts for this pion decaying at the
        # recorded vertex (the importance weight carries the decay
        # probability), so the parent vertex is weighted in Fixed mode.
        weighting=siren.Fixed(),
        kinematics=0.98 * channels.toward_3body("V1_prod", fiducial,
                                                strategy="direct", pair_mass=sx)
                   + 0.02 * channels.physical(),
        expand=(expand.child("V1_prod"),))

    v1 = siren.Vertex(
        "V1_prod", v1_to_chi, position=sv(),
        # daughter by index: the signature is V1 -> chi chi, so "chi" is
        # ambiguous; direct the first chi and expand only that one below
        kinematics=0.98 * channels.toward(0, fiducial)
                   + 0.02 * channels.physical(),
        expand=(expand.child("chi", index=0),))

    chi = siren.Vertex(
        "chi", offshell_xs,
        # bounded: confine the scattering to the active volume. The chi ray
        # starts at the beamline decay vertex, so the bound must cover the
        # full standoff for the fiducial clip to engage.
        position=siren.dist.BoundedVertex(fiducial, max_length),
        # 2->3 scatter: direct V1_sig (pair index 1) at the detector with a
        # Breit-Wigner pair mass on the virtual chi' resonance
        kinematics=0.98 * channels.toward_3body(
            1, fiducial, spectator=2, strategy="recursive",
            pair_mass=channels.PairMass.breit_wigner(M_CHI_PRIME, chip_width))
                   + 0.02 * channels.physical(),
        expand=(expand.child("V1_sig"),))

    v1s = siren.Vertex(
        "V1_sig", visible_decay, position=sv(),
        kinematics=0.5 * channels.physical() + 0.5 * channels.isotropic(0),
        expand=(expand.depth_below(0),))  # terminal: e+ e- are final

    return primary, (v1, chi, v1s)


def main():
    ap = argparse.ArgumentParser(
        description="Vector-portal dark matter at SBND from BNB dk2nu pions")
    ap.add_argument("--events", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tune", action="store_true",
                    help="Kleiss-Pittau channel-weight tuning before production")
    ap.add_argument("--dk2nu-dir", default=os.environ.get(
        "DK2NU_DIR",
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", "..", "..", "..", "..",
                     "g4bnb", "sources", "G4BNB")),
        help="directory holding the G4BNB dk2nu ROOT files")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.dk2nu_dir, "*dk2nu*.root")))
    if not files:
        files = sorted(glob.glob(os.path.join(args.dk2nu_dir, "nubeam*.root")))
    if not files:
        print("No dk2nu files found in %s" % args.dk2nu_dir)
        sys.exit(1)

    print("Reading %d dk2nu file(s) ..." % len(files))
    data = DK2NU.read_dk2nu(files, parent_pdg=[DK2NU.PTYPE_PIPLUS])
    DK2NU.print_summary(data)

    print("Loading the SBN detector model (SBND) ...")
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    external = DK2NU.dk2nu_to_primary_distribution(data, detector_model)
    fiducial = siren.geometry.Box(widths=ACTIVE_WIDTHS, center=ACTIVE_CENTER)

    models = build_models()
    pion_decay, v1_to_chi, offshell_xs, visible_decay = models
    print("  Pion 3-body width: %.4e GeV" % pion_decay._total_width)
    print("  V1->chi chi width: %.4e GeV" % v1_to_chi.total_width())
    print("  V1->e+e- width: %.4e GeV" % visible_decay.total_width())
    print("  chi scattering threshold: %.4f GeV" % offshell_xs._ups.Ethreshold)
    print("  Loaded %d pion rows" % external.GetPhysicalNumEvents())

    primary, secondaries = build_vertices(models, external, fiducial)

    injector = Injector(detector=detector_model, primary=primary,
                        secondaries=secondaries,
                        events=args.events, seed=args.seed)
    weighter = Weighter(injector, primary_physical=primary.physical)

    if args.tune:
        print(siren.tune.tune(injector, weighter, events=200, rounds=3))
        injector.reset()

    results = siren.generate(injector, weighter, events=args.events,
                             on_shortfall="warn")
    results.summary()
    print(results.variance_report())
    # SM secondaries pruned by the expand rules (mu+, nu, and the recoil
    # nucleus) appear in the ledger as UnregisteredSecondaryType; that is
    # the declared pruning, not a misconfiguration.
    print(injector.report())

    # Chain completion counters: a chi scattering vertex and a visible
    # V1 decay mark a complete signal topology.
    n_chi_scatter = 0
    n_v1_decay = 0
    for event in results.events:
        types = [d.record.signature.primary_type for d in event.tree]
        n_chi_scatter += types.count(CHI)
        n_v1_decay += types.count(V1_SIG)
    print("Events with a chi scattering vertex: %d" % n_chi_scatter)
    print("Events with a V1 -> e+ e- signal decay: %d" % n_v1_decay)

    # POT bookkeeping: each dk2nu row is one recorded decay, and the file
    # header carries the protons on target that produced the whole set.
    total_pot = data["pot"]
    n_rows = len(data["E"])
    if total_pot > 0 and n_rows > 0:
        pot_per_row = total_pot / n_rows
        print("POT per dk2nu pion row: %.3e" % pot_per_row)
        print("POT represented by %d sampled rows: %.3e"
              % (results.attempts, results.attempts * pot_per_row))


if __name__ == "__main__":
    main()
