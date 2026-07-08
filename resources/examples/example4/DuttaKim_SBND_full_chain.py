"""
Dutta-Kim vector-portal chain at SBND, spec-form.

On-shell chain (5 vertices):

    pi+ -> mu+ nu V1;   V1 -> chi chi;   chi Ar -> chi' Ar;
    chi' -> chi V1_sig;   V1_sig -> e+ e-

Each vertex is a siren.Vertex carrying its interaction model, its vertex
distributions, a channels-algebra phase space (directed-toward-detector
mixed with the physical channel), and declarative expand rules replacing
a hand-written stopping condition. Runs end-to-end with in-tree
resources: SBND geometry and a monoenergetic 2 GeV pion. A dk2nu-driven
pion flux needs external beam files (see VectorPortal_SBND_dk2nu.py).
The off-shell (4-vertex, virtual chi') variant is assembled in
tests/python/test_dutta_kim_chain.py.

Usage:
    python DuttaKim_SBND_full_chain.py [--events N] [--seed S] [--tune]
"""

import argparse
import os

import numpy as np

import siren
from siren import _util, channels, dataclasses, distributions, expand
from siren.Injector import Injector
from siren.Weighter import Weighter

_PROC = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
_MESON = _util.load_module("DuttaKim_MesonProduction",
                           os.path.join(_PROC, "MesonProduction.py"))
_VP = _util.load_module("DuttaKim_VectorPortal",
                        os.path.join(_PROC, "VectorPortal.py"))

M_PION = 0.13957039
M_CHI, M_CHI_PRIME = 8e-3, 50e-3
M_V1, M_V2 = 17e-3, 200e-3
M_ARGON40 = 37.215
G_D, EPSILON_1, EPSILON_2, G_MU = 1.0, 7e-5, 1e-4, 1e-3

PION = dataclasses.Particle.ParticleType(211)
V1_PROD = siren.particles.define("V1_prod", 5922, M_V1)
CHI = siren.particles.define("chi", 5917, M_CHI)
CHI_PRIME = siren.particles.define("chi_prime", 5918, M_CHI_PRIME)
V1_SIG = siren.particles.define("V1_sig", 5923, M_V1)


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
            V1_SIG: [visible_decay],
        },
        "models": {
            "v1_to_chi": v1_to_chi,
            "upscatter": upscatter,
            "chi_prime_decay": chi_prime_decay,
            "visible_decay": visible_decay,
        },
    }


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


def build_vertices(models, fiducial, pion_energy=2.0):
    """The chain as spec-form Vertex objects with channels algebra."""
    m = models["models"]
    sx = channels.PairMass.tabulated(*build_sX_cdf(models["pion_decay"]))
    sv = distributions.SecondaryPhysicalVertexDistribution

    primary = siren.Vertex(
        PION, models["pion_decay"],
        distributions=[
            distributions.PrimaryMass(M_PION),
            distributions.Monoenergetic(pion_energy),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
            distributions.PointSourcePositionDistribution([0, 0, 0], 300.0),
        ],
        physical=[
            distributions.Monoenergetic(pion_energy),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        ],
        kinematics=0.98 * channels.toward_3body("V1_prod", fiducial,
                                                strategy="direct", pair_mass=sx)
                   + 0.02 * channels.physical(),
        expand=(expand.child("V1_prod"),))

    v1 = siren.Vertex(
        "V1_prod", m["v1_to_chi"], position=sv(),
        # daughter by index: the signature is V1 -> chi chi, so "chi" is
        # ambiguous; direct the first chi and expand only that one below
        kinematics=0.98 * channels.toward(0, fiducial)
                   + 0.02 * channels.physical(),
        expand=(expand.child("chi", index=0),))

    chi = siren.Vertex(
        "chi", m["upscatter"],
        # bounded: confine the upscatter to the fiducial volume
        position=distributions.SecondaryBoundedVertexDistribution(fiducial),
        kinematics=channels.physical(),
        expand=(expand.child("chi_prime"),))

    chip = siren.Vertex(
        "chi_prime", m["chi_prime_decay"], position=sv(),
        kinematics=channels.physical(),
        expand=(expand.child("V1_sig"),))

    v1s = siren.Vertex(
        "V1_sig", m["visible_decay"], position=sv(),
        kinematics=0.5 * channels.physical() + 0.5 * channels.isotropic(0),
        expand=(expand.depth_below(0),))  # terminal: e+ e- are final

    return primary, (v1, chi, chip, v1s)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events", type=int, default=25)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tune", action="store_true",
                    help="Kleiss-Pittau channel-weight tuning before production")
    args = ap.parse_args()

    detector = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(
        widths=(2 * (4.5 + 2.022), 2 * 4.074645, 2 * 5.010))

    models = build_onshell_models()
    primary, secondaries = build_vertices(models, fiducial)

    injector = Injector(detector=detector, primary=primary,
                        secondaries=secondaries,
                        events=args.events, seed=args.seed)
    weighter = Weighter(injector, primary_physical=primary.physical)

    if args.tune:
        print(siren.tune.tune(injector, weighter, events=200, rounds=3))
        injector.reset()

    results = siren.generate(injector, weighter,
                             events=args.events, on_shortfall="warn")
    results.summary()
    print(results.variance_report())
    # SM secondaries pruned by the expand rules (mu+, nu, e+ e-) appear in
    # the ledger as UnregisteredSecondaryType; that is the declared pruning,
    # not a misconfiguration.
    print(injector.report())


if __name__ == "__main__":
    main()
