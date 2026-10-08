"""Dipole-portal heavy neutral lepton (N4) at CCM.

A monoenergetic pi+ decay-at-rest muon neutrino upscatters to N4 through the
DarkNews dipole portal (primary vertex) and the N4 decays to a photon inside
the detector (secondary vertex). The primary's ``expand`` rule recurses only
into the N4; the N4 vertex is terminal. Injection points the neutrino from
the lower tungsten target into a cone covering the detector, and weighting
uses the isotropic physical direction.
"""
import argparse

import numpy as np

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=1)
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default="output/CCM_Dipole")
args = parser.parse_args()

model_kwargs = {
    "m4": 0.0235,
    "mu_tr_mu4": 6e-7,
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

NuMu, N4 = siren.particles.NuMu, siren.particles.N4

detector_model = siren.load_detector("CCM")
fiducial = siren.get_fiducial_volume("CCM")

dn_version = siren.utilities.darknews_version()
table_name = f"DarkNewsTables-v{dn_version}/"
table_name += "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])

bundle = siren.load_processes(
    "DarkNewsTables",
    primary_type=NuMu,
    detector_model=detector_model,
    table_name=table_name,
    **model_kwargs,
)

target_origin = siren.math.Vector3D(0, 0, -0.241)
detector_origin = siren.math.Vector3D(23, 0, -0.65)
beam_dir = detector_origin - target_origin
beam_dir.normalize()

energy = siren.dist.Monoenergetic(0.02965)  # pi+ decay at rest

primary = siren.Vertex(
    NuMu, bundle.primary[NuMu],
    distributions=[
        siren.dist.PrimaryMass(0),
        energy,
        siren.dist.Cone(beam_dir, np.arctan(5 / 23.0)),
        siren.dist.PointSource(target_origin - detector_origin, 25),
    ],
    physical=[energy, siren.dist.IsotropicDirection()],
    expand=(siren.expand.child("N4"),),
)

hnl = siren.Vertex(
    N4, bundle.secondary[N4],
    position=siren.dist.BoundedVertex(fiducial, 25),
    expand=(siren.expand.depth_below(0),),  # terminal: photon and neutrino are final
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector=detector_model,
    primary=primary, secondaries=[hnl])
results = sim.run()
results.summary()
# save_hepmc3=True additionally writes <output>.hepmc3 (HepMC3/NuHepMC).
results.save(args.output, save_hepmc3=True)
