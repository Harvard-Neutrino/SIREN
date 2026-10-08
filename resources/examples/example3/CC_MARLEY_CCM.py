"""Charged-current electron-neutrino interactions at CCM with MARLEY.

Electron neutrinos from pi+ decay at rest are injected from one of the two
tungsten targets into a cone covering the detector; MARLEY provides the
low-energy CC cross section on argon. The primary ``siren.Vertex`` shares the
pi-DAR energy spectrum between injection and weighting and declares the
isotropic physical direction, so the weight corrects the cone sampling.
"""
import argparse

import numpy as np

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--target", type=int, default=0,
    help="0 for upper tungsten target, 1 for lower tungsten target",
)
parser.add_argument("--events", type=int, default=10000)
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default="output/CCM_MARLEY")
args = parser.parse_args()

NuE = siren.particles.NuE

detector_model = siren.load_detector("CCM")
fiducial = siren.get_fiducial_volume("CCM")

interactions = siren.load_processes(
    "MarleyCrossSection", primary_types=[NuE], process_types=["CC"],
).primary[NuE]

if args.target == 0:
    target_origin = siren.math.Vector3D(0, 0, 0.1375)
    target_cylinder = siren.geometry.Cylinder(
        siren.geometry.Placement(target_origin - siren.math.Vector3D(23, 0, -0.65)),
        0.05, 0.0, 0.091,
    )
else:
    target_origin = siren.math.Vector3D(0, 0, -0.241)
    target_cylinder = siren.geometry.Cylinder(
        siren.geometry.Placement(target_origin - siren.math.Vector3D(23, 0, -0.65)),
        0.05, 0.0, 0.298,
    )

detector_origin = siren.math.Vector3D(23, 0, -0.65)
beam_dir = detector_origin - target_origin
beam_dir.normalize()
opening_angle = np.arcsin(1.21 / 23.0)

energy = siren.dist.PiDARNuEDistribution()

primary = siren.Vertex(
    NuE, interactions,
    distributions=[
        siren.dist.PrimaryMass(0),
        energy,
        siren.dist.Cone(beam_dir, opening_angle),
        siren.dist.FixedTargetPositionDistribution(target_cylinder, fiducial, 25),
    ],
    physical=[energy, siren.dist.IsotropicDirection()],
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector=detector_model,
    primary=primary)
results = sim.run()
results.summary()
results.save(args.output)
