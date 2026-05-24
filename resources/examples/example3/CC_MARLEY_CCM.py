import numpy as np
import siren
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--target", type=int, default=0,
    help="0 for upper tungsten target, 1 for lower tungsten target",
)
args = parser.parse_args()

detector_model = siren.load_detector("CCM")
fiducial = siren.get_fiducial_volume("CCM")

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

sim = siren.Simulation(
    n_events=10000,
    detector=detector_model,
    primary="NuE",
    interactions="MarleyCrossSection",
    process="CC",
    energy=siren.distributions.PiDARNuEDistribution(),
    injection_direction=siren.dist.Cone(beam_dir, opening_angle),
    physical_direction=siren.dist.IsotropicDirection(),
    position=siren.distributions.FixedTargetPositionDistribution(
        target_cylinder, fiducial, 25,
    ),
)

results = sim.run()
results.summary()
results.save("output/CCM_MARLEY")
