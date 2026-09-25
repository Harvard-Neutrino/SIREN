"""Muon-neutrino charged-current DIS in the ATLAS tile calorimeter from a
high-energy supernova flux.

Injection samples the tabulated flux shape (the sampler always uses the
unit-normalized table), points the neutrinos along +z, and draws vertices
uniformly in the ``tilecal`` sector volume. The physical target reuses the
same table and direction; ``physically_normalized=False`` on the physical
side weights per unit flux rather than in the table's flux units.
"""
import argparse

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=int(1e5))
parser.add_argument("--seed", type=int, default=99)
parser.add_argument("--output", default="output/ATLAS")
args = parser.parse_args()

NuMu = siren.particles.NuMu

detector_model = siren.load_detector("ATLAS")

interactions = siren.load_processes(
    "CSMSDISSplines",
    primary_types=[NuMu],
    target_types=[siren.particles.Nucleon],
    isoscalar=True,
    process_types=["CC"],
).primary[NuMu]

direction = siren.dist.FixedDirection([0, 0, 1])

primary = siren.Vertex(
    NuMu, interactions,
    distributions=[
        siren.dist.PrimaryMass(0),
        siren.load_flux(
            "HE_SN", tag="numu", min_energy=100, max_energy=1e6,
            physically_normalized=True),
        direction,
        siren.get_volume_position_distribution_from_sector(
            detector_model, "tilecal"),
    ],
    physical=[
        siren.load_flux(
            "HE_SN", tag="numu", min_energy=100, max_energy=1e6,
            physically_normalized=False),
        direction,
    ],
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector=detector_model,
    primary=primary)
results = sim.run()
results.summary()
results.save(args.output)
