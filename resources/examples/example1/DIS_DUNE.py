"""Muon-neutrino charged-current DIS in one DUNE far-detector module.

The DUNE far-detector model requires choosing a module design: "HD"
(horizontal drift) or "VD" (vertical drift). The primary ``siren.Vertex``
samples an E^-1 spectrum for injection and declares the same spectrum and
direction as the physical target, so the weight corrects only the
column-depth position sampling and carries the interaction probability.
"""
import argparse

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=int(1e5))
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default="output/DUNE")
args = parser.parse_args()

NuMu = siren.particles.NuMu

detector = siren.load_detector("DUNEFD", detector="HD")

interactions = siren.load_processes(
    "CSMSDISSplines",
    primary_types=[NuMu],
    target_types=[siren.particles.Nucleon],
    isoscalar=True,
    process_types=["CC"],
).primary[NuMu]

energy = siren.dist.PowerLaw(1, 1e3, 1e6)
direction = siren.dist.IsotropicDirection()

primary = siren.Vertex(
    NuMu, interactions,
    distributions=[
        siren.dist.PrimaryMass(0),
        energy,
        direction,
        siren.dist.ColumnDepth(
            60, 60.0, siren.distributions.LeptonDepthFunction()),
    ],
    physical=[energy, direction],
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector=detector, primary=primary)
results = sim.run()
results.summary()
results.save(args.output)
