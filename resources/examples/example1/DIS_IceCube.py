"""Muon-neutrino charged-current DIS in IceCube.

One ``siren.Vertex`` declares the primary: ``distributions`` sample the
injection (mass, energy, direction, and a column-depth position along the
muon range) and ``physical`` names the flux and direction factors used for
weighting. The energy and direction objects are shared between the two
lists, so they cancel in the weight, which then corrects only the position
sampling and carries the interaction probability.
"""
import argparse

import siren

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=int(1e5))
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default="output/IceCube")
args = parser.parse_args()

NuMu = siren.particles.NuMu

# CSMS DIS cross sections: NuMu CC on an isoscalar nucleon target.
interactions = siren.load_processes(
    "CSMSDISSplines",
    primary_types=[NuMu],
    target_types=[siren.particles.Nucleon],
    isoscalar=True,
    process_types=["CC"],
).primary[NuMu]

energy = siren.dist.PowerLaw(2, 1e3, 1e6)
direction = siren.dist.IsotropicDirection()

primary = siren.Vertex(
    NuMu, interactions,
    distributions=[
        siren.dist.PrimaryMass(0),
        energy,
        direction,
        siren.dist.ColumnDepth(
            600, 600.0, siren.distributions.LeptonDepthFunction()),
    ],
    physical=[energy, direction],
)

sim = siren.Simulation(
    events=args.events, seed=args.seed, detector="IceCube", primary=primary)
results = sim.run()
results.summary()
results.save(args.output)
