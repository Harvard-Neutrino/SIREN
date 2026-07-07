import siren

# The DUNE far-detector model requires choosing a module design:
# "HD" (horizontal drift) or "VD" (vertical drift).
detector = siren.load_detector("DUNEFD", detector="HD")

sim = siren.Simulation(
    events=int(1e5),
    detector=detector,
    primary="NuMu",
    interactions="CSMSDISSplines",
    targets="Nucleon",
    process="CC",
    energy=siren.dist.PowerLaw(1, 1e3, 1e6),
    direction=siren.dist.IsotropicDirection(),
    position=siren.dist.ColumnDepth(
        60, 60.0, siren.distributions.LeptonDepthFunction()
    ),
)

results = sim.run()
results.save("output/DUNE")
