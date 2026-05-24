import siren

sim = siren.Simulation(
    n_events=int(1e5),
    detector="DUNEFD",
    primary="NuMu",
    interactions="CSMSDISSplines",
    target="Nucleon",
    process="CC",
    energy=siren.dist.PowerLaw(1, 1e3, 1e6),
    direction=siren.dist.IsotropicDirection(),
    position=siren.dist.ColumnDepth(
        60, 60.0, siren.distributions.LeptonDepthFunction()
    ),
)

results = sim.run()
results.save("output/DUNE")
