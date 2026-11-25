import siren

sim = siren.Simulation(
    events=int(1e5),
    detector="IceCube",
    primary="NuMu",
    interactions="CSMSDISSplines",
    targets="Nucleon",
    process="CC",
    energy=siren.dist.PowerLaw(2, 1e3, 1e6),
    direction=siren.dist.IsotropicDirection(),
    position=siren.dist.ColumnDepth(
        600, 600.0, siren.distributions.LeptonDepthFunction()
    ),
)

results = sim.run()
results.save("output/IceCube")
