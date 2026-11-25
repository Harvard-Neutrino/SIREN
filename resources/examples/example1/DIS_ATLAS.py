import siren

detector_model = siren.load_detector("ATLAS")
position = siren.get_volume_position_distribution_from_sector(
    detector_model, "tilecal"
)

sim = siren.Simulation(
    events=int(1e5),
    seed=99,
    detector=detector_model,
    primary="NuMu",
    interactions="CSMSDISSplines",
    targets="Nucleon",
    process="CC",
    injection_energy=siren.load_flux(
        "HE_SN", tag="numu", min_energy=100, max_energy=1e6,
        physically_normalized=True,
    ),
    physical_energy=siren.load_flux(
        "HE_SN", tag="numu", min_energy=100, max_energy=1e6,
        physically_normalized=False,
    ),
    direction=siren.dist.FixedDirection([0, 0, 1]),
    position=position,
)

results = sim.run()
results.save("output/ATLAS")
