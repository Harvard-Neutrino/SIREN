import os
import siren

model_kwargs = {
    "m4": 0.47,
    "mu_tr_mu4": 2.50e-6,
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

detector_model = siren.load_detector("MINERvA")

dn_version = siren.utilities.darknews_version()
table_name = f"DarkNewsTables-v{dn_version}/"
table_name += "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])

bundle = siren.load_processes(
    "DarkNewsTables",
    primary_type=siren.particles.NuMu,
    detector_model=detector_model,
    table_name=table_name,
    **model_kwargs,
)

sim = siren.Simulation(
    events=100000,
    detector=detector_model,
    primary="NuMu",
    interactions=bundle.primary[siren.particles.NuMu],
    injection_energy=siren.load_flux(
        "NUMI", tag="FHC_ME_numu",
        min_energy=model_kwargs["m4"], max_energy=20,
        physically_normalized=False,
    ),
    physical_energy=siren.load_flux(
        "NUMI", tag="FHC_ME_numu", physically_normalized=True,
    ),
    direction=siren.dist.FixedDirection([0, 0, 1]),
    position=siren.dist.RangePosition(
        1.24, 5.0,
        siren.dist.DecayRange(
            model_kwargs["m4"],
            min(d.TotalDecayWidth(siren.dataclasses.InteractionRecord())
                for d in bundle.secondary[siren.particles.N4]),
            3, 240,
        ),
        siren.get_detector_model_targets(detector_model),
    ),
    secondary_interactions=bundle.secondary,
    secondary_position=siren.dist.BoundedVertex(
        siren.get_fiducial_volume("MINERvA"), 25
    ),
    stopping_condition=lambda tree, datum, i: (
        datum.record.signature.secondary_types[i] != siren.particles.N4
    ),
)

results = sim.run()

os.makedirs("output", exist_ok=True)
results.save(
    "output/MINERvA_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])
)
