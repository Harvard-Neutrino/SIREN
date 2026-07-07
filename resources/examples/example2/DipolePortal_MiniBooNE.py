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

detector_model = siren.load_detector("MiniBooNE")

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

def _total_decay_width(decay):
    # TotalDecayWidth(record) only returns a nonzero width when the record's
    # signature matches the decay's own signature; a default-constructed
    # InteractionRecord has an "unknown" signature and always yields zero.
    record = siren.dataclasses.InteractionRecord()
    record.signature = decay.GetPossibleSignatures()[0]
    return decay.TotalDecayWidth(record)

sim = siren.Simulation(
    events=100000,
    detector=detector_model,
    primary="NuMu",
    interactions=bundle.primary[siren.particles.NuMu],
    injection_energy=siren.load_flux(
        "BNB", tag="FHC_numu",
        min_energy=model_kwargs["m4"], max_energy=10,
        physically_normalized=False,
    ),
    physical_energy=siren.load_flux(
        "BNB", tag="FHC_numu", physically_normalized=True,
    ),
    direction=siren.dist.FixedDirection([0, 0, 1]),
    position=siren.dist.RangePosition(
        6.2, 6.2,
        siren.dist.DecayRange(
            model_kwargs["m4"],
            min(_total_decay_width(d)
                for d in bundle.secondary[siren.particles.N4]),
            3, 541,
        ),
        siren.get_detector_model_targets(detector_model),
    ),
    secondary_interactions=bundle.secondary,
    secondary_position=siren.dist.BoundedVertex(
        siren.get_fiducial_volume("MiniBooNE"), 25
    ),
    stopping_condition=lambda tree, datum, i: (
        datum.record.signature.secondary_types[i] != siren.particles.N4
    ),
)

results = sim.run()

os.makedirs("output", exist_ok=True)
results.save(
    "output/MiniBooNE_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])
)
