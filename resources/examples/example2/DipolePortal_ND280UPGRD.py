import os
import siren
from siren._util import get_tabulated_flux_file

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

detector_model = siren.load_detector("ND280UPGRD")

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

flux_file = get_tabulated_flux_file("T2K_NEAR", "PLUS_numu")

sim = siren.Simulation(
    n_events=100000,
    detector=detector_model,
    primary="NuMu",
    interactions=bundle.primary[siren.particles.NuMu],
    injection_energy=siren.distributions.TabulatedFluxDistribution(
        model_kwargs["m4"], 20, flux_file, False,
    ),
    physical_energy=siren.distributions.TabulatedFluxDistribution(
        flux_file, True,
    ),
    direction=siren.dist.FixedDirection([0, 0, 1]),
    position=siren.dist.RangePosition(
        5.0, 9.0,
        siren.dist.DecayRange(
            model_kwargs["m4"],
            min(d.TotalDecayWidth(siren.dataclasses.InteractionRecord())
                for d in bundle.secondary[siren.particles.N4]),
            3, 240,
        ),
        set(bundle.primary.keys()),
    ),
    secondary_interactions=bundle.secondary,
    secondary_position=siren.dist.BoundedVertex(
        siren.get_fiducial_volume("ND280UPGRD"), 25
    ),
    stopping_condition=lambda tree, datum, i: (
        datum.record.signature.secondary_types[i] != siren.particles.N4
    ),
)

results = sim.run()

os.makedirs("output", exist_ok=True)
results.save(
    "output/ND280UPGRD_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])
)
