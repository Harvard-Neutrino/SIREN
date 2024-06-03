import os

import siren
from siren.SIREN_Controller import SIREN_Controller

# Define a DarkNews model
model_kwargs = {
    "m4": 0.47,  # 0.140,
    "mu_tr_mu4": 2.50e-6,  # 1e-6, # GeV^-1
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

# Number of events to inject
events_to_inject = 100000

# Expeirment to run
experiment = "MINERvA"

# Define the controller
controller = SIREN_Controller(events_to_inject, experiment)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

xs_path = siren.utilities.get_cross_section_model_path(f"DarkNewsTables-v{siren.utilities.darknews_version()}", must_exist=False)
# Define DarkNews Model
table_dir = os.path.join(
    xs_path,
    "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]),
)
controller.InputDarkNewsModel(primary_type, table_dir, **model_kwargs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
flux_file = siren.utilities.get_tabulated_flux_file("NUMI","FHC_ME_numu")
edist = siren.distributions.TabulatedFluxDistribution(flux_file, True)
edist_gen = siren.distributions.TabulatedFluxDistribution(
    model_kwargs["m4"], 20, flux_file, False
)
primary_injection_distributions["energy"] = edist_gen
primary_physical_distributions["energy"] = edist

# direction distribution
direction_distribution = siren.distributions.FixedDirection(siren.math.Vector3D(0, 0, 1.0))
primary_injection_distributions["direction"] = direction_distribution
primary_physical_distributions["direction"] = direction_distribution

# position distribution
decay_range_func = siren.distributions.DecayRangeFunction(
    model_kwargs["m4"], controller.DN_min_decay_width, 3, 240
)
position_distribution = siren.distributions.RangePositionDistribution(
    1.24, 5.0, decay_range_func, set(controller.GetDetectorModelTargets()[0])
)
primary_injection_distributions["position"] = position_distribution

# SetProcesses
controller.SetProcesses(
    primary_type, primary_injection_distributions, primary_physical_distributions
)

controller.Initialize()

def stop(datum, i):
    secondary_type = datum.record.signature.secondary_types[i]
    return secondary_type != siren.dataclasses.Particle.ParticleType.N4

controller.injector.SetStoppingCondition(stop)

events = controller.GenerateEvents(fill_tables_at_exit=False)

os.makedirs("output", exist_ok=True)

controller.SaveEvents(
    "output/MINERvA_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]),
    fill_tables_at_exit=True
)
