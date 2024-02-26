import os
import numpy as np

import leptoninjector as LI
from leptoninjector.LIController import LIController

# Define a DarkNews model
model_kwargs = {
    "m4": 0.0235,
    "mu_tr_mu4": 3e-7,  # 1e-6, # GeV^-1
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

# Number of events to inject
events_to_inject = 10000

# Expeirment to run
experiment = "CCM"

# Define the controller
controller = LIController(events_to_inject, experiment)

# Particle to inject
primary_type = LI.dataclasses.Particle.ParticleType.NuMu

xs_path = LI.utilities.get_cross_section_model_path(f"DarkNewsTables-v{LI.utilities.darknews_version()}", must_exist=False)
# Define DarkNews Model
table_dir = os.path.join(
    xs_path,
    "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]),
)
controller.InputDarkNewsModel(primary_type, table_dir, model_kwargs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
nu_energy = 0.02965  # from pi+ DAR
edist = LI.distributions.Monoenergetic(nu_energy)
primary_injection_distributions["energy"] = edist
primary_physical_distributions["energy"] = edist

# Flux normalization:
# using the number quoted in 2105.14020, 4.74e9 nu/m^2/s / (6.2e14 POT/s) * 4*pi*20m^2 to get nu/POT
flux_units = LI.distributions.NormalizationConstant(3.76e-2)
primary_physical_distributions["flux_units"] = flux_units

# direction distribution: cone from lower W target
opening_angle = np.arctan(5 / 23.0)
# slightly larger than CCM
lower_target_origin = LI.math.Vector3D(0, 0, -0.241)
detector_origin = LI.math.Vector3D(23, 0, -0.65)
lower_dir = detector_origin - lower_target_origin
lower_dir.normalize()
lower_inj_ddist = LI.distributions.Cone(lower_dir, opening_angle)
phys_ddist = (
    LI.distributions.IsotropicDirection()
)  # truly we are isotropicprimary_injection_distributions['direction'] = direction_distribution
primary_injection_distributions["direction"] = lower_inj_ddist
primary_physical_distributions["direction"] = phys_ddist

# Position distribution: consider neutrinos from a point source
max_dist = 25
lower_pos_dist = LI.distributions.PointSourcePositionDistribution(
    lower_target_origin - detector_origin, max_dist, set(controller.GetDetectorModelTargets()[0])
)
primary_injection_distributions["position"] = lower_pos_dist

# SetProcesses
controller.SetProcesses(
    primary_type, primary_injection_distributions, primary_physical_distributions
)

controller.Initialize()

def stop(datum, i):
    secondary_type = datum.record.signature.secondary_types[i]
    return secondary_type != LI.dataclasses.Particle.ParticleType.N4

controller.injector.SetStoppingCondition(stop)

events = controller.GenerateEvents(fill_tables_at_exit=False)

os.makedirs("output", exist_ok=True)

controller.SaveEvents(
    "output/CCM_Dipole_M%2.2e_mu%2.2e_example"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]),
    fill_tables_at_exit=False
)
