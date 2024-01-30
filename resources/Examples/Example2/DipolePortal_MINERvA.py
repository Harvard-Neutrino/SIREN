import os
import sys
import numpy as np

import leptoninjector as LI
from leptoninjector import _util
from leptoninjector.LIController import LIController

import DarkNews

darknews_version = _util.normalize_version(DarkNews.__version__)

resources_dir = _util.resource_package_dir()

# Define a DarkNews model
model_kwargs = {
    "m4": 0.47,  # 0.140,
    "mu_tr_mu4": 1.25e-6,  # 1e-6, # GeV^-1
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
experiment = "MINERvA"

# Define the controller
controller = LIController(events_to_inject, experiment)

# Particle to inject
primary_type = LI.dataclasses.Particle.ParticleType.NuMu

xs_path = _util.get_cross_section_model_path(f"DarkNewsTables-v{darknews_version}", must_exist=False)
# Define DarkNews Model
table_dir = os.path.join(
    xs_path,
    "Dipole_M%2.2f_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"]),
)
controller.InputDarkNewsModel(primary_type, table_dir, model_kwargs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
flux_file = os.path.join(
    resources_dir,
    "Fluxes",
    "NUMI_Flux_Tables",
    "ME_FHC_numu.txt",
)
edist = LI.distributions.TabulatedFluxDistribution(flux_file, True)
edist_gen = LI.distributions.TabulatedFluxDistribution(
    1.05 * model_kwargs["m4"], 10, flux_file, False
)
primary_injection_distributions["energy"] = edist_gen
primary_physical_distributions["energy"] = edist

# Flux normalization: go from cm^-2 to m^-2
flux_units = LI.distributions.NormalizationConstant(1e4)
primary_physical_distributions["normalization"] = flux_units

# direction distribution
direction_distribution = LI.distributions.FixedDirection(LI.math.Vector3D(0, 0, 1.0))
primary_injection_distributions["direction"] = direction_distribution
primary_physical_distributions["direction"] = direction_distribution

# position distribution
decay_range_func = LI.distributions.DecayRangeFunction(
    model_kwargs["m4"], controller.DN_min_decay_width, 3, 240
)
position_distribution = LI.distributions.RangePositionDistribution(
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
    return secondary_type != LI.dataclasses.Particle.ParticleType.N4

controller.injector.SetStoppingCondition(stop)

events = controller.GenerateEvents()

controller.SaveEvents(
    "output/MINERvA_Dipole_M%2.2f_mu%2.2e_example.hdf5"
    % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])
)
