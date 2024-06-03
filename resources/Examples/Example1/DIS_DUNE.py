import os

import siren
from siren.SIREN_Controller import SIREN_Controller

# Number of events to inject
events_to_inject = int(1e5)

# Expeirment to run
experiment = "DUNEFD"

# Define the controller
controller = SIREN_Controller(events_to_inject, experiment)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

cross_section_model = "CSMSDISSplines"

xsfiledir = siren.utilities.get_cross_section_model_path(cross_section_model)

# Cross Section Model
target_type = siren.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = siren.interactions.DISFromSpline(
    os.path.join(xsfiledir, "dsdxdy_nu_CC_iso.fits"),
    os.path.join(xsfiledir, "sigma_nu_CC_iso.fits"),
    [primary_type],
    [target_type], "m"
)

primary_xs = siren.interactions.InteractionCollection(primary_type, [DIS_xs])
controller.SetInteractions(primary_xs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
edist = siren.distributions.PowerLaw(1, 1e3, 1e6)
primary_injection_distributions["energy"] = edist
primary_physical_distributions["energy"] = edist

# direction distribution
direction_distribution = siren.distributions.IsotropicDirection()
primary_injection_distributions["direction"] = direction_distribution
primary_physical_distributions["direction"] = direction_distribution

# position distribution
muon_range_func = siren.distributions.LeptonDepthFunction()
position_distribution = siren.distributions.ColumnDepthPositionDistribution(
    60, 60.0, muon_range_func, set(controller.GetDetectorModelTargets()[0])
)
primary_injection_distributions["position"] = position_distribution

# SetProcesses
controller.SetProcesses(
    primary_type, primary_injection_distributions, primary_physical_distributions
)

controller.Initialize()

events = controller.GenerateEvents()

os.makedirs("output", exist_ok=True)

controller.SaveEvents("output/DUNE_DIS")