import os

import leptoninjector as LI
from leptoninjector import _util
from leptoninjector.LIController import LIController

resources_dir = _util.resource_package_dir()

# Number of events to inject
events_to_inject = 10000

# Expeirment to run
experiment = "ATLAS"

# Define the controller
controller = LIController(events_to_inject, experiment, seed=99)

# Particle to inject
primary_type = LI.dataclasses.Particle.ParticleType.NuMu

cross_section_model = "CSMSDISSplines"

xsfiledir = _util.get_cross_section_model_path(cross_section_model)

# Cross Section Model
target_type = LI.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = LI.interactions.DISFromSpline(
    os.path.join(xsfiledir, "dsdxdy_nu_CC_iso.fits"),
    os.path.join(xsfiledir, "sigma_nu_CC_iso.fits"),
    [primary_type],
    [target_type], "m"
)

primary_xs = LI.interactions.InteractionCollection(primary_type, [DIS_xs])
controller.SetInteractions(primary_xs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
# HE SN flux from ATLAS paper
flux_file = _util.get_tabulated_flux_file("HE_SN","numu")
edist = LI.distributions.TabulatedFluxDistribution(100, 1e6, flux_file, True) #bool is whether flux is physical
primary_injection_distributions["energy"] = edist
primary_physical_distributions["energy"] = edist

# direction distribution
# let's just inject upwards
injection_dir = LI.math.Vector3D(0, 0, 1)
injection_dir.normalize()
direction_distribution = LI.distributions.FixedDirection(injection_dir)
primary_injection_distributions["direction"] = direction_distribution
primary_physical_distributions["direction"] = direction_distribution

# position distribution
position_distribution = controller.GetCylinderVolumePositionDistributionFromSector("tilecal")
primary_injection_distributions["position"] = position_distribution

# SetProcesses
controller.SetProcesses(
    primary_type, primary_injection_distributions, primary_physical_distributions
)

controller.Initialize()

events = controller.GenerateEvents()

controller.SaveEvents("output/ATLAS_DIS")
