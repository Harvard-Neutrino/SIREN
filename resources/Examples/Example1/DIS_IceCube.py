import numpy as np
import os

import leptoninjector as LI
import sys
from leptoninjector import _util
from leptoninjector.LIController import LIController

resources_dir = _util.resource_package_dir()

# Number of events to inject
events_to_inject = 1000

# Expeirment to run
experiment = "IceCube"

# Define the controller
controller = LIController(events_to_inject, experiment)

# Particle to inject
primary_type = LI.dataclasses.Particle.ParticleType.NuMu

# Cross Section Model
xsfiledir = os.path.join(resources_dir, "CrossSectionTables", "DISSplines")
target_type = LI.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = LI.interactions.DISFromSpline(
    os.path.join(xsfiledir, "test_xs.fits"),
    os.path.join(xsfiledir, "test_xs_total.fits"),
    [primary_type],
    [target_type],
)

primary_xs = LI.interactions.InteractionCollection(primary_type, [DIS_xs])
controller.SetInteractions(primary_xs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
edist = LI.distributions.PowerLaw(2, 1e3, 1e6)
primary_injection_distributions["energy"] = edist
primary_physical_distributions["energy"] = edist

# direction distribution
direction_distribution = LI.distributions.IsotropicDirection()
primary_injection_distributions["direction"] = direction_distribution
primary_physical_distributions["direction"] = direction_distribution

# position distribution
muon_range_func = LI.distributions.LeptonDepthFunction()
position_distribution = LI.distributions.ColumnDepthPositionDistribution(
    600, 600.0, muon_range_func, set(controller.GetDetectorModelTargets()[0])
)
primary_injection_distributions["position"] = position_distribution

# SetProcesses
controller.SetProcesses(
    primary_type, primary_injection_distributions, primary_physical_distributions
)

controller.Initialize()

events = controller.GenerateEvents()

controller.SaveEvents(
    "IceCube_DIS.hdf5"
)
