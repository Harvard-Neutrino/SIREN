import os
import sys
import numpy as np
import functools

import leptoninjector as LI
from leptoninjector import _util
from leptoninjector.LIController import LIController

@functools.wraps(LIController.GenerateEvents)
def GenerateEvents(self, N=None):
    if N is None:
        N = self.events_to_inject
    count = 0
    while (self.injector.InjectedEvents() < self.events_to_inject) and (count < N):
        print("Injecting Event", count, end="\r")
        tree = self.injector.GenerateEvent()
        self.weighter.EventWeight(tree)
        self.events.append(tree)
        count += 1
    #if hasattr(self, "DN_processes"):
    #    self.DN_processes.SaveCrossSectionTables()
    return self.events

LIController.GenerateEvents = GenerateEvents

resources_dir = _util.resource_package_dir()

# Number of events to inject
events_to_inject = int(1e6)

# Expeirment to run
experiment = "IceCube"

# Define the controller
controller = LIController(events_to_inject, experiment)

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
    [target_type],
)

primary_xs = LI.interactions.InteractionCollection(primary_type, [DIS_xs])
controller.SetInteractions(primary_xs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

mass_dist = LI.distributions.PrimaryMass(0)
primary_injection_distributions["mass"] = mass_dist
primary_physical_distributions["mass"] = mass_dist

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

controller.SaveEvents("IceCube_DIS.hdf5")
