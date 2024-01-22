import numpy as np
import os

import leptoninjector as LI
import sys
sys.path.insert(1,'/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/LeptonInjector/python')
from LIController import LIController
 
LI_SRC = os.environ.get('LEPTONINJECTOR_SRC')

# Number of events to inject
events_to_inject = 1000

# Expeirment to run
experiment = 'IceCube'

# Define the controller
controller = LIController(events_to_inject,
                          experiment)

# Particle to inject
primary_type = LI.dataclasses.Particle.ParticleType.NuMu

# Cross Section Model
xsfiledir = LI_SRC+'resources/CrossSectionTables/DISSplines/'
target_type = LI.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = LI.interactions.DISFromSpline(xsfiledir+'test_xs.fits',
                                        xsfiledir+'test_xs_total.fits',
                                        [primary_type],
                                        [target_type])

primary_xs = LI.interactions.CrossSectionCollection(primary_type, [DIS_xs])
controller.SetCrossSections(primary_xs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
edist = LI.distributions.PowerLaw(2,1e3,1e6)
primary_injection_distributions['energy'] = edist
primary_physical_distributions['energy'] = edist

# direction distribution
direction_distribution = LI.distributions.IsotropicDirection()
primary_injection_distributions['direction'] = direction_distribution
primary_physical_distributions['direction'] = direction_distribution

# position distribution
muon_range_func = LI.distributions.LeptonDepthFunction()
position_distribution = LI.distributions.ColumnDepthPositionDistribution(600, 600.,
                                                                         muon_range_func,
                                                                         set(controller.GetEarthModelTargets()[0]))
primary_injection_distributions['position'] = position_distribution

# SetProcesses
controller.SetProcesses(primary_type,
                        primary_injection_distributions,
                        primary_physical_distributions)

controller.Initialize()

events = controller.GenerateEvents()

controller.SaveEvents('output/MINERvA_Dipole_M%2.2f_mu%2.2e_example.hdf5'%(model_kwargs['m4'],model_kwargs['mu_tr_mu4']))