import numpy as np
import os
import matplotlib.pyplot as plt

import leptoninjector as LI
import sys
sys.path.insert(1,'/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/LeptonInjector/python')
from LIController import LIController
 
LI_SRC = os.environ.get('LEPTONINJECTOR_SRC')

# Define a DarkNews model 
model_kwargs = {
    'm4': 0.47,#0.140,
    'mu_tr_mu4': 2.5e-6, #1e-6, # GeV^-1
    'UD4': 0,
    'Umu4': 0,
    'epsilon': 0.0,
    'gD': 0.0,
    'decay_product':'photon',
    'noHC':True,
    'HNLtype':"dirac",
}

# Number of events to inject
events_to_inject = 1000

# Expeirment to run
experiment = 'MiniBooNE'

# Define the controller
controller = LIController(events_to_inject,
                          experiment)

# Particle to inject
primary_type = LI.dataclasses.Particle.ParticleType.NuMu

# Define DarkNews Model
table_dir = LI_SRC+'/resources/CrossSectionTables/DarkNewsTables/Dipole_M%2.2f_mu%2.2e/'%(model_kwargs['m4'],model_kwargs['mu_tr_mu4'])
controller.InputDarkNewsModel(primary_type,
                              table_dir,
                              model_kwargs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

# energy distribution
flux_file = LI_SRC + '/resources/Fluxes/BNB_Flux_Tables/BNB_numu_flux.txt'
edist = LI.distributions.TabulatedFluxDistribution(flux_file, True)
edist_gen = LI.distributions.TabulatedFluxDistribution(1.05*model_kwargs['m4'], 10, flux_file, False)
primary_injection_distributions['energy'] = edist_gen
primary_physical_distributions['energy'] = edist

# Flux normalization: go from cm^-2 to m^-2
flux_units = LI.distributions.NormalizationConstant(1e4)
primary_physical_distributions['flux_units'] = flux_units

# direction distribution
direction_distribution = LI.distributions.FixedDirection(LI.math.Vector3D(0, 0, 1.0))
primary_injection_distributions['direction'] = direction_distribution
primary_physical_distributions['direction'] = direction_distribution

# position distribution
decay_range_func = LI.distributions.DecayRangeFunction(model_kwargs['m4'],
                                                       controller.DN_min_decay_width,
                                                       3,
                                                       541)
position_distribution = LI.distributions.RangePositionDistribution(6.2, 6.2,
                                                                   decay_range_func,
                                                                   set(controller.GetEarthModelTargets()[0]))
primary_injection_distributions['position'] = position_distribution

# SetProcesses
controller.SetProcesses(primary_type,
                        primary_injection_distributions,
                        primary_physical_distributions)

controller.Initialize()

events = controller.GenerateEvents()

controller.SaveEvents('output/MiniBooNE_Dipole_M%2.2f_mu%2.2e_example.hdf5'%(model_kwargs['m4'],model_kwargs['mu_tr_mu4']))

