import os
import numpy as np

import leptoninjector as LI
import sys
sys.path.insert(1,'/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/LeptonInjector/python')
from LIController import LIController
 
LI_SRC = os.environ.get('LEPTONINJECTOR_SRC')

# Define a DarkNews model 
model_kwargs = {
    'm4': 0.02,
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
experiment = 'CCM'

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
nu_energy = 0.02965 # from pi+ DAR
edist = LI.distributions.Monoenergetic(nu_energy)
primary_injection_distributions['energy'] = edist
primary_physical_distributions['energy'] = edist

# Flux normalization: 
# using the number quoted in 2105.14020, 4.74e9 nu/m^2/s / (6.2e14 POT/s) * 4*pi*20m^2 to get nu/POT
flux_units = LI.distributions.NormalizationConstant(3.76e-2)
primary_physical_distributions['flux_units'] = flux_units

# direction distribution: cone from lower W target
opening_angle = np.arctan(12/23.); # slightly larger than CCM
lower_target_origin = LI.math.Vector3D(0, 0, -0.241)
detector_origin = LI.math.Vector3D(23, 0, -0.65)
lower_dir = detector_origin - lower_target_origin
lower_dir.normalize()
lower_inj_ddist = LI.distributions.Cone(lower_dir,opening_angle)
phys_ddist = LI.distributions.IsotropicDirection() # truly we are isotropicprimary_injection_distributions['direction'] = direction_distribution
primary_injection_distributions['direction'] = lower_inj_ddist
primary_physical_distributions['direction'] = phys_ddist

# Position distribution: consider neutrinos from a point source
max_dist = 25
lower_pos_dist = LI.distributions.PointSourcePositionDistribution(lower_target_origin, max_dist, set(controller.GetDetectorModelTargets()[0]))
primary_injection_distributions['position'] = lower_pos_dist

# SetProcesses
controller.SetProcesses(primary_type,
                        primary_injection_distributions,
                        primary_physical_distributions)

controller.Initialize()

events = controller.GenerateEvents()

controller.SaveEvents('output/CCM_Dipole_M%2.2f_mu%2.2e_example.hdf5'%(model_kwargs['m4'],model_kwargs['mu_tr_mu4']))