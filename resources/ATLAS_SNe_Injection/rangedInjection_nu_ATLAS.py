# Import python libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.special import erf
import pandas as pd

# Import the CCM modules
import leptoninjector as LI

# Make an instance of random for use 
random = LI.utilities.LI_random()

# Make the earth model, add ATLAS detector layout and materials file

materials_file = '../earthparams/materials/ATLAS.dat'
earth_model_file = '../earthparams/densities/PREM_ATLAS.dat'

earth_model = LI.detector.EarthModel()
earth_model.LoadMaterialModel(materials_file)
earth_model.LoadEarthModel(earth_model_file)

primary_type = LI.dataclasses.Particle.ParticleType.NuMu

primary_injection_process = LI.injection.InjectionProcess()
primary_injection_process.primary_type = primary_type

primary_physical_process = LI.injection.PhysicalProcess()
primary_physical_process.primary_type = primary_type


# DISFromSpline(std::string differential_filename, 
#               std::string total_filename, 
#               std::vector<LI::dataclasses::Particle::ParticleType> primary_types, 
#               std::vector<LI::dataclasses::Particle::ParticleType> target_types);

xsfiledir = '/n/home01/awen/prometheus/resources/cross_section_splines'
target_type = LI.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = LI.crosssections.DISFromSpline(xsfiledir+'/dsdxdy_nu_CC_iso.fits',
                                        xsfiledir+'/sigma_nu_CC_iso.fits',
                                        [primary_type],
                                        [target_type])

primary_xs = LI.crosssections.CrossSectionCollection(primary_type, [DIS_xs])
primary_injection_process.SetCrossSections(primary_xs)
primary_physical_process.SetCrossSections(primary_xs)



nu_energy = 100 #GeV
edist = LI.distributions.Monoenergetic(nu_energy)

primary_injection_process.AddInjectionDistribution(edist)
primary_physical_process.AddPhysicalDistribution(edist)

flux_units = LI.distributions.NormalizationConstant(3.76e-2)
primary_physical_process.AddPhysicalDistribution(flux_units)




ATLAS_origin = LI.math.Vector3D(0, 0, 6371234) #6371234 #90m below earth surface, ~depth of altas cavern
earth_origin = LI.math.Vector3D(0, 0, 0)

injection_dir = ATLAS_origin - earth_origin
injection_dir.normalize()

inj_ddist = LI.distributions.FixedDirection(injection_dir)
phys_ddist = LI.distributions.FixedDirection(injection_dir)

primary_injection_process.AddInjectionDistribution(inj_ddist)
primary_physical_process.AddPhysicalDistribution(phys_ddist)




target_momentum_distribution = LI.distributions.TargetAtRest()
primary_injection_process.AddInjectionDistribution(target_momentum_distribution)
primary_physical_process.AddPhysicalDistribution(target_momentum_distribution)

helicity_distribution = LI.distributions.PrimaryNeutrinoHelicityDistribution()
primary_injection_process.AddInjectionDistribution(helicity_distribution)
primary_physical_process.AddPhysicalDistribution(helicity_distribution)


# Put it all together! One injector for each W target
events_to_inject = 10000

#ATLAS_injector = LI.injection.InjectorBase(events_to_inject, earth_model, primary_injection_process, [], random)

ATLAS_ColumnDepthInjector = LI.injection.ColumnDepthLeptonInjector(events_to_inject,
                                                                  earth_model,
                                                                  primary_injection_process,
                                                                  [],
                                                                  random, 
                                                                  LI.distributions.LeptonDepthFunction(), 
                                                                  23, 23) #radius, endcap length

print('Finished building injector!')

# stopping condition for interaction tree generation
#
# this function should return true if the input datum should terminate
# the simulation for the current branch of the interaction
#
# for this test, stop after any secondary interaction tree datum is created
def StoppingCondition(datum):
    return True

ATLAS_ColumnDepthInjector.SetStoppingCondition(StoppingCondition)


ATLAS_weighter = LI.injection.LeptonTreeWeighter([ATLAS_ColumnDepthInjector], 
                                                 earth_model, 
                                                 primary_physical_process, 
                                                 []) #empty last argument since no secondary physical process

print('Built weighter!')

tree = ATLAS_ColumnDepthInjector.GenerateEvent()

weight = ATLAS_weighter.EventWeight(tree)
print(weight)

print('Reached end, after weighter!')