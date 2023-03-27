# Import python libraries
import numpy as np

# Import the CCM modules
import Math
import Utilities
import Dataclasses
import Geometry
import Detector
import CrossSections
import Distributions
import Injection

# Make an instance of random for use
random = Utilities.LI_random()

# HNL parameters
hnl_mass_str = "0.02035" # make sure this matches one of the cross section tables
hnl_mass = float(hnl_mass_str) # GeV
d_dipole = 3e-7 # GeV^-1

# Injection parameters
events_to_inject = int(1e4)
primary_type = Dataclasses.Particle.ParticleType.NuMu

# For upscattering cross section
z_samp = False
inv_GeV = True
inelastic = True
target_types = [Dataclasses.Particle.ParticleType.HNucleus,
                Dataclasses.Particle.ParticleType.Be9Nucleus,
                Dataclasses.Particle.ParticleType.C12Nucleus,
                Dataclasses.Particle.ParticleType.N14Nucleus,
                Dataclasses.Particle.ParticleType.O16Nucleus,
                Dataclasses.Particle.ParticleType.Na23Nucleus,
                Dataclasses.Particle.ParticleType.Al27Nucleus,
                Dataclasses.Particle.ParticleType.Si28Nucleus,
                Dataclasses.Particle.ParticleType.Ar40Nucleus,
                Dataclasses.Particle.ParticleType.Ca40Nucleus,
                Dataclasses.Particle.ParticleType.Mn55Nucleus,
                Dataclasses.Particle.ParticleType.Fe56Nucleus,
                Dataclasses.Particle.ParticleType.Cu63Nucleus,
                Dataclasses.Particle.ParticleType.Cu65Nucleus,
                Dataclasses.Particle.ParticleType.W183Nucleus,
                Dataclasses.Particle.ParticleType.Pb208Nucleus]
ZA_combos = [[1,1],
             [4,9],
             [6,12],
             [7,14],
             [8,16],
             [11,23],
             [13,27],
             [14,28],
             [18,40],
             [20,40],
             [25,55],
             [26,56],
             [29,63],
             [29,65],
             [74,183],
             [82,208]]
# file locations
tot_xsec_path = '../../../../Sandbox/Dipole_xsec_tables/tot_xsec_E0_Ne20'
diff_xsec_path = '../../../../Sandbox/Dipole_xsec_tables/diff_xsec_y_Enu_0.055_Ne20'
tot_hf_files = ['%s/xsec_Z_%i_A_%i_mHNL_%s_hf.dat'%(tot_xsec_path,
                                                    Z,A,
                                                    hnl_mass_str) for Z,A in ZA_combos]
tot_hc_files = ['%s/xsec_Z_%i_A_%i_mHNL_%s_hc.dat'%(tot_xsec_path,
                                                    Z,A,
                                                    hnl_mass_str) for Z,A in ZA_combos]
diff_hf_files = ['%s/dxsec_Z_%i_A_%i_mHNL_%s_hf.dat'%(diff_xsec_path,
                                                      Z,A,
                                                      hnl_mass_str) for Z,A in ZA_combos]
diff_hc_files = ['%s/dxsec_Z_%i_A_%i_mHNL_%s_hc.dat'%(diff_xsec_path,
                                                      Z,A,
                                                      hnl_mass_str) for Z,A in ZA_combos]

# Make the earth model, add CCM detector layout and materials file

materials_file = '../earthparams/materials/CCM.dat'
earth_model_file = '../earthparams/densities/PREM_ccm.dat'

earth_model = Detector.EarthModel()
earth_model.LoadMaterialModel(materials_file)
earth_model.LoadEarthModel(earth_model_file)

# Define injection processes for each target
primary_injection_process_upper_target = Injection.InjectionProcess()
primary_injection_process_lower_target = Injection.InjectionProcess()
secondary_injection_processes = []
primary_injection_process_upper_target.primary_type = primary_type
primary_injection_process_lower_target.primary_type = primary_type

# Define physical processes for each target
primary_physical_process_upper_target = Injection.PhysicalProcess()
primary_physical_process_lower_target = Injection.PhysicalProcess()
secondary_physical_processes = []
primary_physical_process_upper_target.primary_type = primary_type
primary_physical_process_lower_target.primary_type = primary_type

# Define upscattering cross section classes

cross_sections = []
hf_xs = CrossSections.DipoleFromTable(hnl_mass,
                                      d_dipole,
                                      CrossSections.DipoleFromTable.HelicityChannel.Flipping,
                                      z_samp,inv_GeV,inelastic)
hc_xs = CrossSections.DipoleFromTable(hnl_mass,
                                      d_dipole,
                                      CrossSections.DipoleFromTable.HelicityChannel.Conserving,
                                      z_samp,inv_GeV,inelastic)
for i in range(len(target_types)):
  hf_xs.AddTotalCrossSectionFile(tot_hf_files[i],target_types[i])
  hf_xs.AddDifferentialCrossSectionFile(diff_hf_files[i],target_types[i])
  hc_xs.AddTotalCrossSectionFile(tot_hc_files[i],target_types[i])
  hc_xs.AddDifferentialCrossSectionFile(diff_hc_files[i],target_types[i])

cross_sections.append(hf_xs)
cross_sections.append(hc_xs)
primary_cross_sections = CrossSections.CrossSectionCollection(primary_type, cross_sections)
primary_injection_process_upper_target.SetCrossSections(primary_cross_sections)
primary_injection_process_lower_target.SetCrossSections(primary_cross_sections)
primary_physical_process_upper_target.SetCrossSections(primary_cross_sections)
primary_physical_process_lower_target.SetCrossSections(primary_cross_sections)

# Energy distribution: monoenergetic neutrino from pion decay at rest

nu_energy = 0.02965 # from pi+ DAR
edist = Distributions.Monoenergetic(nu_energy)
primary_injection_process_upper_target.AddInjectionDistribution(edist)
primary_injection_process_lower_target.AddInjectionDistribution(edist)
primary_physical_process_upper_target.AddPhysicalDistribution(edist)
primary_physical_process_lower_target.AddPhysicalDistribution(edist)

# Flux normalization: 
# using the number quoted in 2105.14020, 4.74e9 nu/m^2/s / (6.2e14 POT/s) * 4*pi*20m^2 to get nu/POT
flux_units = Distributions.NormalizationConstant(3.76e-2)
primary_physical_process_upper_target.AddPhysicalDistribution(flux_units)
primary_physical_process_lower_target.AddPhysicalDistribution(flux_units)

# Primary direction: A cone around CCM

opening_angle = np.arctan(12/23.); # slightly larger than CCM
upper_target_origin = Math.Vector3D(0, 0, 0.1375)
lower_target_origin = Math.Vector3D(0, 0, -0.241)
detector_origin = Math.Vector3D(23, 0, -0.65)
upper_dir = detector_origin - upper_target_origin
upper_dir.normalize()
lower_dir = detector_origin - lower_target_origin
lower_dir.normalize()
upper_inj_ddist = Distributions.Cone(upper_dir,opening_angle)
lower_inj_ddist = Distributions.Cone(lower_dir,opening_angle)
phys_ddist = Distributions.IsotropicDirection() # truly we are isotropic
primary_injection_process_upper_target.AddInjectionDistribution(upper_inj_ddist)
primary_injection_process_lower_target.AddInjectionDistribution(lower_inj_ddist)
primary_physical_process_upper_target.AddPhysicalDistribution(phys_ddist)
primary_physical_process_lower_target.AddPhysicalDistribution(phys_ddist);

# Target momentum distribution: target is at rest

target_momentum_distribution = Distributions.TargetAtRest()
primary_injection_process_upper_target.AddInjectionDistribution(target_momentum_distribution)
primary_injection_process_lower_target.AddInjectionDistribution(target_momentum_distribution)
primary_physical_process_upper_target.AddPhysicalDistribution(target_momentum_distribution)
primary_physical_process_lower_target.AddPhysicalDistribution(target_momentum_distribution)

# Helicity distribution: primary neutrino helicity

helicity_distribution = Distributions.PrimaryNeutrinoHelicityDistribution()
primary_injection_process_upper_target.AddInjectionDistribution(helicity_distribution)
primary_injection_process_lower_target.AddInjectionDistribution(helicity_distribution)
primary_physical_process_upper_target.AddPhysicalDistribution(helicity_distribution)
primary_physical_process_lower_target.AddPhysicalDistribution(helicity_distribution)

# Position distribution: consider neutrinos from a point source

max_dist = 25
upper_pos_dist = Distributions.PointSourcePositionDistribution(upper_target_origin, max_dist, primary_cross_sections.TargetTypes())
lower_pos_dist = Distributions.PointSourcePositionDistribution(lower_target_origin, max_dist, primary_cross_sections.TargetTypes())
primary_injection_process_upper_target.AddInjectionDistribution(upper_pos_dist)
primary_injection_process_lower_target.AddInjectionDistribution(lower_pos_dist)

# Secondary process definition

secondary_decay_injection_process = Injection.InjectionProcess()
secondary_decay_physical_process = Injection.PhysicalProcess()
secondary_decay_injection_process.primary_type = Dataclasses.Particle.ParticleType.NuF4
secondary_decay_physical_process.primary_type = Dataclasses.Particle.ParticleType.NuF4

# Secondary cross section class: HNL decay

sec_decay = CrossSections.NeutrissimoDecay(hnl_mass, d_dipole, CrossSections.NeutrissimoDecay.ChiralNature.Majorana)
secondary_cross_sections = CrossSections.CrossSectionCollection(Dataclasses.Particle.ParticleType.NuF4, [sec_decay])
secondary_decay_injection_process.SetCrossSections(secondary_cross_sections)
secondary_decay_physical_process.SetCrossSections(secondary_cross_sections)

# Secondary position distribution

secondary_pos_dist = Distributions.SecondaryPositionDistribution()
for sector in earth_model.GetSectors():
  if sector.name=='ccm_inner_argon':
    fid_vol = sector.geo
    secondary_pos_dist = Distributions.SecondaryPositionDistribution(sector.geo)

secondary_decay_injection_process.AddInjectionDistribution(secondary_pos_dist)

secondary_injection_processes.append(secondary_decay_injection_process)
secondary_physical_processes.append(secondary_decay_physical_process)

# Put it all together!
upper_injector = Injection.InjectorBase(events_to_inject, earth_model,
                                        primary_injection_process_upper_target,
                                        secondary_injection_processes,
                                        random)

lower_injector = Injection.InjectorBase(events_to_inject, earth_model,
                                        primary_injection_process_lower_target,
                                        secondary_injection_processes,
                                        random)

def StoppingCondition(datum):
  return True

upper_injector.SetStoppingCondition(StoppingCondition)
lower_injector.SetStoppingCondition(StoppingCondition)

# Weighter instances

upper_weighter = Injection.LeptonTreeWeighter([upper_injector],
                                              earth_model,
                                              primary_physical_process_upper_target,
                                              secondary_physical_processes)
lower_weighter = Injection.LeptonTreeWeighter([lower_injector],
                                              earth_model,
                                              primary_physical_process_upper_target,
                                              secondary_physical_processes)


gamma_energy_list = []
gamma_weight_list = []
while lower_injector.InjectedEvents() < events_to_inject:
  print(lower_injector.InjectedEvents(),end='\r')
  tree = lower_injector.GenerateEvent()
  weight = lower_weighter.EventWeight(tree)
  for datum in tree.tree:
    if(datum.record.signature.primary_type == Dataclasses.Particle.ParticleType.NuF4):
      HNL_vtx = Math.Vector3D(datum.record.interaction_vertex)
      HNL_dir = Math.Vector3D(datum.record.primary_momentum[1:])
      HNL_dir.normalize()
      if(fid_vol.IsInside(HNL_vtx,HNL_dir)):
#         print('decay length:',secondary_cross_sections.TotalDecayLength(datum.record))
#         print('upscattering vertex:',datum.parent.record.interaction_vertex)
#         print('decay vertex:',datum.record.interaction_vertex)
#         print('weight:',weight)
        for sID,sP in zip(datum.record.signature.secondary_types,
                          datum.record.secondary_momenta):
          if(sID==Dataclasses.Particle.ParticleType.Gamma):
            gamma_energy_list.append(sP[0])
            gamma_weight_list.append(weight)
#     if(datum.record.signature.primary_type == Dataclasses.Particle.ParticleType.NuMu):
#       print(primary_cross_sections.)

import matplotlib.pyplot as plt
from scipy.stats import poisson

POT = 2.25e22
bkg = 30 # optimistic: 10 evts/year

bins = np.linspace(1e-3,3e-2,20)

weights = POT*np.array(gamma_weight_list)
sig = sum(weights)
poisson_prob = poisson.cdf(sig+bkg,bkg)
print(sig / np.sqrt(bkg))

plt.hist(gamma_energy_list,weights=weights,bins=bins,label='%2.2f total events, %2.2e C.L.'%(sig,poisson_prob))
plt.xlabel('Gamma Energy [GeV]')
plt.ylabel('Number of Events in CCM')
plt.title('HNL Mass: %s GeV; Coupling: %2.2e GeV^-1'%(hnl_mass_str,d_dipole))
plt.legend()
plt.show()

