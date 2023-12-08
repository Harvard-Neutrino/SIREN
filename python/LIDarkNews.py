import numpy as np
import os
import datetime
import json
import ntpath
import pickle
import glob
from collections import OrderedDict

import leptoninjector as LI
from leptoninjector.crosssections import DarkNewsCrossSection,DarkNewsDecay
from leptoninjector.dataclasses import Particle
#from leptoninjector.crosssections import DarkNewsDecay
from DarkNews import phase_space
from DarkNews.ModelContainer import ModelContainer
from DarkNews.processes import *
from DarkNews.nuclear_tools import NuclearTarget

# Class containing all upscattering and decay modes available in DarkNews
class PyDarkNewsCrossSectionCollection:

    def __init__(self, table_dir=None, param_file=None, **kwargs):
        # Defines a series of upscattering and decay objects
        # Each derive from the respective LeptonInjector classes
        
        # Get our model container with all ups_case and dec_case DarkNews objects
        self.models = ModelContainer(param_file,**kwargs)
        self.table_dir = table_dir
        
        # Default table_dir settings
        if self.table_dir is None:
            self.table_dir = os.environ.get('LEPTONINJECTOR_SRC') + '/resources/CrossSectionTables/DarkNewsTables/'
            ct = datetime.datetime.now().strftime('%Y_%m_%d__%H:%M/')
            self.table_dir += ct
        
        # Make the table directory where will we store cross section integrators
        table_dir_exists = False
        if(os.path.exists(self.table_dir)):
            print("Directory '%s' already exists"%self.table_dir)
            table_dir_exists = True
        else:
            try:
                os.makedirs(self.table_dir, exist_ok = False)
                print("Directory '%s' created successfully"%self.table_dir)
            except OSError as error:
                print("Directory '%s' cannot be created"%self.table_dir)
                exit(0)

        if table_dir_exists:
            # Ensure that the model requested matches the model file already in the dictionary
            if param_file is not None:
                # ensure the param filename already exists
                param_filename = ntpath.basename(param_file) # should be OS-independent
                assert(os.path.isfile(self.table_dir + param_filename))
            # Make sure the model arguments agree
            with open(self.table_dir+'model_parameters.json',) as f:
                _model_args_dict = json.load(f)
                assert(self.models.model_args_dict==_model_args_dict)
        else:
            # Write a file to the directory containing infomration on the parameters used to create the model
            if param_file is not None:
                # Copy the param_file to the folder
                command = 'scp ' + param_file + ' ' + self.table_dir
                os.system(command)
            # Dump the model arguments
            with open(self.table_dir+'model_parameters.json','w') as f:
                json.dump(self.models.model_args_dict,f)
        
        self.cross_sections = []
        for ups_key,ups_case in self.models.ups_cases.items():
            table_subdirs = ''
            for x in ups_key:
                if(type(x)==NuclearTarget):
                    x = x.name
                table_subdirs += '%s/'%str(x)
            self.cross_sections.append(PyDarkNewsCrossSection(ups_case, table_dir = self.table_dir + table_subdirs))
        self.decays = []
        for dec_key,dec_case in self.models.dec_cases.items():
            self.decays.append(PyDarkNewsDecay(dec_case))
        


# A class representing a single ups_case DarkNews class
# Only handles methods concerning the upscattering part
class PyDarkNewsCrossSection(DarkNewsCrossSection):

    def __init__(self,
                 ups_case,
                 table_dir=None,
                 tolerance=1e-3,
                 interp_tolerance=5e-3):
        DarkNewsCrossSection.__init__(self) # C++ constructor
        self.ups_case = ups_case
        self.tolerance = tolerance
        self.interp_tolerance = interp_tolerance
        self.table_dir = table_dir

        # objects that will help with interpolation
        # will remain empty if table_dir is not set
        self.cross_section_integrator = OrderedDict() # holds vegas integrator objects
        self.cross_section_norms = OrderedDict() # for the pre-calculated norms of batch_integrand objects
        self.cross_section_energies_float = np.empty(0,dtype=float)
        self.cross_section_energies_string = np.empty(0,dtype=str)
        
        if table_dir is None: 
            print('No table_dir specified; disabling interpolation\nWARNING: this will siginficantly slow down event generation')
            return

        # Make the table directory where will we store cross section integrators
        table_dir_exists = False
        if(os.path.exists(self.table_dir)):
            print("Directory '%s' already exists"%self.table_dir)
            table_dir_exists = True
        else:
            try:
                os.makedirs(self.table_dir, exist_ok = False)
                print("Directory '%s' created successfully"%self.table_dir)
            except OSError as error:
                print("Directory '%s' cannot be created"%self.table_dir)
                exit(0)
        
        # Look in table dir and record filenames of existing integrator dumps
        if table_dir_exists:
            # First find the existing energies and sort an ascending order
            existing_files = glob.glob(self.table_dir + 'cross_section_E*.pkl')
            for file in existing_files:
                energy_str = file[file.rfind('E'):file.find('.pkl')]
                self.cross_section_integrator[energy_str] = file
                self.cross_section_energies_string = np.append(self.cross_section_energies_string,energy_str)
                self.cross_section_energies_float = np.append(self.cross_section_energies_float,float(energy_str))
                norm_file = self.table_dir + 'norm_E%s.json'%energy_str
                assert(os.path.isfile(norm_file)) # require that the norm file exists
                self.cross_section_norms[energy_str] = norm_file
            self._sort_interpolation_objects()

    def _sort_interpolation_objects(self):
        # sort in ascending energy order
        energy_srt_idxs = self.cross_section_energies_float.argsort()
        self.cross_section_energies_float = self.cross_section_energies_float[energy_srt_idxs]
        self.cross_section_energies_string = self.cross_section_energies_string[energy_srt_idxs]
        # sort the OrderedDicts
        for energy_str in self.cross_section_energies_string:
            self.cross_section_integrator.move_to_end(energy_str)
            self.cross_section_norms.move_to_end(energy_str)
    
    ##### START METHODS FOR SERIALIZATION #########
    # def get_initialized_dict(config):
    #     # do the intitialization step
    #     pddn = PyDerivedDarkNews(config)
    #     return pddn.__dict__
    #     # return the conent of __dict__ for PyDerivedDarkNews

    # @staticmethod  
    # def get_config(self):
    #     return self.config
    ##### END METHODS FOR SERIALIZATION #########

    def GetPossiblePrimaries(self):
        return [Particle.ParticleType(self.ups_case.nu_projectile.pdgid)]
    
    def GetPossibleTargetsFromPrimary(self, primary_type):
        if Particle.ParticleType(self.ups_case.nu_projectile.pdgid) == primary_type:
            return [Particle.Particle.ParticleType(self.ups_case.nuclear_target.pdgid)]
        return []

    def GetPossibleTargets(self):
        return [Particle.ParticleType(self.ups_case.nuclear_target.pdgid)]

    def GetPossibleSignatures(self):
        signature = LI.dataclasses.InteractionSignature()
        signature.primary_type = Particle.ParticleType(self.ups_case.nu_projectile.pdgid)
        signature.target_type = Particle.ParticleType(self.ups_case.nuclear_target.pdgid)
        signature.secondary_types = []
        signature.secondary_types.append(Particle.ParticleType(self.ups_case.nu_upscattered.pdgid))
        signature.secondary_types.append(Particle.ParticleType(self.ups_case.nuclear_target.pdgid))
        return [signature]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        if (Particle.ParticleType(self.ups_case.nu_projectile.pdgid) == primary_type) and (Particle.ParticleType(self.ups_case.nuclear_target.pdgid) == target_type):
            signature = LI.dataclasses.InteractionSignature()
            signature.primary_type = Particle.ParticleType(self.ups_case.nu_projectile.pdgid)
            signature.target_type = Particle.ParticleType(self.ups_case.nuclear_target.pdgid)
            secondary_types = []
            secondary_types.append(Particle.ParticleType(self.ups_case.nu_upscattered.pdgid))
            secondary_types.append(Particle.ParticleType(self.ups_case.nuclear_target.pdgid))
            signature.secondary_types = secondary_types
            return [signature]
        return []

    def DifferentialCrossSection(self, primary, target, energy, Q2):
        if primary != self.ups_case.nu_projectile:
            return 0
        interaction = LI.dataclasses.InteractionRecord()
        interaction.signature.primary_type = primary
        interaction.signature.target_type = target
        interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            return 0
        return self.ups_case.diff_xsec_Q2(energy, Q2)
    
    def SetUpscatteringMasses(self, interaction):
        interaction.primary_mass = 0
        interaction.target_mass = self.ups_case.MA
        secondary_masses = []
        secondary_masses.append(self.ups_case.m_ups)
        secondary_masses.append(self.ups_case.MA)
        interaction.secondary_masses = secondary_masses
        self.m_ups = self.ups_case.m_ups
        self.m_target = self.ups_case.MA

    def SetUpscatteringHelicities(self, interaction):
        secondary_helicities = []
        secondary_helicities.append(self.ups_case.h_upscattered * interaction.primary_helicity)
        secondary_helicities.append(interaction.target_helicity)
        interaction.secondary_helicity = secondary_helicities
        self.h_ups = self.ups_case.m_ups
        self.h_target = self.ups_case.MA
    
    def GetXsecFromFiles(self, integrator_file, norm_file):
        # given integrator file and norm file, return the total cross section
        # assumes boeth files exist
        with open(integrator_file, 'rb') as ifile:
            results, integrator = pickle.load(ifile)
        with open(norm_file,) as nfile:
            norm = json.load(nfile)
        return results["diff_xsec"].mean * norm["diff_xsec"]
    
    def TotalCrossSection(self, arg1, energy=None, target=None):
        # Handle overloaded arguments
        if type(arg1==LI.dataclasses.InteractionRecord):
            primary = arg1.signature.primary_type
            energy = arg1.primary_momentum[0]
            target = arg1.signature.target_type
        elif energy is not None and target is not None:
            primary = arg1
        else:
            print('Incorrect function call to TotalCrossSection!')
            exit(0)
        if primary != self.ups_case.nu_projectile:
            return 0
        
        # first check if we have saved close-enough integrator(s) already
        if len(self.cross_section_energies_float) > 0: 
            closest_idx = np.argmin(np.abs(self.cross_section_energies_float - energy))
            diff = self.cross_section_energies_float[closest_idx] - energy
            if np.abs(diff) < self.tolerance:
                # We are close enough to use one existing integrator
                existing_integrator = self.cross_section_integrator[self.cross_section_energies_string[closest_idx]]
                existing_norm = self.cross_section_norms[self.cross_section_energies_string[closest_idx]]
                return self.GetXsecFromFiles(existing_integrator,existing_norm)
            elif np.abs(diff)<self.interp_tolerance:
                # closest existing energy is within interpolation range
                interpolate = True # bool to tell us whether to interpolate
                if diff>0:
                    # closest existing energy is above requested energy
                    diff_above = diff
                    idx_above = closest_idx
                    # check if we are at the boundary
                    if closest_idx == 0: interpolate = False
                    idx_below = closest_idx-1
                    diff_below = energy - self.cross_section_energies_float[idx_below]
                    # check if the node below is also within the interpolation tolerance
                    if diff_below>=self.interp_tolerance: interpolate = False       
                elif diff<0 and -diff<self.interp_tolerance:
                    # closest existing energy is below requested energy
                    diff_below = -diff
                    idx_below = closest_idx
                    # check if we are at boundary
                    if closest_idx >= len(self.cross_section_energies_float): interpolate = False
                    idx_above = closest_idx+1
                    diff_above = self.cross_section_energies_float[idx_above] - energy
                    # check if the node above is also within the interpolation tolerance
                    if diff_above>=self.interp_tolerance: interpolate = False 
                if interpolate:
                    # carry out linear interpolation
                    integrator_below = self.cross_section_integrator[self.cross_section_energies_string[idx_below]]
                    norm_below = self.cross_section_norms[self.cross_section_energies_string[idx_below]]
                    xsec_below = self.GetXsecFromFiles(integrator_below,norm_below)
                    integrator_above = self.cross_section_integrator[self.cross_section_energies_string[idx_above]]
                    norm_above = self.cross_section_norms[self.cross_section_energies_string[idx_above]]
                    xsec_above = self.GetXsecFromFiles(integrator_above,norm_above)
                    return (xsec_below/diff_below + xsec_above/diff_above) / (1./diff_below + 1./diff_above)
        
        # If we have reached this block, we must compute the cross section using vegas
        print('calculating a new xsec')
        interaction = LI.dataclasses.InteractionRecord()
        interaction.signature.primary_type = primary
        interaction.signature.target_type = target
        interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            ret = 0
        
        self.cross_section_energies_float = np.append(self.cross_section_energies_float,energy)
        energy_str = '%3.3e'%energy
        self.cross_section_energies_string = np.append(self.cross_section_energies_string,energy_str)
        self.cross_section_integrator[energy_str] = self.table_dir + 'cross_section_E%s.pkl'%energy_str
        self.cross_section_norms[energy_str] = self.table_dir + 'norm_E%s.json'%energy_str
        return self.ups_case.scalar_total_xsec(energy,
                                               savefile_xsec=self.cross_section_integrator[energy_str],
                                               savefile_norm=self.cross_section_norms[energy_str])
        

    def InteractionThreshold(self, interaction):
        return 1.05 * self.ups_case.Ethreshold

    def Q2Min(self, interaction):
        return phase_space.upscattering_Q2min(interaction.primary_momentum[0], self.ups_case.m_ups, interaction.target_mass)

    def Q2Max(self, interaction):
        return phase_space.upscattering_Q2max(interaction.primary_momentum[0], self.ups_case.m_ups, interaction.target_mass)

# A class representing a single decay_case DarkNews class
# Only handles methods concerning the decay part
class PyDarkNewsDecay(DarkNewsDecay):

    def __init__(self, dec_case, **kwargs):
        DarkNewsDecay.__init__(self) # C++ constructor
        self.dec_case = dec_case 

    ##### START METHODS FOR SERIALIZATION #########
    # def get_initialized_dict(config):
    #     # do the intitialization step
    #     pddn = PyDerivedDarkNews(config)
    #     return pddn.__dict__
    #     # return the conent of __dict__ for PyDerivedDarkNews

    # @staticmethod  
    # def get_config(self):
    #     return self.config
    ##### END METHODS FOR SERIALIZATION #########

    def GetPossibleSignatures(self):
        signature = LI.dataclasses.InteractionSignature()
        signature.primary_type = Particle.ParticleType(self.dec_case.nu_parent.pdgid)
        signature.target_type = Particle.ParticleType.Decay
        secondary_types = []
        secondary_types.append(Particle.ParticleType(self.dec_case.nu_daughter.pdgid))
        for secondary in self.dec_case.secondaries:
            secondary_types.append(Particle.ParticleType(secondary.pdgid))
        signature.secondary_types = secondary_types
        return [signature]

    def GetPossibleSignaturesFromParent(self, primary_type):
        if (Particle.ParticleType(self.dec_case.nu_parent.pdgid) == primary_type):
            signature = LI.dataclasses.InteractionSignature()
            signature.primary_type = Particle.ParticleType(self.dec_case.nu_parent.pdgid)
            signature.target_type = Particle.ParticleType.Decay
            secondary_types = []
            secondary_types.append(Particle.ParticleType(self.dec_case.nu_daughter.pdgid))
            for secondary in self.dec_case.secondaries:
                secondary_types.append(Particle.ParticleType(secondary.pdgid))
            signature.secondary_types = secondary_types
            return [signature]
        return []

    def DifferentialDecayWidth(self, record):
        if type(self.dec_case)==FermionSinglePhotonDecay:
            gamma_idx = 0
            for secondary in record.signature.secondary_types:
                if secondary == LI.dataclasses.Particle.ParticleType.Gamma:
                    break
                gamma_idx += 1
            if gamma_idx >= len(record.secondary_momenta):
                print('No gamma found in the list of secondaries!')
                exit(0)

            E1,p1x,p1y,p1z = record.primary_momentum
            E2,p2x,p2y,p2z = record.secondary_momenta[gamma_idx]
            cost = p2z / E2
            return self.dec_case.differential_width(cost)
        else:
            #TODO: implement dilepton case
            return 0
    
    def TotalDecayWidth(self, arg1):
        if type(arg1)==LI.dataclasses.InteractionRecord:
            primary = arg1.signature.primary_type
        elif type(arg1)==LI.dataclasses.Particle.ParticleType:
            primary = arg1
        else:
            print('Incorrect function call to TotalDecayWidth!')
            exit(0)
        if primary != self.dec_case.nu_parent:
            return 0
        return self.dec_case.total_width()
    
    def TotalDecayWidthForFinalState(self,record):
        sig = self.GetPossibleSignatures()[0]
        if (record.signature.primary_type != sig.primary_type) or \
           (record.signature.target_type != sig.target_type) or \
           (len(record.signature.secondary_types) != len(sig.secondary_types)) or \
           (np.any([record.signature.secondary_types[i] != sig.secondary_types[i] for i in range(len(sig.secondary_types))])):
            return 0
        return self.dec_case.total_width()
    
    def DensityVariables(self):
        if type(self.dec_case)==FermionSinglePhotonDecay:
            return "cost"
        elif type(self.dec_case)==FermionDileptonDecay:
            if self.dec_case.vector_on_shell and self.dec_case.scalar_on_shell:
                print('Can\'t have both the scalar and vector on shell')
                exit(0)
            elif (self.dec_case.vector_on_shell and self.dec_case.scalar_off_shell) or \
                 (self.dec_case.vector_on_shell and self.dec_case.scalar_off_shell):
                return "cost"
            elif self.dec_case.vector_off_shell and self.dec_case.scalar_off_shell:
                return "PS"
        return ""
    
    def SampleFinalState(self,record,random):
        # TODO: implement this after talking to Matheus about how to modify DarkNews
        return