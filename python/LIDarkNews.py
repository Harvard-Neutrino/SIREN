import numpy as np

import leptoninjector as LI
from leptoninjector.crosssections import DarkNewsCrossSection,DarkNewsDecay
from leptoninjector.dataclasses import Particle
#from leptoninjector.crosssections import DarkNewsDecay
from DarkNews import phase_space
from DarkNews.ModelContainer import ModelContainer
from DarkNews.processes import *

# Class containing all upscattering and decay modes available in DarkNews
class PyDarkNewsCrossSectionCollection:

    def __init__(self, param_file=None, **kwargs):
        # This constructor follows closely from the GenLauncher constructor from DarkNews
        # but only enough information to define upscattering/decay objects
        self.models = ModelContainer(param_file,**kwargs)
        self.cross_sections = []
        for ups_key,ups_case in self.models.ups_cases.items():
            self.cross_sections.append(PyDarkNewsCrossSection(ups_case))
        self.decays = []
        for dec_key,dec_case in self.models.dec_cases.items():
            self.decays.append(PyDarkNewsDecay(dec_case))
        


# A class representing a single ups_case DarkNews class
# Only handles methods concerning the upscattering part
class PyDarkNewsCrossSection(DarkNewsCrossSection):

    def __init__(self, ups_case, table_dir=None, **kwargs):
        DarkNewsCrossSection.__init__(self) # C++ constructor
        #CrossSection.__init__(self) # C++ constructor
        self.ups_case = ups_case
        # self.table_dir = table_dir
        # if table_dir is None:
        #     self.table_dir = os.environ.get('LEPTONINJECTOR_SRC') + '/resources/CrossSectionTables/DarkNewsTables/'

        # Define the vegas integration object to be adapted
        # DIM = 3 # TODO: check this for upscattering only
        # self.integrand = integrands.UpscatteringXsec(DIM,
        #                                             )


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
    
    def TotalCrossSection(self, arg1, energy=None, target=None):
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
        interaction = LI.dataclasses.InteractionRecord()
        interaction.signature.primary_type = primary
        interaction.signature.target_type = target
        interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            ret = 0
        ret = self.ups_case.total_xsec(energy)
        return ret

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