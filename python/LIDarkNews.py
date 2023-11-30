import leptoninjector as LI
from leptoninjector.crosssections import DarkNewsCrossSection
from leptoninjector.dataclasses import Particle
#from leptoninjector.crosssections import DarkNewsDecay
from DarkNews import GenLauncher,phase_space

# Class containing all upscattering and decay modes available in DarkNews
class PyDarkNewsCrossSectionCollection:

    def __init__(self, param_file=None, **kwargs):
        self.gen = GenLauncher(param_file,**kwargs)
        self.cross_sections = {}
        self.decays = {}
        for gen_case in self.gen.gen_cases:
            cross_section_key = (gen_case.nu_projectile, # initial state neutrino
                                 gen_case.nuclear_target, # nuclear target
                                 gen_case.nu_upscattered, # final state neutrino
                                 gen_case.helicity, # helicity conserving or flipping
                                 gen_case.scope["scattering_regime"]) # coherent or p/n-elastic
            if cross_section_key not in self.cross_sections.keys():
                self.cross_sections[cross_section_key] = PyDarkNewsCrossSection(gen_case)
            # TODO: add decays


# A class representing a single MC_events DarkNews class
# Only handles methods concerning the upscattering part
class PyDarkNewsCrossSection(DarkNewsCrossSection):

    def __init__(self, gen_case, table_dir=None, **kwargs):
        DarkNewsCrossSection.__init__(self) # C++ constructor
        #CrossSection.__init__(self) # C++ constructor
        self.gen_case = gen_case
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
        return [Particle.ParticleType(self.gen_case.nu_projectile.pdgid)]
        # primaries = []
        # for gen_case in self.gen.gen_cases:
        #     primaries.append(gen_case.nu_projectile)
        # return primaries
    
    def GetPossibleTargetsFromPrimary(self, primary_type):
        if Particle.ParticleType(self.gen_case.nu_projectile.pdgid) == primary_type:
            return [Particle.Particle.ParticleType(self.gen_case.nuclear_target.pdgid)]
        return []
        # targets = []
        # for gen_case in self.gen.gen_cases:
        #     if gen_case.nu_projectile == primary_type:
        #         targets.append(gen_case.nuclear_target)
        # return targets

    def GetPossibleTargets(self):
        return [Particle.ParticleType(self.gen_case.nuclear_target.pdgid)]
        # targets = []
        # for gen_case in self.gen.gen_cases:
        #     targets.append(gen_case.nuclear_target)
        # return targets

    def GetPossibleSignatures(self):
        signature = LI.dataclasses.InteractionSignature()
        signature.primary_type = Particle.ParticleType(self.gen_case.nu_projectile.pdgid)
        signature.target_type = Particle.ParticleType(self.gen_case.nuclear_target.pdgid)
        signature.secondary_types = []
        signature.secondary_types.append(Particle.ParticleType(self.gen_case.nu_upscattered.pdgid))
        signature.secondary_types.append(Particle.ParticleType(self.gen_case.nuclear_target.pdgid))
        return [signature]
        # signatures = []
        # for gen_case in self.gen.gen_cases:
        #     signature = LI.dataclasses.InteractionSignature
        #     signature.primary_type = gen_case.nu_projectile
        #     signature.target_type = gen_case.nuclear_target
        #     signature.secondary_types = []
        #     signature.secondary_types.append(gen_case.nu_upscattered)
        #     signature.secondary_types.append(gen_case.nuclear_target)
        #     signatures.append(signature)
        # return signatures

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        if (Particle.ParticleType(self.gen_case.nu_projectile.pdgid) == primary_type) and (Particle.ParticleType(self.gen_case.nuclear_target.pdgid) == target_type):
            signature = LI.dataclasses.InteractionSignature()
            signature.primary_type = Particle.ParticleType(self.gen_case.nu_projectile.pdgid)
            signature.target_type = Particle.ParticleType(self.gen_case.nuclear_target.pdgid)
            secondary_types = []
            secondary_types.append(Particle.ParticleType(self.gen_case.nu_upscattered.pdgid))
            secondary_types.append(Particle.ParticleType(self.gen_case.nuclear_target.pdgid))
            signature.secondary_types = secondary_types
            return [signature]
        return []
        # signatures = []
        # for gen_case in self.gen.gen_cases:
        #     if (gen_case.nu_projectile == primary_type) and (gen_case.nuclear_target == target_type):
        #         signature = LI.dataclasses.InteractionSignature
        #         signature.primary_type = gen_case.nu_projectile
        #         signature.target_type = gen_case.nuclear_target
        #         signature.secondary_types = []
        #         signature.secondary_types.append(gen_case.nu_upscattered)
        #         signature.secondary_types.append(gen_case.nuclear_target)
        #         signatures.append(signature)
        # return signatures

    def DifferentialCrossSection(self, primary, target, energy, Q2):
        if primary != self.gen_case.nu_projectile:
            return 0
        interaction = LI.dataclasses.InteractionRecord()
        interaction.signature.primary_type = primary
        interaction.signature.target_type = target
        interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            return 0
        return self.gen_case.ups_case.diff_xsec_Q2(energy, Q2)
    
    def SetUpscatteringMasses(self, interaction):
        interaction.primary_mass = 0
        interaction.target_mass = self.gen_case.ups_case.MA
        secondary_masses = []
        secondary_masses.append(self.gen_case.ups_case.m_ups)
        secondary_masses.append(self.gen_case.ups_case.MA)
        interaction.secondary_masses = secondary_masses
        self.m_ups = self.gen_case.ups_case.m_ups
        self.m_target = self.gen_case.ups_case.MA
    
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
        if primary != self.gen_case.nu_projectile:
            return 0
        interaction = LI.dataclasses.InteractionRecord()
        interaction.signature.primary_type = primary
        interaction.signature.target_type = target
        interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            ret = 0
        ret = self.gen_case.ups_case.total_xsec(energy)
        print(ret)
        return ret

    def InteractionThreshold(self, interaction):
        return 1.05 * self.gen_case.ups_case.Ethreshold

    def Q2Min(self, interaction):
        return phase_space.upscattering_Q2min(interaction.primary_momentum[0], self.gen_case.ups_case.m_ups, interaction.target_mass)

    def Q2Max(self, interaction):
        return phase_space.upscattering_Q2max(interaction.primary_momentum[0], self.gen_case.ups_case.m_ups, interaction.target_mass)
