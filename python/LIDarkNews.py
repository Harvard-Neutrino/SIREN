import leptoninjector as LI
from leptoninjector.crosssections import DarkNewsCrossSection as DNCS
#from leptoninjector.crosssections import DarkNewsDecay as DND
import DarkNews
from DarkNews import integrands

# Class containing all upscattering and decay modes available in DarkNews
class PyDarkNewsCrossSectionCollection:

    def __init__(self, param_file=None, **kwargs):
        self.gen = DarkNews.GenLauncher(param_file,**kwargs)
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
class PyDarkNewsCrossSection(DNCS):

    def __init__(self, gen_case, table_dir=None, **kwargs):
        DNCS.__init__(self) # C++ constructor
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
        return [self.gen_case.nu_projectile]
        # primaries = []
        # for gen_case in self.gen.gen_cases:
        #     primaries.append(gen_case.nu_projectile)
        # return primaries
    
    def GetPossibleTargetsFromPrimary(self, primary_type):
        if self.gen_case.nu_projectile == primary_type:
            return [self.gen_case.nuclear_target]
        return []
        # targets = []
        # for gen_case in self.gen.gen_cases:
        #     if gen_case.nu_projectile == primary_type:
        #         targets.append(gen_case.nuclear_target)
        # return targets

    def GetPossibleTargets(self):
        return [self.gen_case.nuclear_target]
        # targets = []
        # for gen_case in self.gen.gen_cases:
        #     targets.append(gen_case.nuclear_target)
        # return targets

    def GetPossibleSignatures(self):
        signature = LI.dataclasses.InteractionSignature
        signature.primary_type = self.gen_case.nu_projectile
        signature.target_type = self.gen_case.nuclear_target
        signature.secondary_types = []
        signature.secondary_types.append(self.gen_case.nu_upscattered)
        signature.secondary_types.append(self.gen_case.nuclear_target)
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
        if (self.gen_case.nu_projectile == primary_type) and (self.gen_case.nuclear_target == target_type):
            signature = LI.dataclasses.InteractionSignature
            signature.primary_type = self.gen_case.nu_projectile
            signature.target_type = self.gen_case.nuclear_target
            signature.secondary_types = []
            signature.secondary_types.append(self.gen_case.nu_upscattered)
            signature.secondary_types.append(self.gen_case.nuclear_target)
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
        #interaction.signature.primary_type = primary
        #interaction.signature.target_type = target
        #interaction.primary_momentum[0] = energy
        if energy < self.InteractionThreshold(interaction):
            return 0
        return self.gen_case.ups_case.diff_xsec_Q2(energy, Q2)



    def InteractionThreshold(self, interaction):
        return 1.05 * self.gen_case.ups_case.Ethreshold

    def Q2Min(self, interaction):
        return 0

    def Q2Max(self, interaction):
        return 0