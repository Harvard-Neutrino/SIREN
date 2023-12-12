import os

import leptoninjector as LI
from LIDarkNews import PyDarkNewsCrossSectionCollection

# Parent python class for handling event generation
class LI:

    def __init__(self,
                 events_to_inject,
                 experiment,
                 primary_type,
                 primary_injection_distributions,
                 primary_physical_distributions,
                 primary_cross_sections,
                 secondary_types,
                 secondary_injection_distributions,
                 secondary_physical_distributions,
                 secondary_cross_sections
                 seed=0):
        """
        LI class constructor.
        :param int event_to_inject: number of events to generate
        :param str experiment: experiment name in string
        :param ParticleType primary_type: The primary particle being generated
        :param list<InjectionDistribution> primary_injection_distributions: The list of injection distributions for the primary process
        :param list<PhysicalDistribution> primary_physical_distributions: The list of physical distributions for the primary process
        :param list<ParticleType> secondary_types: The secondary particles being generated
        :param list<list<InjectionDistribution> secondary_injection_distributions: List of list of injection distributions for each secondary process
        :param list<list<PhysicalDistribution> secondary_physical_distributions: List of list of physical distributions for each secondary process
        :param int seed: Optional random number generator seed
        """

        # Initialize a random number generator
        self.random = LI.utilities.LI_random(seed=seed)

        # Find the density and materials files
        LI_SRC = os.environ.get('LEPTONINJECTOR_SRC')
        materials_file = LI_SRC + '/earthparams/materials/%s.dat'%experiment
        if experiment in ['ATLAS','dune']:
            earth_model_file = LI_SRC + '/earthparams/densities/PREM_%s.dat'%experiment
        else:
            earth_model_file = LI_SRC + '/earthparams/densities/%s.dat'%experiment
        
        self.earth_model = LI.detector.EarthModel()
        self.earth_model.LoadMaterialModel(materials_file)
        self.earth_model.LoadEarthModel(earth_model_file)
        
        # Keep track of available targets in Materials model
        count = 0
        targets = []
        target_strs = []
        while self.earth_modelearth_model.Materials.HasMaterial(count):
            for _target in self.earth_model.Materials.GetMaterialTargets(count):
                targets.append(_target)
                if str(_target).find('Nucleus') == -1: 
                    target_strs.append(str(_target))
                else:
                    target_str = str(_target)[str(_target).find('Type')+5:str(_target).find('Nucleus')]
                    if target_str == 'H': target_str = 'H1'
                    target_strs.append(target_str)
            count += 1

        # Define the primary injection and physical process
        primary_injection_process = LI.injection.InjectionProcess()
        primary_physical_process = LI.injection.PhysicalProcess()
        primary_injection_process.primary_type = primary_type
        primary_physical_process.primary_type = primary_type

        # Add all injection distributions
        for idist in primary_injection_distributions:
            primary_injection_process.AddInjectionDistribution(idist)
        # Add all physical distributions
        for pdist in primary_physical_distributions:
            primary_physical_process.AddPhysicalDistribution(pdist)

        # Define lists for the secondary injection and physical processes
        secondary_injection_processes = []
        secondary_physical_processes = []

        # Loop through possible secondary interactions
        for i_sec,secondary_type in enumerate(secondary_types):
            secondary_injection_process = LI.injection.InjectionProcess()
            secondary_physical_process = LI.injection.PhysicalProcess()
            secondary_injection_process.primary_type = secondary_type
            secondary_physical_process.primary_type = secondary_type

             # Add all injection distributions
            for idist in secondary_injection_distributions[i_sec]:
                secondary_injection_process.AddInjectionDistribution(idist)
            # Add all physical distributions
            for pdist in secondary_physical_distributions[i_sec]:
                secondary_physical_process.AddPhysicalDistribution(pdist)
            
            secondary_injection_processes.append(secondary_injection_process)
            secondary_physical_processes.append(secondary_physical_process)

        # Define stopping condition
        # TODO: make this more general
        def StoppingCondition(datum):
            return True
        
        # Define the injector object
        self.injector = LI.injection.InjectorBase(events_to_inject,
                                                  self.earth_model, 
                                                  primary_injection_process, 
                                                  secondary_injection_processes,
                                                  self.random)

        self.weighter = LI.injection.LeptonTreeWeighter([self.injector],
                                                        self.earth_model, 
                                                        primary_physical_process, 
                                                        secondary_physical_processes)
