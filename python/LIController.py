import os
import h5py
import numpy as np

import leptoninjector as LI
from LIDarkNews import PyDarkNewsInteractionCollection

# For determining fiducial volume of different experiments
fid_vol_dict = {'MiniBooNE':'fid_vol',
                'CCM':'ccm_inner_argon',
                'MINERvA':'fid_vol'}

# Parent python class for handling event generation
class LIController:

    def __init__(self,
                 events_to_inject,
                 experiment,
                 seed=0):
        """
        LI class constructor.
        :param int event_to_inject: number of events to generate
        :param str experiment: experiment name in string
        :param int seed: Optional random number generator seed
        """

        # Initialize a random number generator
        self.random = LI.utilities.LI_random(seed)
        
        # Save number of events to inject
        self.events_to_inject = events_to_inject
        self.experiment = experiment

        # Empty list for our interaction trees
        self.events = []

        # Find the density and materials files
        self.LI_SRC = os.environ.get('LEPTONINJECTOR_SRC')
        materials_file = self.LI_SRC + '/resources/DetectorParams/materials/%s.dat'%experiment
        if experiment in ['ATLAS','dune']:
            detector_model_file = self.LI_SRC + '/resources/DetectorParams/densities/PREM_%s.dat'%experiment
        else:
            detector_model_file = self.LI_SRC + '/resources/DetectorParams/densities/%s.dat'%experiment
        
        self.detector_model = LI.detector.DetectorModel()
        self.detector_model.LoadMaterialModel(materials_file)
        self.detector_model.LoadDetectorModel(detector_model_file)

        # Define the primary injection and physical process
        self.primary_injection_process = LI.injection.InjectionProcess()
        self.primary_physical_process = LI.injection.PhysicalProcess()
        
        # Define lists for the secondary injection and physical processes
        self.secondary_injection_processes = []
        self.secondary_physical_processes = []

    def SetProcesses(self,
                     primary_type,
                     primary_injection_distributions,
                     primary_physical_distributions,
                     secondary_types=[],
                     secondary_injection_distributions=[],
                     secondary_physical_distributions=[]):
        """
        LI process setter.
        :param ParticleType primary_type: The primary particle being generated
        :param dict<str,InjectionDistribution> primary_injection_distributions: The dict of injection distributions for the primary process
        :param dict<str,PhysicalDistribution> primary_physical_distributions: The dict of physical distributions for the primary process
        :param list<ParticleType> secondary_types: The secondary particles being generated
        :param list<dict<str,InjectionDistribution> secondary_injection_distributions: List of dict of injection distributions for each secondary process
        :param list<dict<str,PhysicalDistribution> secondary_physical_distributions: List of dict of physical distributions for each secondary process
        """

        # Define the primary injection and physical process
        self.primary_injection_process.primary_type = primary_type
        self.primary_physical_process.primary_type = primary_type

        # Add all injection distributions
        for _,idist in primary_injection_distributions.items():
            self.primary_injection_process.AddInjectionDistribution(idist)
        # Add all physical distributions
        for _,pdist in primary_physical_distributions.items():
            self.primary_physical_process.AddPhysicalDistribution(pdist)
        
        # Default injection distributions
        if 'target' not in primary_injection_distributions.keys():
            self.primary_injection_process.AddInjectionDistribution(LI.distributions.TargetAtRest())
        if 'helicity' not in primary_injection_distributions.keys():
            self.primary_injection_process.AddInjectionDistribution(LI.distributions.PrimaryNeutrinoHelicityDistribution())
    
        # Default injection distributions
        if 'target' not in primary_physical_distributions.keys():
            self.primary_physical_process.AddPhysicalDistribution(LI.distributions.TargetAtRest())
        if 'helicity' not in primary_physical_distributions.keys():
            self.primary_physical_process.AddPhysicalDistribution(LI.distributions.PrimaryNeutrinoHelicityDistribution())

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
            
            # Add the position distribution
            fid_vol = self.GetFiducialVolume()
            if fid_vol is not None:
                secondary_injection_process.AddInjectionDistribution(LI.distributions.SecondaryPositionDistribution(fid_vol))
            else:
                secondary_injection_process.AddInjectionDistribution(LI.distributions.SecondaryPositionDistribution())
            
            self.secondary_injection_processes.append(secondary_injection_process)
            self.secondary_physical_processes.append(secondary_physical_process)

    def InputDarkNewsModel(self,
                           primary_type,
                           table_dir,
                           model_kwargs):
        """
        Sets up the relevant processes and cross section/decay objects related to a provided DarkNews model dictionary.
        Will handle the primary cross section collection as well as the entire list of secondary processes

        :param LI.dataclasses.Particle.ParticleType primary_type: primary particle to be generated
        :param string table_dir: Directory for storing cross section and decay tables
        :param dict<str,val> model_kwargs: The dict of DarkNews model parameters
        """
        # Add nuclear targets to the model arguments
        model_kwargs['nuclear_targets'] = self.GetDetectorModelTargets()[1]
        # Initialize DarkNews cross sections and decays
        self.DN_processes = PyDarkNewsInteractionCollection(table_dir=table_dir,
                                                             **model_kwargs)
        
        # Initialize primary InteractionCollection
        # Loop over available cross sections and save those which match primary type
        primary_cross_sections = []
        for cross_section in self.DN_processes.cross_sections:
            if primary_type == LI.dataclasses.Particle.ParticleType(cross_section.ups_case.nu_projectile.pdgid):
                primary_cross_sections.append(cross_section)
        primary_interaction_collection = LI.crosssections.InteractionCollection(primary_type, primary_cross_sections)
        
        # Initialize secondary processes and define secondary InteractionCollection objects
        secondary_decays = {}
        # Also keep track of the minimum decay width for defining the position distribution later
        self.DN_min_decay_width = 0
        # Loop over available decays, group by parent type
        for decay in self.DN_processes.decays:
            secondary_type = LI.dataclasses.Particle.ParticleType(decay.dec_case.nu_parent.pdgid)
            if secondary_type not in secondary_decays.keys():
                secondary_decays[secondary_type] = []
            secondary_decays[secondary_type].append(decay)
            total_decay_width = decay.TotalDecayWidth(secondary_type)
            if total_decay_width < self.DN_min_decay_width:
                self.DN_min_decay_width = total_decay_width
        # Now make the list of secondary cross section collections
        # Add new secondary injection and physical processes at the same time
        fid_vol = self.GetFiducialVolume() # find fiducial volume for secondary position distirbutions
        secondary_interaction_collections = []
        for secondary_type,decay_list in secondary_decays.items():
            
            # Define a sedcondary injection distribution
            secondary_injection_process = LI.injection.InjectionProcess()
            secondary_physical_process = LI.injection.PhysicalProcess()
            secondary_injection_process.primary_type = secondary_type
            secondary_physical_process.primary_type = secondary_type
            
            # Add the secondary position distribution
            if fid_vol is not None:
                secondary_injection_process.AddInjectionDistribution(LI.distributions.SecondaryPositionDistribution(fid_vol))
            else:
                secondary_injection_process.AddInjectionDistribution(LI.distributions.SecondaryPositionDistribution())
            
            self.secondary_injection_processes.append(secondary_injection_process)
            self.secondary_physical_processes.append(secondary_physical_process)

            secondary_interaction_collections.append(LI.crosssections.InteractionCollection(secondary_type, decay_list))

        self.SetCrossSections(primary_interaction_collection,secondary_interaction_collections)
    
        

    def GetFiducialVolume(self):
        """
        :return: identified fiducial volume for the experiment, None if not found
        """
        fid_vol = None
        for sector in self.detector_model.Sectors:
            if self.experiment in fid_vol_dict.keys():
                if sector.name==fid_vol_dict[self.experiment]:
                    fid_vol = sector.geo
        return fid_vol

    def GetDetectorModelTargets(self):
        """
        Determines the targets that exist inside the detector model
        :return: lists of targets and strings
        :rtype: (list<ParticleType>, list<str>)
        """
        count = 0
        targets = []
        target_strs = []
        while self.detector_model.Materials.HasMaterial(count):
            for _target in self.detector_model.Materials.GetMaterialTargets(count):
                if _target not in targets:
                    targets.append(_target)
                if str(_target).find('Nucleus') == -1: 
                    continue
                else:
                    target_str = str(_target)[str(_target).find('Type')+5:str(_target).find('Nucleus')]
                    if target_str == 'H': target_str = 'H1'
                    if target_str not in target_strs:
                        target_strs.append(target_str)
            count += 1
        return targets, target_strs
    
    def SetCrossSections(self,
                         primary_interaction_collection,
                         secondary_interaction_collections):
        """
        Set cross sections for the primary and secondary processes
        :param InteractionCollection primary_interaction_collection: The cross section collection for the primary process
        :param list<InteractionCollection> secondary_interaction_collections: The list of cross section collections for the primary process
        """
        # Set primary cross sections
        self.primary_injection_process.interactions = primary_interaction_collection
        self.primary_physical_process.interactions = primary_interaction_collection
        
        # Loop through secondary processes
        for sec_inj,sec_phys in zip(self.secondary_injection_processes,
                                    self.secondary_physical_processes):
            assert(sec_inj.primary_type == sec_phys.primary_type)
            record = LI.dataclasses.InteractionRecord()
            record.signature.primary_type = sec_inj.primary_type
            found_collection = False
            # Loop through possible seconday cross sections
            for sec_xs in secondary_interaction_collections:
                # Match cross section collection on  the primary type
                if sec_xs.MatchesPrimary(record):
                    sec_inj.interactions = sec_xs
                    sec_phys.interactions = sec_xs
                    found_collection = True
            if(not found_collection):
                print('Couldn\'t find cross section collection for secondary particle %s; Exiting'%record.primary_type)
                exit(0)



    
    def Initialize(self):
        
        # Define stopping condition
        # TODO: make this more general
        def StoppingCondition(datum):
            return True
        
        # Define the injector object
        self.injector = LI.injection.Injector(self.events_to_inject,
                                              self.detector_model, 
                                              self.primary_injection_process, 
                                              self.secondary_injection_processes,
                                              self.random)
        
        self.injector.SetStoppingCondition(StoppingCondition)

        # Define the weighter object
        self.weighter = LI.injection.LeptonTreeWeighter([self.injector],
                                                        self.detector_model, 
                                                        self.primary_physical_process, 
                                                        self.secondary_physical_processes)
        
    def GenerateEvents(self,N=None):
        if N is None:
            N = self.events_to_inject
        count = 0
        while (self.injector.InjectedEvents() < self.events_to_inject) \
               and (count < N):
            print('Injecting Event %d/%d  '%(count,N),end='\r')
            tree = self.injector.GenerateEvent()
            self.events.append(tree)
            count += 1
        DN_processes = getattr(self,"DN_processes")
        if DN_processes is not None:
            DN_processes.SaveCrossSectionTables()
        return self.events

    def SaveEvents(self,filename):
        fout = h5py.File(filename,'w')
        fout.attrs['num_events'] = len(self.events)
        for ie,event in enumerate(self.events):
            print('Saving Event %d/%d  '%(ie,len(self.events)),end='\r')
            event_group = fout.require_group("event%d"%ie)
            event_group.attrs['event_weight'] = self.weighter.EventWeight(event)
            event_group.attrs['num_interactions'] = len(event.tree)
            for id,datum in enumerate(event.tree):
                interaction_group = event_group.require_group("interaction%d"%id)
                
                # Add metadata on interaction signature
                interaction_group.attrs['primary_type'] = str(datum.record.signature.primary_type)
                interaction_group.attrs['target_type'] = str(datum.record.signature.target_type)
                for isec,secondary in enumerate(datum.record.signature.secondary_types):
                    interaction_group.attrs['secondary_type%d'%isec] = str(secondary)

                # Save vertex as dataset
                interaction_group.create_dataset('vertex',data=np.array(datum.record.interaction_vertex,dtype=float))
                
                # Save each four-momenta as a dataset
                interaction_group.create_dataset('primary_momentum',data=np.array(datum.record.primary_momentum,dtype=float))
                interaction_group.create_dataset('target_momentum',data=np.array(datum.record.primary_momentum,dtype=float))
                for isec_momenta,sec_momenta in enumerate(datum.record.secondary_momenta):
                    interaction_group.create_dataset('secondary_momentum%d'%isec_momenta,data=np.array(sec_momenta,dtype=float))

        fout.close()
        DN_processes = getattr(self,"DN_processes")
        if DN_processes is not None:
            DN_processes.SaveCrossSectionTables()

