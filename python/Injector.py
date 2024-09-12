from . import utilities as _utilities
from . import math as _math
from . import dataclasses as _dataclasses
from . import geometry as _geometry
from . import detector as _detector
from . import interactions as _interactions
from . import distributions as _distributions
from . import injection as _injection

import collections
from functools import wraps

from typing import Tuple, List, Dict, Optional, Union
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import siren

_Injector = _injection.Injector

ParticleType = _dataclasses.Particle.ParticleType
CrossSection = _interactions.CrossSection
Decay = _interactions.Decay
Interaction = _interactions.Interaction
DetectorModel = _detector.DetectorModel
SIREN_random = _utilities.SIREN_random
PrimaryInjectionDistribution = _distributions.PrimaryInjectionDistribution
SecondaryInjectionDistribution = _distributions.SecondaryInjectionDistribution
SecondaryInjectionProcess = _injection.SecondaryInjectionProcess

class Injector:
    def __init__(
        self,
        number_of_events: int,
        detector_model: "DetectorModel",
        random: "SIREN_random",
        primary_interactions: Dict["ParticleType", List[Union["CrossSection", "Decay", "Interaction"]]],
        primary_injection_distributions: List["PrimaryInjectionDistribution"],
        secondary_interactions: Optional[Dict["ParticleType", List[Union["CrossSection", "Decay", "Interaction"]]]] = None,
        secondary_injection_distributions: Optional[Dict["ParticleType", List["SecondaryInjectionDistribution"]]] = None,
    ):
        self.number_of_events = number_of_events

        self.detector_model = detector_model

        if len(primary_interactions) != 1:
            raise ValueError(f"len(primary_interactions) != 1")

        if (secondary_interactions is None) != (secondary_injection_distributions is None):
            raise ValueError("Neither or both of secondary_interactions and secondary_injection_distributions must be provided")

        if secondary_interactions is None:
            secondary_interactions = dict()
            secondary_injection_distributions = dict()

        self.primary_interactions = primary_interactions
        self.primary_injection_distributions = primary_injection_distributions

        primary_type, primary_interactions = list(primary_interactions.items())[0]

        self.primary_interaction_collection = _interactions.InteractionCollection(
            primary_type, primary_interactions
        )
        self.primary_process = _injection.PrimaryInjectionProcess(
            primary_type, self.primary_interaction_collection
        )
        for dist in primary_injection_distributions:
            self.primary_process.AddPrimaryInjectionDistribution(dist)


        self.secondary_interactions = secondary_interactions
        self.secondary_injection_distributions = secondary_injection_distributions

        self.secondary_interaction_collections = []
        self.secondary_processes = []
        for secondary_type, secondary_interactions in secondary_interactions.items():
            secondary_distributions = self.secondary_injection_distributions[secondary_type]
            secondary_process = SecondaryInjectionProcess(secondary_type, secondary_interactions)
            for dist in secondary_distributions:
                secondary_process.AddSecondaryInjectionDistribution(dist)
            self.secondary_processes.append(secondary_process)

        self.injector = _injection.Injector(
            self.number_of_events,
            self.detector_model,
            self.primary_process,
            self.secondary_processes,
            self.random,
        )

    # TODO define wrapper functions that modify the internal state of the python object
    @wraps(_Injector.SetPrimaryProcess)
    def SetPrimaryProcess(self, primary_process):
        # Get the internals first
        primary_injection_distributions = primary_process.GetPrimaryInjectionDistributions()
        primary_interaction_collection = primary_process.GetInteractionCollection()
        primary_interactions = list(primary_interaction_collection.GetCrossSections()) + list(primary_interaction_collection.GetDecays())

        # Now we can overwite things
        self.injector.SetPrimaryProcess(primary_process)
        self.primary_process = primary_process
        self.primary_injection_distributions = primary_injection_distributions
        self.primary_interaction_collection = primary_interaction_collection
        self.primary_interactions = {primary_process.primary_type: primary_interactions}

    @wraps(_Injector.SetStoppingCondition)
    def SetStoppingCondition(self, stopping_condition):
        self.stopping_condition = stopping_condition
        self.injector.SetStoppingCondition(stopping_condition)

    @wraps(_Injector.AddSecondaryProcess)
    def AddSecondaryProcess(self, secondary_process):
        # Update internal state
        secondary_type = secondary_process.secondary_type
        secondary_distributions = secondary_process.GetSecondaryInjectionDistributions()
        secondary_interaction_collection = secondary_process.GetInteractionCollection()
        secondary_interactions = list(secondary_interaction_collection.GetCrossSections()) + list(secondary_interaction_collection.GetDecays())

        # Update class attributes
        self.secondary_processes.append(secondary_process)
        if secondary_type not in self.secondary_interactions:
            self.secondary_interactions[secondary_type] = []
        self.secondary_interactions[secondary_type].extend(secondary_interactions)
        if secondary_type not in self.secondary_injection_distributions:
            self.secondary_injection_distributions[secondary_type] = []
        self.secondary_injection_distributions[secondary_type].extend(secondary_distributions)

        # Update the underlying C++ object
        self.injector.AddSecondaryProcess(secondary_process)
    
    @wraps(_Injector.GetPrimaryProcess)
    def GetPrimaryProcess(self):
        return self.primary_process

    @wraps(_Injector.GetSecondaryProcessMap)
    def GetSecondaryProcesses(self):
        return self.secondary_processes

    @wraps(_Injector.NewRecord)
    def NewRecord(self):
        return self.injector.NewRecord()

    @wraps(_Injector.SetRandom)
    def SetRandom(self, random):
        self.injector.SetRandom(random)

    @wraps(_Injector.GenerateEvent)
    def GenerateEvent(self):
        return self.injector.GenerateEvent()

    @wraps(_Injector.DensityVariables)
    def DensityVariables(self):
        return self.injector.DensityVariables()

    @wraps(_Injector.Name)
    def Name(self):
        return self.injector.Name()
    
    @wraps(_Injector.GetPrimaryInjectionDistributions)
    def GetPrimaryInjectionDistributions(self):
        return self.primary_injection_distributions

    @wraps(_Injector.GetDetectorModel)
    def GetDetectorModel(self):
        return self.detector_model

    @wraps(_Injector.GetInteractions)
    def GetInteractions(self):
        return self.injector.GetInteractions()

    @wraps(_Injector.InjectedEvents)
    def InjectedEvents(self):
        return self.injector.InjectedEvents()

    @wraps(_Injector.EventsToInject)
    def EventsToInject(self):
        return self.injector.EventsToInject()

    @wraps(_Injector.ResetInjectedEvents)
    def ResetInjectedEvents(self):
        self.injector.ResetInjectedEvents()

    @wraps(_Injector.SaveInjector)
    def SaveInjector(self, filename):
        self.injector.SaveInjector(filename)

    @wraps(_Injector.LoadInjector)
    def LoadInjector(self, filename):
        self.injector.LoadInjector(filename)
        # Update Python object state after loading
        self.primary_process = self.injector.GetPrimaryProcess()
        self.secondary_processes = self.injector.GetSecondaryProcesses()
        self.primary_injection_distributions = self.primary_process.GetPrimaryInjectionDistributions()
        self.primary_interaction_collection = self.primary_process.GetInteractionCollection()
        self.primary_interactions = {self.primary_process.primary_type: list(self.primary_interaction_collection.GetCrossSections()) + list(self.primary_interaction_collection.GetDecays())}
        # Update secondary interactions and distributions
        self.secondary_interactions = {}
        self.secondary_injection_distributions = {}
        for process in self.secondary_processes:
            secondary_type = process.secondary_type
            self.secondary_interactions[secondary_type] = list(process.GetInteractionCollection().GetCrossSections()) + list(process.GetInteractionCollection().GetDecays())
            self.secondary_injection_distributions[secondary_type] = process.GetSecondaryInjectionDistributions()

