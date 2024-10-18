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

from typing import Tuple, List, Dict, Optional, Union, Callable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import siren

_Injector = _injection._Injector

ParticleType = _dataclasses.ParticleType
CrossSection = _interactions.CrossSection
Decay = _interactions.Decay
DetectorModel = _detector.DetectorModel
SIREN_random = _utilities.SIREN_random
PrimaryInjectionDistribution = _distributions.PrimaryInjectionDistribution
SecondaryInjectionDistribution = _distributions.SecondaryInjectionDistribution
SecondaryInjectionProcess = _injection.SecondaryInjectionProcess
InteractionTreeDatum = _dataclasses.InteractionTreeDatum

class Injector:
    def __init__(
        self,
        number_of_events: Optional[int] = None,
        detector_model: Optional[_detector.DetectorModel] = None,
        seed: Optional[int] = None,
        primary_type: Optional[_dataclasses.ParticleType] = None,
        primary_interactions: Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]] = None,
        primary_injection_distributions: List[_distributions.PrimaryInjectionDistribution] = None,
        secondary_interactions: Optional[Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]]] = None,
        secondary_injection_distributions: Optional[Dict[_dataclasses.ParticleType, List[_distributions.SecondaryInjectionDistribution]]] = None,
        stopping_condition: Optional[Callable[[_dataclasses.InteractionTreeDatum, int], bool]] = None,
    ):
        self.__seed = None
        self.__number_of_events = 0
        self.__detector_model = None

        self.__primary_type = None
        self.__primary_interactions = []
        self.__primary_injection_distributions = []

        self.__secondary_interactions = {}
        self.__secondary_injection_distributions = {}
        self.__stopping_condition = None

        self.__injector = None

        if seed is not None:
            self.__seed = seed
        if number_of_events is not None:
            self.__number_of_events = number_of_events
        if detector_model is not None:
            self.__detector_model = detector_model
        if primary_type is not None:
            self.__primary_type = primary_type
        if primary_interactions is not None:
            self.__primary_interactions = primary_interactions
        if primary_injection_distributions is not None:
            self.__primary_injection_distributions = primary_injection_distributions
        if secondary_interactions is not None:
            self.__secondary_interactions = secondary_interactions
        if secondary_injection_distributions is not None:
            self.__secondary_injection_distributions = secondary_injection_distributions
        if stopping_condition is not None:
            self.__stopping_condition = stopping_condition


    def __initialize_injector(self):
        if self.__seed is None:
            random = _utilities.SIREN_random()
            self.__seed = random.get_seed()
        else:
            random = _utilities.SIREN_random(self.__seed)

        if self.__number_of_events is None:
            raise ValueError("number_of_events must be provided")
        elif self.__number_of_events <= 0:
            raise ValueError("number_of_events must be positive")

        if self.__detector_model is None:
            raise ValueError("detector_model must be provided")

        if self.__primary_type is None:
            raise ValueError("primary_type must be provided")

        if len(self.__primary_interactions) == 0:
            raise ValueError("primary_interactions must be provided")

        if len(self.__primary_injection_distributions) == 0:
            raise ValueError("primary_injection_distributions must be provided")

        if list(sorted(self.__secondary_interactions.keys())) != list(sorted(self.__secondary_injection_distributions.keys())):
            raise ValueError("secondary_interactions and secondary_injection_distributions must have the same keys")

        primary_type = self.primary_type
        primary_interaction_collection = _interactions.InteractionCollection(
            primary_type, self.primary_interactions
        )
        primary_process = _injection.PrimaryInjectionProcess(
            primary_type, primary_interaction_collection
        )
        primary_process.distributions = self.primary_injection_distributions

        secondary_interactions = self.secondary_interactions
        secondary_injection_distributions = self.secondary_injection_distributions

        secondary_interaction_collections = []
        secondary_processes = []
        for secondary_type, secondary_interactions in secondary_interactions.items():
            secondary_interaction_collection = _interactions.InteractionCollection(
                secondary_type, secondary_interactions
            )
            secondary_process = SecondaryInjectionProcess(
                secondary_type, secondary_interaction_collection
            )
            secondary_process.distributions = secondary_injection_distributions[secondary_type]
            secondary_processes.append(secondary_process)

        self.__injector = _Injector(
            self.number_of_events,
            self.detector_model,
            primary_process,
            secondary_processes,
            random,
        )

        if self.__stopping_condition is not None:
            self.__injector.SetStoppingCondition(self.__stopping_condition)

    # Custom method to retrieve the internal state for pickling
    def __getstate__(self):
        # The seed and stopping condition are the only things we cannot serialize
        state = (self.__seed, self.__stopping_condition, self.__injector.__getstate__())
        return state

    # Custom method to restore the state from a pickle
    def __setstate__(self, state):

        self.__seed, self.__stopping_condition, injector_state = state

        # Create a new instance of the C++ class and restore its state
        self.__injector = _Injector.__new__(_Injector)  # Create an empty instance
        if self.__injector is None:
            raise TypeError("Failed to create C++ Injector object")
        print(self.__injector)
        self.__injector.__setstate__(injector_state)
        print(self.__injector)
        self.__number_of_events = self.__injector.EventsToInject()
        self.__detector_model = self.__injector.GetDetectorModel()
        primary_process = self.__injector.GetPrimaryProcess()
        self.__primary_type = primary_process.primary_type
        self.__primary_interactions = list(primary_process.interactions.GetCrossSections()) + list(primary_process.interactions.GetDecays())
        self.__primary_injection_distributions = list(primary_process.distributions)

        self.__secondary_interactions = {}
        self.__secondary_injection_distributions = {}
        for secondary_type, secondary_process in self.__injector.GetSecondaryProcessMap():
            self.__secondary_interactions[secondary_type] = list(secondary_process.interactions.GetCrossSections()) + list(secondary_process.interactions.GetDecays())
            self.__secondary_injection_distributions[secondary_type] = list(secondary_process.distributions)

    @property
    def seed(self):
        return self.__seed

    @seed.setter
    def seed(self, seed):
        self.__seed = seed
        if self.__injector is not None:
            self.__injector.GetRandom().set_seed(seed)

    @property
    def number_of_events(self):
        if self.__injector is not None:
            return self.__injector.EventsToInject()
        return self.__number_of_events

    @number_of_events.setter
    def number_of_events(self, number_of_events):
        if self.__injector is not None:
            raise ValueError("Cannot change the number of events after initialization")
        self.__number_of_events = number_of_events

    @property
    def detector_model(self):
        if self.__injector is not None:
            return self.__injector.GetDetectorModel()
        return self.__detector_model

    @detector_model.setter
    def detector_model(self, detector_model):
        if self.__injector is not None:
            self.__injector.SetDetectorModel(detector_model)
        self.__detector_model = detector_model

    @property
    def primary_type(self):
        return self.__primary_type

    @primary_type.setter
    def primary_type(self, primary_type):
        if self.__injector is not None:
            primary_process = self.__injector.GetPrimaryProcess()
            primary_process.primary_type = primary_type
        self.__primary_type = primary_type

    @property
    def primary_interactions(self):
        return self.__primary_interactions

    @primary_interactions.setter
    def primary_interactions(self, primary_interactions):
        if self.__injector is not None:
            primary_process = self.__injector.GetPrimaryProcess()
            primary_interaction_collection = _interactions.InteractionCollection(
                self.primary_type, primary_interactions
            )
            primary_process.interactions = primary_interaction_collection
        self.__primary_interactions = primary_interactions

    @property
    def primary_injection_distributions(self):
        return self.__primary_injection_distributions

    @primary_injection_distributions.setter
    def primary_injection_distributions(self, primary_injection_distributions):
        if self.__injector is not None:
            primary_process = self.__injector.GetPrimaryProcess()
            primary_process.distributions = primary_injection_distributions
        self.__primary_injection_distributions = primary_injection_distributions

    @property
    def secondary_interactions(self):
        return self.__secondary_interactions

    @secondary_interactions.setter
    def secondary_interactions(self, secondary_interactions):
        if self.__injector is not None:
            secondary_processes = self.__injector.GetSecondaryProcessMap()
            current_secondary_types = sorted(list(secondary_processes.keys()))
            new_secondary_types = sorted(list(secondary_interactions.keys()))
            if current_secondary_types != new_secondary_types:
                raise ValueError("Cannot change the secondary types after initialization")
            for secondary_type, secondary_process in secondary_processes.items():
                secondary_process.interactions = secondary_interactions[secondary_type]
        self.__secondary_interactions = secondary_interactions

    @property
    def secondary_injection_distributions(self):
        return self.__secondary_injection_distributions

    @secondary_injection_distributions.setter
    def secondary_injection_distributions(self, secondary_injection_distributions):
        if self.__injector is not None:
            secondary_processes = self.__injector.GetSecondaryProcesses()
            current_secondary_types = sorted(list(secondary_processes.keys()))
            new_secondary_types = sorted(list(secondary_injection_distributions.keys()))
            if current_secondary_types != new_secondary_types:
                raise ValueError("Cannot change the secondary types after initialization")
            for secondary_type, secondary_process in secondary_injection_distributions.items():
                secondary_process.distributions = secondary_distributions[secondary_type]
        self.__secondary_injection_distributions = secondary_injection_distributions

    @property
    def stopping_condition(self):
        return self.__stopping_condition

    @stopping_condition.setter
    def stopping_condition(self, stopping_condition):
        if self.__injector is not None:
            self.__injector.SetStoppingCondition(stopping_condition)
        self.__stopping_condition = stopping_condition

    @wraps(_Injector.NewRecord)
    def new_record(self):
        return self.__injector.NewRecord()
    new_record.__name__ = "new_record"
    new_record.__doc__ = _Injector.NewRecord.__doc__.replace("NewRecord", "new_record")

    @wraps(_Injector.GenerateEvent)
    def generate_event(self):
        if self.__injector is None:
            self.__initialize_injector()
        return self.__injector.GenerateEvent()
    generate_event.__name__ = "generate_event"
    generate_event.__doc__ = _Injector.GenerateEvent.__doc__.replace("GenerateEvent", "generate_event")

    @property
    def density_variables(self):
        if self.__injector is not None:
            return self.__injector.DensityVariables()
        return None

    @property
    def injected_events(self):
        if self.__injector is not None:
            return self.__injector.InjectedEvents()
        return 0

    @wraps(_Injector.ResetInjectedEvents)
    def reset_injected_events(self):
        if self.__injector is not None:
            self.__injector.ResetInjectedEvents()
    reset_injected_events.__name__ = "reset_injected_events"
    reset_injected_events.__doc__ = _Injector.ResetInjectedEvents.__doc__.replace("ResetInjectedEvents", "reset_injected_events")

    @wraps(_Injector.SaveInjector)
    def save(self, filename):
        self.__injector.SaveInjector(filename)
    save.__name__ = "save"
    save.__doc__ = _Injector.SaveInjector.__doc__.replace("SaveInjector", "save")

    @wraps(_Injector.LoadInjector)
    def load(self, filename):
        self.__injector.LoadInjector(filename)
        # Update Python object state after loading
        self.__number_of_events = self.__injector.EventsToInject()
        self.__detector_model = self.__injector.GetDetectorModel()
        primary_process = self.__injector.GetPrimaryProcess()
        self.__primary_type = primary_process.primary_type
        self.__primary_interactions = list(primary_process.interactions.GetCrossSections()) + list(primary_process.interactions.GetDecays())
        self.__primary_injection_distributions = list(primary_process.distributions)

        self.__secondary_interactions = {}
        self.__secondary_injection_distributions = {}
        for secondary_type, secondary_process in self.__injector.GetSecondaryProcessMap():
            self.__secondary_interactions[secondary_type] = list(secondary_process.interactions.GetCrossSections()) + list(secondary_process.interactions.GetDecays())
            self.__secondary_injection_distributions[secondary_type] = list(secondary_process.distributions)

        self.__stopping_condition = self.__injector.GetStoppingCondition()

