from . import utilities as _utilities
from . import math as _math
from . import dataclasses as _dataclasses
from . import geometry as _geometry
from . import detector as _detector
from . import interactions as _interactions
from . import distributions as _distributions
from . import injection as _injection
from . import Injector as _Injector_module

from typing import Tuple, List, Dict, Optional, Union, Callable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import siren

_Injector = _injection._Injector
_Weighter = _injection._Weighter

_PyInjector = _Injector_module.Injector

ParticleType = _dataclasses.ParticleType
CrossSection = _interactions.CrossSection
Decay = _interactions.Decay
DetectorModel = _detector.DetectorModel
InteractionTree = _dataclasses.InteractionTree

class Weighter:
    """
    A wrapper for the C++ Weighter class, handling event weight calculations.

    Besides the pooled central-value weight (``__call__`` / ``event_weight``),
    two factorized per-vertex quantities are exposed for a single chosen
    injector: ``interaction_probabilities`` and ``survival_probabilities`` (see
    those methods). Still finer per-vertex factors -- the individual generation
    and physical probability terms -- are reachable through the bound
    ``siren.injection.PrimaryProcessWeighter`` /
    ``siren.injection.SecondaryProcessWeighter`` classes, which expose
    ``InteractionProbability``, ``NormalizedPositionProbability``,
    ``PhysicalProbability``, ``GenerationProbability`` and ``EventWeight`` per
    interaction datum.
    """

    def __init__(self, 
        injectors: Optional[List[_Injector]] = None,
        detector_model: Optional[DetectorModel] = None,
        primary_type: Optional[_dataclasses.ParticleType] = None,
        primary_interactions: Optional[Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]]] = None,
        primary_physical_distributions: Optional[List[_distributions.WeightableDistribution]] = None,
        secondary_interactions: Optional[Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]]] = None,
        secondary_physical_distributions: Optional[Dict[_dataclasses.ParticleType, List[_distributions.WeightableDistribution]]] = None,
    ):
        """
        Initialize the Weighter with interactions and physical processes.

        Args:
            injectors: List of injector objects.
            detector_model: The detector model.
            primary_type: The primary particle type.
            primary_interactions: Dictionary of primary particle interactions.
            primary_physical_distributions: List of primary physical distributions.
            secondary_interactions: Dictionary of secondary particle interactions.
            secondary_physical_distributions: Dictionary of secondary physical distributions.

        Note:
            All parameters are optional and can be set later using property setters.
        """

        self.__injectors = None
        self.__detector_model = None

        self.__primary_type = None
        self.__primary_interactions = []
        self.__primary_physical_distributions = []

        self.__secondary_interactions = {}
        self.__secondary_physical_distributions = {}

        self.__weighter = None

        if injectors is not None:
            self.injectors = injectors
        if detector_model is not None:
            self.__detector_model = detector_model
        if primary_type is not None:
            self.__primary_type = primary_type
        if primary_interactions is not None:
            self.__primary_interactions = primary_interactions
        if primary_physical_distributions is not None:
            self.__primary_physical_distributions = primary_physical_distributions
        if secondary_interactions is not None:
            self.__secondary_interactions = secondary_interactions
        if secondary_physical_distributions is not None:
            self.__secondary_physical_distributions = secondary_physical_distributions

    def __initialize_weighter(self):
        """
        Initialize the internal C++ Weighter object.

        This method creates the C++ Weighter object using the configured parameters.
        It is called automatically when needed and should not be called directly.

        Raises:
            ValueError: If any required attributes are not set.
        """

        if self.__injectors is None:
            raise ValueError("Injectors have not been set.")
        if self.__detector_model is None:
            raise ValueError("Detector model has not been set.")
        if self.__primary_type is None:
            raise ValueError("Primary type has not been set.")
        if len(self.__primary_interactions) == 0:
            raise ValueError("Primary interactions have not been set.")
        if len(self.__primary_physical_distributions) == 0:
            raise ValueError("Primary physical distributions have not been set.")
    
        injectors = [injector._Injector__injector if isinstance(injector, _PyInjector) else injector for injector in self.__injectors]

        primary_type = self.primary_type
        primary_interaction_collection = _interactions.InteractionCollection(
            primary_type, self.primary_interactions
        )
        primary_process = _injection.PhysicalProcess(
            primary_type, primary_interaction_collection
        )
        primary_process.distributions = self.primary_physical_distributions

        secondary_interactions = self.secondary_interactions
        secondary_physical_distributions = self.secondary_physical_distributions

        secondary_processes = []
        for secondary_type, secondary_interactions in secondary_interactions.items():
            secondary_interaction_collection = _interactions.InteractionCollection(
                secondary_type, secondary_interactions
            )
            secondary_process = _injection.PhysicalProcess(
                secondary_type, secondary_interaction_collection
            )
            if secondary_type in secondary_physical_distributions:
                secondary_process.distributions = secondary_physical_distributions[secondary_type]
            else:
                secondary_process.distributions = []
            secondary_processes.append(secondary_process)

        self.__weighter = _Weighter(
            injectors,
            self.detector_model,
            primary_process,
            secondary_processes,
        )

    @property
    def injectors(self) -> List[_Injector]:
        """
        Get the list of injectors.

        Returns:
            List[_Injector]: The current list of injector objects.
        """
        return self.__injectors

    @injectors.setter
    def injectors(self, injectors: List[_Injector]):
        """
        Set the list of injectors.

        Args:
            injectors: A list of Injector objects.

        Raises:
            ValueError: If the weighter has already been initialized.
            TypeError: If the input is not a list of Injector objects.
            ValueError: If any of the injectors are not initialized.
        """

        if self.__weighter is not None:
            raise ValueError("Cannot set injectors after weighter has been initialized.")
        if not isinstance(injectors, list):
            raise TypeError("Injectors must be a list.")
        if not all(isinstance(injector, (_Injector, _PyInjector)) for injector in injectors):
            raise TypeError("All injectors must be of type Injector.")
        if not all(injector._Injector__injector is not None for injector in injectors if isinstance(injector, _PyInjector)):
            raise ValueError("All injectors must be initialized.")
        self.__injectors = injectors

    @property
    def detector_model(self) -> DetectorModel:
        """
        Get the detector model.

        Returns:
            DetectorModel: The current detector model.
        """
        return self.__detector_model

    @detector_model.setter
    def detector_model(self, detector_model: DetectorModel):
        """
        Set the detector model.

        Args:
            detector_model: The DetectorModel object to set.

        Raises:
            ValueError: If the weighter has already been initialized.
            TypeError: If the input is not a DetectorModel object.
        """

        if self.__weighter is not None:
            raise ValueError("Cannot set detector model after weighter has been initialized.")
        if not isinstance(detector_model, DetectorModel):
            raise TypeError("Detector model must be of type DetectorModel.")
        self.__detector_model = detector_model

    @property
    def primary_type(self) -> ParticleType:
        return self.__primary_type

    @primary_type.setter
    def primary_type(self, primary_type: ParticleType):
        if self.__weighter is not None:
            raise ValueError("Cannot set primary type after weighter has been initialized.")
        if not isinstance(primary_type, ParticleType):
            raise TypeError("Primary type must be of type ParticleType.")
        self.__primary_type = primary_type

    @property
    def primary_interactions(self) -> Dict[ParticleType, List[Union[CrossSection, Decay]]]:
        return self.__primary_interactions

    @primary_interactions.setter
    def primary_interactions(self, primary_interactions: List[Union[CrossSection, Decay]]):
        if self.__weighter is not None:
            raise ValueError("Cannot set primary interactions after weighter has been initialized.")
        if not isinstance(primary_interactions, list):
            raise TypeError("Primary interactions must be a list.")
        if not all(isinstance(interaction, (CrossSection, Decay)) for interaction in primary_interactions):
            raise TypeError("All interactions in primary interactions must be of type CrossSection or Decay.")
        self.__primary_interactions = primary_interactions

    @property
    def primary_physical_distributions(self) -> List[_distributions.WeightableDistribution]:
        return self.__primary_physical_distributions

    @primary_physical_distributions.setter
    def primary_physical_distributions(self, primary_physical_distributions: List[_distributions.WeightableDistribution]):
        if self.__weighter is not None:
            raise ValueError("Cannot set primary physical distributions after weighter has been initialized.")
        if not isinstance(primary_physical_distributions, list):
            raise TypeError("Primary physical distributions must be a list.")
        if not all(isinstance(distribution, _distributions.WeightableDistribution) for distribution in primary_physical_distributions):
            raise TypeError("All distributions in primary physical distributions must be of type WeightableDistribution.")
        self.__primary_physical_distributions = primary_physical_distributions

    @property
    def secondary_interactions(self) -> Dict[ParticleType, List[Union[CrossSection, Decay]]]:
        return self.__secondary_interactions

    @secondary_interactions.setter
    def secondary_interactions(self, secondary_interactions: Dict[ParticleType, List[Union[CrossSection, Decay]]]):
        if self.__weighter is not None:
            raise ValueError("Cannot set secondary interactions after weighter has been initialized.")
        if not isinstance(secondary_interactions, dict):
            raise TypeError("Secondary interactions must be a dictionary.")
        if not all(isinstance(particle_type, ParticleType) for particle_type in secondary_interactions.keys()):
            raise TypeError("All keys in secondary interactions must be of type ParticleType.")
        if not all(isinstance(interactions, list) for interactions in secondary_interactions.values()):
            raise TypeError("All values in secondary interactions must be lists.")
        if not all(isinstance(interaction, (CrossSection, Decay)) for interactions in secondary_interactions.values() for interaction in interactions):
            raise TypeError("All interactions in secondary interactions must be of type CrossSection or Decay.")
        self.__secondary_interactions = secondary_interactions

    @property
    def secondary_physical_distributions(self) -> Dict[ParticleType, List[_distributions.WeightableDistribution]]:
        return self.__secondary_physical_distributions

    @secondary_physical_distributions.setter
    def secondary_physical_distributions(self, secondary_physical_distributions: Dict[ParticleType, List[_distributions.WeightableDistribution]]):
        if self.__weighter is not None:
            raise ValueError("Cannot set secondary physical distributions after weighter has been initialized.")
        if not isinstance(secondary_physical_distributions, dict):
            raise TypeError("Secondary physical distributions must be a dictionary.")
        if not all(isinstance(particle_type, ParticleType) for particle_type in secondary_physical_distributions.keys()):
            raise TypeError("All keys in secondary physical distributions must be of type ParticleType.")
        if not all(isinstance(distributions, list) for distributions in secondary_physical_distributions.values()):
            raise TypeError("All values in secondary physical distributions must be lists.")
        if not all(isinstance(distribution, _distributions.WeightableDistribution) for distributions in secondary_physical_distributions.values() for distribution in distributions):
            raise TypeError("All distributions in secondary physical distributions must be of type WeightableDistribution.")
        self.__secondary_physical_distributions = secondary_physical_distributions

    def __call__(self, interaction_tree: InteractionTree) -> float:
        """
        Calculate the event weight for a given interaction tree.

        This method initializes the weighter if necessary and then calculates the event weight.

        Args:
            interaction_tree: The interaction tree to weight.

        Returns:
            float: The calculated event weight.
        """

        if self.__weighter is None:
            self.__initialize_weighter()
        return self.__weighter.EventWeight(interaction_tree)

    def event_weight(self, interaction_tree: InteractionTree) -> float:
        """
        Calculate the event weight for a given interaction tree.

        This method is an alias for __call__ and provides the same functionality.

        Args:
            interaction_tree: The interaction tree to weight.

        Returns:
            float: The calculated event weight.
        """
        return self(interaction_tree)

    def interaction_probabilities(self, interaction_tree: InteractionTree, i_inj: int = 0) -> List[float]:
        """
        Per-vertex physical interaction probability for one chosen injector.

        Returns one value per interaction datum in ``interaction_tree`` (in tree
        order): the probability that the primary interacts within injector
        ``i_inj``'s injection bounds -- i.e. over the segment spanning that
        injector's PrimaryInjectionBounds (or SecondaryInjectionBounds for
        non-root vertices). This is the interaction-region segment only.

        It is DISJOINT from ``survival_probabilities`` (below): survival covers
        the pre-injection segment leading up to the injection region, while this
        covers the injection region itself. The two are NOT complementary
        probabilities of one another (they integrate the column density over
        different, non-overlapping segments), so do not expect them to sum to 1.

        Args:
            interaction_tree: The interaction tree to evaluate.
            i_inj: Index into this weighter's injector list (default 0).

        Returns:
            List[float]: One interaction probability per interaction datum.
        """
        if self.__weighter is None:
            self.__initialize_weighter()
        return list(self.__weighter.GetInteractionProbabilities(interaction_tree, i_inj))

    def survival_probabilities(self, interaction_tree: InteractionTree, i_inj: int = 0) -> List[float]:
        """
        Per-vertex survival probability over the pre-injection segment.

        Returns one value per interaction datum in ``interaction_tree`` (in tree
        order): the probability that the primary survives (does not interact)
        over the pre-injection segment, from ``primary_initial_position`` up to
        the near edge of injector ``i_inj``'s injection region (the first element
        of its PrimaryInjectionBounds / SecondaryInjectionBounds).

        This segment is DISJOINT from the one measured by
        ``interaction_probabilities``: survival covers everything before the
        injection region, the interaction probability covers the injection region
        itself. They are therefore not complementary probabilities and need not
        sum to 1.

        Args:
            interaction_tree: The interaction tree to evaluate.
            i_inj: Index into this weighter's injector list (default 0).

        Returns:
            List[float]: One survival probability per interaction datum.
        """
        if self.__weighter is None:
            self.__initialize_weighter()
        return list(self.__weighter.GetSurvivalProbabilities(interaction_tree, i_inj))

    def save(self, filename: str):
        """
        Serialize the weighter to ``<filename>.siren_weighter``.

        Args:
            filename: Base path; the ".siren_weighter" suffix is added.
        """
        if self.__weighter is None:
            self.__initialize_weighter()
        self.__weighter.SaveWeighter(filename)

    def load(self, filename: str):
        """
        Restore the weighter from ``<filename>.siren_weighter``.

        Constructs the underlying C++ weighter via its ``(injectors, filename)``
        constructor, which loads the detector model and physical processes from
        the file. The detector / interactions / distributions therefore do NOT
        need to be configured first (unlike a freshly built weighter). Only the
        injectors are used -- to bind the weighter to your live injection
        processes for the generation-probability cancellation; if they were not
        set, the injectors serialized in the file are used instead.

        Args:
            filename: Base path; the ".siren_weighter" suffix is added.
        """
        if self.__injectors is not None:
            injectors = [injector._Injector__injector if isinstance(injector, _PyInjector) else injector
                         for injector in self.__injectors]
        else:
            injectors = []
        self.__weighter = _Weighter(injectors, filename)
