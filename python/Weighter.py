from . import utilities as _utilities
from . import math as _math
from . import dataclasses as _dataclasses
from . import geometry as _geometry
from . import detector as _detector
from . import interactions as _interactions
from . import distributions as _distributions
from . import injection as _injection
from . import Injector as _Injector_module
from .Injector import _is_trampoline
from ._validation import validate_reweighting_compatibility

from typing import Tuple, List, Dict, Optional, Union, Callable
from typing import TYPE_CHECKING
import math
from numbers import Real
import warnings

import numpy as np

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


def _checked_weight(value, label="event weight", *label_args):
    """Validate a real scalar, formatting diagnostic labels only on failure."""
    if not isinstance(value, Real):
        raise _utilities.WeightCalculationError(
            "{} must be a real scalar, got {!r}".format(
                label.format(*label_args), value))
    # A negative Real (e.g. Fraction) can underflow to -0.0 during float().
    if value < 0:
        raise _utilities.WeightCalculationError(
            "{} must be finite and nonnegative, got {!r}".format(
                label.format(*label_args), value))
    try:
        weight = float(value)
    except OverflowError as exc:
        raise _utilities.WeightCalculationError(
            "{} is not finite: {!r}".format(
                label.format(*label_args), value)) from exc
    if not math.isfinite(weight) or weight < 0:
        raise _utilities.WeightCalculationError(
            "{} must be finite and nonnegative, got {!r}".format(
                label.format(*label_args), value))
    return weight


class Weighter:
    """
    A wrapper for the C++ Weighter class, handling event weight calculations.

    The pooled weight (``__call__`` / ``event_weight``) includes the optional
    ``event_factor``. The per-vertex ``interaction_probabilities`` and
    ``survival_probabilities`` expose native quantities for one chosen injector
    without this factor. Still finer per-vertex factors -- the individual generation
    and physical probability terms -- are reachable through the bound
    ``siren.injection.PrimaryProcessWeighter`` /
    ``siren.injection.SecondaryProcessWeighter`` classes, which expose
    ``InteractionProbability``, ``NormalizedPositionProbability``,
    ``PhysicalProbability``, ``GenerationProbability`` and ``EventWeight`` per
    interaction datum.
    """

    # Default for pickles written before event_factor existed.
    __event_factor = None

    def __init__(self,
        *args,
        injectors: Optional[List[_Injector]] = None,
        detector_model: Optional[DetectorModel] = None,
        primary_type: Optional[_dataclasses.ParticleType] = None,
        primary_interactions: Optional[Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]]] = None,
        primary_physical_distributions: Optional[List[_distributions.WeightableDistribution]] = None,
        secondary_interactions: Optional[Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]]] = None,
        secondary_physical_distributions: Optional[Dict[_dataclasses.ParticleType, List[_distributions.WeightableDistribution]]] = None,
        primary_physical: Optional[List[_distributions.WeightableDistribution]] = None,
        secondary_physical: Optional[Dict[_dataclasses.ParticleType, List[_distributions.WeightableDistribution]]] = None,
        overrides: Optional[Dict[str, object]] = None,
        event_factor: Optional[Callable[[InteractionTree], float]] = None,
    ):
        """
        Initialize the Weighter with interactions and physical processes.

        Two calling conventions are supported.

        Legacy keyword form::

            Weighter(injectors=[...], detector_model=..., primary_type=...,
                     primary_interactions=..., primary_physical_distributions=...,
                     secondary_interactions=..., secondary_physical_distributions=...)

        Spec form inherits the first injector's detector, types, and vertex
        physical models/distributions for the shared physical target::

            Weighter(injector)

        ``event_factor(tree)`` multiplies the physical event weight once, after
        combining injectors. It must be deterministic for a fixed tree and
        model, leave the tree unchanged, and return a finite nonnegative scalar.
        For correlations, supply the joint/reference physical density ratio;
        the sampler must cover the joint density's support. Generation densities
        and per-vertex factors are unchanged. ``None`` applies no correction.

        ``primary_physical`` and ``secondary_physical`` replace the inherited
        distributions, including when explicitly empty.

        Injectors without Vertex specifications default to their sampling
        models and no physical distributions.

        ``overrides`` (spec form only) is a dict of legacy field names
        (``detector_model``, ``primary_type``, ``primary_interactions``,
        ``secondary_interactions``) to override what would otherwise be
        inherited from the injectors.

        Args:
            injectors: List of injector objects (legacy keyword form).
            detector_model: The detector model.
            primary_type: The primary particle type.
            primary_interactions: Dictionary of primary particle interactions.
            primary_physical_distributions: List of primary physical distributions.
            secondary_interactions: Dictionary of secondary particle interactions.
            secondary_physical_distributions: Dictionary of secondary physical distributions.
            primary_physical: Primary physical distributions (spec form).
            secondary_physical: Secondary physical distributions (spec form).
            overrides: Legacy-field overrides for the spec form.
            event_factor: Optional physical factor evaluated on the whole tree.

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
        self.event_factor = event_factor

        spec_injectors = self.__detect_spec_injectors(args, injectors)

        if spec_injectors is not None:
            self.__init_from_injectors(
                spec_injectors, primary_physical, secondary_physical, overrides)
            return

        if len(args) == 1 and injectors is None:
            # Legacy form with injectors passed positionally as a list.
            injectors = args[0]
        elif len(args) > 1:
            raise TypeError(
                "Weighter() accepts either legacy keyword arguments or "
                "one-or-more positional Injector arguments, not {} "
                "positional arguments".format(len(args)))

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

    @staticmethod
    def __detect_spec_injectors(args, injectors_kwarg):
        """Return the positional injectors list if this is the spec-form call.

        Spec form is detected by one-or-more positional args that are each an
        Injector (python wrapper or raw C++ engine), with no ``injectors=``
        keyword given. A single positional list (the legacy form) is left for
        the caller to handle.
        """
        if injectors_kwarg is not None:
            return None
        if len(args) == 0:
            return None
        if len(args) == 1 and isinstance(args[0], list):
            return None
        if all(isinstance(a, (_Injector, _PyInjector)) for a in args):
            return list(args)
        return None

    def __init_from_injectors(self, injectors, primary_physical,
                               secondary_physical, overrides):
        """Inherit the first injector's physical target, then apply overrides."""
        overrides = overrides or {}

        self.injectors = injectors

        first = injectors[0]
        engine0 = first.engine if isinstance(first, _PyInjector) else first
        primary_proc = engine0.GetPrimaryProcess()
        primary_vertex = first.primary if isinstance(first, _PyInjector) else None

        self.__detector_model = overrides.get(
            "detector_model",
            first.detector_model if isinstance(first, _PyInjector)
            else engine0.GetDetectorModel())
        self.__primary_type = overrides.get(
            "primary_type", primary_proc.primary_type)
        primary_interactions = (
            list(primary_proc.interactions.GetCrossSections())
            + list(primary_proc.interactions.GetDecays()))
        if primary_vertex is not None and primary_vertex.physical_interactions is not None:
            primary_interactions = list(primary_vertex.physical_interactions)
        self.__primary_interactions = overrides.get(
            "primary_interactions", primary_interactions)

        secondary_interactions = overrides.get("secondary_interactions", None)
        if secondary_interactions is None:
            secondary_interactions = {}
            for ptype, sproc in engine0.GetSecondaryProcessMap().items():
                secondary_interactions[ptype] = (
                    list(sproc.interactions.GetCrossSections())
                    + list(sproc.interactions.GetDecays()))
            if isinstance(first, _PyInjector):
                for vertex in first.secondaries:
                    if vertex.physical_interactions is not None:
                        secondary_interactions[vertex._resolved_particle] = list(
                            vertex.physical_interactions)
        self.__secondary_interactions = secondary_interactions

        if primary_physical is None:
            primary_physical = primary_vertex.physical if primary_vertex is not None else []
        if secondary_physical is None:
            secondary_physical = {
                vertex._resolved_particle: list(vertex.physical)
                for vertex in first.secondaries
            } if isinstance(first, _PyInjector) else {}
        self.__primary_physical_distributions = list(primary_physical)
        self.__secondary_physical_distributions = {
            ptype: list(dists) for ptype, dists in secondary_physical.items()}

    @property
    def event_factor(self) -> Optional[Callable[[InteractionTree], float]]:
        """Whole-tree physical multiplier; ``None`` disables it.

        Saving or pickling a weighter with this callback is unsupported.
        Shallow copies retain the callback; deep copies copy it with the
        remaining state, subject to the contained objects' copy support.
        """
        return self.__event_factor

    @event_factor.setter
    def event_factor(self, factor):
        if factor is not None and not callable(factor):
            raise TypeError("event_factor must be callable or None")
        self.__event_factor = factor

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

        Note:
            A python Injector wrapper need not already have built its
            underlying engine -- it is unwrapped via its ``engine`` property
            (which builds it lazily) at weighter-initialize time.
        """

        if self.__weighter is not None:
            raise ValueError("Cannot set injectors after weighter has been initialized.")
        if not isinstance(injectors, list):
            raise TypeError("Injectors must be a list.")
        if not all(isinstance(injector, (_Injector, _PyInjector)) for injector in injectors):
            raise TypeError("All injectors must be of type Injector.")
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
        Calculate the event weight, including ``event_factor`` when configured.

        Args:
            interaction_tree: The interaction tree to weight.

        Returns:
            float: The calculated event weight.
        """

        base = _checked_weight(self.engine.EventWeight(interaction_tree),
                               "base event weight")
        if self.event_factor is None:
            return base
        factor = _checked_weight(self.event_factor(interaction_tree), "event factor")
        return _checked_weight(base * factor,
                               "corrected event weight (base={!r}, event_factor={!r})",
                               base, factor)

    def event_weight(self, interaction_tree: InteractionTree) -> float:
        """Alias for ``__call__``."""
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

        Phase space maps are archived with their processes. Configurations that
        still cannot survive the round-trip -- Python trampoline-derived
        interactions or distributions, a Python injector's stopping condition,
        and ``event_factor`` -- raise ``NotSerializableError`` instead.

        Args:
            filename: Base path; the ".siren_weighter" suffix is added.
        """
        from . import errors as _errors
        self._guard_event_factor_serializable()
        if self.__weighter is None:
            self.__initialize_weighter()
        # Building spec injectors resolves their models and expansion callbacks.
        self._guard_serializable(_errors)
        self.__weighter.SaveWeighter(filename)

    def _guard_event_factor_serializable(self):
        if self.event_factor is not None:
            from .errors import NotSerializableError
            raise NotSerializableError(
                "Weighter serialization does not support event_factor.",
                offenders=["event_factor"])

    def __reduce_ex__(self, protocol):
        self._guard_event_factor_serializable()
        return super().__reduce_ex__(protocol)

    def __copy__(self):
        return self._copy_state()

    def __deepcopy__(self, memo):
        return self._copy_state(memo)

    def _copy_state(self, memo=None):
        """Copy in-memory state without invoking the weighter's pickle guard."""
        from copy import deepcopy
        from types import MemberDescriptorType

        cls = type(self)
        result = object.__new__(cls)
        if memo is not None:
            memo[id(self)] = result
            result.__dict__ = deepcopy(self.__dict__, memo)
        else:
            result.__dict__ = self.__dict__.copy()
        # Slots can be inherited or name-mangled; their descriptors handle both.
        for base in cls.__mro__:
            for member in vars(base).values():
                if isinstance(member, MemberDescriptorType):
                    try:
                        value = member.__get__(self, cls)
                    except AttributeError:
                        continue
                    member.__set__(result, value if memo is None
                                   else deepcopy(value, memo))
        return result

    def _guard_serializable(self, _errors):
        """Reject configurations the native archive cannot preserve."""
        self._guard_event_factor_serializable()
        offenders = []
        for injector in (self.__injectors or []):
            if isinstance(injector, _PyInjector):
                try:
                    injector._guard_serializable(_errors, for_weighter=True)
                except _errors.NotSerializableError as exc:
                    offenders.extend("injector: " + o for o in exc.offenders)
        for interaction in self.__primary_interactions:
            if _is_trampoline(interaction):
                offenders.append(
                    "primary interaction {!r} is a Python subclass (not "
                    "serializable)".format(type(interaction).__name__))
        for dist in self.__primary_physical_distributions:
            if _is_trampoline(dist):
                offenders.append(
                    "primary distribution {!r} is a Python subclass (not "
                    "serializable)".format(type(dist).__name__))
        for stype, interactions in self.__secondary_interactions.items():
            for interaction in interactions:
                if _is_trampoline(interaction):
                    offenders.append(
                        "secondary interaction {!r} for type {} is a Python "
                        "subclass (not serializable)".format(
                            type(interaction).__name__, str(stype)))
        for stype, dists in self.__secondary_physical_distributions.items():
            for dist in dists:
                if _is_trampoline(dist):
                    offenders.append(
                        "secondary distribution {!r} for type {} is a Python "
                        "subclass (not serializable)".format(
                            type(dist).__name__, str(stype)))
        if offenders:
            raise _errors.NotSerializableError(
                "this weighter cannot be saved without silently changing "
                "physics on reload:\n  - " + "\n  - ".join(offenders),
                offenders=offenders)

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

        Load with ``event_factor=None``; a factor can be attached after loading.

        Args:
            filename: Base path; the ".siren_weighter" suffix is added.
        """
        self._guard_event_factor_serializable()
        if self.__injectors is not None:
            injectors = [injector.engine if isinstance(injector, _PyInjector) else injector
                         for injector in self.__injectors]
        else:
            injectors = []
        self.__weighter = _Weighter(injectors, filename)
        self.__detector_model = self.__weighter.GetDetectorModel()
        primary_process = self.__weighter.GetPrimaryPhysicalProcess()
        self.__primary_type = primary_process.primary_type
        self.__primary_interactions = list(
            primary_process.interactions.GetCrossSections()) + list(
            primary_process.interactions.GetDecays())
        self.__primary_physical_distributions = list(primary_process.distributions)
        self.__secondary_interactions = {}
        self.__secondary_physical_distributions = {}
        for secondary_process in self.__weighter.GetSecondaryPhysicalProcesses():
            stype = secondary_process.primary_type
            self.__secondary_interactions[stype] = list(
                secondary_process.interactions.GetCrossSections()) + list(
                secondary_process.interactions.GetDecays())
            self.__secondary_physical_distributions[stype] = list(
                secondary_process.distributions)

    def weight_all(self, events) -> "np.ndarray":
        """
        Calculate weights by calling ``self(event)`` for each event.

        Args:
            events: A list of InteractionTree objects.

        Returns:
            numpy.ndarray: The calculated event weights.

        Invalid results, including subclass corrections, raise
        WeightCalculationError with the event and its zero-based event_index.
        """
        weights = []
        for index, event in enumerate(events):
            try:
                weights.append(_checked_weight(self(event)))
            except _utilities.WeightCalculationError as exc:
                error = _utilities.WeightCalculationError(
                    "Event {}: {}".format(index, exc))
                error.event_index = index
                error.event = event
                raise error from exc
        return np.array(weights, dtype=float)

    @property
    def engine(self) -> _Weighter:
        """The raw C++ weighter, excluding the Python ``event_factor``."""
        if self.__weighter is None:
            self.__initialize_weighter()
        return self.__weighter

    def explain(self, interaction_tree: InteractionTree) -> "report.WeightBreakdown":
        """Explain the final weight and its native vertex factors.

        Returns a ``report.WeightBreakdown`` whose ``total`` includes the
        optional ``event_factor``; ``base_total`` and vertex factors remain
        native. Invalid returned factors or corrected weights produce
        event-level flags and a NaN total. Callback exceptions propagate.
        A valid zero base still evaluates the factor; an invalid native
        weight does not.
        """
        from .report import WeightBreakdown
        breakdown = WeightBreakdown.from_engine(
            self.engine.EventWeightWithBreakdown(interaction_tree))
        if self.event_factor is None:
            return breakdown
        try:
            _checked_weight(breakdown.base_total, "base event weight")
        except _utilities.WeightCalculationError:
            return breakdown
        breakdown.event_factor = self.event_factor(interaction_tree)
        try:
            factor = _checked_weight(breakdown.event_factor, "event factor")
            breakdown.total = _checked_weight(
                breakdown.base_total * factor, "corrected event weight")
            breakdown.event_factor = factor
        except _utilities.WeightCalculationError as exc:
            breakdown.flags.append(str(exc))
            breakdown.total = float("nan")
        return breakdown

    def breakdown(self, interaction_tree: InteractionTree) -> "report.WeightBreakdown":
        """Deprecated alias of explain().

        Retained for backward compatibility; emits a DeprecationWarning and
        returns the same ``report.WeightBreakdown`` as explain().
        """
        warnings.warn(
            "Weighter.breakdown() is deprecated, use Weighter.explain() instead",
            DeprecationWarning, stacklevel=2)
        return self.explain(interaction_tree)

    def __initialize_weighter(self):
        return self.__do_initialize_weighter()

    def __do_initialize_weighter(self):
        if self.__injectors is None:
            raise ValueError("Injectors have not been set.")
        if self.__detector_model is None:
            raise ValueError("Detector model has not been set.")
        if self.__primary_type is None:
            raise ValueError("Primary type has not been set.")
        if len(self.__primary_interactions) == 0:
            raise ValueError("Primary interactions have not been set.")

        injectors = [
            injector.engine
            if isinstance(injector, _PyInjector) else injector
            for injector in self.__injectors
        ]

        primary_type = self.primary_type
        primary_interaction_collection = _interactions.InteractionCollection(
            primary_type, self.primary_interactions
        )
        primary_process = _injection.PhysicalProcess(
            primary_type, primary_interaction_collection
        )
        primary_process.distributions = self.primary_physical_distributions

        # Copy weighting mode from the injector's primary process
        inj0_cpp = injectors[0]
        if inj0_cpp is not None:
            primary_process.weighting_mode = inj0_cpp.GetPrimaryProcess().GetWeightingMode()

        # The downstream mode contributes implicit physical factors (most
        # importantly the normalized position density), so compatibility can
        # only be decided here, after the mode has been copied from the
        # injector.  Every injector in a pooled weighter must cover the shared
        # physical measure.
        for injector in injectors:
            validate_reweighting_compatibility(
                injector.GetPrimaryProcess().distributions,
                self.primary_physical_distributions,
                compute_interaction_probability=(
                    primary_process.weighting_mode
                    .compute_interaction_probability
                ),
                compute_position_probability=(
                    primary_process.weighting_mode
                    .compute_position_probability
                ),
            )

        self.__weighter_primary_phys = primary_process

        secondary_interactions = self.secondary_interactions
        secondary_physical_distributions = self.secondary_physical_distributions

        # Copy weighting modes from injector's secondary processes
        sec_modes = {}
        if inj0_cpp is not None:
            for ptype, proc in inj0_cpp.GetSecondaryProcessMap().items():
                sec_modes[ptype] = proc.GetWeightingMode()

        secondary_processes = []
        self.__weighter_secondary_phys = {}
        for secondary_type, sec_ints in secondary_interactions.items():
            secondary_interaction_collection = _interactions.InteractionCollection(
                secondary_type, sec_ints
            )
            secondary_process = _injection.PhysicalProcess(
                secondary_type, secondary_interaction_collection
            )
            if secondary_type in secondary_physical_distributions:
                secondary_process.distributions = secondary_physical_distributions[secondary_type]
            else:
                secondary_process.distributions = []
            if secondary_type in sec_modes:
                secondary_process.weighting_mode = sec_modes[secondary_type]
            secondary_processes.append(secondary_process)
            self.__weighter_secondary_phys[secondary_type] = secondary_process

        self.__weighter = _Weighter(
            injectors,
            self.detector_model,
            primary_process,
            secondary_processes,
        )
