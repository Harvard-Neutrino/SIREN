"""Simulation forwards Vertex declarations to Injector and Weighter."""
import difflib
import warnings

from .Injector import Injector
from .Weighter import Weighter
from .vertex import Vertex
from ._util import load_detector
from .generate import generate


class Simulation:
    """Run a declared primary/secondary chain and return weighted Results.

    Each run constructs a fresh Injector. The seed and generation policy are
    forwarded unchanged; validation and compilation belong to Injector/Vertex.
    Deprecated convenience slots translate to Vertex once at construction.
    Model-loader arguments belong in an explicit call or process_options.
    """

    def __init__(self, *, events=None, detector, primary, secondaries=(), seed=None,
                 event_factor=None, max_attempts=None, stopping_condition=None,
                 process_options=None, **legacy):
        from ._simulation_legacy import SLOTS, translate
        unknown = set(legacy) - SLOTS - {'n_events'}
        if unknown:
            key = sorted(unknown)[0]
            suggestions = difflib.get_close_matches(
                key, ['events', 'detector', 'primary', 'secondaries', 'seed',
                      'event_factor', 'max_attempts', 'process_options'] + list(SLOTS), n=1)
            hint = '; did you mean ' + suggestions[0] + '?' if suggestions else ''
            raise TypeError('Unknown Simulation argument: ' + key + hint)
        if 'n_events' in legacy:
            if events is not None:
                raise TypeError("Specify only events, not both events and n_events")
            warnings.warn('n_events is deprecated; use events', DeprecationWarning, stacklevel=2)
            events = legacy.pop('n_events')
        if events is None:
            raise TypeError('Simulation requires events=')
        self._detector_model = load_detector(detector) if isinstance(detector, str) else detector
        self._process_metadata = ()
        if legacy or not isinstance(primary, Vertex):
            primary, secondaries, self._process_metadata = translate(
                primary, secondaries, self._detector_model, process_options or {}, legacy)
        elif process_options:
            raise TypeError('process_options belongs to a named interaction loader')
        self._primary = primary
        self._secondaries = tuple(secondaries)
        self._events = events
        self._seed = seed
        self._max_attempts = max_attempts
        self._stopping_condition = stopping_condition
        if event_factor is not None and not callable(event_factor):
            raise TypeError('event_factor must be callable or None')
        self._event_factor = event_factor
        self._injector = self._weighter = self._last_results = None

    def _build_injector(self):
        return Injector(events=self._events, detector=self._detector_model,
                        primary=self._primary, secondaries=self._secondaries,
                        seed=self._seed, max_attempts=self._max_attempts,
                        stopping_condition=self._stopping_condition)

    def _build_weighter(self, injector):
        return Weighter(injector, event_factor=self._event_factor)

    def run(self, *, on_shortfall='warn', on_failure='retry', progress=None,
            min_efficiency=None):
        """Generate a fresh run; strict production uses both policies='raise'."""
        injector = self._build_injector()
        weighter = self._build_weighter(injector)
        results = generate(injector, weighter, events=self._events,
                           on_shortfall=on_shortfall, on_failure=on_failure,
                           progress=progress, min_efficiency=min_efficiency)
        self._injector, self._weighter, self._last_results = injector, weighter, results
        return results

    def reweight(self, *, physical_energy=None, physical_direction=None,
                 physical_distributions=None, interactions=None,
                 secondary_interactions=None, secondary_physical_distributions=None,
                 event_factor=...):
        """Reweight the last Results using its retained generation snapshot."""
        if self._last_results is None:
            raise RuntimeError('Must call run() before reweight()')
        from . import distributions as d
        if physical_distributions is not None and (
                physical_energy is not None or physical_direction is not None):
            raise TypeError('physical_distributions cannot be mixed with named physical slots')
        physical = physical_distributions
        if physical is None and (physical_energy is not None or physical_direction is not None):
            physical = list(self._last_results._snapshot_weighter().primary_physical_distributions)
            for replacement, base in ((physical_energy, d.PrimaryEnergyDistribution),
                                      (physical_direction, d.PrimaryDirectionDistribution)):
                if replacement is not None:
                    physical = [x for x in physical if not isinstance(x, base)]
                    physical.append(replacement)
        return self._last_results.reweight(primary_physical=physical,
            secondary_physical=secondary_physical_distributions,
            primary_interactions=interactions, secondary_interactions=secondary_interactions,
            event_factor=event_factor)

    @property
    def injector(self):
        if self._injector is None:
            self._injector = self._build_injector()
        return self._injector

    @property
    def weighter(self):
        if self._weighter is None:
            self._weighter = self._build_weighter(self.injector)
        return self._weighter

    @property
    def detector_model(self):
        return self._detector_model

    @property
    def primary_type(self):
        return self._primary._resolved_particle

    @property
    def injection_distributions(self):
        return list(self._primary.distributions)

    @property
    def physical_distributions(self):
        return list(self._primary.physical)

    @property
    def process_metadata(self):
        return self._process_metadata

    @property
    def phase_spaces(self):
        return self.injector.engine.GetPrimaryProcess().GetPhaseSpaceMap()

    def describe(self):
        print('Simulation: {} events, primary {}, {} secondary declarations'.format(
            self._events, self.primary_type, len(self._secondaries)))

    def __repr__(self):
        return 'Simulation(events={}, primary={!r}, seed={!r})'.format(
            self._events, self.primary_type, self._seed)

    def __getstate__(self):
        from .errors import NotSerializableError
        if self._event_factor is not None:
            raise NotSerializableError('Simulation serialization does not support event_factor',
                                       offenders=['event_factor'])
        # The declaration may own live Python physics: ordinary pickle must
        # preserve it or fail, never replace it with a partial native archive.
        return self.__dict__
