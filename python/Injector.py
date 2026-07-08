from . import utilities as _utilities
from . import math as _math
from . import dataclasses as _dataclasses
from . import geometry as _geometry
from . import detector as _detector
from . import interactions as _interactions
from . import distributions as _distributions
from . import injection as _injection

import json
import warnings
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
    """Biased-sampling injector.

    Two authoring surfaces compile to one engine configuration:

    * the spec form ``Injector(detector=..., primary=Vertex(...),
      secondaries=(...), events=N)`` builds each process from a Vertex; and
    * the legacy keyword form (``number_of_events``, ``primary_type``,
      ``primary_interactions``, the four ``secondary_*`` dicts,
      ``stopping_condition``) is compiled through the same validation and
      assembly path -- a translation shim, not a parallel code path.

    ``generate(events)`` counts successful trees; ``report()`` renders the
    engine failure ledger; ``engine`` exposes the raw C++ injector.
    """

    def __init__(
        self,
        number_of_events: Optional[int] = None,
        detector_model: Optional[_detector.DetectorModel] = None,
        seed: Optional[int] = None,
        primary_type: Optional[_dataclasses.ParticleType] = None,
        primary_interactions: List[Union[_interactions.CrossSection, _interactions.Decay]] = None,
        primary_injection_distributions: List[_distributions.PrimaryInjectionDistribution] = None,
        primary_phase_spaces: Optional[Dict[_dataclasses.InteractionSignature, _injection.MultiChannelPhaseSpace]] = None,
        primary_weighting_mode=None,
        secondary_interactions: Optional[Dict[_dataclasses.ParticleType, List[Union[_interactions.CrossSection, _interactions.Decay]]]] = None,
        secondary_injection_distributions: Optional[Dict[_dataclasses.ParticleType, List[_distributions.SecondaryInjectionDistribution]]] = None,
        secondary_phase_spaces: Optional[Dict[_dataclasses.ParticleType, Dict[_dataclasses.InteractionSignature, _injection.MultiChannelPhaseSpace]]] = None,
        secondary_weighting_modes: Optional[Dict[_dataclasses.ParticleType, object]] = None,
        stopping_condition: Optional[Callable[[_dataclasses.InteractionTree, _dataclasses.InteractionTreeDatum, int], bool]] = None,
        *,
        detector: Optional[_detector.DetectorModel] = None,
        primary=None,
        secondaries=(),
        events: Optional[int] = None,
        max_attempts: Optional[int] = None,
    ):
        self.__seed = None
        self.__number_of_events = 0
        self.__detector_model = None

        self.__primary_type = None
        self.__primary_interactions = []
        self.__primary_injection_distributions = []
        self.__primary_phase_spaces = {}

        self.__secondary_interactions = {}
        self.__secondary_injection_distributions = {}
        self.__secondary_phase_spaces = {}
        self.__primary_weighting_mode = None
        self.__secondary_weighting_modes = {}
        self.__stopping_condition = None

        self.__injector = None

        # Spec-form fields (Vertex-based authoring).
        self.__primary_vertex = None
        self.__secondary_vertices = []
        self.__primary_compiled_process = None
        self.__secondary_compiled_processes = {}
        self.__max_attempts = max_attempts

        # Strong references to Python models/distributions/mixtures so their
        # pybind trampolines survive gc while the compiled processes hold them.
        self.__keepalive = []

        detector_model = detector_model if detector_model is not None else detector
        if events is not None:
            number_of_events = events

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
        if primary_phase_spaces is not None:
            self.__primary_phase_spaces = primary_phase_spaces
        if secondary_interactions is not None:
            self.__secondary_interactions = secondary_interactions
        if secondary_injection_distributions is not None:
            self.__secondary_injection_distributions = secondary_injection_distributions
        if secondary_phase_spaces is not None:
            self.__secondary_phase_spaces = secondary_phase_spaces
        if primary_weighting_mode is not None:
            self.__primary_weighting_mode = primary_weighting_mode
        if secondary_weighting_modes is not None:
            self.__secondary_weighting_modes = secondary_weighting_modes
        if stopping_condition is not None:
            self.__stopping_condition = stopping_condition

        if primary is not None:
            self.__primary_vertex = primary
            self.__secondary_vertices = list(secondaries)

    # ------------------------------------------------------------------ #
    #  Build (single validation choke point)                              #
    # ------------------------------------------------------------------ #

    def _build(self):
        """Validate the configuration and construct the engine injector.

        Every path -- spec Vertices or legacy kwargs -- lands here. It resolves
        the RNG streams, compiles processes, checks the injection
        distributions, expand wiring, per-signature name resolution and mixture
        measures, runs the detector-dependent channel-density probe, and builds
        the C++ injector.
        """
        from . import _validation
        from . import errors as _errors
        from . import expand as _expand

        if self.__seed is None:
            random = _utilities.SIREN_random()
            self.__seed = random.get_seed()
        else:
            random = _utilities.SIREN_random(self.__seed)
        # An independent stream for config-time validation probes, so a probe
        # draw never perturbs the generation stream. The seed is masked to the
        # 32-bit range the RNG constructor accepts.
        validation_random = _utilities.SIREN_random(int(self.__seed) & 0x7FFFFFFF)

        if self.__number_of_events is None:
            raise ValueError("number_of_events must be provided")
        elif self.__number_of_events <= 0:
            raise ValueError("number_of_events must be positive")

        if self.__detector_model is None:
            raise ValueError("detector_model must be provided")

        if self.__primary_vertex is not None:
            self._compile_from_vertices(_validation, _expand)

        if self.__primary_type is None:
            raise ValueError("primary_type must be provided")
        if len(self.__primary_interactions) == 0:
            raise ValueError("primary_interactions must be provided")
        if len(self.__primary_injection_distributions) == 0:
            raise ValueError("primary_injection_distributions must be provided")

        _validation.validate_injection_distributions(
            self.__primary_injection_distributions)

        # Cross-key validation of every secondary-keyed dict (a typo'd key in
        # any of the four dicts raises here rather than silently vanishing).
        _validation.validate_secondary_keys(
            ("secondary_interactions", self.__secondary_interactions),
            ("secondary_injection_distributions", self.__secondary_injection_distributions),
            ("secondary_phase_spaces", self.__secondary_phase_spaces),
            ("secondary_weighting_modes", self.__secondary_weighting_modes),
        )
        if sorted(self.__secondary_interactions.keys()) != sorted(
                self.__secondary_injection_distributions.keys()):
            raise ValueError(
                "secondary_interactions and secondary_injection_distributions "
                "must have the same keys")

        # A chain with secondary processes but no expansion control -- neither
        # a legacy stopping_condition nor Vertex expand/continue_if wiring (the
        # latter is already checked above, in _compile_from_vertices, via
        # validate_expansion_wiring) -- would silently prune every secondary:
        # the engine's default stopping condition accepts primary-only trees
        # unconditionally. This check runs regardless of how the secondaries
        # were assembled, so it also covers the legacy secondary_interactions/
        # secondary_injection_distributions dict form, which never reaches
        # _compile_from_vertices.
        if self.__secondary_interactions and self.__stopping_condition is None:
            raise _errors.ConfigurationError(
                "secondary processes are configured ({}) but no stopping "
                "condition governs their expansion; the engine's default "
                "stopping condition would silently prune every secondary. "
                "Provide stopping_condition=..., or (spec form) wire "
                "Vertex(expand=[...]).".format(
                    sorted(str(k) for k in self.__secondary_interactions)))

        primary_process = self._assemble_primary_process()
        secondary_processes = self._assemble_secondary_processes()

        self._probe_phase_spaces(
            _validation, primary_process, secondary_processes, validation_random)

        self.__injector = _Injector(
            self.__number_of_events,
            self.__detector_model,
            primary_process,
            secondary_processes,
            random,
        )

        if self.__stopping_condition is not None:
            self.__injector.SetStoppingCondition(self.__stopping_condition)

    def _compile_from_vertices(self, _validation, _expand):
        """Translate spec Vertices into the legacy field storage.

        The compiled processes and their distributions/phase spaces are read
        back into the same private fields the legacy path populates, so a
        single assembly path serves both surfaces.
        """
        primary = self.__primary_vertex
        all_vertices = [primary] + list(self.__secondary_vertices)
        has_expand = any(v.expand or v.continue_if is not None
                         for v in all_vertices)
        _validation.check_expand_vs_legacy_stopping(
            self.__stopping_condition is not None, has_expand)

        self.__keepalive.append(primary)
        self.__primary_type = primary._resolved_particle
        self.__primary_interactions = list(primary.interactions)
        self.__primary_injection_distributions = list(primary.distributions)

        specs = [primary.as_vertex_spec()]
        for sv in self.__secondary_vertices:
            self.__keepalive.append(sv)
            ptype = sv._resolved_particle
            # Two secondary vertices resolving to the same particle type would
            # each overwrite the other's interactions/distributions/process in
            # the by-type maps below; require a single vertex per type.
            if ptype in self.__secondary_interactions:
                from . import errors as _errors
                raise _errors.ConfigurationError(
                    "two secondary vertices resolve to the same particle type "
                    "{!r}; a chain admits one vertex per secondary type".format(
                        str(ptype)))
            self.__secondary_interactions[ptype] = list(sv.interactions)
            self.__secondary_injection_distributions[ptype] = list(sv.distributions)
            specs.append(sv.as_vertex_spec())

        # Validate the expansion wiring whenever the chain carries secondaries
        # and no legacy stopping_condition governs it. This is the only place
        # the "secondaries declared but no expand rule anywhere" misconfiguration
        # is caught: without it the engine's default stopping condition prunes
        # every secondary, silently collapsing the chain to primary-only trees.
        # A legacy stopping_condition takes control of expansion, so it passes
        # through untouched. Compiling the expansion callable stays gated on
        # has_expand; a legacy stopping_condition alone is left as supplied.
        if self.__secondary_vertices and self.__stopping_condition is None:
            _validation.validate_expansion_wiring(specs)
        if self.__secondary_vertices and has_expand:
            self.__stopping_condition = _expand.compile_expansion(specs)

        # Compile per-signature phase spaces from each Vertex's kinematics.
        primary_process = primary.compile(
            is_primary=True, detector=self.__detector_model)
        self.__primary_compiled_process = primary_process
        self.__secondary_compiled_processes = {}
        for sv in self.__secondary_vertices:
            proc = sv.compile(is_primary=False, detector=self.__detector_model)
            self.__secondary_compiled_processes[sv._resolved_particle] = proc

    def _assemble_primary_process(self):
        if self.__primary_compiled_process is not None:
            return self.__primary_compiled_process
        primary_interaction_collection = _interactions.InteractionCollection(
            self.__primary_type, self.__primary_interactions)
        primary_process = _injection.PrimaryInjectionProcess(
            self.__primary_type, primary_interaction_collection)
        primary_process.distributions = self.__primary_injection_distributions
        for sig, ps in self.__primary_phase_spaces.items():
            primary_process.SetPhaseSpace(sig, ps)
        if self.__primary_weighting_mode is not None:
            primary_process.weighting_mode = self.__primary_weighting_mode
        return primary_process

    def _assemble_secondary_processes(self):
        if self.__secondary_compiled_processes:
            return list(self.__secondary_compiled_processes.values())
        secondary_processes = []
        for secondary_type, sec_ints in self.__secondary_interactions.items():
            secondary_interaction_collection = _interactions.InteractionCollection(
                secondary_type, sec_ints)
            secondary_process = SecondaryInjectionProcess(
                secondary_type, secondary_interaction_collection)
            secondary_process.distributions = (
                self.__secondary_injection_distributions[secondary_type])
            if secondary_type in self.__secondary_phase_spaces:
                for sig, ps in self.__secondary_phase_spaces[secondary_type].items():
                    secondary_process.SetPhaseSpace(sig, ps)
            if secondary_type in self.__secondary_weighting_modes:
                secondary_process.weighting_mode = (
                    self.__secondary_weighting_modes[secondary_type])
            secondary_processes.append(secondary_process)
        return secondary_processes

    def _probe_phase_spaces(self, _validation, primary_process,
                            secondary_processes, validation_random):
        """Run the detector-dependent ValidateChannelDensities probe.

        Only mixtures actually registered on a process are probed. A typed
        SIREN error (a measure incompatibility, a configuration error, a weight
        calculation error, and every other error in ``errors.py``) is a real
        fault and propagates; a probe that simply cannot run on the synthetic
        template (unsolvable kinematics, missing target mass) surfaces as a bare
        ValueError/TypeError/RuntimeError and is skipped -- the config-time
        Mixture.validate() already screened measures without a detector.
        """
        typed_errors = _typed_siren_errors()
        processes = [primary_process] + list(secondary_processes)
        for process in processes:
            ps_map = process.GetPhaseSpaceMap()
            for sig, mcps in ps_map.items():
                try:
                    _validation.probe_channel_densities(
                        mcps, sig, self.__detector_model, validation_random)
                except typed_errors:
                    raise
                except (ValueError, TypeError, RuntimeError):
                    continue

    def __initialize_injector(self):
        self._build()

    # ------------------------------------------------------------------ #
    #  Generation                                                          #
    # ------------------------------------------------------------------ #

    def generate(self, events=None, *, on_shortfall="warn", progress=None,
                 min_efficiency=None):
        """Generate `events` successful trees.

        Retries failed ``GenerateEvent`` calls (never yields an empty tree) up
        to ``max_attempts`` (default ``events * 1000``). Counts SUCCESSES, not
        attempts. ``min_efficiency`` aborts early if the success rate falls
        below it. On shortfall, ``on_shortfall`` chooses ``'warn'`` (default),
        ``'raise'``, or ``'ignore'``, carrying an InjectionReport.
        """
        from . import errors as _errors

        if self.__injector is None:
            self._build()
        if events is None:
            events = self.__number_of_events
        if events <= 0:
            raise ValueError("events must be positive")

        max_attempts = self.__max_attempts
        if max_attempts is None:
            max_attempts = events * 1000

        # The engine throws once its own attempt quota (events_to_inject) is
        # reached, so raise that quota to the retry budget; otherwise a
        # sub-unity efficiency could never reach `events` successes.
        self.__injector.ResetInjectedEvents(max_attempts)

        # Give the efficiency estimate a minimum sample before it can abort, so
        # an unlucky early run is not cut short.
        min_efficiency_warmup = 50

        trees = []
        attempts = 0
        while len(trees) < events and attempts < max_attempts:
            attempts += 1
            try:
                tree = self.__injector.GenerateEvent()
            except RuntimeError as err:
                # The engine raises once its own attempt quota is exhausted;
                # any other error is a real fault and must not be swallowed.
                if "maximum number of injection attempts" in str(err):
                    break
                raise
            if len(tree.tree) == 0:
                continue
            trees.append(tree)
            if progress is not None:
                progress(len(trees), events)
            if (min_efficiency is not None and attempts >= min_efficiency_warmup
                    and (len(trees) / attempts) < min_efficiency):
                report = self.report()
                message = (
                    "injection efficiency {:.2%} fell below min_efficiency "
                    "{:.2%} after {} attempts ({} succeeded)".format(
                        len(trees) / attempts, min_efficiency, attempts,
                        len(trees)))
                self._handle_shortfall(on_shortfall, message, report, _errors)
                return trees

        if len(trees) < events:
            report = self.report()
            message = (
                "requested {} events but only {} succeeded in {} attempts"
                .format(events, len(trees), attempts))
            self._handle_shortfall(on_shortfall, message, report, _errors)
        return trees

    def _handle_shortfall(self, on_shortfall, message, report, _errors):
        if on_shortfall == "ignore":
            return
        if on_shortfall == "raise":
            raise _errors.InjectionShortfall(message, report=report)
        if on_shortfall == "warn":
            warnings.warn(_errors.InjectionShortfall(message, report=report))
            return
        raise ValueError(
            "on_shortfall must be 'warn', 'raise', or 'ignore', got {!r}"
            .format(on_shortfall))

    def reset(self, events=None):
        """Reset the injected/attempt counters and failure ledger.

        With ``events`` given, also change the target count.
        """
        if self.__injector is None:
            if events is not None:
                self.__number_of_events = events
            return
        if events is not None:
            self.__injector.ResetInjectedEvents(events)
            self.__number_of_events = events
        else:
            self.__injector.ResetInjectedEvents()

    def report(self):
        """Render the engine failure ledger as an InjectionReport."""
        from .report import InjectionReport
        if self.__injector is None:
            return InjectionReport(0, 0, [], None)
        ledger = self.__injector.GetFailureLedger()
        attempts = self.__injector.InjectionAttempts()
        successes = self.__injector.InjectedEvents()
        last_tree = self.__injector.GetLastFailedTree()
        return InjectionReport.from_ledger(
            ledger, attempts=attempts, successes=successes,
            last_failed_tree=last_tree)

    # ------------------------------------------------------------------ #
    #  Optimizer tuning I/O                                                #
    # ------------------------------------------------------------------ #

    def export_tuning(self, path=None):
        """Return (and optionally write) each mixture's channel weights.

        The result maps a stable mixture index to its weight list, so a tuned
        run's channel weights can be reloaded via apply_tuning(). Reloading them
        retunes the target run's generation density, so the same ordering rule
        applies there as for apply_tuning(): weight only events generated after
        the weights are applied.
        """
        if self.__injector is None:
            self._build()
        tuning = {}
        for i, mcps in enumerate(self.__injector.GetPhaseSpaces()):
            tuning[str(i)] = list(mcps.weights)
        if path is not None:
            with open(path, "w") as f:
                json.dump(tuning, f)
        return tuning

    def apply_tuning(self, dict_or_path):
        """Set each mixture's channel weights from a dict or a JSON path.

        Retuning the channel weights changes the generation density g. Events
        generated before this call carry no snapshot of the density that
        produced them, so weighting them afterwards silently applies the new g
        where the old one generated, biasing their weights. Always regenerate
        after retuning and before weighting; the shipped drivers do (Simulation
        .run tunes, resets, then generates before it weights).
        """
        if self.__injector is None:
            self._build()
        if isinstance(dict_or_path, str):
            with open(dict_or_path) as f:
                tuning = json.load(f)
        else:
            tuning = dict_or_path
        mixtures = self.__injector.GetPhaseSpaces()
        for key, weights in tuning.items():
            idx = int(key)
            if idx < len(mixtures):
                mixtures[idx].weights = list(weights)
                mixtures[idx].Normalize()

    # ------------------------------------------------------------------ #
    #  Introspection                                                       #
    # ------------------------------------------------------------------ #

    @property
    def engine(self):
        """The raw C++ injector, building it on first access."""
        if self.__injector is None:
            self._build()
        return self.__injector

    def describe(self):
        """A short human-readable description of the configuration."""
        lines = ["Injector(events={}, primary={})".format(
            self.__number_of_events, self.__primary_type)]
        if self.__secondary_interactions:
            lines.append("  secondaries: {}".format(
                sorted(str(k) for k in self.__secondary_interactions)))
        return "\n".join(lines)

    def to_dict(self):
        """A plain-dict snapshot of the configuration counts."""
        return {
            "events": self.__number_of_events,
            "seed": self.__seed,
            "primary_type": str(self.__primary_type),
            "secondary_types": [str(k) for k in self.__secondary_interactions],
        }

    def save(self, filename):
        """Serialize the injector, refusing configurations that cannot survive.

        The engine archives processes, including weighting modes and phase-space
        maps, but not a stopping condition, and it cannot serialize Python
        trampoline-derived interactions/distributions. Unsupported Python state
        raises NotSerializableError instead. Pickle preserves the stopping
        condition, so pickle a chain that needs one.
        """
        from . import errors as _errors
        if self.__injector is None:
            self._build()
        self._guard_serializable(_errors)
        self.__injector.SaveInjector(filename)

    def _guard_serializable(self, _errors, for_pickle=False):
        offenders = []
        # The stopping condition rides along in the pickle state tuple, so
        # pickle preserves it; the standalone C++ archive cannot.
        if self.__stopping_condition is not None and not for_pickle:
            offenders.append(
                "a stopping condition is set (not archived; use pickle, which "
                "preserves it)")
        for interaction in self.__primary_interactions:
            if _is_trampoline(interaction):
                offenders.append(
                    "primary interaction {!r} is a Python subclass (not "
                    "serializable)".format(type(interaction).__name__))
        for dist in self.__primary_injection_distributions:
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
        for stype, dists in self.__secondary_injection_distributions.items():
            for dist in dists:
                if _is_trampoline(dist):
                    offenders.append(
                        "secondary distribution {!r} for type {} is a Python "
                        "subclass (not serializable)".format(
                            type(dist).__name__, str(stype)))
        if offenders:
            raise _errors.NotSerializableError(
                "this injector cannot be saved without silently changing "
                "physics on reload:\n  - " + "\n  - ".join(offenders),
                offenders=offenders)

    @classmethod
    def load(cls, filename):
        """Construct an Injector from a saved archive.

        The archive omits the random engine, so the reloaded injector is given
        a FRESH, unseeded generation stream. RNG streams do not resume across
        save/load: a saved seed is not restored, and generation begins anew.
        (save() drops any stopping condition; pickle a chain that needs one.)
        """
        obj = cls()
        # Construct the engine through the archive-loading constructor rather
        # than loading into a __new__-allocated shell: __new__ allocates the
        # C++ storage without running any constructor, so the shell has no
        # vtable and the first virtual call (GetDetectorModel below) crashes.
        # The constructor also takes the random engine the archive omits, so
        # generate() does not dereference a null RNG; the event count is a
        # placeholder LoadInjector immediately overwrites from the archive.
        # load() starts a fresh unseeded stream and records the drawn seed so
        # obj.seed reflects the live stream.
        _random = _utilities.SIREN_random()
        obj.__injector = _Injector(1, filename, _random)
        obj.__seed = _random.get_seed()
        obj.__number_of_events = obj.__injector.EventsToInject()
        obj.__detector_model = obj.__injector.GetDetectorModel()
        primary_process = obj.__injector.GetPrimaryProcess()
        obj.__primary_type = primary_process.primary_type
        obj.__primary_interactions = list(
            primary_process.interactions.GetCrossSections()) + list(
            primary_process.interactions.GetDecays())
        obj.__primary_injection_distributions = list(primary_process.distributions)
        obj.__primary_phase_spaces = primary_process.GetPhaseSpaceMap()
        obj.__secondary_interactions = {}
        obj.__secondary_injection_distributions = {}
        obj.__secondary_phase_spaces = {}
        for stype, sproc in obj.__injector.GetSecondaryProcessMap().items():
            obj.__secondary_interactions[stype] = list(
                sproc.interactions.GetCrossSections()) + list(
                sproc.interactions.GetDecays())
            obj.__secondary_injection_distributions[stype] = list(sproc.distributions)
            obj.__secondary_phase_spaces[stype] = sproc.GetPhaseSpaceMap()
        return obj

    # ------------------------------------------------------------------ #
    #  Pickling                                                            #
    # ------------------------------------------------------------------ #

    def __getstate__(self):
        from . import errors as _errors
        # Build first: an unbuilt injector has no engine to pickle, and building
        # resolves the generation seed so the state tuple carries a concrete
        # seed (a fresh injector re-seeded from it reproduces this stream from
        # the start). The guard then refuses exactly what save() refuses, minus
        # the stopping condition -- the state tuple preserves that.
        if self.__injector is None:
            self._build()
        self._guard_serializable(_errors, for_pickle=True)
        return (self.__seed, self.__stopping_condition,
                self.__injector.__getstate__())

    def __setstate__(self, state):
        self.__seed, self.__stopping_condition, injector_state = state
        self.__injector = _Injector.__new__(_Injector)
        if self.__injector is None:
            raise TypeError("Failed to create C++ Injector object")
        self.__injector.__setstate__(injector_state)
        self.__number_of_events = self.__injector.EventsToInject()
        self.__detector_model = self.__injector.GetDetectorModel()
        primary_process = self.__injector.GetPrimaryProcess()
        self.__primary_type = primary_process.primary_type
        self.__primary_interactions = list(primary_process.interactions.GetCrossSections()) + list(primary_process.interactions.GetDecays())
        self.__primary_injection_distributions = list(primary_process.distributions)
        self.__primary_phase_spaces = primary_process.GetPhaseSpaceMap()
        self.__secondary_interactions = {}
        self.__secondary_injection_distributions = {}
        self.__secondary_phase_spaces = {}
        self.__primary_vertex = None
        self.__secondary_vertices = []
        self.__primary_compiled_process = None
        self.__secondary_compiled_processes = {}
        self.__keepalive = []
        self.__max_attempts = None
        for secondary_type, secondary_process in self.__injector.GetSecondaryProcessMap().items():
            self.__secondary_interactions[secondary_type] = list(secondary_process.interactions.GetCrossSections()) + list(secondary_process.interactions.GetDecays())
            self.__secondary_injection_distributions[secondary_type] = list(secondary_process.distributions)
            self.__secondary_phase_spaces[secondary_type] = secondary_process.GetPhaseSpaceMap()
        # The C++ archive omits the random engine, so install a fresh one:
        # generate() skips _build once the engine exists and would otherwise
        # dereference a null RNG. RNG streams do NOT resume across a pickle
        # round-trip -- the stored seed re-seeds a fresh stream, so a fresh
        # injector and this unpickled one produce the same stream from the
        # start (the seed was resolved during __getstate__'s build).
        if self.__seed is not None:
            self.__injector.SetRandom(_utilities.SIREN_random(self.__seed))
        else:
            self.__injector.SetRandom(_utilities.SIREN_random())
        # The stopping condition rode along in the state tuple but must be
        # re-attached to the engine; otherwise the engine default prunes every
        # secondary and silently truncates chains.
        if self.__stopping_condition is not None:
            self.__injector.SetStoppingCondition(self.__stopping_condition)

    # ------------------------------------------------------------------ #
    #  Properties (legacy surface, retained)                              #
    # ------------------------------------------------------------------ #

    @property
    def seed(self):
        return self.__seed

    @seed.setter
    def seed(self, seed):
        self.__seed = seed
        if self.__injector is not None:
            # The engine exposes only SetRandom; install a fresh engine
            # seeded from the new value (streams restart, not resume).
            self.__injector.SetRandom(_utilities.SIREN_random(seed))

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
    def primary_phase_spaces(self):
        return self.__primary_phase_spaces

    @primary_phase_spaces.setter
    def primary_phase_spaces(self, phase_spaces):
        if self.__injector is not None:
            primary_process = self.__injector.GetPrimaryProcess()
            for sig, ps in phase_spaces.items():
                primary_process.SetPhaseSpace(sig, ps)
        self.__primary_phase_spaces = phase_spaces

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
            secondary_processes = self.__injector.GetSecondaryProcessMap()
            current_secondary_types = sorted(list(secondary_processes.keys()))
            new_secondary_types = sorted(list(secondary_injection_distributions.keys()))
            if current_secondary_types != new_secondary_types:
                raise ValueError("Cannot change the secondary types after initialization")
            for secondary_type, secondary_process in secondary_processes.items():
                secondary_process.distributions = secondary_injection_distributions[secondary_type]
        self.__secondary_injection_distributions = secondary_injection_distributions

    @property
    def secondary_phase_spaces(self):
        return self.__secondary_phase_spaces

    @secondary_phase_spaces.setter
    def secondary_phase_spaces(self, phase_spaces):
        if self.__injector is not None:
            secondary_processes = self.__injector.GetSecondaryProcessMap()
            for secondary_type, phase_space in phase_spaces.items():
                if secondary_type not in secondary_processes:
                    raise ValueError("Cannot set a phase space for an unknown secondary type")
                for sig, ps in phase_space.items():
                    secondary_processes[secondary_type].SetPhaseSpace(sig, ps)
        self.__secondary_phase_spaces = phase_spaces

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
            self._build()
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

    @property
    def injection_attempts(self):
        if self.__injector is not None:
            return self.__injector.InjectionAttempts()
        return 0

    @wraps(_Injector.ResetInjectedEvents)
    def reset_injected_events(self):
        if self.__injector is not None:
            self.__injector.ResetInjectedEvents()
    reset_injected_events.__name__ = "reset_injected_events"

    def __iter__(self):
        """Iterate raw generation attempts (yields empty trees on failure).

        Deprecated: kept for backward compatibility. Prefer generate(), which
        counts successes and never yields an empty tree.
        """
        warnings.warn(
            "iterating an Injector yields raw attempts including empty trees; "
            "use generate(events) which counts successes",
            DeprecationWarning, stacklevel=2)
        if self.__injector is None:
            self._build()
        for _ in range(self.number_of_events):
            yield self.__injector.GenerateEvent()

    def __len__(self):
        return self.number_of_events


def _typed_siren_errors():
    """The tuple of typed SIREN exception classes exported by ``errors.py``.

    Every name in ``errors.py`` derives from ``RuntimeError``; a probe raising
    any of them signals a real fault that must propagate rather than be treated
    as probe-inapplicability.
    """
    from . import errors as _errors
    classes = []
    for name in dir(_errors):
        obj = getattr(_errors, name)
        if isinstance(obj, type) and issubclass(obj, BaseException) \
                and issubclass(obj, RuntimeError):
            classes.append(obj)
    return tuple(classes)


def _is_trampoline(obj):
    """Whether `obj` is a Python subclass of a pybind base (a trampoline).

    A C++-native pybind object's type is defined in an extension module; a
    Python subclass is defined in a .py module, so its class __module__ does
    not start with 'siren.'.
    """
    module = type(obj).__module__ or ""
    return not module.startswith("siren.") and module not in ("builtins",)
