"""Container for simulation results (events + weights).

A Results is an immutable snapshot of one generation run: the event trees,
their weights (stored write-protected), the per-event generation times, and
the injector/weighter that produced them. The run-level counts (attempts,
requested, injected) are captured at construction so a Results never drifts
when the injector is reused for a later run.
"""

import math

import numpy as np

from .errors import ConfigurationError
from ._util import SaveEvents as _SaveEvents


class Results:
    """Weighted event collection returned by :meth:`Simulation.run`.

    Attributes
    ----------
    events : list
        List of ``InteractionTree`` objects (a shallow snapshot copy).
    weights : numpy.ndarray
        Per-event physics weights, write-protected (immutable).
    gen_times : list[float]
        Per-event generation times. :meth:`Simulation.run` fills these with
        the uniform run-average -- the total generation wall-clock time
        divided evenly across the events -- not individually measured
        per-event times.
    attempts : int
        Injection attempts made during the generation run.
    requested : int
        Number of events the run was asked to produce. :meth:`Simulation.run`
        passes this explicitly (see the ``requested`` constructor argument)
        because an injector's live ``number_of_events`` is not a reliable
        post-generate snapshot of the request.
    injected : int
        Number of events actually injected (realized successes).
    """

    def __init__(self, events, weights, gen_times, weighter, injector,
                 *, requested=None):
        self.events = list(events)
        self.gen_times = list(gen_times)
        w = np.array(weights, dtype=float, copy=True)
        if w.ndim != 1 or len(w) != len(self.events) or len(self.gen_times) != len(self.events):
            raise ValueError("Results requires one weight and generation time per event")
        self._weights = np.frombuffer(w.tobytes(), dtype=float)
        self._weighter = None
        self._injector = None
        self._config_key = None
        self._snapshot_error = None
        self._factor_serializable = True
        self._run_ids = frozenset()
        self.attempts = int(getattr(injector, "injection_attempts", 0) or 0)
        self.injected = int(getattr(injector, "injected_events", 0) or 0)
        self.requested = int(getattr(injector, "number_of_events", 0) or 0
                             if requested is None else requested)
        self._capture_configuration(weighter, injector)

    def _capture_configuration(self, weighter, injector):
        from .Weighter import Weighter
        from .Injector import Injector
        import uuid
        if not isinstance(weighter, Weighter) or not isinstance(injector, Injector):
            self._snapshot_error = "a public Injector and Weighter are required"
            return
        sources = weighter.injectors or []
        if len(sources) != 1 or (sources[0] is not injector and sources[0] is not injector.engine):
            raise ConfigurationError("Results requires a weighter for the same single injector")
        if not hasattr(injector, "_results_source_id"):
            injector._results_source_id = uuid.uuid4().hex
        self._run_ids = frozenset([injector._results_source_id])
        self._freeze_weighter(weighter)

    def _freeze_weighter(self, weighter):
        from . import errors
        import copy
        import pickle
        factor = weighter.event_factor
        factor_key = None
        if factor is not None:
            try:
                factor_key = pickle.dumps(factor)
                factor = pickle.loads(factor_key)
            except (TypeError, ValueError, RuntimeError, AttributeError, pickle.PickleError):
                self._factor_serializable = False
        try:
            source = copy.copy(weighter)
            source.event_factor = None
            source.engine
            source._guard_serializable(errors)
            frozen = pickle.loads(pickle.dumps(source))
        except (TypeError, ValueError, RuntimeError, AttributeError, pickle.PickleError) as exc:
            self._snapshot_error = str(exc)
            return
        frozen.event_factor = factor
        self._weighter = frozen
        self._injector = frozen.injectors[0]
        self._set_config_key(factor_key)

    def _set_config_key(self, factor_key=None):
        from . import injection
        eng, phys = self._injector, self._weighter.engine
        expansion = eng.GetSecondaryExpansion()
        self._config_key = (
            injection._detector_state(eng.GetDetectorModel()),
            injection._detector_state(phys.GetDetectorModel()),
            expansion.__getstate__() if expansion is not None else None,
            _process_config(eng.GetPrimaryProcess()),
            {p.primary_type: _process_config(p) for p in eng.GetSecondaryProcesses()},
            _process_config(phys.GetPrimaryPhysicalProcess()),
            {p.primary_type: _process_config(p) for p in phys.GetSecondaryPhysicalProcesses()},
            factor_key,
        )

    def _snapshot_weighter(self):
        if self._weighter is None:
            raise ConfigurationError("Results cannot freeze this live configuration for "
                                     "reweighting/pooling: " + str(self._snapshot_error))
        return self._weighter


    @property
    def weights(self):
        """Per-event weights as a write-protected numpy array."""
        return self._weights.view()


    def __len__(self):
        return len(self.events)

    def __iter__(self):
        """Iterate over (event, weight) pairs."""
        return zip(self.events, self._weights)

    def __getitem__(self, idx):
        """Index or slice into the results.

        A slice returns a new Results over the selected subset that carries the
        same weighter/injector and the parent's run-level counts (attempts /
        requested / injected describe the generation run, not the subset). An
        integer returns an (event, weight) pair.
        """
        if isinstance(idx, slice):
            return self._view(self.events[idx], self._weights[idx], self.gen_times[idx])
        return self.events[idx], self._weights[idx]

    def _view(self, events, weights, times):
        import copy
        child = copy.copy(self)
        child.events = list(events)
        child.gen_times = list(times)
        child._weights = np.frombuffer(np.asarray(weights, dtype=float).tobytes(), dtype=float)
        return child

    @property
    def n_events(self):
        """Number of generated events."""
        return len(self.events)


    def __copy__(self):
        from ._copy import copy_state
        return copy_state(self)

    def __getstate__(self):
        from .errors import NotSerializableError
        if self._weighter is None or self._weighter.event_factor is not None:
            raise NotSerializableError("Results configuration cannot be archived; "
                                       "export stored events/weights with save() instead",
                                       offenders=["live physics or event_factor"])
        state = self.__dict__.copy()
        state.pop("_config_key", None)
        state.pop("_injector", None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._weights = np.frombuffer(np.asarray(self._weights, dtype=float).tobytes(), dtype=float)
        self._injector = self._weighter.injectors[0]
        self._set_config_key()

    @classmethod
    def load(cls, path):
        """Load a trusted Results snapshot written with save(snapshot=True)."""
        import pickle
        with open(str(path) + ".siren_results", "rb") as stream:
            result = pickle.load(stream)
        if not isinstance(result, cls):
            raise TypeError("File does not contain Results")
        return result

    def save(self, path, *, pot=None, snapshot=False, **kwargs):
        """Save the stored weights and run counts to HDF5/Parquet/native files.

        Saving never reevaluates the weighter or ``event_factor``. Weight and
        injector overrides are not accepted; use reweight() to create new Results.

        Parameters
        ----------
        path : str
            Output path (without extension).
        pot : float, optional
            Protons-on-target normalization, recorded as run metadata in the
            saved HDF5 output (the ``pot`` attribute of its ``Events`` group).
        **kwargs
            Additional keyword arguments forwarded to
            :func:`siren._util.SaveEvents`.
        """
        reserved = {"weights", "weighter", "injector", "run_counts", "gen_times", "output_filename"}
        if reserved.intersection(kwargs):
            raise TypeError("Results.save uses stored weights and counts; overrides are not accepted")
        payload = None
        if snapshot:
            import pickle
            payload = pickle.dumps(self)
        headers = [(list(event.header.weights), dict(event.header.provenance))
                   for event in self.events]
        try:
            _SaveEvents(
                self.events, gen_times=self.gen_times, weights=self._weights,
                run_counts={"attempted_events": self.attempts,
                            "accepted_events": self.injected,
                            "events_to_inject": self.requested},
                output_filename=str(path), pot=pot, **kwargs)
        finally:
            for event, (weights, provenance) in zip(self.events, headers):
                event.header.weights = weights
                event.header.provenance = provenance

        if payload is not None:
            with open(str(path) + ".siren_results", "wb") as stream:
                stream.write(payload)

    def reweight(self, *, primary_physical=None, secondary_physical=None,
                 primary_interactions=None, secondary_interactions=None, event_factor=...):
        """Apply a new physical target to the retained generation proposal."""
        from .Weighter import Weighter
        old = self._snapshot_weighter()
        if event_factor is ... and not self._factor_serializable:
            raise ConfigurationError("Provide event_factor explicitly; the original callable "
                                     "cannot be frozen")
        new = Weighter(*old.injectors,
            primary_physical=(old.primary_physical_distributions if primary_physical is None
                              else primary_physical),
            secondary_physical=(old.secondary_physical_distributions if secondary_physical is None
                                else secondary_physical),
            overrides={"detector_model": old.detector_model,
                       "primary_interactions": old.primary_interactions if primary_interactions is None
                                               else primary_interactions,
                       "secondary_interactions": old.secondary_interactions if secondary_interactions is None
                                                 else secondary_interactions},
            event_factor=old.event_factor if event_factor is ... else event_factor)
        child = self._view(self.events, new.weight_all(self.events), self.gen_times)
        child._weighter = child._injector = child._config_key = None
        child._snapshot_error = None
        child._factor_serializable = True
        child._freeze_weighter(new)
        return child


    def bad(self):
        """Indices of pathological weights, grouped by kind.

        Returns a dict with keys ``'inf'``, ``'nan'``, and ``'zero'``, each
        mapping to a sorted list of integer indices. A weight of +/-inf lands
        in ``'inf'``, NaN in ``'nan'``, and exactly-zero in ``'zero'``; the
        three lists are disjoint.
        """
        w = self._weights
        inf_idx = [int(i) for i in np.nonzero(np.isinf(w))[0]]
        nan_idx = [int(i) for i in np.nonzero(np.isnan(w))[0]]
        zero_idx = [int(i) for i in np.nonzero(w == 0.0)[0]]
        return {"inf": inf_idx, "nan": nan_idx, "zero": zero_idx}

    def ess(self):
        """Effective sample size over finite weights.

        ESS = (sum w)^2 / sum(w^2), computed over the finite (non-inf, non-nan)
        weights. Returns 0.0 when there is no finite positive support.
        """
        w = self._weights
        finite = w[np.isfinite(w)]
        scale = float(np.max(np.abs(finite))) if finite.size else 0.0
        if scale == 0.0:
            return 0.0
        finite = finite / scale
        denom = float(np.sum(finite * finite))
        if denom <= 0.0:
            return 0.0
        num = float(np.sum(finite))
        return (num * num) / denom

    def where(self, pred):
        """A new Results over events satisfying ``pred(event, weight)``.

        The filtered Results carries the same weighter/injector and the parent's
        run-level counts (the filter is a view of the same run). ``pred`` is
        called with each (event, weight) pair and kept when it returns truthy.
        """
        keep = [k for k, (e, w) in enumerate(zip(self.events, self._weights))
                if pred(e, w)]
        return self._view([self.events[k] for k in keep], self._weights[keep],
                          [self.gen_times[k] for k in keep])

    def summary(self):
        """Print a multi-line summary of the run and its weight quality."""
        n = len(self._weights)
        if n == 0:
            print("Results: 0 events")
            return

        w = self._weights
        bad = self.bad()
        n_inf = len(bad["inf"])
        n_nan = len(bad["nan"])
        n_zero = len(bad["zero"])
        n_nonzero = n - n_zero - n_inf - n_nan

        finite = w[np.isfinite(w)]
        w_sum = float(np.sum(finite))
        pos = finite[finite > 0]
        w_max = float(np.max(pos)) if pos.size else 0.0
        w_min_nonzero = float(np.min(pos)) if pos.size else 0.0

        ess = self.ess()
        ess_frac = ess / n if n > 0 else 0.0
        eff = (self.injected / self.attempts) if self.attempts > 0 else float("nan")
        eff_str = "{:.1%}".format(eff) if math.isfinite(eff) else "n/a"

        print("Results: {} events".format(n))
        print("  Requested / injected: {} / {}".format(self.requested, self.injected))
        print("  Attempts:             {}".format(self.attempts))
        print("  Efficiency:           {}".format(eff_str))
        print("  ESS:                  {:.1f} ({:.1%} of {})".format(
            ess, ess_frac, n))
        print("  Non-zero weights:     {} / {} ({:.1f}%)".format(
            n_nonzero, n, 100 * n_nonzero / n))
        print("  Weight sum:           {:.6e}".format(w_sum))
        print("  Weight range:         [{:.3e}, {:.3e}]".format(
            w_min_nonzero, w_max))
        print("  Bad weights:          inf={} nan={} zero={}".format(
            n_inf, n_nan, n_zero))


    @classmethod
    def merge(cls, results_list):
        """Pool same-config runs into a single Results.

        Concatenates events/weights/gen_times across runs and rescales each
        run's weights by ``n_i / n_total``, where n_i is that run's attempted
        count and n_total is the sum of attempted counts. This makes
        pooled weights agree with a single run of size n_total. Each run must
        have a positive attempted count.

        Configured models, distributions, detector, phase spaces, expansion,
        and event_factor must match and remain unchanged across runs. Requested
        counts may differ.
        Unknown configurations are rejected. Stored weights are rescaled without
        reevaluating event_factor.
        """
        runs = list(results_list)
        if not runs:
            raise ConfigurationError(
                "Results.merge: empty results list; nothing to pool")

        if any(not r._factor_serializable for r in runs):
            raise ConfigurationError("Results.merge: event_factor cannot be frozen")
        key = runs[0]._config_key
        if key is None or any(r._config_key is None for r in runs):
            raise ConfigurationError("Results.merge: cannot verify an unknown configuration")
        run_ids = set()
        for r in runs:
            if run_ids.intersection(r._run_ids):
                raise ConfigurationError("Results.merge: runs must use independent injectors")
            run_ids.update(r._run_ids)
        seeds = set()
        for r in runs:
            current = {i.GetRandom().get_seed() for i in r._snapshot_weighter().injectors}
            if seeds.intersection(current):
                raise ConfigurationError("Results.merge: runs must use independent seeds")
            seeds.update(current)
        for r in runs[1:]:
            try:
                matches = r._config_key == key
            except TypeError as exc:
                raise ConfigurationError(
                    "Results.merge: cannot compare configurations; model and "
                    "distribution equality must accept the compared types") from exc
            if not matches:
                raise ConfigurationError(
                    "Results.merge: runs have differing config; only "
                    "same-config runs may be pooled")

        if any(r.attempts <= 0 for r in runs):
            raise ConfigurationError("Results.merge: each run needs a positive attempted count")
        n_total = sum(r.attempts for r in runs)
        events = []
        gen_times = []
        weights = []
        for r in runs:
            scale = r.attempts / n_total
            events.extend(r.events)
            gen_times.extend(r.gen_times)
            weights.extend(w * scale for w in r._weights)

        from .Weighter import Weighter
        first = runs[0]._snapshot_weighter()
        injectors = [inj for run in runs for inj in run._snapshot_weighter().injectors]
        weighter = Weighter(*injectors,
            primary_physical=first.primary_physical_distributions,
            secondary_physical=first.secondary_physical_distributions,
            overrides={"detector_model": first.detector_model,
                       "primary_interactions": first.primary_interactions,
                       "secondary_interactions": first.secondary_interactions},
            event_factor=first.event_factor)
        merged = runs[0]._view(events, weights, gen_times)
        merged._weighter = weighter
        merged._config_key = key
        merged._run_ids = frozenset(run_ids)
        merged.attempts = sum(r.attempts for r in runs)
        merged.requested = sum(r.requested for r in runs)
        merged.injected = sum(r.injected for r in runs)
        return merged


def _process_config(process):
    """Compare full process declarations, including the sampled mixture weights."""
    models = process.interactions
    mode = process.GetWeightingMode()
    return (process.primary_type, models.GetDecayChannels(),
            tuple(_model_config(m) for m in
                  list(models.GetCrossSections()) + list(models.GetDecays())),
            tuple(_distribution_config(d) for d in process.distributions),
            (mode.compute_interaction_probability, mode.compute_position_probability,
             int(mode.bound_source), mode.survival_from_creation),
            {sig: _mixture_config(mc)
             for sig, mc in process.GetPhaseSpaceMap().items()})


def _model_config(model):
    # Typed pybind equality methods reject objects from unrelated model bases.
    # Compare model families first, then honor declared value equality, even
    # when two subclasses of the same native base implement it across types.
    from .interactions import CrossSection
    return ("cross_section" if isinstance(model, CrossSection) else "decay", model)


def _distribution_config(distribution):
    from .distributions import PhysicallyNormalizedDistribution
    normalization = None
    if isinstance(distribution, PhysicallyNormalizedDistribution):
        # Native energy-distribution equality compares the shape, not its
        # physical normalization. Snapshot the latter separately.
        normalization = (distribution.IsNormalizationSet(), distribution.normalization)
    return (distribution, normalization)


def _mixture_config(mixture):
    # Accumulators and diagnostic labels do not affect the sampled density.
    return (tuple(_channel_config(c) for c in mixture.channels), tuple(mixture.weights))


def _channel_config(channel):
    from . import injection as inj
    cls = type(channel)
    if cls is inj.NestedMixtureChannel:
        return (cls, _mixture_config(channel.mixture))
    if cls in (inj.PhysicalDecayChannel, inj.PhysicalCrossSectionChannel):
        model = (channel.GetDecay() if cls is inj.PhysicalDecayChannel
                 else channel.GetCrossSection())
        # These adapters store only their model, topology, and measure. Do not
        # serialize the model: it may be a Python authoring class.
        return (cls, _model_config(model), channel.Topology(), channel.Measure())
    if cls in (inj.Isotropic2BodyChannel,):
        try:
            # Native leaf-channel archives contain all proposal parameters
            # (geometry, indices, mappings, etc.), without runtime caches.
            return (cls, "native", channel.__getstate__())
        except (TypeError, RuntimeError):
            # A custom dependency may not be serializable. Require that opaque
            # channel to be shared instead of guessing equality from its name.
            pass
    return (cls, "opaque", channel)
