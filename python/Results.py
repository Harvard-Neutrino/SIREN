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


class VertexVarianceShare:
    """Per-channel variance shares for one vertex mixture.

    W is the per-channel Kleiss-Pittau statistic W_i = mean(w^2 g_i/g).
    shares is W normalized to sum to 1 over the mixture's channels; it is the
    fraction of the vertex's weight-variance attributable to each channel.
    labels names each channel. degenerate is True when no finite-weight event
    reached this mixture (shares is then a uniform placeholder).
    """

    __slots__ = ("index", "labels", "W", "shares", "count", "degenerate")

    def __init__(self, index, labels, W, shares, count, degenerate):
        self.index = index
        self.labels = list(labels)
        self.W = list(W)
        self.shares = list(shares)
        self.count = count
        self.degenerate = degenerate

    def __repr__(self):
        return ("VertexVarianceShare(index={} channels={} degenerate={})"
                .format(self.index, len(self.shares), self.degenerate))


class VarianceReport:
    """Per-vertex decomposition of weight variance across channels.

    vertices is a list of VertexVarianceShare, one per multi-channel mixture in
    the injector. Each entry's shares sum to 1 (a within-vertex fraction, not a
    global one). An injector with no mixtures yields vertices == [].
    """

    __slots__ = ("vertices",)

    def __init__(self, vertices):
        self.vertices = list(vertices)

    def __len__(self):
        return len(self.vertices)

    def __iter__(self):
        return iter(self.vertices)

    def __str__(self):
        if not self.vertices:
            return "VarianceReport: no multi-channel mixtures"
        lines = ["VarianceReport: {} vertex mixture(s)".format(len(self.vertices))]
        for v in self.vertices:
            tag = " (degenerate)" if v.degenerate else ""
            lines.append("  vertex {}{}: {} channel(s)".format(
                v.index, tag, len(v.shares)))
            for label, share in zip(v.labels, v.shares):
                lines.append("    {:<40} {:6.1%}".format(label, share))
        return "\n".join(lines)


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
        w = np.asarray(weights, dtype=float)
        w.setflags(write=False)
        self._weights = w
        self.gen_times = list(gen_times)
        self._weighter = weighter
        self._injector = injector
        # Snapshot the run-level counts so this Results reflects the run at
        # creation even if the injector is reused afterwards. A raw engine (or
        # None) may lack these attributes; fall back to 0.
        self.attempts = int(getattr(injector, "injection_attempts", 0) or 0)
        self.injected = int(getattr(injector, "injected_events", 0) or 0)
        # `requested` is the count the run was asked to produce. Prefer an
        # explicit override: an injector's live number_of_events is not a
        # reliable snapshot of the request, because siren.Injector.generate
        # raises the engine attempt quota to events * 1000, so a post-generate
        # read reports the retry budget rather than the request. When no
        # override is given, fall back to the injector attribute for callers
        # that keep number_of_events stable.
        if requested is None:
            self.requested = int(getattr(injector, "number_of_events", 0) or 0)
        else:
            self.requested = int(requested)

    # ------------------------------------------------------------------ #
    #  Immutable weights                                                   #
    # ------------------------------------------------------------------ #

    @property
    def weights(self):
        """Per-event weights as a write-protected numpy array."""
        return self._weights

    # ------------------------------------------------------------------ #
    #  Sequence protocol                                                   #
    # ------------------------------------------------------------------ #

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
            child = Results(
                self.events[idx],
                self._weights[idx],
                self.gen_times[idx],
                self._weighter,
                self._injector,
            )
            child._inherit_counts(self)
            return child
        return self.events[idx], self._weights[idx]

    def _inherit_counts(self, parent):
        """Copy run-level counts from a parent Results (view semantics)."""
        self.attempts = parent.attempts
        self.requested = parent.requested
        self.injected = parent.injected

    @property
    def n_events(self):
        """Number of generated events."""
        return len(self.events)

    # ------------------------------------------------------------------ #
    #  Persistence                                                         #
    # ------------------------------------------------------------------ #

    def save(self, path, *, pot=None, **kwargs):
        """Save events and weights to HDF5/Parquet/native files.

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
        _SaveEvents(
            self.events,
            self._weighter,
            self.gen_times,
            output_filename=path,
            pot=pot,
            **kwargs,
        )

    # ------------------------------------------------------------------ #
    #  Weight diagnostics                                                  #
    # ------------------------------------------------------------------ #

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
        denom = float(np.sum(finite * finite))
        if denom <= 0.0:
            return 0.0
        num = float(np.sum(finite))
        return (num * num) / denom

    def explain(self, i):
        """Per-vertex weight breakdown for event ``i``.

        Returns a :class:`report.WeightBreakdown` from the retained weighter,
        decomposing the event's weight into per-vertex generation and physical
        factors. Useful for diagnosing inf/0/NaN weights.
        """
        return self._weighter.explain(self.events[i])

    def where(self, pred):
        """A new Results over events satisfying ``pred(event, weight)``.

        The filtered Results carries the same weighter/injector and the parent's
        run-level counts (the filter is a view of the same run). ``pred`` is
        called with each (event, weight) pair and kept when it returns truthy.
        """
        keep = [k for k, (e, w) in enumerate(zip(self.events, self._weights))
                if pred(e, w)]
        child = Results(
            [self.events[k] for k in keep],
            self._weights[keep],
            [self.gen_times[k] for k in keep],
            self._weighter,
            self._injector,
        )
        child._inherit_counts(self)
        return child

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

    # ------------------------------------------------------------------ #
    #  Variance decomposition                                             #
    # ------------------------------------------------------------------ #

    def variance_report(self):
        """Per-vertex, per-channel weight-variance decomposition.

        Feeds the snapshot events through the injector's multi-channel mixtures
        to accumulate the bare Kleiss-Pittau statistic W_i = mean(w^2 g_i/g) per
        channel, then reports each channel's share of its vertex's variance
        (W_i / sum_j W_j, summing to 1 per vertex). Reading the accumulators
        does not call UpdateWeights, so the injector's tuning is untouched.

        Returns a :class:`VarianceReport`; an injector with no mixtures yields
        an empty report (vertices == []).
        """
        eng = getattr(self._injector, "engine", self._injector)
        if eng is None or not hasattr(eng, "GetPhaseSpaces"):
            return VarianceReport([])

        mixtures = eng.GetPhaseSpaces()
        if not mixtures:
            return VarianceReport([])

        # Start clean, credit each finite-weight event to the mixtures it
        # touched, then read the per-channel statistic. discount_fallback and
        # recurse credit grouped/nested channels against the same outer g.
        for mc in mixtures:
            mc.ResetAccumulators(True)
        for event, weight in zip(self.events, self._weights):
            w = float(weight)
            if not math.isfinite(w):
                continue
            eng.AccumulateEventToMixtures(event, w, True, True)

        vertices = []
        for k, mc in enumerate(mixtures):
            acc = list(mc.kp_accumulator)
            cnt = int(mc.kp_count)
            labels = _channel_labels(mc)
            if cnt <= 0 or sum(acc) <= 0.0:
                # No finite-weight event reached this mixture; report uniform
                # shares over its channels so the per-vertex shares still sum
                # to 1. The accumulator may be unsized (empty) until the first
                # Accumulate, so size from the channel list.
                m = len(acc) if acc else len(mc.channels)
                uniform = [1.0 / m] * m if m > 0 else []
                vertices.append(VertexVarianceShare(
                    k, labels, [0.0] * m, uniform, cnt, degenerate=True))
                continue
            W = [a / cnt for a in acc]
            total = sum(W)
            shares = [wi / total for wi in W]
            vertices.append(VertexVarianceShare(
                k, labels, W, shares, cnt, degenerate=False))

        # Leave no residue in the accumulators (idempotent for callers that
        # optimize the injector afterwards).
        for mc in mixtures:
            mc.ResetAccumulators(True)

        return VarianceReport(vertices)

    # ------------------------------------------------------------------ #
    #  Pooling                                                             #
    # ------------------------------------------------------------------ #

    @property
    def _config_key(self):
        """Best-effort hashable signature of the generation configuration.

        Used by :meth:`merge` to reject pooling across differing configs. Keys
        on the requested count and, when cheaply available, the injector's
        primary type; permissive by design so same-config runs always match
        (a fresh DetectorModel object per run rules out identity comparison).
        """
        primary = None
        inj = self._injector
        if inj is not None:
            try:
                primary = str(getattr(inj, "primary_type", None))
            except Exception:
                primary = None
        return (self.requested, primary)

    @classmethod
    def merge(cls, results_list):
        """Pool same-config runs into a single Results.

        Concatenates events/weights/gen_times across runs and rescales each
        run's weights by ``n_i / n_total``, where n_i is that run's injected
        (realized) count and n_total is the sum of injected counts. This makes
        the pooled per-event mean equal a single run of size n_total.

        All runs must share a configuration (see :attr:`_config_key`); a
        mismatch raises :class:`ConfigurationError`. The merged Results carries
        the first run's weighter/injector so explain/save still work.
        """
        runs = list(results_list)
        if not runs:
            raise ConfigurationError(
                "Results.merge: empty results list; nothing to pool")

        key = runs[0]._config_key
        for r in runs[1:]:
            if r._config_key != key:
                raise ConfigurationError(
                    "Results.merge: runs have differing config; only "
                    "same-config runs may be pooled")

        n_total = sum(r.injected for r in runs)
        if n_total <= 0:
            raise ConfigurationError(
                "Results.merge: total injected count is zero; cannot pool")

        events = []
        gen_times = []
        weights = []
        for r in runs:
            scale = r.injected / n_total
            events.extend(r.events)
            gen_times.extend(r.gen_times)
            weights.extend(w * scale for w in r._weights)

        merged = cls(events, weights, gen_times,
                     runs[0]._weighter, runs[0]._injector)
        merged.attempts = sum(r.attempts for r in runs)
        merged.requested = sum(r.requested for r in runs)
        merged.injected = n_total
        return merged


def _channel_labels(mixture):
    """Human labels for a mixture's channels, falling back to an index tag."""
    labels = []
    for j, ch in enumerate(mixture.channels):
        label = None
        name = getattr(ch, "Name", None)
        if callable(name):
            try:
                label = name()
            except Exception:
                label = None
        if not label:
            label = getattr(ch, "label", None)
        labels.append(str(label) if label else "channel[{}]".format(j))
    return labels
