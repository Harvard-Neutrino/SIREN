"""
Container for simulation results (events + weights).
"""

from ._util import SaveEvents as _SaveEvents


class Results:
    """Weighted event collection returned by :meth:`Simulation.run`.

    Attributes
    ----------
    events : list
        List of ``InteractionTree`` objects.
    weights : list[float]
        Per-event physics weights.
    gen_times : list[float]
        Per-event generation wall-clock times.
    """

    def __init__(self, events, weights, gen_times, weighter, injector):
        self.events = events
        self.weights = weights
        self.gen_times = gen_times
        self._weighter = weighter
        self._injector = injector

    def __len__(self):
        return len(self.events)

    def __iter__(self):
        """Iterate over (event, weight) pairs."""
        return zip(self.events, self.weights)

    def __getitem__(self, idx):
        """Index or slice into the results."""
        if isinstance(idx, slice):
            return Results(
                self.events[idx],
                self.weights[idx],
                self.gen_times[idx],
                self._weighter,
                self._injector,
            )
        return self.events[idx], self.weights[idx]

    def save(self, filename, **kwargs):
        """Save events and weights to HDF5/Parquet.

        Parameters
        ----------
        filename : str
            Output path (without extension).
        **kwargs
            Additional keyword arguments forwarded to
            :func:`siren._util.SaveEvents`.
        """
        _SaveEvents(
            self.events,
            self._weighter,
            self.gen_times,
            output_filename=filename,
            **kwargs,
        )

    @property
    def n_events(self):
        """Number of generated events."""
        return len(self.events)

    def summary(self):
        """Print a brief summary of the results."""
        import math

        n = len(self.weights)
        if n == 0:
            print("Results: 0 events")
            return

        w_sum = sum(self.weights)
        w_nonzero = sum(1 for w in self.weights if w > 0)
        w_max = max(self.weights)
        w_min_nonzero = min((w for w in self.weights if w > 0), default=0)

        print(f"Results: {n} events")
        print(f"  Non-zero weights: {w_nonzero} / {n} ({100*w_nonzero/n:.1f}%)")
        print(f"  Weight sum:       {w_sum:.6e}")
        print(f"  Weight range:     [{w_min_nonzero:.3e}, {w_max:.3e}]")
