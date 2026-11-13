"""Generate and weight one fresh run through the public Layer-2 objects."""
import time

from .Results import Results
from .Weighter import _checked_weight
from .errors import ConfigurationError, WeightCalculationError


def generate(injector, weighter, *, events, on_shortfall='warn', progress=None,
             min_efficiency=None, on_failure='retry'):
    """Return a Results snapshot; use Injector.generate for cumulative batches.

    A Results run starts with fresh counters. This prevents returning only a
    later batch with weights normalized by attempts from unrepresented batches.
    """
    from .Injector import Injector
    from .Weighter import Weighter
    if not isinstance(injector, Injector) or not isinstance(weighter, Weighter):
        raise TypeError('siren.generate requires public Injector and Weighter objects')
    sources = weighter.injectors or []
    if len(sources) != 1 or (sources[0] is not injector and sources[0] is not injector.engine):
        raise ConfigurationError('siren.generate requires a weighter for the same injector')
    if injector.injection_attempts or injector.injected_events:
        raise ConfigurationError('siren.generate requires fresh counters; use a new injector '
                                 'or explicitly reset before starting a new Results run')
    started = time.perf_counter()
    trees = injector.generate(events, on_shortfall=on_shortfall, progress=progress,
                              min_efficiency=min_efficiency, on_failure=on_failure)
    elapsed = time.perf_counter() - started
    weights = list(weighter.weight_all(trees))
    if len(weights) != len(trees):
        raise WeightCalculationError('Expected {} weights, got {}'.format(len(trees), len(weights)))
    weights = [_checked_weight(w, 'Event {} weight', i) for i, w in enumerate(weights)]
    times = [elapsed / len(trees)] * len(trees) if trees else []
    return Results(trees, weights, times, weighter, injector, requested=events)
