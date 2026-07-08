"""One-call generation + weighting entry point.

generate(injector, weighter, *, events, ...) drives the injector to `events`
successful trees and weights them through the weighter, returning a Results
snapshot of the trees and their weights.
"""

from __future__ import annotations

from .Results import Results
from .Weighter import _checked_weight
from .errors import WeightCalculationError


def generate(injector, weighter, *, events, on_shortfall="warn",
             progress=None, min_efficiency=None, on_failure="retry"):
    """Generate `events` weighted trees.

    Delegates generation to ``injector.generate`` (which counts successes and
    honours ``on_shortfall``/``min_efficiency``/``on_failure``), weights via
    ``weighter.weight_all``, and returns a Results over the trees and weights.
    Every tree must have one finite nonnegative weight, including when a custom
    weighter supplies the batch. Invalid weights raise WeightCalculationError.
    """
    trees = injector.generate(
        events, on_shortfall=on_shortfall, progress=progress,
        min_efficiency=min_efficiency, on_failure=on_failure)
    weights = list(weighter.weight_all(trees))
    if len(weights) != len(trees):
        raise WeightCalculationError(
            "Expected {} event weights, got {}".format(len(trees), len(weights)))
    weights = [_checked_weight(weight, "Event {} weight", index)
               for index, weight in enumerate(weights)]
    gen_times = [0.0] * len(trees)
    return Results(list(trees), list(weights), gen_times, weighter, injector,
                   requested=events)
