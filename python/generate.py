"""One-call generation + weighting entry point.

generate(injector, weighter, *, events, ...) drives the injector to `events`
successful trees and weights them through the weighter, returning a Results
snapshot of the trees and their weights.
"""

from __future__ import annotations

from .Results import Results


def generate(injector, weighter, *, events, on_shortfall="warn",
             progress=None, min_efficiency=None):
    """Generate `events` weighted trees.

    Delegates generation to ``injector.generate`` (which counts successes and
    honours ``on_shortfall``/``min_efficiency``), weights via
    ``weighter.weight_all``, and returns a Results over the trees and weights.
    """
    trees = injector.generate(
        events, on_shortfall=on_shortfall, progress=progress,
        min_efficiency=min_efficiency)
    weights = weighter.weight_all(trees)
    gen_times = [0.0] * len(trees)
    return Results(list(trees), list(weights), gen_times, weighter, injector,
                   requested=events)
