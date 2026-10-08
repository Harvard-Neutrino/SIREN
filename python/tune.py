"""
Multi-channel weight optimization for phase-space biasing.

Provides iterative optimization of the per-channel weights in a
MultiChannelPhaseSpace, using the conditional Kleiss-Pittau update
formula with total event weights.
"""

from __future__ import annotations

import math
from typing import Callable, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from . import injection as _injection
    from . import detector as _detector
    from . import dataclasses as _dataclasses


def _kp_update(
    old_weights: List[float],
    W: List[float],
    update_rule: str = "sqrt_W",
    min_weight: float = 0.0,
) -> Optional[List[float]]:
    """One Kleiss-Pittau channel-weight update step.

    Given the current weights ``old_weights`` and the per-channel variance
    contributions ``W[i] = <(g_i/g) w^2>`` (already averaged over the
    batch), return the new normalized weights (before any damping blend).

    update_rule:
      ``"sqrt_W"``        -- new_i ~ sqrt(W_i)             (memoryless;
                             fixed point alpha_i ~ sqrt(W_i))
      ``"alpha_sqrt_W"``  -- new_i ~ alpha_i * sqrt(W_i)   (canonical
                             Kleiss-Pittau; fixed point W_i = const, the
                             variance minimum -- every channel contributes
                             equal variance)

    Both rules leave the integral estimate unbiased (Kleiss-Pittau holds
    for any valid alpha); they differ only in which weight set they
    converge to.  Returns ``None`` when the update is degenerate (all
    contributions zero), so the caller can keep the previous weights.
    """
    n = len(W)
    roots = [math.sqrt(v) if v > 0.0 else 0.0 for v in W]
    if update_rule == "alpha_sqrt_W":
        new_alpha = [old_weights[i] * roots[i] for i in range(n)]
    elif update_rule == "sqrt_W":
        new_alpha = list(roots)
    else:
        raise ValueError(
            f"unknown update_rule {update_rule!r} "
            "(expected 'sqrt_W' or 'alpha_sqrt_W')")

    total = sum(new_alpha)
    if total <= 0.0:
        return None
    new_alpha = [a / total for a in new_alpha]

    if min_weight > 0.0:
        new_alpha = [max(a, min_weight) for a in new_alpha]
        total = sum(new_alpha)
        new_alpha = [a / total for a in new_alpha]

    return new_alpha


def group_directed_channels(mc, label="DirectedGroup"):
    """Wrap a flat mixture's directed channels into one NestedMixtureChannel.

    A detector-directed channel exposes ``DirectingActive``; an
    already-nested group exposes ``mixture``.  This collects the (non-nested)
    directed channels of ``mc`` into a single ``NestedMixtureChannel`` and
    returns a new ``MultiChannelPhaseSpace`` of ``[non-directed channels...,
    group]``, leaving the overall density g(x) EXACTLY unchanged (the group's
    outer weight is the sum of the directed weights and its inner weights are
    those weights renormalized).

    Why: the directed channels share an isotropic 1/4pi fallback, so when they
    are inert their variance contributions are degenerate.  As N separate
    channels the optimizer drives each to ``min_weight`` independently, leaving
    ``N * min_weight`` of wasted fallback admixture; as one group it sees a
    single "direct vs not" weight that floors at one ``min_weight`` (N times
    smaller) and converges faster, while the inner "which target" weights are
    tuned recursively only on directing-active events.  The group reports
    ``DirectingActive`` = "any member directs", so ``discount_fallback`` drives
    the whole group down when every member is in fallback.

    Returns ``mc`` unchanged if it has fewer than two groupable directed
    channels.
    """
    from . import injection as inj

    channels = list(mc.channels)
    weights = list(mc.weights)

    def groupable(ch):
        return hasattr(ch, "DirectingActive") and not hasattr(ch, "mixture")

    dir_idx = [i for i, ch in enumerate(channels) if groupable(channels[i])]
    if len(dir_idx) < 2:
        return mc

    dir_set = set(dir_idx)
    dir_weight_total = sum(weights[i] for i in dir_idx)

    inner = inj.MultiChannelPhaseSpace()
    inner.channels = [channels[i] for i in dir_idx]
    if dir_weight_total > 0:
        inner.weights = [weights[i] / dir_weight_total for i in dir_idx]
    else:
        inner.weights = [1.0 / len(dir_idx)] * len(dir_idx)
    group = inj.NestedMixtureChannel(inner)
    group.label = label

    out = inj.MultiChannelPhaseSpace()
    out.channels = [channels[i] for i in range(len(channels)) if i not in dir_set] + [group]
    out.weights = [weights[i] for i in range(len(channels)) if i not in dir_set] + [dir_weight_total]
    return out


def optimize_multichannel_weights(
    mc: "_injection.MultiChannelPhaseSpace",
    physical_density_fn,
    random: object,
    detector_model: "_detector.DetectorModel",
    template_record: "_dataclasses.InteractionRecord",
    n_iterations: int = 5,
    batch_size: int = 1000,
    damping: float = 0.5,
    min_weight: float = 0.005,
    update_rule: str = "sqrt_W",
    discount_fallback: bool = True,
) -> List[float]:
    """Optimize channel weights to minimize variance of f/g.

    Uses a Kleiss-Pittau update of the per-channel variance contribution
    W_i = mean(w^2 * g_i / g), where w = f(x) / g(x).  See ``_kp_update``
    for the two rules selected by ``update_rule`` ("sqrt_W" or the
    canonical "alpha_sqrt_W").  Channels that wrap a sub-mixture (e.g.
    NestedMixtureChannel) have their inner weights optimized recursively as
    a nested level, with the same outer denominator g.

    Parameters
    ----------
    mc : MultiChannelPhaseSpace
        The multi-channel to optimize.  Weights are modified in-place.
    physical_density_fn : callable(record) -> float
        Evaluates the physical density f(x) at the phase-space point
        described by the record (e.g. model.FinalStateProbability).
    random : SIREN_random
        Random number generator.
    detector_model : DetectorModel or None
        Detector model (passed to Sample/Density).
    template_record : InteractionRecord
        Template record with signature, primary momentum, masses, etc.
        Secondary momenta will be overwritten by sampling.
    n_iterations : int
        Number of optimization iterations.
    batch_size : int
        Events per iteration.
    damping : float
        Blending factor: new_alpha = damping * update + (1-damping) * old.
        Use 1.0 for pure update, 0.5 for damped.
    min_weight : float
        Minimum weight for any channel (prevents channels from being
        turned off entirely).

    Returns
    -------
    list[float]
        The optimized weights (also written to mc.weights).
    """
    import copy

    n_channels = len(mc.channels)
    if n_channels < 2:
        return list(mc.weights)

    # The mixture owns its Kleiss-Pittau statistics.  Each warm-up point is fed
    # to mc.Accumulate (which credits every channel's bare contribution
    # w^2 * g_i/g, honors discount_fallback, and recurses into nested groups
    # against the same outer g); mc.UpdateWeights then applies the chosen rule,
    # damping, the min_weight floor and the degenerate-keep, and resets the
    # accumulators.
    for _ in range(n_iterations):
        mc.ResetAccumulators()

        for _ in range(batch_size):
            record = copy.deepcopy(template_record)
            mc.Sample(random, detector_model, record)

            f = physical_density_fn(record)
            g = mc.Density(detector_model, record)
            if g <= 0 or not math.isfinite(g):
                continue
            if f <= 0 or not math.isfinite(f):
                continue

            mc.Accumulate(detector_model, record, f / g, discount_fallback)

        mc.UpdateWeights(update_rule, damping, min_weight)

    return list(mc.weights)


def _resolve_engine(injector):
    """Return the raw C++ injector for either a wrapper or raw Injector.

    Forces the wrapper to build (which requires one round-trip event
    generation to warm up the mixture accumulators) via its public
    ``engine`` property; a raw C++ injector is returned unchanged.
    """
    from .Injector import Injector as _PyInjector
    if isinstance(injector, _PyInjector):
        return injector.engine
    return injector


def _run_accumulation_round(cpp_inj, mixtures, weighter, batch_size, metric,
                            discount_fallback, recurse_nested):
    """Run one warm-up batch: generate events, accumulate their weights and
    selection status into every vertex mixture, and return the batch's
    non-zero finite weights (for the ESS estimate) plus event/failure counts.

    Shared by ``optimize_chain_weights`` and ``tune`` so both drive identical
    per-round accumulation semantics.
    """
    cpp_inj.ResetInjectedEvents(batch_size)
    for mc in mixtures:
        mc.ResetAccumulators(recurse_nested)

    n_events = 0
    n_failures = 0
    weights_seen = []
    for _ in range(batch_size):
        event = cpp_inj.GenerateEvent()
        if not event.tree:
            # A partial/failed tree: feed its per-vertex selection so the
            # failure penalty inflates channels that tend to fail.
            ft = cpp_inj.GetLastFailedTree()
            if ft.tree:
                cpp_inj.AccumulateSelectionToMixtures(ft, True)
                n_failures += 1
            continue
        w = weighter(event)
        if metric is not None:
            w = metric(event, w)
        if not math.isfinite(w):
            w = 0.0
        # Feed the TOTAL event weight to every vertex mixture this event
        # touched (the C++ accumulator credits each channel w^2*g_i/g), plus
        # the success selection for the failure penalty.  The signature ->
        # mixture routing happens in C++.  With recurse_nested, the inner
        # per-target channels of grouped directed channels are credited too,
        # against the same outer vertex g.
        cpp_inj.AccumulateEventToMixtures(event, w, discount_fallback, recurse_nested)
        cpp_inj.AccumulateSelectionToMixtures(event, False)
        n_events += 1
        weights_seen.append(w)

    return weights_seen, n_events, n_failures


def _ess(weights_seen: List[float]) -> float:
    """Effective sample size ``(sum w)^2 / sum(w^2)`` over non-zero finite
    weights; 0.0 if the batch has no usable weight."""
    w_nz = [w for w in weights_seen if w != 0.0 and math.isfinite(w)]
    if not w_nz:
        return 0.0
    s1 = sum(w_nz)
    s2 = sum(w * w for w in w_nz)
    if s2 <= 0.0:
        return 0.0
    return (s1 * s1) / s2


def optimize_chain_weights(
    injector,
    weighter,
    n_iterations: int = 3,
    batch_size: int = 500,
    damping: float = 0.5,
    min_weight: float = 0.005,
    update_rule: str = "alpha_sqrt_W",
    discount_fallback: bool = True,
    metric=None,
    verbose: bool = False,
    recurse_nested: bool = True,
    failure_mode: str = "throughput",
) -> None:
    """Optimize all multi-channel weights across a full injection chain.

    Uses conditional optimization: each vertex's weights are updated
    using the total event weight (not just the per-vertex weight).
    This captures cross-vertex correlations.

    With ``recurse_nested`` (default), the optimizer also descends into any
    grouped directed channels (``NestedMixtureChannel``, e.g. produced by
    ``group_directed_channels``) and tunes their inner per-target weights with
    the same total-event-weight statistic -- so the global optimization balances
    the distribution among targets, not just the outer "direct vs physical"
    weight.

    Parameters
    ----------
    injector : siren.Injector
        A configured injector (will be used to generate events).
    weighter : siren.Weighter
        A configured weighter (will be used to compute total weights).
    n_iterations : int
        Number of full-chain optimization iterations.
    batch_size : int
        Events per iteration.
    damping : float
        Blending factor for weight updates.
    min_weight : float
        Minimum per-channel weight.
    update_rule : str
        Kleiss-Pittau update variant.  Default is the canonical
        "alpha_sqrt_W" (alpha_i * sqrt(W_i); fixed point W_i = const, the
        variance minimum).  The memoryless "sqrt_W" is also available but
        CANNOT turn off a channel whose density already covers the support
        (e.g. a directed channel in its isotropic fallback competing with the
        physical channel): near the optimum W_i ~ 1 for every covering
        channel, so alpha_i ~ sqrt(W_i) parks the redundant channel at a
        finite weight and dilutes the physical channel.  The multiplicative
        canonical rule decays it to ``min_weight`` instead.  See ``_kp_update``.
    metric : callable(event, weight) -> float, optional
        Transforms the event weight before it enters the variance
        computation.  The optimizer minimizes variance of metric(event, w)
        rather than variance of w.

        Common patterns:
          - ``None`` (default): minimize variance of the raw weight.
          - Selection: ``lambda ev, w: w if passes_cut(ev) else 0``
            concentrates sampling on events that pass the cut.
          - Observable: ``lambda ev, w: w * observable(ev)``
            minimizes variance of a weighted observable.
          - Region: ``lambda ev, w: w if in_fiducial(ev) else 0``
            optimizes for events in a specific detector region.
    verbose : bool
        Print progress and weight updates.
    recurse_nested : bool
        If True (default), also tune the inner weights of grouped directed
        channels (NestedMixtureChannel) -- i.e. the distribution among targets --
        using the total event weight.  Set False to tune only the outer vertex
        weights (the pre-existing behavior), leaving each group's internal
        target split at its initial value.
    failure_mode : str
        How injection failures feed back into the per-channel weight (acts on the
        outer per-vertex channels, where failed/successful selection is recorded):
          - "throughput" (default): W_i *= (1 - f_i) -- down-weight channels that
            disproportionately feed failed trees, so the sampling tracks the
            successful contribution to the integral.  Converges to the variance
            optimum fastest and is insensitive to the warm-up batch size.
          - "ignore": no failure adjustment (the success-weighted statistic
            already discounts failures implicitly); reaches the same fixed point
            as "throughput" but more slowly, and a bit faster with fewer samples.
          - "coverage": W_i /= (1 - f_i) -- up-weight lossy channels to keep their
            region sampled (the original behavior; tends to retain lossy channels
            and gives the worst ESS where directing is lossy).
        where f_i is channel i's fraction of selection mass on failed trees.
    """
    import numpy as np

    # Access the C++ injector (the Python wrapper holds it lazily).
    from .Injector import Injector as _PyInjector
    if isinstance(injector, _PyInjector):
        # Force initialization if needed
        if injector._Injector__injector is None:
            # Generate one event to trigger init
            for _ in injector:
                break
            injector._Injector__injector.ResetInjectedEvents(batch_size)
        cpp_inj = injector._Injector__injector
    else:
        cpp_inj = injector

    prev_events_to_inject = cpp_inj.EventsToInject()

    # The injector enumerates its own multi-channel mixtures and (in C++, via the
    # signature map) routes each generated record to the mixture that sampled it.
    # So the Python side no longer walks the process structure or re-matches
    # signatures by hand.
    mixtures = cpp_inj.GetPhaseSpaces()
    if not mixtures:
        if verbose:
            print("No multi-channel phase spaces found to optimize.")
        return

    if verbose:
        print(f"Found {len(mixtures)} multi-channel phase spaces to optimize.")

    for iteration in range(n_iterations):
        weights_seen, n_events, n_failures = _run_accumulation_round(
            cpp_inj, mixtures, weighter, batch_size, metric,
            discount_fallback, recurse_nested)

        if verbose:
            w_nz = np.array([x for x in weights_seen if x != 0.0])
            n_total = n_events + n_failures
            if len(w_nz) > 0 and n_total > 0:
                eff = ((w_nz.sum() ** 2) / (len(w_nz) * (w_nz ** 2).sum())
                       * 100 * len(w_nz) / n_total)
            else:
                eff = 0.0
            print(f"\n  Iteration {iteration}: {n_events} events, "
                  f"{n_failures} failures, eff = {eff:.1f}%")

        # One Kleiss-Pittau update per vertex mixture from the accumulated
        # total-weight statistics (the failure penalty is folded in by
        # UpdateWeights).  With recurse_nested, the inner per-target weights of
        # grouped directed channels are updated too (the failure penalty stays
        # outer-only; inner groups have no selection data).
        for mc in mixtures:
            old = list(mc.weights)
            mc.UpdateWeights(update_rule, damping, min_weight, recurse_nested,
                             failure_mode)
            if verbose:
                print(f"    {[f'{w:.3f}' for w in old]} -> "
                      f"{[f'{w:.3f}' for w in mc.weights]}")

    cpp_inj.ResetInjectedEvents(prev_events_to_inject)

    if verbose:
        print("\nOptimization complete. Final weights:")
        for mc in mixtures:
            print(f"  {[f'{w:.3f}' for w in mc.weights]}")
