"""
Multi-channel weight optimization for phase-space biasing.

Provides iterative optimization of the per-channel weights in a
MultiChannelPhaseSpace, using the conditional Kleiss-Pittau update
formula with total event weights.
"""

from __future__ import annotations

import math
from typing import List, Optional, Tuple, Dict, TYPE_CHECKING

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


def _submixture(channel):
    """Return a channel's inner MultiChannelPhaseSpace if it wraps one
    (e.g. a NestedMixtureChannel), else None.  Used to optimize a nested
    channel's inner weights as a sub-level."""
    inner = getattr(channel, "mixture", None)
    if inner is not None and len(getattr(inner, "channels", [])) >= 2:
        return inner
    return None


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

    # Channels that wrap a sub-mixture: their inner weights are optimized as a
    # nested level using the SAME outer denominator g.  The inner sub-channel
    # densities g_j contribute to g weighted by this channel's outer alpha, so
    # inner channel j's variance contribution is <w^2 g_j / g> -- the same KP
    # form, just with g_j the inner density and g the full outer density.
    nested = [(i, sub) for i, sub in
              ((k, _submixture(ch)) for k, ch in enumerate(mc.channels))
              if sub is not None]

    for iteration in range(n_iterations):
        var_contrib = [0.0] * n_channels
        inner_contrib = {i: [0.0] * len(sub.channels) for i, sub in nested}
        n_nonzero = 0

        for _ in range(batch_size):
            record = copy.deepcopy(template_record)
            mc.Sample(random, detector_model, record)

            f = physical_density_fn(record)
            g = mc.Density(detector_model, record)

            if g <= 0 or not math.isfinite(g):
                continue
            if f <= 0 or not math.isfinite(f):
                continue

            w2 = (f / g) ** 2
            n_nonzero += 1

            for i in range(n_channels):
                gi = mc.channels[i].Density(detector_model, record)
                if gi > 0 and math.isfinite(gi):
                    var_contrib[i] += w2 * gi / g

            for i, sub in nested:
                acc = inner_contrib[i]
                for j in range(len(sub.channels)):
                    gj = sub.channels[j].Density(detector_model, record)
                    if gj > 0 and math.isfinite(gj):
                        acc[j] += w2 * gj / g

        if n_nonzero == 0:
            continue

        old = list(mc.weights)
        new_alpha = _kp_update(
            old, [v / n_nonzero for v in var_contrib], update_rule, min_weight)
        if new_alpha is not None:
            mc.weights = [
                damping * new_alpha[i] + (1 - damping) * old[i]
                for i in range(n_channels)
            ]

        for i, sub in nested:
            old_in = list(sub.weights)
            new_in = _kp_update(
                old_in, [v / n_nonzero for v in inner_contrib[i]],
                update_rule, min_weight)
            if new_in is not None:
                sub.weights = [
                    damping * new_in[k] + (1 - damping) * old_in[k]
                    for k in range(len(sub.channels))
                ]

    return list(mc.weights)


def optimize_chain_weights(
    injector,
    weighter,
    n_iterations: int = 3,
    batch_size: int = 500,
    damping: float = 0.5,
    min_weight: float = 0.005,
    update_rule: str = "sqrt_W",
    metric=None,
    verbose: bool = False,
) -> None:
    """Optimize all multi-channel weights across a full injection chain.

    Uses conditional optimization: each vertex's weights are updated
    using the total event weight (not just the per-vertex weight).
    This captures cross-vertex correlations.

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
        Kleiss-Pittau update variant, "sqrt_W" (default) or the canonical
        "alpha_sqrt_W" (alpha_i * sqrt(W_i)).  See ``_kp_update``.
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
    """
    import copy
    import numpy as np

    def _sig_match(a, b):
        """Compare signatures by fields (workaround for pybind __eq__)."""
        if int(a.primary_type) != int(b.primary_type):
            return False
        if int(a.target_type) != int(b.target_type):
            return False
        a_sec = [int(s) for s in a.secondary_types]
        b_sec = [int(s) for s in b.secondary_types]
        return a_sec == b_sec

    # Access the C++ injector to get phase space maps
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

    def _collect_phase_spaces(process, ptype):
        """Collect multi-channel phase spaces from a process."""
        results = []
        if not process.HasAnyPhaseSpace():
            return results
        interactions = process.interactions
        for sig_list_getter in [interactions.GetDecays, interactions.GetCrossSections]:
            try:
                models = sig_list_getter()
            except Exception:
                continue
            for model in models:
                try:
                    sigs = model.GetPossibleSignatures()
                except Exception:
                    try:
                        sigs = model.GetPossibleSignaturesFromParent(ptype)
                    except Exception:
                        continue
                for sig in sigs:
                    if process.HasPhaseSpace(sig):
                        mc = process.GetPhaseSpace(sig)
                        if len(mc.channels) >= 2:
                            results.append({
                                'ptype': ptype,
                                'sig': sig,
                                'mc': mc,
                                'n_channels': len(mc.channels),
                            })
        return results

    vertex_info = []

    # Collect primary process phase spaces
    primary_proc = cpp_inj.GetPrimaryProcess()
    primary_ptype = primary_proc.primary_type
    for vi in _collect_phase_spaces(primary_proc, primary_ptype):
        vi['label'] = 'primary'
        vertex_info.append(vi)

    # Collect secondary process phase spaces
    sec_map = cpp_inj.GetSecondaryProcessMap()
    for ptype, process in sec_map.items():
        for vi in _collect_phase_spaces(process, ptype):
            vi['label'] = 'secondary'
            vertex_info.append(vi)

    if not vertex_info:
        if verbose:
            print("No multi-channel phase spaces found to optimize.")
        return

    if verbose:
        print(f"Found {len(vertex_info)} vertices with multi-channel phase spaces:")
        for vi in vertex_info:
            print(f"  {int(vi['ptype'])} ({vi['label']}): "
                  f"{vi['n_channels']} channels, "
                  f"weights = {[f'{w:.3f}' for w in vi['mc'].weights]}")

    for iteration in range(n_iterations):
        # Generate full-chain events, collecting both successes and
        # partial trees from failures
        cpp_inj.ResetInjectedEvents(batch_size)
        events = []
        total_weights = []
        failed_trees = []

        for attempt in range(batch_size):
            event = cpp_inj.GenerateEvent()
            if not event.tree:
                ft = cpp_inj.GetLastFailedTree()
                if ft.tree:
                    failed_trees.append(ft)
                continue
            w = weighter(event)
            if metric is not None:
                w = metric(event, w)
            if not math.isfinite(w):
                w = 0.0
            events.append(event)
            total_weights.append(w)

        if len(events) == 0:
            if verbose:
                print(f"  Iteration {iteration}: no valid events generated")
            continue

        n_total = len(events) + len(failed_trees)
        n_nonzero = np.count_nonzero(total_weights)
        if verbose:
            w_nz = np.array([x for x in total_weights if x != 0])
            if len(w_nz) > 0:
                eff = (w_nz.sum()**2) / (len(w_nz) * (w_nz**2).sum()) * 100 * n_nonzero / n_total
            else:
                eff = 0.0
            print(f"\n  Iteration {iteration}: {len(events)} events, "
                  f"{len(failed_trees)} failures "
                  f"({n_nonzero} pass metric), eff = {eff:.1f}%")

        # Update weights for each vertex using total event weights
        for vi in vertex_info:
            mc = vi['mc']
            n_ch = vi['n_channels']
            ptype_int = int(vi['ptype'])
            sig = vi['sig']

            var_contrib = [0.0] * n_ch
            n_matched = 0

            for event, w_total in zip(events, total_weights):
                for datum in event.tree:
                    r = datum.record
                    if (int(r.signature.primary_type) == ptype_int
                            and _sig_match(r.signature, sig)):
                        g = mc.Density(None, r)
                        if g <= 0 or not math.isfinite(g):
                            continue
                        w2 = w_total * w_total
                        for i in range(n_ch):
                            gi = mc.channels[i].Density(None, r)
                            if gi > 0 and math.isfinite(gi):
                                var_contrib[i] += w2 * gi / g
                        n_matched += 1

            if n_matched == 0:
                continue

            # Estimate per-channel failure rates from partial trees.
            # For each failed tree, find this vertex's record and
            # compute the channel selection probability p_i = a_i*g_i/g.
            # Channels that are more likely to have been selected for
            # failed events get penalized.
            fail_select = [0.0] * n_ch
            succ_select = [0.0] * n_ch
            n_fail_matched = 0
            n_succ_matched = 0

            for ft in failed_trees:
                for datum in ft.tree:
                    r = datum.record
                    if (int(r.signature.primary_type) == ptype_int
                            and _sig_match(r.signature, sig)):
                        g = mc.Density(None, r)
                        if g <= 0 or not math.isfinite(g):
                            break
                        for i in range(n_ch):
                            gi = mc.channels[i].Density(None, r)
                            if gi > 0 and math.isfinite(gi):
                                fail_select[i] += mc.weights[i] * gi / g
                        n_fail_matched += 1
                        break

            for event in events:
                for datum in event.tree:
                    r = datum.record
                    if (int(r.signature.primary_type) == ptype_int
                            and _sig_match(r.signature, sig)):
                        g = mc.Density(None, r)
                        if g <= 0 or not math.isfinite(g):
                            break
                        for i in range(n_ch):
                            gi = mc.channels[i].Density(None, r)
                            if gi > 0 and math.isfinite(gi):
                                succ_select[i] += mc.weights[i] * gi / g
                        n_succ_matched += 1
                        break

            # Apply failure penalty: inflate variance contribution
            # for channels with high failure rates.
            # Effective variance scales as 1/(1-f_i) where f_i is
            # the channel's failure rate.
            if n_fail_matched > 0 and n_succ_matched > 0:
                for i in range(n_ch):
                    total_i = succ_select[i] + fail_select[i]
                    if total_i > 0:
                        f_i = fail_select[i] / total_i
                        if f_i < 1.0:
                            var_contrib[i] /= (1.0 - f_i)

            W = [v / n_matched for v in var_contrib]
            old = list(mc.weights)
            new_alpha = _kp_update(old, W, update_rule, min_weight)
            if new_alpha is None:
                continue

            mc.weights = [
                damping * new_alpha[i] + (1 - damping) * old[i]
                for i in range(n_ch)
            ]

            if verbose:
                print(f"    {ptype_int} ({vi.get('label', '?')}): "
                      f"{[f'{w:.3f}' for w in old]} -> "
                      f"{[f'{w:.3f}' for w in mc.weights]}")

    cpp_inj.ResetInjectedEvents(prev_events_to_inject)

    if verbose:
        print("\nOptimization complete. Final weights:")
        for vi in vertex_info:
            print(f"  {int(vi['ptype'])} ({vi.get('label', '?')}): "
                  f"{[f'{w:.3f}' for w in vi['mc'].weights]}")
