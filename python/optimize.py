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
) -> List[float]:
    """Optimize channel weights to minimize variance of f/g.

    Uses the Kleiss-Pittau update: alpha_i^new proportional to
    sqrt(mean(w^2 * g_i / g)) where w = f(x) / g(x).

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

    for iteration in range(n_iterations):
        var_contrib = [0.0] * n_channels
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

            w = f / g
            w2 = w * w
            n_nonzero += 1

            for i in range(n_channels):
                gi = mc.channels[i].Density(detector_model, record)
                if gi > 0 and math.isfinite(gi):
                    var_contrib[i] += w2 * gi / g

        if n_nonzero == 0:
            continue

        new_alpha = [math.sqrt(max(v / n_nonzero, 0.0)) for v in var_contrib]
        total = sum(new_alpha)
        if total <= 0:
            continue

        new_alpha = [a / total for a in new_alpha]

        # Apply minimum weight floor
        for i in range(n_channels):
            new_alpha[i] = max(new_alpha[i], min_weight)
        total = sum(new_alpha)
        new_alpha = [a / total for a in new_alpha]

        # Damped update
        old = list(mc.weights)
        mc.weights = [
            damping * new_alpha[i] + (1 - damping) * old[i]
            for i in range(n_channels)
        ]

    return list(mc.weights)


def optimize_chain_weights(
    injector,
    weighter,
    n_iterations: int = 3,
    batch_size: int = 500,
    damping: float = 0.5,
    min_weight: float = 0.005,
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
        # Generate full-chain events
        cpp_inj.ResetInjectedEvents(batch_size)
        events = []
        total_weights = []

        for attempt in range(batch_size):
            event = cpp_inj.GenerateEvent()
            if not event.tree:
                continue
            w = weighter(event)
            # Apply metric to transform the weight
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

        n_nonzero = np.count_nonzero(total_weights)
        if verbose:
            w_nz = np.array([x for x in total_weights if x != 0])
            if len(w_nz) > 0:
                eff = (w_nz.sum()**2) / (len(w_nz) * (w_nz**2).sum()) * 100 * n_nonzero / len(events)
            else:
                eff = 0.0
            print(f"\n  Iteration {iteration}: {len(events)} events "
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
                # Find the record at this vertex
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

            new_alpha = [math.sqrt(max(v / n_matched, 0)) for v in var_contrib]
            total = sum(new_alpha)
            if total <= 0:
                continue
            new_alpha = [a / total for a in new_alpha]

            for i in range(n_ch):
                new_alpha[i] = max(new_alpha[i], min_weight)
            total = sum(new_alpha)
            new_alpha = [a / total for a in new_alpha]

            old = list(mc.weights)
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
