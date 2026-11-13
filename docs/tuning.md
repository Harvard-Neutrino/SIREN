# Diagnostics and proposal tuning

`weighter.explain(tree)` returns the same total as scalar weighting for valid
weights. It includes per-vertex generation and physical densities, the injector
index in a pool, cancelled factors, channel densities and their measure, and
creation survival. A skipped injector has unevaluated physical/channel values;
those NaNs describe an omitted contribution. Invalid inputs produce a NaN total
with diagnostic flags, while scalar weighting raises `WeightCalculationError`.
The event-factor callback's own exceptions propagate. Native creation-survival
underflow, previously thrown from the breakdown method, now appears as a flagged
NaN total there; scalar weighting still raises. A positive-factor product that underflows
is an error, distinct from an exact zero physical factor. This checks the
framework's probability composition; a numerical zero already returned by a model
is outside this guarantee. Unrepresentable generation densities fail explicitly.
A negligible positive inverse-weight contribution may round to zero when another
pooled injector supplies a representable total.

`results.explain(i)` uses the saved configuration and verifies the explanation
against the stored weight. `results.bad()` and `results.ess()` summarize stored
weights. Generation failures retain their typed reason in `injector.report()`.
Terminal daughters encountered across all attempts, including discarded trees,
are counted separately in `report.unregistered_secondaries`. They do not discard the parent event or
enter failure counts. Explicit invalid expansion wiring is rejected by the compiler.

## Separate pilot and production runs

```python
from siren.tune import Plan

result = simulation.run(
    optimize=Plan(events=2000, rounds=5, seed=1729, min_weight=0.001),
    on_failure="raise", on_shortfall="raise",
)
report = simulation.optimization_report
print(report.attempts, report.ess_trajectory)
```

The Simulation must have an explicit production seed different from the pilot
seed. Each pilot round makes exactly `events` attempts, including geometric
misses. Only production events enter Results and its attempt normalization.
Pilot weights are evaluated after the whole batch, then discarded. Production
uses a pristine clone of the original RNG and sampler state with fresh counters;
only the finalized proposal coefficients are transferred from the pilot.

For direct use, `report = siren.tune.tune(injector, weighter, seed=1729)` returns
`report.injector` and `report.weighter`. Generate with that pair using
`siren.generate`. The input objects are unchanged. Use fresh input counters;
previously generated events cannot be assigned a changed proposal law.

Tuning requires a serializable native configuration. Separately declared
physical models/distributions are retained in the same graph archive as the
samplers. Picklable stopping conditions and event factors are copied independently. Python model implementations
and unsupported live callbacks fail explicitly. A picklable `event_factor` is
copied for the pilot/production pair; Results retains its existing archive guard.
No pilot is silently omitted when cloning fails.

## Update and acceptance

The default Kleiss–Pittau update is proportional to
`alpha_i * sqrt(mean(w**2 * g_i/g))`. `sqrt_W` omits the old coefficient.
Damping blends with the old proposal. Each component has a strictly positive
floor, including a physical fallback, so tuning cannot delete its support.
With `recurse=True`, nested channels, including U07 tiling coefficients, are
tunable. `recurse=False` keeps their interior coefficients fixed. Nested channels
use the outer mixture's density in the statistic, with all
bare densities converted to the same measure. Repeated references and cycles
in nested mixture graphs are refused. The native accumulator uses log sums to
retain finite extreme weights without squaring overflow.

`metric(event, weight)` can select a region or target a signed observable. It
must return a finite real scalar. `siren.metrics.vertex_in(volume)` builds a
selection predicate; use `lambda event, w: w if predicate(event) else 0`.
Failures and exact zero contributions remain in attempt counts. Failed draws
receive no positive second-moment contribution. No heuristic failure penalty
or fallback discount modifies the canonical statistic.

ESS is `(sum w)**2 / sum(w**2)` over the pilot contributions. Dividing it by
attempts includes failed and zero-contribution draws. `converged` only means
the last two coefficient updates differed by less than 0.001; it does not
establish a converged physical integral. Validate rates with independent
production streams and an analytic or converged numerical reference.

`MultiChannelPhaseSpace.Accumulate`, `UpdateWeights`, and `ResetAccumulators`
are the native interface. Intermediate accumulators cannot be serialized;
finish or reset the update first. The resulting proposal weights use the
existing archive format. Results snapshots remain independent of later tuning.
