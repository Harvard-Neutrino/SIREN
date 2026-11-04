# Results and simulation runs

`Simulation(primary=Vertex(...), secondaries=[Vertex(...)], ...)` delegates to
Injector's existing compiler and generation policies. Every `run()` constructs
a fresh injector; repeating a fixed seed repeats the stream. Use independent
seeds for independent samples. Production examples set both `on_failure` and
`on_shortfall` to `"raise"`; defaults retain retry/warn compatibility.

`siren.generate(injector, weighter, events=...)` returns Results for one fresh
run. It refuses used counters. Use `Injector.generate()` for cumulative batches
and weight all contributing events after generating the last batch. Explicit
`reset()` starts a new run; it cannot retroactively change earlier Results.

Results stores read-only weights, generation times (the uniform run average),
and attempted/requested/accepted counts. Event trees are shallow references;
modifying a tree changes the event and invalidates its stored weight. A slice
or `where(predicate)` retains the full run counts, even if no events survive.

For serializable native configurations, Results retains a separate copy of the
proposal and physical target, including counters, detector, selected decay
channels, phase-space mixtures, expansion and creation-survival mode. Changes
to the original injector or physical distributions do not affect later
`result.reweight(...)`. Reweighting changes the physical target while retaining
the original proposal; it cannot create missing support.

Live Python physics that cannot be archived remains usable for generation and
event export. Reweighting and pooling refuse configurations they cannot freeze.
A live `event_factor` that cannot be copied must be supplied explicitly for
reweighting. Callables must be deterministic for fixed model state; a pickleable
function can still refer to mutable globals and does not make those globals
part of a snapshot.

`Results.merge([a, b])` pools matching, frozen configurations from independent
injectors. Weights from run i are scaled by `attempts_i / total_attempts`,
including failed attempts and zero-success runs. Requested and accepted counts
are metadata, not the normalization denominator. Results from the same source
injector (including across resets), overlapping subsets, unknown configurations,
and differing physical targets are refused. Distinct injectors with identical
seeds are refused as correlated samples. This check does not prove independence
of externally supplied sources; callers retain that responsibility.

`result.save("output/run")` exports stored weights without evaluating physics
again. HDF5 includes the run counts and optional `pot` metadata; POT does not
multiply weights. Parquet and native event files retain event weights.
`save(..., snapshot=True)` additionally writes `.siren_results`, carrying the
configuration needed by `Results.load("output/run").reweight(...)`. Native event
files alone are not Results snapshots. Live unarchivable physics and all non-None event factors (including pickleable
functions) are refused before writing any snapshot output, matching Weighter
archive policy. Load pickle snapshots only
from trusted sources and use a compatible SIREN version.

Legacy Simulation energy/direction/interaction slots issue a deprecation warning
and translate once to Vertex. Loader options must be in `process_options={...}`
or an explicit `load_processes` call. Unknown constructor arguments fail before
loading data, including a suggested spelling for close matches. `n_events` is a
deprecated alias for `events`. `SIREN_Controller` and `GenerateEvents` remain
available with deprecation warnings. Tuning, weight explanations, beam readers
and directed-channel convenience APIs are separate stack units.
