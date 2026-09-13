# Simulation and Results

`Simulation` resolves convenience arguments into an Injector and a Weighter.
A primary or secondary `Vertex` contributes its sampling interactions and
injection distributions, plus its `physical_interactions` and `physical`
declarations for weighting. An empty `physical` list stays empty. Named
physical energy/direction arguments replace those roles; explicit physical
distribution lists append to the declaration. Vertex weighting modes,
kinematics, and expansion rules are retained.

Migration from the earlier Simulation facade: a primary `Vertex` no longer
infers physical energy or direction factors from its injection distributions.
This also applies to `injection_energy=` and `injection_direction=` with a
Vertex primary. To retain the previous normalization, declare those factors
explicitly, for example:

```python
energy = siren.dist.PowerLaw(2, 0.5, 5.0)
direction = siren.dist.IsotropicDirection()
primary = siren.Vertex(
    "NuMu", interactions,
    distributions=[mass, energy, direction, position],
    physical=[energy, direction],
)
```

Here `interactions`, `mass`, and `position` are the models and injection
distributions from the existing configuration. Alternatively, supply
`physical_energy=energy` and `physical_direction=direction` to Simulation.
Leave `physical=[]` when the intended physical target has neither factor.
Avoid adding a factor twice when it is already declared on the Vertex.

Pass `event_factor=callable` to multiply the whole-tree physical weight once.
`sim.run(on_failure="raise")` stops on diagnosed sampling or numerical failures
and retries geometric misses. The policy also applies to warm-up tuning when
`optimize=True` or a `tune.Plan` is supplied. The default remains `"retry"`.

`sim.reweight(...)` uses the last run's events and the original physical target
for omitted arguments. A new `event_factor` replaces the original callback;
`event_factor=None` disables it. `physical_distributions=[]` clears the primary
physical factors and `secondary_physical_distributions={}` clears the secondary
ones. Reweighting does not modify earlier Results or accumulate changes from
previous reweight calls. Changed injector counts are rejected. Keep the
injection distributions and channel weights unchanged when reweighting old
events: event records do not snapshot the generation density.

Results store their weights and run counts. Slicing and `where(...)` retain the
original run counts, so a selection does not change the run's normalization.
`Results.merge(...)` combines independent runs with matching sampling and physical
configurations, scaling each run by its fraction of the pooled injected count.
Each run needs a positive injected count. Models and distributions must compare
equal; detector objects and opaque callbacks must be shared. Keep model,
distribution, and callback state unchanged across runs. The merged explanation
uses a pooled Weighter and includes the event factor.
Both the injection detector and the physical detector must be shared across
runs. Distribution normalizations are compared separately from their shapes.
Built-in channels compare by proposal configuration, including nested mixture
weights; their tuning statistics are ignored. Opaque channels whose configuration
cannot be inspected must be shared.

Only one kinematics declaration may configure a given vertex/signature:
combining Vertex `kinematics` with overlapping `biasing` or `bias_targets`
raises a configuration error. Each secondary compiles only signatures for its
configured particle type, so a shared model or biasing declaration can cover
multiple secondary types without replacing another type's channels.
`sim.phase_spaces` includes primary and secondary
registrations. For a signature present at both, inspect the injector's individual
processes to distinguish the two mixtures.

`results.save(...)` writes the stored weights consistently to HDF5, Parquet,
native events, and HepMC3 when available. It does not call the weighter again.
HDF5 and HepMC3 carry the captured run counts; HDF5 also accepts `pot=...`.
Weight-policy and injector overrides are rejected to avoid replacing the
snapshot. Saving restores the in-memory event headers afterward.

`explain(i)` evaluates the retained weighter and checks that its total still
matches the stored weight. It raises if the calculation has drifted. Output
files can preserve corrected weights without serializing the callback;
Simulation, Results, and Weighter archives reject an active `event_factor`.
Shallow copies retain the in-memory callback and shared configuration objects.
A drift error includes any diagnostic flags and exposes the evaluated breakdown
as `error.breakdown`.

Tuning restores the original event quota and clears warm-up counters even if it
raises. A strict generation error retains its failure report. Channel weights
are not rolled back; generate a fresh sample before weighting after tuning.
