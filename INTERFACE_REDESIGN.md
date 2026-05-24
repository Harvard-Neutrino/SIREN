# SIREN Interface Redesign Plan

## Design Principles

1. **Say it once.** Shared configuration (detector, particle type, interactions)
   should be specified exactly once.
2. **Progressive disclosure.** A DIS-in-IceCube script should be ~15 lines.
   BSM with secondaries should be ~30. The full Injector/Weighter API remains
   available underneath for power users.
3. **Validate early, fail clearly.** Missing distributions, type mismatches, and
   injection/physical inconsistencies should be caught at construction time
   with human-readable errors, not deep in C++.
4. **One idiomatic path.** Replace the competing Controller/Injector-Weighter
   APIs with a single `Simulation` class. Deprecate `SIREN_Controller`.
5. **No hidden essentials.** `GenerateEvents`, `SaveEvents`, particle type
   aliases, and distribution constructors should all be importable from
   top-level `siren.*`.
6. **No manual plumbing.** New C++ distributions, particle types, and enum
   values are reflected in the Python interface automatically via
   introspection of the C++ bindings.

---

## Completed Work

### Phase 1: Top-Level Convenience

- `siren.particles` sub-namespace: auto-generated from C++ `ParticleType`
  enum at import time. Human-friendly aliases (Electron, Muon, Proton)
  added on top. String resolution via `particles.resolve("NuMu")`.
- `siren.dist` sub-namespace: auto-generated from `siren.distributions`
  module. Short aliases (ColumnDepth, BoundedVertex, etc.) mapped manually.
- Top-level exports: `siren.Simulation`, `siren.Results`,
  `siren.GenerateEvents`, `siren.SaveEvents`, `siren.LoadEvents`,
  `siren.load_detector`, `siren.load_flux`, `siren.load_processes`,
  `siren.get_fiducial_volume`.

### Phase 2: Simulation Class

- `siren.Simulation`: model-agnostic high-level API. Accepts strings for
  detector, particle type, and interaction model. Extra model kwargs flow
  through via `**process_kwargs`. Owns both Injector and Weighter.
- Distribution resolution: `energy=X` for both injection/physical;
  `injection_energy`/`physical_energy` for split; `flux` alias for
  physical energy; `position` injection-only; mass auto-inferred.
- `siren.Results`: container with `events`, `weights`, `gen_times`.
  Supports iteration, slicing, indexing, `save()`, `summary()`.

### Phase 3: Ergonomics

- `Injector.__iter__` / `__len__` (Python wrapper + C++ pybind `__len__`).
- `Weighter.weight_all()`.
- C++ pybind overloads: `FixedDirection([0,0,1])`, `Cone([0,0,1], 0.1)`,
  `PointSourcePositionDistribution([x,y,z], d)` accept lists/tuples
  natively via `std::array<double,3>` constructors.

### C++ Distribution Variable System

- `DistributionVariable` enum: `PrimaryMass`, `PrimaryEnergy`,
  `PrimaryDirection`, `PrimaryHelicity`, `PrimaryArea`,
  `InitialPosition`, `InteractionVertex`, `InteractionParameters`.
- `SetVariables()`: pure virtual on `PrimaryInjectionDistribution`.
  Returns what variables the distribution determines.
- `RequiredVariables()`: virtual with default empty set. Returns what
  must be in the record before `Sample()` runs.
- Delta functions identified by `SetVariables() \ DensityVariables()`.
- `Monoenergetic::DensityVariables()` fixed to return `{}`.
- `PrimaryExternalDistribution` computes `SetVariables` dynamically
  from CSV column headers at construction time.
- Python validation uses the C++ methods for completeness, ordering,
  and reweighting compatibility checks.

---

## Phase 4: Cleanup and Deprecation (current)

### 4a. Deprecate `SIREN_Controller`

Add `DeprecationWarning` to `SIREN_Controller.__init__`. Keep it
functional for at least one version cycle.

### 4b. Deduplicate code

`SaveEvents`, `GetDetectorModelTargets`, `GetFiducialVolume` are
duplicated between `_util.py` and `SIREN_Controller.py`. Consolidate
into canonical locations; have both modules import from the same source.

### 4c. Standardize `load_processes` return type

DarkNews returns 4 values; everything else returns 2. Wrap in a
`ProcessBundle` NamedTuple so callers get a consistent interface.
Existing `primary, secondary = load_processes(...)` unpacking still
works. Extra metadata accessible as attributes.

### 4d. Update examples

Rewrite examples in `resources/examples/` to use the `Simulation` API.
Move old versions to `resources/examples/legacy/`.

### 4e. Clean up resources import machinery

The 230-line `_SIRENResourcesMetaPathImporter` in `resources.py` is
used for one call site (DarkNews table saving). Replace with an
explicit function and remove the importer.

---

## Phase 5: Reweight API (future)

```python
results = sim.run()
new_weights = sim.reweight(mu_tr_mu4=1e-6)
```

Reweight existing events to new physical parameters without
re-generating. May require C++ Weighter changes to accept new
physical distributions after construction.
