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

---

## Phase 1: Top-Level Convenience (non-breaking, additive)

**Goal:** Reduce import friction and verbosity without changing any existing code.

### 1a. Particle type aliases

Add to `siren/__init__.py`:

```python
# Particle shorthands -- the 6-deep access path is still valid
NuE      = dataclasses.ParticleType.NuE
NuMu     = dataclasses.ParticleType.NuMu
NuTau    = dataclasses.ParticleType.NuTau
NuEBar   = dataclasses.ParticleType.NuEBar
NuMuBar  = dataclasses.ParticleType.NuMuBar
NuTauBar = dataclasses.ParticleType.NuTauBar
Electron = dataclasses.ParticleType.EMinus
Muon     = dataclasses.ParticleType.MuMinus
Nucleon  = dataclasses.ParticleType.Nucleon
N4       = dataclasses.ParticleType.N4
# ... etc for all commonly used types
```

**Before:** `siren.dataclasses.Particle.ParticleType.NuMu`
**After:** `siren.NuMu`

### 1b. Distribution aliases

```python
# Distribution shorthands
PowerLaw         = distributions.PowerLaw
Monoenergetic    = distributions.Monoenergetic
IsotropicDirection = distributions.IsotropicDirection
FixedDirection   = distributions.FixedDirection
Cone             = distributions.Cone
PrimaryMass      = distributions.PrimaryMass
ColumnDepthPosition = distributions.ColumnDepthPositionDistribution
CylinderVolumePosition = distributions.CylinderVolumePositionDistribution
PointSourcePosition = distributions.PointSourcePositionDistribution
# ... etc
```

**Before:** `siren.distributions.ColumnDepthPositionDistribution`
**After:** `siren.ColumnDepthPosition`

### 1c. Expose public functions

```python
# Top-level public API
from ._util import GenerateEvents, SaveEvents, LoadEvents
load_detector  = utilities.load_detector
load_flux      = utilities.load_flux
load_processes = utilities.load_processes
get_fiducial_volume = utilities.get_fiducial_volume
```

**Before:** `from siren._util import GenerateEvents, SaveEvents`
**After:** `siren.GenerateEvents(injector)` or `from siren import GenerateEvents`

### 1d. Normalize `load_processes` return type

Wrap all process loaders to return a consistent `ProcessBundle` named tuple:

```python
from collections import namedtuple

ProcessBundle = namedtuple("ProcessBundle", [
    "primary",       # dict: ParticleType -> [CrossSection/Decay]
    "secondary",     # dict: ParticleType -> [CrossSection/Decay]  (empty dict if none)
])
```

The `DarkNewsTables` loader currently returns 4 values. Fold the extra keys
into the bundle or make them accessible as attributes:

```python
class ProcessBundle(NamedTuple):
    primary: dict
    secondary: dict
    # optional metadata
    primary_keys: list = []
    secondary_keys: list = []
```

Existing 2-value unpacking `primary, secondary = load_processes(...)` still
works because NamedTuple supports positional unpacking.

### Files changed

- `python/__init__.py` -- add aliases, expose functions
- `python/_util.py` -- wrap `load_processes` return in `ProcessBundle`
- All existing code continues to work unchanged.

---

## Phase 2: The `Simulation` Class

**Goal:** A single entry point that owns the full inject-weight pipeline and
eliminates all redundancy.

### API Design

```python
class Simulation:
    def __init__(
        self,
        *,  # keyword-only after this
        n_events: int,
        detector: Union[str, DetectorModel],
        primary: ParticleType,
        interactions: Union[str, list],

        # Distributions -- used for both injection and physical by default
        energy: Distribution = None,
        direction: Distribution = None,
        position: Distribution = None,  # injection only (auto-excluded from physical)
        flux: Distribution = None,      # physical energy (if different from injection energy)

        # Overrides when injection != physical
        injection_energy: Distribution = None,
        physical_energy: Distribution = None,
        injection_direction: Distribution = None,
        physical_direction: Distribution = None,

        # Optional
        seed: int = None,
        target: ParticleType = None,
        process: str = None,            # "CC", "NC", etc
        isoscalar: bool = True,

        # BSM / secondary processes
        secondary_interactions: dict = None,
        secondary_position: dict = None,
        stopping_condition: callable = None,

        # DarkNews shorthand
        darknews_model: dict = None,
    ):
```

### Distribution Resolution Rules

The constructor resolves injection and physical distributions with these rules:

1. **energy/direction**: If only `energy` is set, it is used for both injection
   and physical. If `injection_energy` and `physical_energy` are both set, they
   take precedence. If `flux` is set, it becomes `physical_energy`; if `energy`
   is also set, it becomes `injection_energy`; if `energy` is not set,
   `injection_energy` is derived from `flux` (same shape, unnormalized).

2. **position**: Always injection-only. Never included in physical
   distributions.

3. **mass**: Auto-inferred from `primary` particle type. User never specifies
   `PrimaryMass(0)` manually.

4. If any required distribution is missing after resolution, raise a
   `ValueError` with a clear message listing what is needed.

### `Results` Object

```python
class Results:
    def __init__(self, events, weights, gen_times, weighter, injector):
        self.events = events
        self.weights = weights
        self.gen_times = gen_times
        self._weighter = weighter
        self._injector = injector

    def __len__(self):
        return len(self.events)

    def __iter__(self):
        return zip(self.events, self.weights)

    def save(self, filename, **kwargs):
        """Save to HDF5/Parquet."""
        SaveEvents(self.events, self._weighter, self.gen_times,
                   output_filename=filename, **kwargs)
```

### `run()` Method

```python
def run(self) -> Results:
    """Generate all events and compute weights."""
    injector = self._build_injector()
    events, gen_times = GenerateEvents(injector)
    weighter = self._build_weighter(injector)
    weights = [weighter(event) for event in events]
    return Results(events, weights, gen_times, weighter, injector)
```

### Before/After: IceCube DIS

**Before (62 lines):**
```python
import os
import siren
from siren._util import GenerateEvents, SaveEvents

events_to_inject = int(1e5)
experiment = "IceCube"
detector_model = siren.utilities.load_detector(experiment)
primary_type = siren.dataclasses.Particle.ParticleType.NuMu
cross_section_model = "CSMSDISSplines"
primary_processes, _ = siren.utilities.load_processes(
    cross_section_model,
    primary_types=[primary_type],
    target_types=[siren.dataclasses.Particle.ParticleType.Nucleon],
    isoscalar=True,
    process_types=["CC"]
)
primary_cross_sections = primary_processes[primary_type]
injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_cross_sections
injector.primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),
    siren.distributions.PowerLaw(2, 1e3, 1e6),
    siren.distributions.IsotropicDirection(),
    siren.distributions.ColumnDepthPositionDistribution(
        600, 600.0, siren.distributions.LeptonDepthFunction())
]
events, gen_times = GenerateEvents(injector)
weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_cross_sections
weighter.primary_physical_distributions = [
    siren.distributions.PowerLaw(2, 1e3, 1e6),
    siren.distributions.IsotropicDirection()
]
os.makedirs("output", exist_ok=True)
SaveEvents(events, weighter, gen_times, output_filename="output/IceCube")
```

**After (14 lines):**
```python
import siren

sim = siren.Simulation(
    n_events=100_000,
    detector="IceCube",
    primary=siren.NuMu,
    interactions="CSMSDISSplines",
    target=siren.Nucleon,
    process="CC",
    energy=siren.PowerLaw(2, 1e3, 1e6),
    direction=siren.IsotropicDirection(),
    position=siren.ColumnDepthPosition(600, 600.0),
)

results = sim.run()
results.save("output/IceCube")
```

### Before/After: ATLAS with Flux

**Before (99 lines):** see `DIS_ATLAS.py`

**After (~20 lines):**
```python
import siren

detector_model = siren.load_detector("ATLAS")
fiducial = siren.get_fiducial_volume("ATLAS", sector="tilecal")

sim = siren.Simulation(
    n_events=100_000,
    seed=99,
    detector=detector_model,
    primary=siren.NuMu,
    interactions="CSMSDISSplines",
    target=siren.Nucleon,
    process="CC",
    flux=siren.load_flux("HE_SN", tag="numu", min_energy=100, max_energy=1e6),
    direction=siren.FixedDirection([0, 0, 1]),
    position=siren.CylinderVolumePosition(fiducial),
)

results = sim.run()
results.save("output/ATLAS")
```

### Before/After: CCM Dipole Portal (BSM)

**Before (136 lines):** see `DipolePortal_CCM.py`

**After (~35 lines):**
```python
import numpy as np
import siren

detector_model = siren.load_detector("CCM")
fiducial = siren.get_fiducial_volume("CCM")

# Source geometry
target_origin = siren.math.Vector3D(0, 0, -0.241)
detector_origin = siren.math.Vector3D(23, 0, -0.65)
beam_dir = detector_origin - target_origin
beam_dir.normalize()

sim = siren.Simulation(
    n_events=1000,
    detector=detector_model,
    primary=siren.NuMu,
    darknews_model={
        "m4": 0.0235,
        "mu_tr_mu4": 6e-7,
        "UD4": 0, "Umu4": 0,
        "epsilon": 0.0, "gD": 0.0,
        "decay_product": "photon",
        "noHC": True, "HNLtype": "dirac",
    },
    energy=siren.Monoenergetic(0.02965),
    injection_direction=siren.Cone(beam_dir, np.arctan(5 / 23.0)),
    physical_direction=siren.IsotropicDirection(),
    position=siren.PointSourcePosition(target_origin - detector_origin, 25),
    secondary_position=siren.SecondaryBoundedVertex(fiducial, 25),
    stopping_condition=lambda datum, i: (
        datum.record.signature.secondary_types[i] != siren.N4
    ),
)

results = sim.run()
results.save("output/CCM_Dipole")
```

### Files changed

- `python/Simulation.py` (new) -- ~300 lines
- `python/Results.py` (new) -- ~80 lines
- `python/__init__.py` -- import and expose `Simulation`, `Results`

---

## Phase 3: Ergonomic Improvements

### 3a. Iterator protocol on Injector

```python
# In Injector.py
def __iter__(self):
    if self.__injector is None:
        self.__initialize_injector()
    for _ in range(self.number_of_events):
        yield self.generate_event()

def __len__(self):
    return self.number_of_events
```

**Before:** `events = [injector.generate_event() for _ in range(N)]`
**After:** `events = list(injector)`

### 3b. Batch weighting on Weighter

```python
# In Weighter.py
def weight_all(self, events):
    """Compute weights for a list of events."""
    return [self(event) for event in events]
```

**Before:** `weights = [weighter(event) for event in events]`
**After:** `weights = weighter.weight_all(events)`

### 3c. `FixedDirection` from list/tuple

Currently requires `siren.math.Vector3D`. Add convenience:

```python
# In a wrapper or __init__.py
_orig_FixedDirection = distributions.FixedDirection
def FixedDirection(direction):
    if isinstance(direction, (list, tuple)):
        direction = math.Vector3D(*direction)
    return _orig_FixedDirection(direction)
```

**Before:** `siren.distributions.FixedDirection(siren.math.Vector3D(0, 0, 1))`
**After:** `siren.FixedDirection([0, 0, 1])`

### 3d. `CylinderVolumePosition` from sector

Add a helper that extracts sector geometry automatically:

```python
def CylinderVolumePosition(source):
    """Accept a Cylinder, a DetectorSector, or a fiducial volume geometry."""
    if isinstance(source, geometry.Cylinder):
        return distributions.CylinderVolumePositionDistribution(source)
    # ... handle sector extraction (the 10-line boilerplate from DIS_ATLAS.py)
```

### 3e. `get_fiducial_volume` improvements

Add optional `sector` kwarg to `get_fiducial_volume`:

```python
def get_fiducial_volume(experiment, sector=None):
    """Return the fiducial volume geometry.

    If sector is specified, extract the geometry from that named sector
    of the detector model instead of using the default fiducial volume.
    """
```

### 3f. Improve validation messages

Add a `_validate_distributions` function called by both `Simulation.__init__`
and `Injector.__initialize_injector`:

```python
def _validate_distributions(distributions, context="injection"):
    """Check that required distribution types are present."""
    types_present = {type(d).__name__ for d in distributions}
    has_energy = any(isinstance(d, (PowerLaw, Monoenergetic, ...)) for d in distributions)
    has_direction = any(isinstance(d, (IsotropicDirection, ...)) for d in distributions)
    has_position = any(isinstance(d, (...)) for d in distributions)

    missing = []
    if not has_energy:
        missing.append("energy (e.g., PowerLaw, Monoenergetic, TabulatedFluxDistribution)")
    if not has_direction:
        missing.append("direction (e.g., IsotropicDirection, FixedDirection, Cone)")
    if context == "injection" and not has_position:
        missing.append("position (e.g., ColumnDepthPosition, CylinderVolumePosition)")

    if missing:
        raise ValueError(
            f"Missing {context} distributions:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )
```

### Files changed

- `python/Injector.py` -- add `__iter__`, `__len__`
- `python/Weighter.py` -- add `weight_all`
- `python/__init__.py` -- add wrapper constructors
- `python/_validation.py` (new) -- validation helpers

---

## Phase 4: Deprecation and Cleanup

### 4a. Deprecate `SIREN_Controller`

Add deprecation warning:

```python
class SIREN_Controller:
    def __init__(self, *args, **kwargs):
        import warnings
        warnings.warn(
            "SIREN_Controller is deprecated. Use siren.Simulation instead. "
            "See examples/ for migration guide.",
            DeprecationWarning,
            stacklevel=2,
        )
        # ... existing code
```

Keep it working for at least one major version cycle.

### 4b. Deduplicate code

- Move `SaveEvents` to a single canonical location (e.g., `python/io.py`)
- `_util.py` and `SIREN_Controller.py` both import from `io.py`
- Same for `GetDetectorModelTargets`, `GetFiducialVolume`

### 4c. Clean up resources import machinery

The 230-line meta path importer in `resources.py` can be simplified. The
dynamic attribute access (`siren.resources.processes.DarkNewsTables`) is used
in exactly one place (the DarkNews table saving). Replace with an explicit
function:

```python
siren.save_darknews_tables(table_dir, processes)
```

### 4d. Update all examples

Rewrite every example in `resources/examples/` to use the new `Simulation` API.
Keep old examples in `resources/examples/legacy/` for reference during the
deprecation period.

---

## Phase 5: Weighter Auto-Construction (Future)

Enable `Simulation.reweight(**new_params)` for parameter scans:

```python
# Original simulation
sim = siren.Simulation(
    detector="CCM",
    primary=siren.NuMu,
    darknews_model={"m4": 0.0235, "mu_tr_mu4": 6e-7, ...},
    ...
)
results = sim.run()

# Reweight to a different coupling without re-generating
new_weights = sim.reweight(mu_tr_mu4=1e-6)
```

This requires the Weighter to accept new physical distributions after
construction, which may need C++ changes to support.

---

## Implementation Priority

| Phase | Effort | Breaking? | Impact |
|-------|--------|-----------|--------|
| 1: Top-level aliases | Small (1-2 days) | No | Medium -- reduces verbosity immediately |
| 2: Simulation class | Medium (3-5 days) | No | High -- the main user-facing improvement |
| 3: Ergonomics | Small (1-2 days) | No | Medium -- quality-of-life |
| 4: Deprecation | Medium (2-3 days) | Soft | Low immediate, high long-term (code health) |
| 5: Reweighting | Large (1+ weeks) | No | High for BSM users doing scans |

Phases 1-3 are purely additive Python changes. No C++ modifications required.
Phase 4 is soft-breaking (deprecation warnings only). Phase 5 may require C++
changes.

---

## Backward Compatibility

- All existing imports (`siren.injection.Injector`, `siren.distributions.*`,
  `siren._util.GenerateEvents`) continue to work.
- `SIREN_Controller` continues to work with a deprecation warning.
- The `Simulation` class is purely additive -- it wraps the existing
  Injector/Weighter internally.
- New aliases (`siren.NuMu`, `siren.PowerLaw`) do not shadow any existing
  names.
