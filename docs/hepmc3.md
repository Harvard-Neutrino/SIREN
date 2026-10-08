# HepMC3 and NuHepMC output

SIREN can write events as [HepMC3](https://gitlab.cern.ch/hepmc/HepMC3) ASCII, following the [NuHepMC](https://github.com/NuHepMC/Spec) conventions for neutrino generators, so events can be read by any HepMC3-aware tool. This requires HepMC3 >= 3.3 (the writer uses attribute types absent from 3.2.x); with an older or missing HepMC3 the package still builds and every other output format works, but the HepMC3 path raises at call time.

## Enabling it

Pass `save_hepmc3=True` to `Results.save` or to `SaveEvents` (or `hepmc3=True` to the deprecated `SIREN_Controller.SaveEvents`). This writes `<output_filename>.hepmc3` alongside the usual HDF5/Parquet/`.siren_events` files. Add `hepmc3_gzip=True` to gzip the output (the `.gz` suffix is added automatically).

```python
from siren._util import SaveEvents
SaveEvents(events, weighter, gen_times, output_filename="my_output", save_hepmc3=True)
```

## Weight policy

The central-value weight written into each event's header is controlled by `hepmc3_weights` (default `"auto"`), and the resolved policy is recorded in the run-level `siren.weights_state` attribute:

| `hepmc3_weights` | `siren.weights_state` | Meaning |
|------------------|-----------------------|---------|
| `"auto"` (default) | `"computed"` | If a weighter is available, every event weight is computed once and written to the headers. |
| `"auto"` (no weighter) | `"header"` or `"unweighted"` | Existing header weights are trusted if present, else headers are left untouched. |
| a `Weighter` or callable | `"computed"` | Weights are computed the same as `"auto"`-with-weighter. |
| a sequence of floats | `"header"` | One central value per event, written to the headers (length must equal `len(events)`). |
| `"header"` | `"header"` | Existing header weights are trusted as-is. |
| `"none"` / `None` | `"unweighted"` | Headers are left untouched; no weights are claimed. |

An `"unweighted"` export is a **plain HepMC3 file**: it carries the `siren.*` provenance attributes but **no `NuHepMC.*` attributes** (no flux-averaged cross section, no status registries), because those metadata are only meaningful for weighted output.

## Workflows

**1. Weight at save time.** Give `SaveEvents` a weighter, as in the example above, to compute and write central-value weights in one pass.

**2. Generate now, weight later.** Save without a weighter, then convert the native `.siren_events` file with the converter, passing the (in-memory) `weighter` object to recompute the CV. Reweighting reads the native file, **not** a HepMC3 round trip (the HepMC3 reader is run-level-lossy by design). `convert_siren_events_to_hepmc3` does not read a `.siren_weighter` companion file itself -- you build or load the `Weighter` and pass it in directly:

```python
from siren._util import convert_siren_events_to_hepmc3
convert_siren_events_to_hepmc3("my_output.siren_events", weighter=weighter)
```

or from the shell, without weights (the CLI has no `weighter` option, so the CVs are left as whatever the native file already carries):

```bash
python -m siren.hepmc3_convert my_output.siren_events -o my_output.hepmc3
```

**3. Factorized analysis.** The `Weighter` exposes the per-vertex probability factors behind the central value — `interaction_probabilities`, `survival_probabilities`, and the finer generation/physical terms on the bound `PrimaryProcessWeighter` / `SecondaryProcessWeighter` — so downstream analyses can recombine the partial factors themselves instead of consuming a single scalar.

**4. Combining simulation sets.** To pool several independent runs into one correctly-normalized HepMC3 file, use `combine_and_export_hepmc3` on the native `.siren_events` + `.siren_weighter` pairs. Pooling rebuilds one `Weighter` over the union of injectors and recomputes every weight; there is no valid shortcut through per-set stored weights, and combining via HepMC3 round-trips is **not** supported:

```python
from siren._util import combine_and_export_hepmc3
combine_and_export_hepmc3(["run_a", "run_b"], out_path="combined.hepmc3")
```

Each set entry is a shared base path (`<base>.siren_events` + `<base>.siren_weighter`), an `(events, weighter)` path pair, or a dict; every set needs a usable `.siren_weighter` and must share the same physical configuration.

## Reading the output

Any HepMC3 reader works. To read back into SIREN interaction trees:

```python
from siren._util import LoadEventsFromHepMC3
events = LoadEventsFromHepMC3("my_output.hepmc3")
```

Or with [pyhepmc](https://github.com/scikit-hep/pyhepmc) for a generator-agnostic view:

```python
import pyhepmc
with pyhepmc.open("my_output.hepmc3") as f:
    events = [evt for evt in f]
weights_state = str(events[0].run_info.attributes["siren.weights_state"])
central_value = events[0].weights[0]
```

See the [packaging guide](packaging.md) for builds that require HepMC3 support
and for wheel acceptance checks.
