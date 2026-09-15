[![Build](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml)

# SIREN

SIREN (**S**ampling and **I**njection for **R**are **E**ve**N**ts) generates and
weights particle interactions in detector geometries. It supports neutrino
and beyond-Standard-Model processes, including interaction chains with
secondary decays, through Python and C++ interfaces.

SIREN grew out of [LeptonInjector](https://github.com/icecube/LeptonInjector)
and supports atmospheric, accelerator, and decay-at-rest neutrino experiments.
Users can supply their own processes, fluxes, and detector geometries alongside
the models included below.

## Quick start

```bash
python -m pip install siren
```

SIREN supports Linux and macOS and requires Python >= 3.8. Python dependencies
are installed automatically. Download the [datasets](#dataset-download) needed
by your simulation before running the example. For source builds and C++ use,
see [Installation](#installation).

The following example injects 1e4 muon-neutrino DIS events in IceCube and computes their physical weights:

```python
import siren

# Load detector geometry and interaction cross sections
detector_model = siren.utilities.load_detector("IceCube")
primary_type = siren.dataclasses.Particle.ParticleType.NuMu
primary_processes, _ = siren.utilities.load_processes(
    "CSMSDISSplines",
    primary_types=[primary_type],
    target_types=[siren.dataclasses.Particle.ParticleType.Nucleon],
    isoscalar=True,
    process_types=["CC"],
)

# Configure the injector: energy spectrum, direction, and position
injector = siren.injection.Injector()
injector.number_of_events = int(1e4)
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_processes[primary_type]
injector.primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),
    siren.distributions.PowerLaw(2, 1e3, 1e6),
    siren.distributions.IsotropicDirection(),
    siren.distributions.ColumnDepthPositionDistribution(
        600, 600.0, siren.distributions.LeptonDepthFunction()
    ),
]

# Generate events
from siren._util import GenerateEvents, SaveEvents
events, gen_times = GenerateEvents(injector)

# Weight events using physical distributions
weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_processes[primary_type]
weighter.primary_physical_distributions = [
    siren.distributions.PowerLaw(2, 1e3, 1e6),
    siren.distributions.IsotropicDirection(),
]

weights = [weighter(event) for event in events]

# Save results to HDF5 and Parquet
SaveEvents(events, weighter, gen_times, output_filename="my_output")
```

More examples — including BSM dipole-portal injection and MARLEY low-energy interactions — are in [`resources/examples/`](https://github.com/Harvard-Neutrino/SIREN/tree/main/resources/examples/).

## How it works

SIREN separates **injection** from **weighting**. Injection samples vertices,
energies, directions, and secondary processes to produce an interaction tree.
Weighting corrects for the difference between those sampling distributions and
the physical flux and interaction probabilities. A generated sample can then
be reweighted for different physical models without regenerating events.

For beam-parent injection, see [beam tables](https://github.com/Harvard-Neutrino/SIREN/blob/main/docs/beam_tables.md) for dk2nu
mass resolution, coordinate transforms, source metadata, and importance weights.

For forced final states, see [decay channels and propagation](https://github.com/Harvard-Neutrino/SIREN/blob/main/docs/decay_channels.md)
for partial widths, branching fractions, lifetime weighting, and archive support.

## Installation

For optional BSM support via DarkNews:

```bash
python -m pip install "siren[DarkNews]"
```

### From source

Source builds require a C++17 compiler, CMake >= 3.20, Python development
headers/libraries, [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/), and
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html).

```bash
git clone --recurse-submodules https://github.com/Harvard-Neutrino/SIREN.git
cd SIREN
python -m pip install . --config-settings='build-dir=build-pip'
```

### Building wheels

From the source checkout, build a wheel with pip:

```bash
python -m pip wheel . --no-deps --wheel-dir dist
```

Or use CMake, selecting the interpreter that will import SIREN:

```bash
cmake -S . -B build -DSIREN_PYTHON_PACKAGE=ON \
  -DPython_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build build --target python_package --parallel
```

#### Installing the wheel

Install the CMake-built wheel explicitly with that interpreter:

```bash
python -m pip install build/dist_wheels/siren-*.whl
```

Repeating this command leaves an unchanged installation in place. After
rebuilding a wheel with the same version, reinstall it with:

```bash
python -m pip install --force-reinstall --no-deps build/dist_wheels/siren-*.whl
```

This assumes the dependencies from the first install are still satisfied.
For the pip-built wheel, use `dist/siren-*.whl` instead. The `install_wheel`
target runs that same reinstall command with the configured interpreter after
rebuilding the wheel:

```bash
cmake --build build --target install_wheel
```

`cmake --install` installs native artifacts only; it never invokes pip.
The Python interpreter and pip determine where the wheel is installed.
Before distributing locally built wheels, repair their native dependencies
with `delocate-wheel` on macOS or `auditwheel repair` on Linux, as release CI
does. See the [packaging guide](https://github.com/Harvard-Neutrino/SIREN/blob/main/docs/packaging.md) for offline builds,
platform settings, staging, and validation.

### C++ library

From the same source checkout, build and install the standalone library:

```bash
cmake -S . -B build-native -DSIREN_PYTHON_PACKAGE=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build-native --parallel
cmake --install build-native
```

Choose your destination with `CMAKE_INSTALL_PREFIX`; `DESTDIR` supports staged
native installations. Set `SIREN_PYTHON_PACKAGE=ON` if this build should also
produce a wheel, then install that wheel separately with pip.

**Use the standalone native library and Python wheel in separate processes.**
Importing SIREN checks for an already loaded standalone core where library
inspection is available; the restriction applies in either loading order.
See [native library layout and loading](https://github.com/Harvard-Neutrino/SIREN/blob/main/docs/packaging.md#native-library-layout-and-loading)
for details and CMake consumer requirements.

## HepMC3 / NuHepMC output

With HepMC3 >= 3.3 available at build time, SIREN can export weighted events
using the NuHepMC conventions:

```python
from siren._util import SaveEvents
SaveEvents(events, weighter, gen_times, output_filename="my_output", save_hepmc3=True)
```

Add `hepmc3_gzip=True` for compressed output. The [HepMC3 guide](https://github.com/Harvard-Neutrino/SIREN/blob/main/docs/hepmc3.md)
explains weight policies, deferred weighting, combining simulation sets, and
reading the output. To require this capability in a source build, configure
CMake with `SIREN_REQUIRE_HEPMC3=ON`.

## Built-in detectors

SIREN includes detector geometry definitions for the following experiments:

| Experiment | Model name |
|------------|------------|
| [ATLAS](https://atlas.cern/) | `ATLAS` |
| [CCM](https://ccm.mit.edu/) | `CCM` |
| [DUNE Far Detector](https://www.dunescience.org/) | `DUNEFD` |
| [Hyper-Kamiokande](https://www.hyperk.org/) | `HyperK` |
| [IceCube](https://icecube.wisc.edu/) | `IceCube` |
| [KM3NeT/ORCA](https://www.km3net.org/) | `KM3NeTORCA` |
| [MINERvA](https://minerva.fnal.gov/) | `MINERvA` |
| [MiniBooNE](https://www-boone.fnal.gov/) | `MiniBooNE` |
| [ND280](https://t2k-experiment.org/) | `ND280` |
| [ND280 Upgrade](https://t2k-experiment.org/) | `ND280UPGRD` |
| [SINE](https://journals.aps.org/prd/abstract/10.1103/z4f4-wdc3) | `SINE` |
| [UNDINE](https://journals.aps.org/prd/abstract/10.1103/z4f4-wdc3) | `UNDINE` |

Each detector is defined by a materials file and a density profile. To load one:

```python
detector_model = siren.utilities.load_detector("IceCube")
```

Contributions of new detector geometries are welcome.

## Built-in process models

For Python model definitions and sampler validation, see
[authoring interaction models](https://github.com/Harvard-Neutrino/SIREN/blob/main/docs/authoring.md).

| Model | Description |
|-------|-------------|
| `CSMSDISSplines` | Deep inelastic scattering (CC and NC) on nucleons, using photospline cross-section tables |
| `MarleyCrossSection` | Low-energy neutrino interactions via [MARLEY](https://www.marleygen.org/) |
| `DarkNewsTables` | BSM processes (dark photons, dipole portal, HNLs) via [DarkNews](https://github.com/mhostert/DarkNews-generator) — see [example2](https://github.com/Harvard-Neutrino/SIREN/tree/main/resources/examples/example2/) |
| `HNLDISSplines` | Heavy neutral lepton (HNL) production via neutrino neutral current deep inelastic scattering on nucleons, using photospline cross-section tables |
| `DipoleHNLDISSplines` | Heavy neutral lepton (HNL) production via neutrino dipole-portal (i.e., via a transition magnetic moment) deep inelastic scattering on nucleons, using photospline cross-section tables |

## Built-in flux models

| Model | Description |
|-------|-------------|
| `BNB` | Booster Neutrino Beam (FHC and RHC modes) |
| `NUMI` | NuMI beamline (low-energy and medium-energy) |
| `T2K_NEAR` | T2K near detector flux |
| `HE_SN` | Supernova neutrino flux |
| `Atmospheric` |  A suite of atmospheric neutrino flux models |

To load a flux model:

```python
flux = siren.utilities.load_flux("BNB", tag="FHC_numu")
```

## Dataset download

Download process tables and flux models after installing SIREN:

```bash
siren-download --processes
siren-download --flux
```

The data are stored under `resources/processes` and `resources/fluxes` in the
SIREN installation.

## Project structure

The [Python API](https://github.com/Harvard-Neutrino/SIREN/tree/main/python/), [C++ modules](https://github.com/Harvard-Neutrino/SIREN/tree/main/projects/), and
[detectors, fluxes, processes, and examples](https://github.com/Harvard-Neutrino/SIREN/tree/main/resources/) live in separate
directories. [CMake support](https://github.com/Harvard-Neutrino/SIREN/tree/main/cmake/) and [wheel tooling](https://github.com/Harvard-Neutrino/SIREN/tree/main/tools/wheels/) handle
builds and packaging.

## Contributing

Create a branch off of `main` named `$GitHubUsername/$YourSubProject`. When your changes are stable, pull from main and open a pull request.

## Citing SIREN

If you use SIREN in your research, please cite the repository:

```
https://github.com/Harvard-Neutrino/SIREN
```

## License

SIREN is licensed under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html).
