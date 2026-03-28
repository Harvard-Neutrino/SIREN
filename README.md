[![Build](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml)

# SIREN

SIREN (**S**ampling and **I**njection for **R**are **E**ve**N**ts) is a framework for injecting and weighting interaction final states of complex topology, with specific concern for the detector geometry. SIREN is designed to support a wide variety of neutrino experimental setups, including atmospheric neutrinos, accelerator beam decay-in-flight neutrinos, and neutrinos from decay-at-rest sources. SIREN grew out of [LeptonInjector](https://github.com/icecube/LeptonInjector), a neutrino injection code developed within the IceCube collaboration to study atmospheric and astrophysical neutrino interactions in the IceCube detector.

SIREN provides a generic interface for user-defined BSM processes (and includes several pre-defined processes). It also supports generation of any number of secondary processes, e.g. the decay of a BSM particle after it has been created by an initial process.

## Quick start

```bash
pip install siren
```

See [Installation](#installation) for other options, including building from source or as a C++ library.

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

More examples — including BSM dipole-portal injection and MARLEY low-energy interactions — are in [`resources/examples/`](resources/examples/).

## How it works

A SIREN workflow has two phases: **injection** and **weighting**.

During **injection**, SIREN samples interaction vertices inside a detector geometry according to user-specified distributions (energy spectrum, direction, position). Each call to `injector.generate_event()` produces an interaction tree — a primary interaction and any secondary processes (e.g. the decay of a BSM particle produced in the primary interaction).

During **weighting**, SIREN computes a physical weight for each injected event. The `Weighter` takes the injection configuration and a set of *physical* distributions (the true flux and cross sections) and returns a weight that corrects for the difference between the injection and physical distributions. This importance-sampling approach allows a single injection run to be reweighted against different physical models without regenerating events.

## Supported detectors

SIREN includes detector geometry definitions for the following experiments:

| Experiment | Model name |
|------------|------------|
| [ATLAS](https://atlas.cern/) | `ATLAS` |
| [CCM](https://ccm.mit.edu/) | `CCM` |
| [DUNE Far Detector](https://www.dunescience.org/) | `DUNEFD` |
| [Hyper-Kamiokande](https://www.hyperk.org/) | `HyperK` |
| [IceCube](https://icecube.wisc.edu/) | `IceCube` |
| [MINERvA](https://minerva.fnal.gov/) | `MINERvA` |
| [MiniBooNE](https://www-boone.fnal.gov/) | `MiniBooNE` |
| [ND280 Upgrade](https://t2k-experiment.org/) | `ND280UPGRD` |

Each detector is defined by a materials file and a density profile. To load one:

```python
detector_model = siren.utilities.load_detector("IceCube")
```

Contributions of new detector geometries are welcome.

## Supported process models

| Model | Description |
|-------|-------------|
| `CSMSDISSplines` | Deep inelastic scattering (CC and NC) on nucleons, using photospline cross-section tables |
| `MarleyCrossSection` | Low-energy neutrino interactions via [MARLEY](https://www.marleygen.org/) |
| `DarkNewsTables` | BSM processes (dark photons, dipole portal, HNLs) via [DarkNews](https://github.com/mhostert/DarkNews-generator) — see [example2](resources/examples/example2/) |

## Supported flux models

| Model | Description |
|-------|-------------|
| `BNB` | Booster Neutrino Beam (FHC and RHC modes) |
| `NUMI` | NuMI beamline (low-energy and medium-energy) |
| `T2K_NEAR` | T2K near detector flux |
| `HE_SN` | Supernova neutrino flux |

To load a flux model:

```python
flux = siren.utilities.load_flux("BNB", tag="FHC_numu")
```

## Installation

### Python (pip)

```bash
pip install siren
```

Requires Python >= 3.8. Dependencies (`numpy`, `scipy`, `awkward`, `pyarrow`, `h5py`) are installed automatically.

For optional BSM support via DarkNews:

```bash
pip install siren DarkNews>=0.4.2
```

### Python (from source)

```bash
git clone https://github.com/Harvard-Neutrino/SIREN.git
cd SIREN
pip install . --config-settings='build-dir=build'
```

### C++ library

SIREN can also be built and installed as a standalone C++ shared library using CMake. This is useful when integrating SIREN into a larger C++ project.

#### Prerequisites

* A C++ compiler with C++17 support
* CMake >= 3.20
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
* [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)

CFITSIO and SuiteSparse are required by photospline for reading cross-section spline tables.

#### Workspace layout

We recommend keeping source, build, and install directories separate. A typical workspace looks like:

```
workspace/
├── env.sh              # environment setup script
├── local/              # install prefix
│   ├── bin/
│   ├── include/
│   └── lib/
└── sources/
    ├── SIREN/
    │   └── build/      # out-of-source build directory
    └── (other projects that depend on SIREN)
```

#### Environment setup

Create an `env.sh` script that configures your compiler and paths. This script should be sourced before building or running anything in the workspace.

```bash
#!/bin/bash

export WORKSPACE=/path/to/workspace
export PREFIX=$WORKSPACE/local

export CC=gcc
export CXX=g++

export PATH=$PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$PREFIX/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$PREFIX/include:$CPLUS_INCLUDE_PATH
export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
PY_SITE_DIR="$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
export PYTHONPATH="$PY_SITE_DIR${PYTHONPATH:+:$PYTHONPATH}"
```

On macOS, use `DYLD_FALLBACK_LIBRARY_PATH` instead of `LD_LIBRARY_PATH`.

#### Build and install

```bash
source env.sh

git clone https://github.com/Harvard-Neutrino/SIREN.git $WORKSPACE/sources/SIREN
cd $WORKSPACE/sources/SIREN
git submodule update --init

mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX
cmake --build . --parallel && cmake --install .
```

After the initial build, only the last line (`cmake --build . --parallel && cmake --install .`) needs to be rerun when source files change. Rerun the first `cmake` step if `CMakeLists.txt` files change or new source files are added.

#### CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_INSTALL_PREFIX` | — | Install destination for libraries, headers, and binaries |
| `SIREN_PYTHON_PACKAGE` | `ON` | Build and install the Python package |
| `SIREN_WITH_MARLEY` | `ON` | Enable MARLEY support (used if found) |
| `SIREN_REQUIRE_MARLEY` | `OFF` | Fail the build if MARLEY is not found |

## Project structure

```
SIREN/
├── projects/              # C++ source modules
│   ├── utilities/         # Interpolation, random numbers, string utilities
│   ├── serialization/     # Serialization support
│   ├── math/              # Vector3D, Matrix3D, Quaternion, interpolation
│   ├── dataclasses/       # Particle, InteractionRecord, InteractionTree
│   ├── geometry/          # Geometric primitives and operations
│   ├── detector/          # Detector geometry and density profiles
│   ├── interactions/      # Cross section and decay implementations
│   ├── distributions/     # Probability distributions for sampling
│   └── injection/         # Injector and Weighter
├── python/                # Python API (Injector, Weighter, utilities)
├── vendor/                # Vendored dependencies (git submodules)
│   ├── cereal             # Serialization
│   ├── delabella          # Delaunay triangulation
│   ├── googletest         # Unit testing
│   ├── pybind11           # Python bindings
│   ├── rk                 # Relativistic kinematics
│   ├── photospline        # B-spline cross-section tables
│   └── NamedType          # Type-safe wrappers
├── resources/
│   ├── detectors/         # Detector geometry definitions
│   ├── fluxes/            # Neutrino flux models
│   ├── processes/         # Cross-section and decay data
│   └── examples/          # Example scripts and notebooks
└── cmake/Packages/        # CMake find-modules for dependencies
```

## Contributing

Create a branch off of `main` named `$GitHubUsername/$YourSubProject`. When your changes are stable, pull from main and open a pull request.

## Citing SIREN

If you use SIREN in your research, please cite the repository:

```
https://github.com/Harvard-Neutrino/SIREN
```

## License

SIREN is licensed under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html).
