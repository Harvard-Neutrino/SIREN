[![Build](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml)

# SIREN

SIREN (**S**ampling and **I**njection for **R**are **E**ve**N**ts) is a framework for injecting and weighting interaction final states of complex topology, with specific concern for the detector geometry. SIREN is designed to support a wide variety of neutrino experimental setups, including atmospheric neutrinos, accelerator beam decay-in-flight neutrinos, and neutrinos from decay-at-rest sources. SIREN grew out of [LeptonInjector](https://github.com/icecube/LeptonInjector), a neutrino injection code developed within the IceCube collaboration to study atmospheric and astrophysical neutrino interactions in the IceCube detector.

SIREN provides a generic interface for user-defined BSM processes (and includes several pre-defined processes). It also supports generation of any number of secondary processes, e.g. the decay of a BSM particle after it has been created by an initial process. SIREN also includes detector geometry definitions for a number of existing HEP experiments, although contributions are always appreciated!

## Python installation

SIREN is distributed on PyPI as `siren`, and can be installed via pip:

```bash
pip install siren
```

For development of SIREN as a Python project, clone the repository and install:

```bash
git clone https://github.com/Harvard-Neutrino/SIREN.git
cd SIREN
pip install . --config-settings='build-dir=build'
```

After installing, verify that the `siren` Python library is available:

```python
import siren
```

## Usage overview

To use SIREN, you will:

1. Define a primary process and a list of secondary processes by specifying a particle type and which interactions (cross sections or decays) each particle can undergo

2. For each (primary or secondary) process, define a set of distributions from which to sample when injecting that particle (e.g. energy, position, direction)

3. Combine this information to define an `Injector` object

4. Generate interaction trees using the `Injector` object

5. Create a `Weighter` object using a list of primary and secondary physical processes

6. Calculate the event weight for each interaction tree using the `Weighter` object

For examples, see the `resources/examples/` directory.

## Dependencies

### Required

* A C++ compiler with C++17 support
* CMake >= 3.20
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) — required by photospline for reading cross-section spline tables
* [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) — required by photospline

### Optional

* [MARLEY](https://www.marleygen.org/) — Monte Carlo event generator for low-energy neutrino interactions. Controlled by the `SIREN_WITH_MARLEY` and `SIREN_REQUIRE_MARLEY` CMake options.

### For Python bindings

* Python >= 3.8
* numpy, scipy, awkward, pyarrow, h5py (installed automatically via pip)
* [DarkNews](https://github.com/mhostert/DarkNews-generator) >= 0.4.2 (optional, for BSM processes)

### Vendored dependencies

These are included automatically as git submodules and do not need to be installed separately:

* [cereal](https://github.com/USCiLab/cereal) — serialization
* [delabella](https://github.com/msokalski/delabella) — Delaunay triangulation for interpolation
* [googletest](https://github.com/google/googletest) — unit testing
* [pybind11](https://github.com/pybind/pybind11) — Python bindings
* [rk](https://rk.hepforge.org/) — relativistic kinematics
* [photospline](https://github.com/icecube/photospline) — B-spline representations of multi-dimensional tables
* [NamedType](https://github.com/joboccara/NamedType) — type-safe wrappers

## C++ installation

SIREN can be built and installed as a C++ library using CMake. The recommended approach keeps source, build, and install directories separate.

### Workspace layout

A typical workspace looks like this:

```
workspace/
├── env.sh           # environment setup script
├── local/           # install prefix (libraries, headers, binaries)
│   ├── bin/
│   ├── include/
│   └── lib/
└── sources/
    └── SIREN/
        └── build/   # out-of-source build directory
```

### Clone and initialize submodules

```bash
git clone https://github.com/Harvard-Neutrino/SIREN.git
cd SIREN
git submodule update --init
```

### Environment setup

Before building, set up your environment so that the compiler and linker can find your installed dependencies and your install prefix. A typical `env.sh` script looks like:

```bash
#!/bin/bash

# Set your workspace root
export WORKSPACE=/path/to/workspace
export PREFIX=$WORKSPACE/local

# Compiler
export CC=gcc
export CXX=g++

# Paths
export PATH=$PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$PREFIX/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$PREFIX/include:$CPLUS_INCLUDE_PATH
export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

# Python (if using Python bindings)
export PYTHONPATH=$PREFIX/lib/python3.11/site-packages:$PYTHONPATH
```

On macOS, use `DYLD_FALLBACK_LIBRARY_PATH` instead of `LD_LIBRARY_PATH`.

Source this script before building or using SIREN:

```bash
source env.sh
```

### Build and install

Create a build directory, run CMake, and build:

```bash
mkdir -p SIREN/build && cd SIREN/build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX
make -j$(nproc) && make install
```

### CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_INSTALL_PREFIX` | — | Install destination for libraries, headers, and binaries |
| `SIREN_PYTHON_PACKAGE` | `ON` | Build and install the Python package |
| `SIREN_WITH_MARLEY` | `ON` | Enable MARLEY support (uses MARLEY if found) |
| `SIREN_REQUIRE_MARLEY` | `OFF` | Fail the build if MARLEY is not found |

### Verifying the installation

After building and installing, verify that the shared library is in place:

```bash
ls $PREFIX/lib/libSIREN.so      # Linux
ls $PREFIX/lib/libSIREN.dylib   # macOS
```

If Python bindings were built, verify the Python package:

```python
import siren
```

## Project structure

```
SIREN/
├── cmake/Packages/    # CMake find-modules for external dependencies
├── projects/          # C++ source modules
│   ├── utilities/     # Interpolation, random numbers, string utilities
│   ├── serialization/ # Serialization support
│   ├── math/          # Vector3D, Matrix3D, Quaternion, interpolation
│   ├── dataclasses/   # Core data structures (Particle, InteractionRecord, etc.)
│   ├── geometry/      # Geometric primitives and operations
│   ├── detector/      # Detector geometry definitions
│   ├── interactions/  # Cross sections and decay implementations
│   ├── distributions/ # Probability distributions for sampling
│   └── injection/     # Injection and weighting engine
├── python/            # Python wrappers and utilities
├── vendor/            # Vendored dependencies (git submodules)
└── resources/
    ├── detectors/     # Detector geometry files (CCM, IceCube, DUNE, ATLAS, etc.)
    ├── fluxes/        # Neutrino flux definitions
    ├── processes/     # Cross-section splines and decay tables
    └── examples/      # Example scripts and notebooks
```

## Making contributions

If you would like to make contributions to this project, please create a branch off of the `main` branch and name it something following the template: `$YourLastName/$YourSubProject`.
Work on this branch until you have made the changes you wished to see and your branch is stable.
Then, pull from main, and create a pull request to merge your branch back into main.
