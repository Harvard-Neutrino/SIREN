# Building and packaging SIREN

The [README](../README.md#installation) covers first installs and basic builds.
This guide covers CMake configuration, staging, native library behavior, and
wheel release checks. Linux and macOS are supported; Windows is unsupported.

## Build requirements

Source builds need a C++17 compiler, CMake >= 3.20, Python >= 3.8 with development
headers/libraries, CFITSIO, and SuiteSparse. Initialize the repository's
submodules with `git submodule update --init --recursive`. For dependencies
installed outside the default search paths, set `CMAKE_PREFIX_PATH` to their
installation prefixes.

HepMC3 output requires HepMC3 >= 3.3. Set `SIREN_REQUIRE_HEPMC3=ON` to fail
configuration if that capability is unavailable; otherwise the package can
build without it, and HepMC3 export raises when called.

## Build and install with CMake

Use the interpreter that will import the wheel:

```bash
cmake -S . -B build -DSIREN_PYTHON_PACKAGE=ON \
  -DPython_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build build --parallel
python -m pip install build/dist_wheels/siren-*.whl
```

The `python_package` target is part of this build and can also be built on its
own. It stages the `PythonWheel` install component into a fresh tree, then
rebuilds the wheel only when the staged contents or packaging inputs change.
The CMake install rules define the payload, including the bytes behind file
symlinks. Both this target and `python -m pip wheel . --no-deps --wheel-dir dist`
use scikit-build-core.

### Repeat installs and rebuilds

A repeated pip install of the same wheel leaves the existing installation in
place. Rebuilding source without changing the package version requires an
explicit reinstall:

```bash
cmake --build build --parallel
python -m pip install --force-reinstall --no-deps build/dist_wheels/siren-*.whl
```

This replaces the previous installation, including files absent from the new
wheel. `--no-deps` assumes the dependencies are already installed and satisfied;
use the ordinary install command for the first install. The wheel's Python ABI
and platform tags must match the selected interpreter.

The `install_wheel` build target runs exactly that reinstall command with the
configured `Python_EXECUTABLE`, after bringing the wheel up to date:

```bash
cmake --build build --target install_wheel
```

It is a build target, not part of `cmake --install`, so `CMAKE_INSTALL_PREFIX`,
`DESTDIR`, and root installs never reach pip, and it never installs
dependencies.

### Native installation and staging

`cmake --install` installs the native library, headers, and CMake exports. It
never invokes pip or installs the SIREN Python package. The vendored
photospline installation retains its own Python wrapper/stubs.

```bash
cmake --install build --prefix /path/to/native-prefix
DESTDIR=/path/to/staging-root cmake --install build --prefix /usr/local
```

The second command stages native artifacts below
`/path/to/staging-root/usr/local`. `CMAKE_INSTALL_PREFIX` sets the default
prefix; `--prefix` overrides it at install time. These native destinations are
independent of the Python environment. Packagers staging a wheel should use
pip's own `--root` and `--prefix` options with the appropriate interpreter.
Installing the wheel through `cmake --install` is unsupported; use pip or the
`install_wheel` build target.

### CMake options

| Option | Default | Purpose |
| --- | --- | --- |
| `CMAKE_INSTALL_PREFIX` | CMake platform default | Native install destination |
| `Python_EXECUTABLE` | Discovered | Interpreter used by the CMake wheel target |
| `SIREN_PYTHON_PACKAGE` | `ON` | Build the wheel; install it separately with pip |
| `SIREN_WHEEL_BUILD_ISOLATION` | `ON` | Provision the wheel backend in an isolated environment |
| `SIREN_WITH_HEPMC3` | `ON` | Enable HepMC3 output if available |
| `SIREN_REQUIRE_HEPMC3` | `OFF` | Fail configuration if HepMC3 support is unavailable |
| `SIREN_WITH_MARLEY` | `ON` | Enable MARLEY support if available |
| `SIREN_REQUIRE_MARLEY` | `OFF` | Fail configuration if MARLEY is unavailable |

For standalone C++ development, `SIREN_PYTHON_PACKAGE=OFF` omits SIREN's Python
extensions and Python core from a plain CMake build. A source-wheel build still
builds those targets, without starting a nested wheel build.

### macOS configuration

Set `CMAKE_OSX_DEPLOYMENT_TARGET` and, when needed, `CMAKE_OSX_ARCHITECTURES` in
the outer CMake build. The wheel target forwards them to its metadata backend.
Use an interpreter and architecture compatible with the intended runtime;
repairing a wheel does not lower the deployment targets of its binaries or
dependencies.

### Offline builds

Provision the configured interpreter with the requirements in
[`pyproject.toml`'s `[build-system]`](../pyproject.toml) and their dependencies
before configuring `SIREN_WHEEL_BUILD_ISOLATION=OFF`. For a direct source-wheel
build, use:

```bash
python -m pip wheel . --no-build-isolation --no-deps --wheel-dir dist
```

Native dependencies must also be available locally.

## Native library layout and loading

Wheels contain the Python package, resources, extensions, and native runtime
libraries, with native ABI/platform tags and relative library lookup paths.
C++ headers and CMake exports come from the standalone native install. CMake
consumers must provide the dependency targets required by the exported targets.

The wheel core is `libSIREN_python.dylib` on macOS or `libSIREN_python.so` on
Linux. The standalone core is `libSIREN`, which links Python for native
applications, including plugin hosts using `dlopen`. The wheel resolves Python
symbols through its host interpreter and must not bundle another Python runtime.
The two cores share compiled components but have distinct names to prevent a
native installation on the loader search path from replacing the wheel core.

Use the standalone native library and Python wheel in **separate processes**.
Loading both duplicates internal state, including particle-ID allocation, even
when they come from the same build. Importing `siren` raises `ImportError` when
it detects a loaded standalone core. This guard is best effort: the restriction
also applies when inspection is unavailable, such as Linux without `/proc`,
and when the native core would be loaded after the wheel. Importing extension
files directly does not make mixed-core use safe.

## Repair and validate release wheels

Local wheels may depend on external native libraries such as CFITSIO and
HepMC3. Before distribution, bundle and repair those dependencies with
`delocate-wheel` on macOS or `auditwheel repair` on Linux. The
[cibuildwheel configuration](../pyproject.toml) performs this for release wheels.

Validate from a fresh virtual environment outside the checkout. Copy
[`test_hepmc3_wheel.py`](../tools/wheels/test_hepmc3_wheel.py) outside the source
tree before making the original source/build directories unavailable. Install
the repaired wheel and run the copied script:

```bash
python -m pip install /path/to/repaired/siren-*.whl packaging
python /path/to/test_hepmc3_wheel.py --clean-environment \
  --foreign-prefix /original/source/path \
  --foreign-prefix /original/dependency/prefix
```

Build with `SIREN_REQUIRE_HEPMC3=ON` for this check. `--clean-environment` runs a
child interpreter with Python and library-search overrides cleared.
`--foreign-prefix` is repeatable and fails if any loaded library comes from
under a named directory; release CI names the checkout and dependency prefix.

The script requires loaded-library inspection. It checks wheel tags, rejects
multiple core files and bundled Python runtimes, and requires the loaded SIREN
core and photospline/spglam to belong to the installed wheel. Other loaded
libraries sharing a bundled dependency's name but resolving outside the wheel
are listed for inspection. That list should be empty in the clean acceptance
environment. The script also exercises native sampling and plain/gzip HepMC3
round trips.

Five Linux release configurations skip runtime import smoke because suitable
NumPy/SciPy binaries are unavailable: `cp38-musllinux*`, `cp39-musllinux*`, and
`cp38-manylinux_aarch64`. Their wheels still build. See `test-skip` in
[`pyproject.toml`](../pyproject.toml) for the current list.

### Diagnosing a competing native installation

After clean acceptance passes, the optional `--standalone-library` diagnostic
checks a wheel with a native installation on its library search path:

```bash
python /path/to/test_hepmc3_wheel.py --standalone-library /native/prefix/lib/libSIREN.dylib
```

On Linux, use the `.so` path, including `lib64` if appropriate. The diagnostic
sets `DYLD_LIBRARY_PATH` or `LD_LIBRARY_PATH` in a fresh interpreter and verifies
that the override arrived. It does not load the standalone core. A competing
prefix can substitute native photospline or another bundled dependency and
cause the diagnostic to fail; this is optional troubleshooting, not the clean
acceptance gate. Clear such overrides before using the wheel normally.
