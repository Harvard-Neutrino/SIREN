[![Build](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/Harvard-Neutrino/SIREN/actions/workflows/build_wheels.yml)

# SIREN

SIREN is a framework for injecting and weighting interaction final states of complex topology, with specific concern for the detector geometry. SIREN is designed to support a wide variety of neutrino experimental setups, including atmospheric neutrinos, accelerator beam decay-in-flight neutrinos, and neutrinos from decay-at-rest sources. The original [SIREN-v1](https://github.com/icecube/SIREN) was developed within the IceCube collaboration to study atmospheric and astrophysical neutrino interactions in the IceCube detector.

This iteration of SIREN provides a generic interface for user-defined BSM processes (with support for a few pre-defined processes). It also supports generation of any number of secondary processes, e.g. the decay of a BSM particle after it has been created by an initial process. SIREN also includes the geometric setup for a number of existing HEP experiments out of the box, although contributions are always appreciated!

# Python installation of SIREN

SIREN is distributed on pypi as `siren`, and can be installed via pip with:

```
pip install siren
```

For development of SIREN as a python project, simply clone the repository:

```
git clone https://github.com/Harvard-Neutrino/SIREN.git
cd SIREN
```

and run the following command to build and install SIREN:

```
pip install . --config-settings='build-dir=build'
```

After the python bindings are installed, you should be able to import the `siren` python library. Open a python interpreter by running `python`, and then run

```
import siren
```

To use SIREN, you will

1. Define a primary process and a list of secondary processes by specifying a particle type and which interactions (cross sections or decays) each particle can undergo

2. For each (primary or secondary) process, define a set of distributions from which to sample when injecting that particle (e.g. energy, position, direction)

3. Combine this information to define an InjectorBase object

4. Generate interaction trees using the InjectorBase object

5. Create a LeptonTreeWeighter object using a list of primary and secondary physical processes

6. Calculate the event weight for each interaction tree using the LeptonTreeWeighter object

For an example of this in action, see `resources/DipoleInjection/inject_HNLs_CCM.{py,ipynb}`

# Dependencies

For local installations, you need the following:

* A C++ compiler with C++14 support.

* Some classes also require Photospline to create and to read cross sections. Read more about it, and its installation at [github.com/icecube/photospline](https://github.com/icecube/photospline). Note that Photospline has dependencies that you will need that are not listed here.

* SIREN requires Photospline's `SuiteSparse` capabilities, whose dependencies are available [here](http://faculty.cse.tamu.edu/davis/suitesparse.html).

For building py-bindings,

* Python > 3.7

* That's it! We use pybind11 to generate our pybindings, which is automatically included in SIREN as a submodule


# Included External Dependencies

These are not ostensibly a part of SIREN, but are included automatically as submodules for its functionality.

* [cereal](https://github.com/USCiLab/cereal): for serialization

* [delabella](https://github.com/msokalski/delabella): for Delaunay triangulation in our interpolation classes

* [googletest](https://github.com/google/googletest): for constructing our tests

* [pybind11](https://github.com/pybind/pybind11): for compiling our python bindings

* [rk](https://rk.hepforge.org/): a relativistic kinematics library used mostly in the CrossSection and Decay subclasses

* [photospline](https://github.com/icecube/photospline): a library that uses the penalized spline technique to efficiently compute, store, and evaluate B-spline representations of such large multi-dimensional tables

# C++ installation of SIREN

To use SIREN in a C++ project, there are a few more steps.

We will be trying to keep our source, build, and install directories separate. To this end, these instructions will assume the following directory structure:

```
| local (for installing built libraries, headers, and binaries)
|  --lib
|  --include
|  --bin
| source (source code for SIREN and other dependencies)
|  --SIREN
|     --build (for building SIREN)
|  --(SIREN dependencies...)
```

`git clone git@github.com:Harvard-Neutrino/SIREN.git`

or

`git clone https://github.com/Harvard-Neutrino/SIREN.git`

to download the source code. To download the submodules, run

`git submodule update --init`

Now `cd SIREN/build` to get to the build directory. We call cmake

`cmake ../ -DCMAKE_INSTALL_PREFIX=../../local`

This tells cmake to install the shared objects in the `local` directory. CMake prepares a `Makefile` which calls the `g++` compiler with the necessary instructions to compile. So now you'll call

`make -j4 && make install`

to build the project and install the project. Now you need to set all the environmental variables so this actually works. We recommend putting the followig commands into a `env.sh` script that can load the environment.

```
export PROJECTSPACE=/path/to/parent/directory
export PROJECTBUILDPATH=$PROJECTSPACE/local
export PROJECTSOURCEPATH=$PROJECTSPACE/source
export PREFIX=$PROJECTBUILDPATH
# On linux:
export LD_LIBRARY_PATH=$PROJECTBUILDPATH/lib/:$LD_LIBRARY_PATH
# On mac:
export DYLD_FALLBACK_LIBRARY_PATH=$PROJECTBUILDPATH/lib/:$DYLD_FALLBACK_LIBRARY_PATH
```

Now you should be good to go!

# Making Contributions
If you would like to make contributions to this project, please create a branch off of the `main` branch and name it something following the template: `$YourLastName/$YourSubProject`.
Work on this branch until you have made the changes you wished to see and your branch is stable.
Then, pull from main, and create a pull request to merge your branch back into main.
