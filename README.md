# LeptonInjector

LeptonInjector is a group of modules used to create events in IceCube. This code represents a standalone version of the original LeptonInjector which has been trimmed of the proprietary Icetray dependencies. It is currently a work in progress, and is *not* in a functional state! 

When finished, this code is intented to be used as a series of modules put together in a python script. 

To use it, you will

    1. Prepare a Injector object (or list of Injector objects). These are represented as `MinimalInjectionConfiguration` objects in the c++ code. 

    2. Use these, along with several generation parameters, to create a Controller object. These were called MultiLeptonInjector in the original code. 

    3. Add/remove any injectors to the controller object. Verify that the controller is properly configured.

    4. Call Controller.Execute(). This will run the simulation. 

# Dependencies

All of the dependencies are already installed on the CVMFS environments on the IceCube Cobalt testbeds. For local installations, you may need to install the below. You only 'may' need to, because you might already have these. All three are definitely necessary. 

This code requires the `HDF5` C libraries. Read more about it here: https://portal.hdfgroup.org/display/support.  

It requires `BOOST`, which can be installed as easily as typing `sudo apt-get install libboost-all-dev` on linux machines. There's probably a homebrew version for mac. 

It also requires Photospline to create and to read cross sections. Read more about it, and its installation at https://github.com/IceCubeOpenSource/photospline

# Included Dependencies

* I3CrossSections: 

* Earthmodel Services: 

# Compilation and Installation

This project uses cmake to compile. First, go to a folder where you would like to build and compile lepton injector. Run 

`git clone git@github.com:IceCubeOpenSource/LeptonInjector.git`

to download the source code. Then, `mv` the folder that is created to rename it to `source`. We will be trying to keep our source, build, and install directories separate. Now

`mkdir build` and `mkdir install`

to prepare. `cd ./build` to move into the build directory. Now, we call cmake

`cmake -DCMAKE_INSTALL_PREFIX=../install ../source`

This tells cmake to install the shared objects in the `install` directory we just made. Cmake prepares a `Makefile` which calls the `g++` compiler with the necessary instructions to... compile. So now you'll call

`make -j4 && make install`

to build the project and install the project. Now you need to set all the environmental variables so this actually works. 

# Structure
The code base is divided into several files. 
* Constants: a header defining various constants 
* Controller: implements a class for managing the simulation
* DataWriter: writes event properties and MCTrees to an HDF5 file
* EventProps: implements an event properties struct for both ranged and volume injectors
* h5write: may be renamed soon. This will be used to write the configurations onto a file
* Helper: defines a few helper functions
* LeptonInjector (the file): defines the Injector objects described above in addition to several configuration objects and event generators 
* Particle: simple implementation of particles 
* Random: object for random number sampling. 

# Cross Sections
More to come on this front later. You will need cross sections, saved as fits files. Maybe we should upload examples XS to the repo?
