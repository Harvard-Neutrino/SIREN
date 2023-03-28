# LeptonInjector

LeptonInjector is a group of modules used to efficiently inject events in a given detector setup. It was originally developed within the IceCube collaboration to study astrophysical neutrino interactions in the IceCube detector. However, this project has evolved significantly since; it is now capable of simulating any interactions the user provides (with support for a few pre-defined processes). It also supports generation of any number of secondary processes, e.g. the decay of a BSM particle after it has been creaetd by an initial process. LeptonInjector also includes the geometric setup for a number of existing HEP experiments out of the box, although contributions are always appreciated!

To use LeptonInjector, you will

    1. Prepare a Injector object by defining a primary process and a list of secondary processes

    2.  

    3. Add more injectors to the controller object using the add injector function. Verify that the controller is properly configured.
    
    4. Specify the full paths and names of the destination output file and LeptonInjector Configuration (LIC) file.

    5. Call Controller.Execute(). This will run the simulation. 

For an example of this in action, see $root/resources/example/inject_muons.py

To learn about the LIC files and weighting, see https://github.com/IceCubeOpenSource/LeptonWeighter

# Dependencies

All of the dependencies are already installed on the CVMFS environments on the IceCube Cobalt testbeds. 

For local installations, you need the following:

* A C++ compiler with C++11 support.

* The `HDF5` C libraries. Read more about it here: https://portal.hdfgroup.org/display/support. These libraries are, of course, used to save the data files. 

* It also requires Photospline to create and to read cross sections. Read more about it, and its installation at https://github.com/IceCubeOpenSource/photospline. Note that Photospline has dependencies that you will need that are not listed here. 

* LeptonInjector requires Photospline's `SuiteSparse` capabilities, whose dependencies are available here http://faculty.cse.tamu.edu/davis/suitesparse.html

For building py-bindings, 

* Python, but you should've known that if you're building pybindings. 

* `BOOST`, which can be installed as easily as typing `sudo apt-get install libboost-all-dev` on linux machines (so long as you use bash and not something like tcsh). There's probably a homebrew version for mac. Boost is primarily needed to compile the python bindings or this software. 


# Included Dependencies

These are not ostensibly a part of LeptonInjector, but are included for its functionality. They were developed by the IceCube Collaboration and modified slightly to use the LeptonInjector datatypes instead of the IceCube proprietary ones. 

* I3CrossSections: provies the tools for sampling DIS and GR cross sections. 

* Earthmodel Services: provides the PREM for column depth calculations. 

# Download, Compilation and Installation

First, go to a folder where you would like to build and compile lepton injector, and run 

`git clone git@github.com:IceCubeOpenSource/LeptonInjector.git`

to download the source code. Then, `mv` the folder that is created to rename it to `source`. We will be trying to keep our source, build, and install directories separate. Now

`mkdir build` and `mkdir install`

so we have target locations. `cd ./build` to move into the build directory. Now, we call cmake

`cmake -DCMAKE_INSTALL_PREFIX=../install ../source`

This tells cmake to install the shared objects in the `install` directory we just made. Cmake prepares a `Makefile` which calls the `g++` compiler with the necessary instructions to... compile. So now you'll call

`make -j4 && make install`

to build the project and install the project. Now you need to set all the environmental variables so this actually works. You will be adding this to your .bashrc or .bash_profile. 

To allow python to find your install directory: 
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/your/install/path/`

To allow the EarthModel to find details about the Earth:
`export EARTH_PARAMS=/your/source/path/resources/earthparams/`

# Structure
The code base is divided into several files. 
* Constants: a header defining various constants. 
* Controller: implements a class for managing the simulation
* DataWriter: writes event properties and MCTrees to an HDF5 file
* EventProps: implements a few structures used to write events in the hdf5 file. 
* h5write: may be renamed soon. This will be used to write the configurations onto a file
* LeptonInjector (the file): defines the Injector objects described above in addition to several configuration objects and event generators 
* Particle: simple implementation of particles. Includes a big enum. 
* Random: object for random number sampling.

# Cross Sections
For generating events you will need fits files of splines specifying the cross sections (total and differential cross sections). These should be made with photospline. 

# Making Contributions
If you would like to make contributions to this project, please create a branch off of the `master` branch and name it something following the template: `$YourLastName/$YourSubProject`. 
Work on this branch until you have made the changes you wished to see and your branch is _stable._ 
Then, pull from master, and create a pull request to merge your branch back into master. 
