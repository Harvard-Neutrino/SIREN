# LeptonInjector

LeptonInjector is a group of modules used to create events in IceCube. This code represents a standalone version of the original LeptonInjector which has been trimmed of the proprietary Icetray dependencies. It is currently a work in progress, and is *not* in a functional state! 

When finished, this code is intented to be used as a series of modules put together in a python script. 

To use it, you will

    1. Prepare a Injector object (or list of Injector objects). These are represented as `MinimalInjectionConfiguration` objects in the c++ code. 

    2. Use these, along with several generation parameters, to create a Controller object. These were called MultiLeptonInjector in the original code. 

    3. Add/remove any injectors to the controller object. Verify that the controller is properly configured.

    4. Call Controller.Execute(). This will run the simulation. 

# Dependencies

This code requires the `HDF5` C libraries. Read more about it here: https://portal.hdfgroup.org/display/support. 

It more than likely require `BOOST` once complete as well. 

What's up with splines? Right? https://github.com/IceCubeOpenSource/photospline

# Compilation and Installation

Will use a configure script, like what's in SQuIDS. So you'll do something like

`./Configure --help`

first, just to know what to do. Then 

`./Configure $YOUR_CHOSEN_OPTIONS`

to configure everything (include paths, etc.). And finally,

`make -j4 && make install`.

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
