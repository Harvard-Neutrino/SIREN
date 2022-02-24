#pragma once
#ifndef LI_Controller_H
#define LI_Controller_H

#include <memory>
#include <iostream>

#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Constants.h> // pi, GeV, degrees, meters, ...
#include <LeptonInjector/DataWriter.h> // adds DataWriter class

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

// This defines the 'Controller' object that manages the overall process
// I tried to mimic the UX of working with IceTrays

namespace LeptonInjector{

class Controller{
private:
    // a deque (basically an array) to hold all of the base injectors. Will be filled with volume
    //      and ranged mode injectors
    std::deque<std::shared_ptr<InjectorBase>> generators;

    // seeds the random number generator
    uint seed = 100;

    // overall generation parameters
    double minimumEnergy, maximumEnergy, powerlawIndex,
           minimumAzimuth, maximumAzimuth, minimumZenith, maximumZenith;

    // parameters used by ranged mode
    double injectionRadius, endcapLength;
    // parameters used by volume mode
    double cylinderRadius, cylinderHeight;

    // Earth's name, left here for historical reasons
    std::string earthmodelname = "Earth";

    // Shared pointer refering to the Earth model used.
    std::shared_ptr<earthmodel::EarthModel> earthModel;
    // Shared pointer refering to the random number generator.
    const std::shared_ptr<LI_random> random = std::make_shared<LI_random>();

    // data file
    std::string out_file="./outfile.h";
    // configuration file
    std::string lic_file="./config.lic";

    // default ranged mode configuration
    BasicInjectionConfiguration rangedConfig = BasicInjectionConfiguration();
    // default volume mode configuration
    BasicInjectionConfiguration volumeConfig = BasicInjectionConfiguration();

    // shared pointer to the object charged with interacting with the OS
    std::shared_ptr<DataWriter> datawriter = std::make_shared<DataWriter>();

public:
    // default constructor will just use some default minimal injection setup
    Controller();

    // sending one will make a single little list...
    Controller(std::shared_ptr<InjectorBase> injector);
    // multilepton injector equivalent

    // changes the Earth model to be used with the injectors
    void SetEarthModel(std::shared_ptr<earthmodel::EarthModel> earthModel);
    // adds a new injector to be used in the process
    void AddInjector(std::shared_ptr<InjectorBase> generator);
    // changes the name of the data file
    void NameOutfile( std::string out_file );
    // changes the name of the configuration file
    void NameLicFile( std::string lic_file );
    // changes whether or not to overwrite any preexisting LIC files or append to them
    void Overwrite( bool overwrite );
    void setSeed( uint seedno ){this->seed = seedno ;}
    // executes the injector process with the configurated parameters
    void Execute();
};

} // namespace LeptonInjector

#endif // LI_Controller_H

