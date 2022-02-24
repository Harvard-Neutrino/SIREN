#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/DataWriter.h"
#include "LeptonInjector/Controller.h"

namespace LeptonInjector {

Controller::Controller(){
}

// constructor if the user only provides one injector and no parameters
Controller::Controller(std::shared_ptr<InjectorBase> generator){
    generators.push_back(std::move(generator));
}

void Controller::AddInjector(std::shared_ptr<InjectorBase> generator){
    generators.push_back(std::move(generator));
}

void Controller::SetEarthModel(std::shared_ptr<earthmodel::EarthModel> earth_model){
    this->earthModel = earth_model;
}

void Controller::NameOutfile( std::string out_file_){
    out_file = out_file_;
}

void Controller::NameLicFile( std::string lic_file_ ){
    this->lic_file = lic_file_;
}

void Controller::Overwrite( bool overwrite_){
    this->datawriter->SetOverwrite( overwrite_ );
}

void Controller::Execute(){
    // setup the injectors!
    std::cout << "Executing Injectors" << std::endl;

    (*this->random).set_seed(seed);

    this->datawriter->OpenLICFile( this->lic_file);
    this->datawriter->OpenFile(this->out_file);

    uint8_t n_gen = 0;
    while(true) {
        std::cout << "Starting up "<< generators.back()->Name() << std::endl;
        //generators.back()->writer_link = this->datawriter;
        //this->datawriter->AddInjector( generators.back()->Name(), generators.back()->isRanged() );
        //this->datawriter->WriteConfig( generators.back()->getConfig(), generators.back()->isRanged() );

        // enters a generating loop. Keep calling generate until it returns FALSE
        while(*generators.back()) {
            generators.back()->GenerateEvent();
        }

        // pop the generator, it's done!
        generators.pop_back();
        if (generators.empty()){ break; } // check if there's another generator, if there isn't, give up
    }

    // deconstruct the hdf5 writer. It closes all the hdf5 stuff.
    // this->datawriter->~DataWriter();
}

} // namespace LeptonInjector

