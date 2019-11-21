#include "Controller.h"

// This is the "controller" of the whole thing
// It is a reimplementation of the multileptoninjector 

namespace LeptonInjector {

    Controller::Controller(){
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = 1200*Constants::m;
        this->endcapLength = 1200*Constants::m;
        this->cylinderRadius = 1200*Constants::m;
        this->cylinderHeight = 1200*Constants::m;

    }

    // constructor if the user only provides one injector and no parameters
    Controller::Controller(MinimalInjectionConfiguration configs_received){
        std::cout << "Creating Controller" << std::endl;
        this->configs.push_back( configs_received );
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = 1200*Constants::m;
        this->endcapLength = 1200*Constants::m;
        this->cylinderRadius = 1200*Constants::m;
        this->cylinderHeight = 1200*Constants::m;
    }

   

    Controller::Controller(MinimalInjectionConfiguration configs_received, double minimumEnergy,
            double maximumEnergy, double powerlawIndex, double minimumAzimuth, 
            double maximumAzimuth, double minimumZenith, double maximumZenith,
            double injectionRadius, double endcapLength,
            double cylinderRadius, double cylinderHeight){
        std::cout << "Creating Controller" << std::endl;
        this->configs.push_back( configs_received );
        this->minimumEnergy   = minimumEnergy;
        this->maximumEnergy   = maximumEnergy;
        this->powerlawIndex   = powerlawIndex;
        this->minimumAzimuth  = minimumAzimuth;
        this->maximumAzimuth  = maximumAzimuth;
        this->minimumZenith   = minimumZenith;
        this->maximumZenith   = maximumZenith;

        this->injectionRadius = injectionRadius;
        this->endcapLength    = endcapLength;
        this->cylinderRadius  = cylinderRadius;
        this->cylinderHeight  = cylinderHeight;

    }


    void Controller::AddInjector( MinimalInjectionConfiguration configs_received ){
        this->configs.push_back( configs_received );
    }

    void Controller::SetEarthModel( std::string new_name ){
        this->earthmodelname = new_name;
    }

    void Controller::NameOutfile( std::string out_file_){
        out_file = out_file_;
    }

    void Controller::Execute(){
        // setup the injectors! 

        bool hasRanged=false, hasVolume=false;
		for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=this->configs.begin(), end=this->configs.end(); genSet!=end; genSet++){
			hasRanged |= genSet->ranged;
			hasVolume |= !genSet->ranged;
		}

        (*this->random).set_seed(seed);

        // sanity check! 
        if (this->minimumEnergy <= 0 ){ std::cout<< "minimum energy must be positive" << std::endl; throw; }
        if (this->maximumEnergy <= 0 ){ std::cout<<  "maximum energy must be positive"<< std::endl; throw; }
        if (this->minimumEnergy > this->maximumEnergy ){ std::cout<<  "Max energy must be greater or equal to minimum energy"<< std::endl; throw; }
        if (this->minimumAzimuth < 0 ){ std::cout<<  "minimum azimuth must be positive"<< std::endl; throw; }
        if (this->maximumAzimuth > 2*Constants::pi ){ std::cout<<  "maximum azimuth must be less than 2pi"<< std::endl; throw; }
        if (this->minimumAzimuth > this->maximumAzimuth ){ std::cout<<  "Max azimuth must be greater or equal to min."<< std::endl; throw; }
        if (this->minimumZenith < 0.0 ){ std::cout<<  "minimum zenith must be positive"<< std::endl; throw;  }
        if (this->minimumZenith > Constants::pi ){ std::cout<<  "maximum zenith must be less than or equal to pi"<< std::endl; throw;  }
        if (this->minimumZenith >this->maximumZenith){std::cout<<  "Max zenith must be greater or equal to min."<< std::endl; throw; }

        // first, construct the template injector configuration objects
        // with only those criteria shared between Configurations 
        std::cout << "min e "<<this->minimumEnergy << std::endl;
        this->rangedConfig.energyMinimum = this->minimumEnergy; 
        this->rangedConfig.energyMaximum = this->maximumEnergy; 
        this->rangedConfig.powerlawIndex = this->powerlawIndex; 
        this->rangedConfig.azimuthMinimum = this->minimumAzimuth; 
        this->rangedConfig.azimuthMaximum = this->maximumAzimuth; 
        this->rangedConfig.zenithMinimum = this->minimumZenith; 
        this->rangedConfig.zenithMaximum = this->maximumZenith; 

        this->volumeConfig.energyMinimum = this->minimumEnergy; 
        this->volumeConfig.energyMaximum = this->maximumEnergy; 
        this->volumeConfig.powerlawIndex = this->powerlawIndex; 
        this->volumeConfig.azimuthMinimum = this->minimumAzimuth; 
        this->volumeConfig.azimuthMaximum = this->maximumAzimuth; 
        this->volumeConfig.zenithMinimum = this->minimumZenith; 
        this->volumeConfig.zenithMaximum = this->maximumZenith; 

        //  SETUP EARTHMODEL

        if(hasRanged){
            this->rangedConfig.injectionRadius = this->injectionRadius; 
            this->rangedConfig.endcapLength = this->endcapLength; 
            // set pointer to earthmodel -- GetParameter("EarthModel",earthModelName);
            
            if(this->rangedConfig.injectionRadius<0){throw": InjectionRadius must be non-negative"; }
            if(this->rangedConfig.endcapLength<0){ throw ": EndcapLength must be non-negative"; }
            //context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthmodelname);
            //earthmodel::EarthModelService actual_model();
            
        }
        
        //get the properties for volume injectors
        if(hasVolume){
            this->volumeConfig.cylinderRadius = this->cylinderRadius; 
            this->volumeConfig.cylinderHeight = this->cylinderHeight;
            
            if(this->volumeConfig.cylinderRadius<0){throw": InjectionRadius must be non-negative"; }
            if(this->volumeConfig.cylinderHeight<0){ throw ": EndcapLength must be non-negative"; }

        }

        //construct all generators
        unsigned int i=0;
        for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=this->configs.begin(), end=this->configs.end(); genSet!=end; genSet++){
            std::cout << std::endl;
            //log_debug_stream("Configuring injector " << i << ":");
            //LeptonInjectorBase* generator;

            if(genSet->ranged){
                //log_debug_stream(" this is a ranged injector");
                RangedLeptonInjector* generator  = new RangedLeptonInjector(this->rangedConfig, this->earthModel, this->random);
                generator->earthModel = this->earthModel;
                generator->Configure( *genSet );//, this->random );
                generators.push_back(generator);
            }
            else{ //volume
                //log_debug_stream(" this is a volume injector");
                VolumeLeptonInjector* generator= new VolumeLeptonInjector(volumeConfig, random);
                generator->Configure( *genSet );//, this->random );
                generators.push_back(generator);
            }
                            
            
            
        } // end for loop constructing generators 
        

        // open the hdf5 file

        this->datawriter->OpenFile(this->out_file);

        uint8_t n_gen = 0;
        while(true){

            // grab the first genereator, get ready to generate! 
            //LeptonInjectorBase* active = generators.front();

            generators.back()->writer_link = this->datawriter;
            this->datawriter->AddInjector( generators.back()->Name(), generators.back()->isRanged() );
            generators.back()->Print_Configuration();

            /*
            active->writer_link = this->datawriter;
            this->datawriter->AddInjector(active->Name(), active->isRanged() );
            active->Print_Configuration();
            
            */
            // enters a generating loop. Keep calling generate until it returns FALSE 
            bool generating = true;
            // std::cout << "running generator for "<< active->config.events << " events" <<std::endl;
            while( generating ){
                generating = generators.back()->Generate();
            }
            
            // pop the generator, it's done! 
            //active = nullptr; // clean that pointer 
            generators.pop_back();
            if (generators.empty()){ break; } // check if there's another generator, if there isn't, give up
        }
        
        // deconstruct the hdf5 writer. It closes all the hdf5 stuff. 
        this->datawriter->~DataWriter();

    } // end execute



} // end namespace LeptonInjector
