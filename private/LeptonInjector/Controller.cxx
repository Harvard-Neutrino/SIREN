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

    //constructor if the user provides a list of injectors and no parameters 
    Controller::Controller(	std::vector<MinimalInjectionConfiguration> configs_received ){
        this->configs = configs_received;
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

    Controller::Controller(std::vector<MinimalInjectionConfiguration> configs_received, double minimumEnergy,
            double maximumEnergy, double powerlawIndex, double minimumAzimuth, 
            double maximumAzimuth, double minimumZenith, double maximumZenith,
            double injectionRadius, double endcapLength,
            double cylinderRadius, double cylinderHeight){

        this->configs = configs_received;
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

        this->injectionRadius = injectionRadius;
        this->endcapLength    = endcapLength;
        this->cylinderRadius  = cylinderRadius;
        this->cylinderHeight  = cylinderHeight;

    }

    Controller::~Controller(){
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

        LI_random random(this->seed);

        // sanity check! 
        if (this->minimumEnergy <= 0 ){ throw "minimum energy must be positive"; }
        if (this->maximumEnergy <= 0 ){ throw "maximum energy must be positive"; }
        if (this->minimumEnergy < this->maximumEnergy ){ throw "Max energy must be greater or equal to minimum energy"; }
        if (this->minimumAzimuth < 0 ){ throw "minimum azimuth must be positive"; }
        if (this->maximumAzimuth > 2*Constants::pi ){ throw "maximum azimuth must be less than 2pi"; }
        if (this->minimumAzimuth < this->maximumAzimuth ){ throw "Max azimuth must be greater or equal to min."; }
        if (this->minimumZenith < 0.0 ){ throw "minimum zenith must be positive"; }
        if (this->minimumZenith > Constants::pi ){ throw "maximum zenith must be less than or equal to pi"; }
        if (this->minimumZenith >this->maximumZenith){throw "Max zenith must be greater or equal to min.";}

        // first, construct the template injector configuration objects
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
//            log_debug_stream("Configuring injector " << i << ":");
            LeptonInjectorBase* generator=NULL;
            try{
                if(genSet->ranged){
                    //log_debug_stream(" this is a ranged injector");
                    RangedLeptonInjector* generator = new RangedLeptonInjector(this->rangedConfig, this->earthModel, this->random);
                    generator->earthModel = this->earthModel;
                }
                else{ //volume
                    //log_debug_stream(" this is a volume injector");
                    VolumeLeptonInjector* generator= new VolumeLeptonInjector(this->volumeConfig, this->random);
                }
                
                
                //set properties not shared with other injectors, or which are not part of the config object


                /*
                generator->GetConfiguration().Set("NEvents",boost::python::object(genSet->events));
                generator->GetConfiguration().Set("FinalType1",boost::python::object(genSet->finalType1));
                generator->GetConfiguration().Set("FinalType2",boost::python::object(genSet->finalType2));
                generator->GetConfiguration().Set("RandomService",boost::python::object(randomServiceName));
                generator->GetConfiguration().Set("DoublyDifferentialCrossSectionFile",boost::python::object(genSet->crossSectionPath));
                generator->GetConfiguration().Set("TotalCrossSectionFile",boost::python::object(genSet->totalCrossSectionPath));
                generator->GetConfiguration().Set("SuspendOnCompletion",boost::python::object(false));
                
                generator->SetName(GetName()+"_Generator_"+boost::lexical_cast<std::string>(i++));
                generator->Configure(); */
                
            }catch(...){
                delete generator;
                throw;
            } // end try/catch
            generators.push_back(generator);
            
        } // end for loop constructing generators 
        
        // open the hdf5 file

        this->datawriter.OpenFile(this->out_file);

        bool generating = true;
        uint8_t n_gen = 0;
        while(generating){

            // grab the first genereator, get ready to generate! 
            LeptonInjectorBase* active = generators.front();
            active->writer_link = this->datawriter;
            this->datawriter.AddInjector(active->Name(), active->isRanged() );

            // enters a generating loop. Keep calling generate until it returns FALSE 
            while( active->Generate() ); 
            
            // pop the generator, it's done! 
            generators.pop_front();
            active = nullptr; // clean that pointer 
            if (generators.empty()){ break; } // check if there's another generator, if there isn't, give up
        }

        // call the destructor on datawriter 
        delete( &this->datawriter );

    } // end execute



}
