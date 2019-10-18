#include <LeptonInjector/Controller.h>

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

    }

    // constructor if the user only provides one injector and no parameters
    Controller::Controller(MinimalInjectionConfiguration configs_received){
        this->configs.push_back( conconfigs_receivedfigs );
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;
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
    }

    Controller(std::vector<MinimalInjectionConfiguration> configs_received, double minimumEnergy
            double maximumEnergy, double powerlawIndex, double minimumAzimuth, 
            double maximumAzimuth, double minimumZenith, double maximumZenith){

        this->configs = configs_received;
        this->minimumEnergy   = 100*Constants::GeV ;
        this->maximumEnergy   = 10000*Constants::GeV;
        this->powerlawIndex   = 2.0;
        this->minimumAzimuth  = 0.0*Constants::degrees;
        this->maximumAzimuth  = 360*Constants::degrees;
        this->minimumZenith   = 80*Constants::degrees;
        this->maximumZenith   = 180*Constants::degrees;

    }

    Controller::~Controller(){
    }

    void Controller::AddInjector( MinimalInjectionConfiguration configs_received ){
        this->configs.push_back( configs_received );
    }

    void Controller::Execute(){
        // setup the injectors! 

        bool hasRanged=false, hasVolume=false;
		for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=configs.begin(), end=configs.end(); genSet!=end; genSet++){
			hasRanged |= genSet->ranged;
			hasVolume |= !genSet->ranged;
		}

        LI_random random(this->seed);

        // sanity check! 
        if (this->energyMinium <= 0 ){ throw "minimum energy must be positive"; }
        if (this->energyMaximum <= 0 ){ throw "maximum energy must be positive"; }
        if (this->energyMaximum < this->energyMinimum ){ throw "Max energy must be greater or equal to minimum energy"; }
        if (this->azimuthMinimum < 0 ){ throw "minimum azimuth must be positive"; }
        if (this->azimuthMaximum > 2*Constants::pi ){ throw "maximum azimuth must be less than 2pi"; }
        if (this->azimuthMinimum < this->azimuthMaximum ){ throw "Max azimuth must be greater or equal to min."; }
        if (this->zenithMinimum < 0.0 ){ throw "minimum zenith must be positive"; }
        if (this->zenithMaximum > Constants::pi ){ throw "maximum zenith must be less than or equal to pi"; }
        if (this->zenithMinimum >this->zenithMaximum){throw "Max zenith must be greater or equal to min."}


			if(rangedConfig.zenithMinimum<0.0)
				log_fatal_stream(GetName() << ": minimum zenith angle must be greater than or equal to zero");
			if(rangedConfig.zenithMaximum>constants::pi<double>())
				log_fatal_stream(GetName() << ": maximum zenith angle must be less than or equal to pi");
			if(rangedConfig.zenithMinimum>rangedConfig.zenithMaximum)
				log_fatal_stream(GetName() << ": minimum zenith angle must be less than or equal to maximum zenith angle");

        // first, construct the template injector configuration objects
        this->rangedConfig.minimumEnergy = this->minimumEnergy; 
        this->rangedConfig.maximumEnergy = this->maximumEnergy; 
        this->rangedConfig.powerlawIndex = this->powerlawIndex; 
        this->rangedConfig.minimumAzimuth = this->minimumAzimuth; 
        this->rangedConfig.maximumAzimuth = this->maximumAzimuth; 
        this->rangedConfig.minimumZenith = this->minimumZenith; 
        this->rangedConfig.maximumZenith = this->maximumZenith; 

        this->volumeConfig.minimumEnergy = this->minimumEnergy; 
        this->volumeConfig.maximumEnergy = this->maximumEnergy; 
        this->volumeConfig.powerlawIndex = this->powerlawIndex; 
        this->volumeConfig.minimumAzimuth = this->minimumAzimuth; 
        this->volumeConfig.maximumAzimuth = this->maximumAzimuth; 
        this->volumeConfig.minimumZenith = this->minimumZenith; 
        this->volumeConfig.maximumZenith = this->maximumZenith; 
    }


}
