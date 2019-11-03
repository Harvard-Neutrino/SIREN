#ifndef LI_CONTROLLER
#define LI_CONTROLLER

#include <LeptonInjector.h>
#include <Constants.h>
#include <Coordinates.h>
#include <DataWriter.h>

namespace LeptonInjector{

class Controller{
    private:
        std::deque<LeptonInjectorBase*> generators;
        std::vector<MinimalInjectionConfiguration> configs; 

        uint seed;

        // overall generation parameters
        double minimumEnergy, maximumEnergy, powerlawIndex,
                minimumAzimuth, maximumAzimuth, minimumZenith, maximumZenith;
        
        // ones used by ranged and/or volume mode
        double injectionRadius, endcapLength;
        double cylinderRadius, cylinderHeight;

        earthmodel::EarthModelService earth();
        LI_random random_obj();

        const std::shared_ptr<earthmodel::EarthModelService> earthModel = std::make_shared<earthmodel::EarthModelService>(earth);
        const std::shared_ptr<LI_random> random = std::make_shared<LI_random>(random_obj);
        std::string earthmodelname;

        std::string out_file;

        RangedInjectionConfiguration rangedConfig;
		VolumeInjectionConfiguration volumeConfig;

        std::shared_ptr<DataWriter> datawriter = std::make_shared<DataWriter>( DataWriter() );

    public:
        // default constructor will just use some default minimal injection setup 
        Controller();
        ~Controller();
        // sending one will make a single little list... 
        Controller(MinimalInjectionConfiguration configs_received );
        // multilepton injector equivalent 
        Controller(	std::vector<MinimalInjectionConfiguration> configs_received );

        // The BEST constructor
        Controller(std::vector<MinimalInjectionConfiguration> configs_received, double minimumEnergy,
            double maximumEnergy, double powerlawIndex, double minimumAzimuth, 
            double maximumAzimuth, double minimumZenith, double maximumZenith,
            double injectionRadius=1200*Constants::m, double endcapLength=1200*Constants::m, 
            double cylinderRadius=1200*Constants::m, double cylinderHeight= 1200*Constants::m);
        
        
        void SetEarthModel( std::string new_name );
        void AddInjector(MinimalInjectionConfiguration configs_received);
        void NameOutfile( std::string out_file );
        void Execute(); 



};

} // end namespace LeptonInjector

#endif
