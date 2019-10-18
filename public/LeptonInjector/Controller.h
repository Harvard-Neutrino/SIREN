#ifndef LI_CONTROLLER
#define LI_CONTROLLER

#include <LeptonInjector.h>
#include <Constants.h>

namespace LeptonInjector{

class Controller{
    private:
        void Generate();
        LeptonInjectorBase& ActiveGenerator;
        std::deque<LeptonInjectorBase*> generators;
        std::vector<MinimalInjectionConfiguration> configs; 

        // overall generation parameters
        double minimumEnergy, maximumEnergy, powerlawIndex,
                minimumAzimuth, maximumAzimuth, minimumZenith, maximumZenith;
        
        RangedInjectionConfiguration rangedConfig;
		VolumeInjectionConfiguration volumeConfig;

    public:
        // default constructor will just use some default minimal injection setup 
        Controller();
        ~Controller();
        // sending one will make a single little list... 
        Controller(MinimalInjectionConfiguration configs_received );
        // multilepton injector equivalent 
        Controller(	std::vector<MinimalInjectionConfiguration> configs_received );

        // The BEST constructor
        Controller(std::vector<MinimalInjectionConfiguration> configs_received, double minimumEnergy
            double maximumEnergy, double powerlawIndex, double minimumAzimuth, 
            double maximumAzimuth, double minimumZenith, double maximumZenith);
        

        void AddInjector(MinimalInjectionConfiguration configs_received);
        void Execute(); 

};

} // end namespace LeptonInjector

#endif
