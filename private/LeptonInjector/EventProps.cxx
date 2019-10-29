#include <EventProps.h>
#include <Particle.h>

namespace LeptonInjector {
	//Event propery objects
	

    // default constructors leave these empty 
    BasicEventProperties::BasicEventProperties(){
    }
    RangedEventProperties::RangedEventProperties(){
    }
    VolumeEventProperties::VolumeEventProperties(){
    }

    // non-default constructors to add everything
    void BasicEventProperties::fill_BasicEventProperties(double totalEnergy, double zenith, double azimuth, double finalStateX, double finalStateY,  int32_t finalType1, int32_t finalType2, int32_t initialType ){
        totalEnergy = totalEnergy;
    }




	BasicEventProperties::~BasicEventProperties(){}
	
	RangedEventProperties::~RangedEventProperties(){}
	
	VolumeEventProperties::~VolumeEventProperties(){}



}// end namespace LeptonInjector