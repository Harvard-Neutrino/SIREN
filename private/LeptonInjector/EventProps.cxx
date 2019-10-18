#include <LeptonInjector/EventProps.h>
#include <LeptonInjector/Particle.h>

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
    void BasicEventProperties::fill_BasicEventProperties(double totalEnergy, double zenith, double azimuth, double finalStateX, double finalStateY,  ParticleType finalType1, ParticleType finalType2, ParticleType initialType ){
        totalEnergy = totalEnergy;
    }




	BasicEventProperties::~BasicEventProperties(){}
	
	RangedEventProperties::~RangedEventProperties(){}
	
	VolumeEventProperties::~VolumeEventProperties(){}

}// end namespace LeptonInjector