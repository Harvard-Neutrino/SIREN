#include <EventProps.h>
#include <Particle.h>

namespace LeptonInjector {
	//Event propery objects
	
    h5Particle::h5Particle(){
    }
	h5Particle::h5Particle(bool initial_, int32_t ptype_, LI_Position pos_, LI_Direction dir_, double energy_){
        initial = initial_;
        ptype = ptype_;
        for (uint8_t iter=0; iter<n_dimensions;iter++){
            pos[iter] = pos_.at(iter);
        }
        dir[0] = dir_.zenith;
        dir[1] = dir_.azimuth;
        energy = energy_;
    }



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





}// end namespace LeptonInjector