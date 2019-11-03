#ifndef LI_EVENT
#define LI_EVENT

#include <Coordinates.h>

// TODO: add support to construct these from their related, complicated, data types! 

namespace LeptonInjector {
    // Event property structures
	
	struct h5Particle{
		h5Particle();
		h5Particle(bool initial_, int32_t ptype_, LI_Position pos_, LI_Direction dir_, double energy_);
		virtual ~h5Particle();

		bool initial;
		int32_t ptype;
		double pos[3];
		double dir[2];
		double energy;
		
	};


	// one of these 
	///Parameters common to events injected in all modes
	struct BasicEventProperties{
		BasicEventProperties();
		~BasicEventProperties();
		
		///Total energy in the final state (lab frame)
		double totalEnergy;
		///Sampled zenith angle (of final state particle 1)
		double zenith;
		///Sampled azimuth angle (of final state particle 1)
		double azimuth;
		///Bjorken x for the interaction
		double finalStateX;
		///Bjorken y for the interaction
		///p1.energy = (1-finalStateY)*totalEnergy
		///p2.energy = finalStateY*totalEnergy
		double finalStateY;
		///Type of first particle which was injected in the final state
		int32_t finalType1;
		///Type of second particle which was injected in the final state
		int32_t finalType2;
		///Type of the neutrino which interacted to produce this event
		int32_t initialType;

		// interaction type [0 - CC, 1 - NC, 2 - GR ]
		uint8_t interaction;

        void fill_BasicEventProperties(double totalEnergy, double zenith, double azimuth, double finalStateX, double finalStateY,  int32_t finalType1, int32_t finalType2, int32_t initialType );
		
	};
	
	///Parameters for events produced in ranged injection mode
	struct RangedEventProperties : public BasicEventProperties{
		RangedEventProperties();
		~RangedEventProperties();
		
		///Sampled distance of the closest approach of the particle path to the
		///origin of the coordinate system
		double impactParameter;
		///The total column depth along the particle path within which the
		///interaction is sampled
		double totalColumnDepth;
		
	};
	
	///Parameters for events produced in volume injection mode
	struct VolumeEventProperties : public BasicEventProperties{
		VolumeEventProperties();
		~VolumeEventProperties();
		
		///Sampled radial cylindrical coordinate of the interaction point
		double radius;
		///Sampled vertical cylindrical coordinate of the interaction point
		double z;
		
	};

    void fill_BasicEventProperties(double totalEnergy, double zenith, double azimuth, double finalStateX, double finalStateY,  int32_t finalType1, int32_t finalType2, int32_t initialType );


}// end namespace LeptonInjector

#endif
