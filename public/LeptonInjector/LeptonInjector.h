#ifndef LEPTONINJECTOR_H_INCLUDED
#define LEPTONINJECTOR_H_INCLUDED

#include <queue>

#include "EarthModelService.h"
#include "LICrossSection.h"

#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>
#include <photospline/bspline.h>

#include <iostream>

#include <Coordinates.h>
#include <Constants.h>
#include <Particle.h>
#include <Random.h>
#include <EventProps.h>
#include <DataWriter.h>

namespace LeptonInjector{
	
	// Generator configuration structures
	
	///Configuration parameters needed for all injection modes
	struct BasicInjectionConfiguration{
		BasicInjectionConfiguration();
		
		///Number of events the generator should/did generate
		uint32_t events;
		///Minimum total event energy to inject
		double energyMinimum;
		///Maximum total event energy to inject
		double energyMaximum;
		///Powerlaw index of the energy spectrum to inject
		double powerlawIndex;
		///Minimum azimuth angle at which to inject events
		double azimuthMinimum;
		///Maximum azimuth angle at which to inject events
		double azimuthMaximum;
		///Minimum zenith angle at which to inject events
		double zenithMinimum;
		///Maximum zenith angle at which to inject events
		double zenithMaximum;
		///Type of first particle to be injected in the final state
		ParticleType finalType1;
		///Type of second particle to be injected in the final state
		ParticleType finalType2;
		
		std::vector<char> crossSectionBlob;
		std::vector<char> totalCrossSectionBlob;

		void setCrossSection(const photospline::splinetable<>& crossSection, const photospline::splinetable<>& totalCrossSection);
	};
	
	///Configuration parameters for ranged injections mode
	struct RangedInjectionConfiguration : BasicInjectionConfiguration{
		RangedInjectionConfiguration();
		
		///Radius around the origin within which to target events
		double injectionRadius;
		///Length of the fixed endcaps add to the distance along which to sample interactions
		double endcapLength;
		
	};
	
	///Configuration parameters for volume injection mode
	struct VolumeInjectionConfiguration : BasicInjectionConfiguration{
		VolumeInjectionConfiguration();
		
		///Radius of the origin-centered vertical cylinder within which to inject events
		double cylinderRadius;
		///Height of the origin-centered vertical cylinder within which to inject events
		double cylinderHeight;
		
	};
	
	///Parameters for injectors placed within a MultiLeptonInjector
	struct MinimalInjectionConfiguration{
		// The below  puts in a constructor for the MinimalInjectionConfiguration. 
		// you just pass the things in order and BAM
		MinimalInjectionConfiguration(unsigned int events,
		  ParticleType finalType1, ParticleType finalType2,
		  const std::string& crossSectionPath, const std::string& totalCrossSectionPath,
		  bool ranged):
		events(events),finalType1(finalType1),finalType2(finalType2),
		crossSectionPath(crossSectionPath),totalCrossSectionPath(totalCrossSectionPath),
		ranged(ranged){}
		
		
		///Number of events the generator should/did generate
		unsigned int events;
		///Type of first particle to be injected in the final state
		ParticleType finalType1;
		///Type of second particle to be injected in the final state
		ParticleType finalType2;
		///
		std::string crossSectionPath;
		///
		std::string totalCrossSectionPath;
		///
		bool ranged;
	};
	bool operator == (const MinimalInjectionConfiguration& one , const MinimalInjectionConfiguration& two);
	//----
	
	
	
	//----
	
	class LeptonInjectorBase {
	public:
		LeptonInjectorBase();
		LeptonInjectorBase(BasicInjectionConfiguration& config);
		//No implementation of DAQ; this base class should be pure virtual
		bool Generate(){ return(false); }
        void Finish();
		//Whether this module has generated as many events already as it was configured to
		bool DoneGenerating() const{ return(eventsGenerated>=config.events); }
		std::string Name(){return("BasicInjector");}
		bool isRanged(){ return(false);}

		void Configure(const MinimalInjectionConfiguration basic, std::shared_ptr<LI_random> pass);

		std::shared_ptr<DataWriter> writer_link;

	protected:
		///Add common I3Module parameters
		//void AddBaseParameters();
		
		///Get common I3Module parameter values
		
		///Sample a random position on a disk with a given size and orientation.
		///The disk is always centered on the origin of the coordinate system.
		///\param radius the radius of the disk
		///\param zenith the zenith angle of the disk normal
		///\param azimuth the azimuth angle of the disk normal
		LI_Position SampleFromDisk(double radius, double zenith=0., double azimuth=0.);
		
		///Sample one energy value from the energy spectrum
		double SampleEnergy();
		
		///Determine the angles of the final state particles with respect to the
		///initial state neutrino direction, in the lab frame
		///\param E_total the energy of the initial neutrino
		///\param x Bjorken x for the interaction
		///\param y Bjorken y for the interaction
		///\return the relative zenith angles for the first and second final state particles
		std::pair<double,double> computeFinalStateAngles(double E_total, double x, double y);
		
		///\brief Construct an I3MCTree representing an interaction
		///
		///Samples a suitable final state and computes all resulting directions
		///and energies.
		///\param vertex the point at which the interaction occurs
		///\param dir the direction of the interacting neutrino
		///\param energy the energy of the interacting neutrino
		///\param properties the associated structure where the event properties should be recorded
		void FillTree(LI_Position vertex, LI_Direction dir, double energy, BasicEventProperties& properties, std::array<h5Particle,3>& particle_tree);
		
		///Random number source
		std::shared_ptr<LI_random> random = nullptr;
		///Configuration structure in which to store parameters
		BasicInjectionConfiguration config;
		///Number of events produced so far
		unsigned int eventsGenerated;
		///Whether an S frame has been written
		bool wroteConfigFrame;
		///Whether to suspend the tray after all events have been generated
		bool suspendOnCompletion;
		///The type of interacting neutrino this instance will produce.
		///Note that in the presence of oscillations this may not be the type of
		///the neutrino which arrived at the surface of the Earth.
		ParticleType initialType;
		
		const photospline::splinetable<>& getCrossSection() const{ return(crossSection.getCrossSection()); }
		const photospline::splinetable<>& getTotalCrossSection() const{ return(crossSection.getTotalCrossSection()); }
	private:
		
		I3CrossSection crossSection;
		
	};
	
	class RangedLeptonInjector : public LeptonInjectorBase{
	public:
		RangedLeptonInjector();
		RangedLeptonInjector(RangedInjectionConfiguration config, std::shared_ptr<earthmodel::EarthModelService> earth);
		bool Generate();
		std::string Name(){return("RangedInjector");}
		bool isRanged(){return(true);}

		// the earthmodel will just be a null poitner at instantiation
		std::shared_ptr<earthmodel::EarthModelService> earthModel = nullptr;

	private:
		RangedInjectionConfiguration config;
		///Model to use for calculating lepton range due to matter		
	};
	
	class VolumeLeptonInjector : public LeptonInjectorBase{
	public:
		VolumeLeptonInjector();
		VolumeLeptonInjector(VolumeInjectionConfiguration config);
		bool Generate();
		std::string Name(){return("VolumeInjector");}
		bool isRanged(){return(false);}
	private:
		VolumeInjectionConfiguration config;
		
	};
	
	//----
	
	
	///Construct a new direction with the given relative angles with respect to
	///an existing direction.
	///\param base the existing base direction
	///\param zenith the angle of the new direction with respect to the base
	///\param azimuth the rotation of the new direction about the base
	std::pair<double,double> rotateRelative(std::pair<double,double> base, double zenith, double azimuth);
	
	
		
	
} //namespace LeptonInjector

#endif
