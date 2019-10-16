#ifndef LEPTONINJECTOR_H_INCLUDED
#define LEPTONINJECTOR_H_INCLUDED

#include <queue>

#include <icetray/I3ConditionalModule.h>
#include <dataclasses/physics/I3Particle.h>
#include <phys-services/I3RandomService.h>
#include <earthmodel-service/EarthModelService.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTree.h>
#include <phys-services/I3CrossSection.h>

namespace LeptonInjector{
	
	// Generator configuration structures
	
	///Configuration parameters needed for all injection modes
	struct BasicInjectionConfiguration : public I3FrameObject{
		BasicInjectionConfiguration();
		~BasicInjectionConfiguration();
		
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
		I3Particle::ParticleType finalType1;
		///Type of second particle to be injected in the final state
		I3Particle::ParticleType finalType2;
		
		std::vector<char> crossSectionBlob;
		std::vector<char> totalCrossSectionBlob;
		
		friend class icecube::serialization::access;
		template <typename Archive>
		void serialize(Archive&, unsigned);
		
		void setCrossSection(const splinetable& crossSection, const splinetable& totalCrossSection);
	};
	
	///Configuration parameters for ranged injections mode
	struct RangedInjectionConfiguration : public BasicInjectionConfiguration{
		RangedInjectionConfiguration();
		~RangedInjectionConfiguration();
		
		///Radius around the origin within which to target events
		double injectionRadius;
		///Length of the fixed endcaps add to the distance along which to sample interactions
		double endcapLength;
		
		friend class icecube::serialization::access;
		template <typename Archive>
		void serialize(Archive&, unsigned);
	};
	
	///Configuration parameters for volume injection mode
	struct VolumeInjectionConfiguration : public BasicInjectionConfiguration{
		VolumeInjectionConfiguration();
		~VolumeInjectionConfiguration();
		
		///Radius of the origin-centered vertical cylinder within which to inject events
		double cylinderRadius;
		///Height of the origin-centered vertical cylinder within which to inject events
		double cylinderHeight;
		
		friend class icecube::serialization::access;
		template <typename Archive>
		void serialize(Archive&, unsigned);
	};
	
	///Parameters for injectors placed within a MultiLeptonInjector
	struct MinimalInjectionConfiguration{
		MinimalInjectionConfiguration(unsigned int events,
		  I3Particle::ParticleType finalType1, I3Particle::ParticleType finalType2,
		  const std::string& crossSectionPath, const std::string& totalCrossSectionPath,
		  bool ranged):
		events(events),finalType1(finalType1),finalType2(finalType2),
		crossSectionPath(crossSectionPath),totalCrossSectionPath(totalCrossSectionPath),
		ranged(ranged){}
		
		///Number of events the generator should/did generate
		unsigned int events;
		///Type of first particle to be injected in the final state
		I3Particle::ParticleType finalType1;
		///Type of second particle to be injected in the final state
		I3Particle::ParticleType finalType2;
		///
		std::string crossSectionPath;
		///
		std::string totalCrossSectionPath;
		///
		bool ranged;
	};
	bool operator==(const MinimalInjectionConfiguration&, const MinimalInjectionConfiguration&);
	
	//----
	
	// Event property structures
	
	///Parameters common to events injected in all modes
	struct BasicEventProperties : public I3FrameObject{
		BasicEventProperties(){}
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
		I3Particle::ParticleType finalType1;
		///Type of second particle which was injected in the final state
		I3Particle::ParticleType finalType2;
		///Type of the neutrino which interacted to produce this event
		I3Particle::ParticleType initialType;
		
		friend class icecube::serialization::access;
		template <typename Archive>
		void serialize(Archive&, unsigned);
	};
	
	///Parameters for events produced in ranged injection mode
	struct RangedEventProperties : public BasicEventProperties{
		RangedEventProperties(){}
		~RangedEventProperties();
		
		///Sampled distance of the closest approach of the particle path to the
		///origin of the coordinate system
		double impactParameter;
		///The total column depth along the particle path within which the
		///interaction is sampled
		double totalColumnDepth;
		
		friend class icecube::serialization::access;
		template <typename Archive>
		void serialize(Archive&, unsigned);
	};
	
	///Parameters for events produced in volume injection mode
	struct VolumeEventProperties : public BasicEventProperties{
		VolumeEventProperties(){}
		~VolumeEventProperties();
		
		///Sampled radial cylindrical coordinate of the interaction point
		double radius;
		///Sampled vertical cylindrical coordinate of the interaction point
		double z;
		
		friend class icecube::serialization::access;
		template <typename Archive>
		void serialize(Archive&, unsigned);
	};
	
	//----
	
	class LeptonInjectorBase : public I3ConditionalModule{
	public:
		LeptonInjectorBase(const I3Context& context, BasicInjectionConfiguration& config);
		virtual ~LeptonInjectorBase();
		//No implementation of DAQ; this base class should be pure virtual
		virtual void DAQ(boost::shared_ptr<I3Frame>)=0;
		void Finish();
		//Whether this module has generated as many events already as it was configured to
		bool DoneGenerating() const{ return(eventsGenerated>=config.events); }
	protected:
		///Add common I3Module parameters
		void AddBaseParameters();
		
		///Get common I3Module parameter values
		void BaseConfigure();
		
		///Sample a random position on a disk with a given size and orientation.
		///The disk is always centered on the origin of the coordinate system.
		///\param radius the radius of the disk
		///\param zenith the zenith angle of the disk normal
		///\param azimuth the azimuth angle of the disk normal
		I3Position SampleFromDisk(double radius, double zenith=0., double azimuth=0.);
		
		///Sample one energy value from the energy spectrum
		double SampleEnergy();
		
		///Sample either baseType or its antiparticle depending on config.toggleAntiparticles
		I3Particle::ParticleType SampleParticleType(I3Particle::ParticleType baseType);
		
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
		boost::shared_ptr<I3MCTree> FillTree(I3Position vertex, I3Direction dir, double energy, BasicEventProperties& properties);
		
		///Random number source
		boost::shared_ptr<I3RandomService> random;
		///Configuration structure in which to store parameters
		BasicInjectionConfiguration& config;
		///Number of events produced so far
		unsigned int eventsGenerated;
		///Whether an S frame has been written
		bool wroteConfigFrame;
		///Whether to suspend the tray after all events have been generated
		bool suspendOnCompletion;
		///The type of interacting neutrino this instance will produce.
		///Note that in the presence of oscillations this may not be the type of
		///the neutrino which arrived at the surface of the Earth.
		I3Particle::ParticleType initialType;
		
		const splinetable& getCrossSection() const{ return(crossSection.getCrossSection()); }
		const splinetable& getTotalCrossSection() const{ return(crossSection.getTotalCrossSection()); }
	private:
		
		I3CrossSection crossSection;
		
		SET_LOGGER("LeptonInjectorBase");
	};
	
	class RangedLeptonInjector : public LeptonInjectorBase{
	public:
		RangedLeptonInjector(const I3Context& context);
		RangedLeptonInjector(const I3Context& context, RangedInjectionConfiguration config);
		void Configure();
		void DAQ(boost::shared_ptr<I3Frame> frame);
	private:
		void init();
		RangedInjectionConfiguration config;
		///Model to use for calculating lepton range due to matter
		boost::shared_ptr<earthmodel::EarthModelService> earthModel;
		
		SET_LOGGER("RangedLeptonInjector");
	};
	
	class VolumeLeptonInjector : public LeptonInjectorBase{
	public:
		VolumeLeptonInjector(const I3Context& context);
		VolumeLeptonInjector(const I3Context& context, VolumeInjectionConfiguration config);
		void Configure();
		void DAQ(boost::shared_ptr<I3Frame> frame);
	private:
		void init();
		VolumeInjectionConfiguration config;
		
		SET_LOGGER("VolumeLeptonInjector");
	};
	
	//----
	
	// Utility functions
	
	///Determine whether a given particle type is a lepton.
	///\note Only knows about the types used by LeptonInjector.
	bool isLepton(I3Particle::ParticleType p);
	
	///Determine whether a given particle type is charged.
	///\note Only knows about the types used by LeptonInjector.
	bool isCharged(I3Particle::ParticleType p);
	
	///Extract the name for a given particle type.
	std::string particleName(I3Particle::ParticleType p);
	
	///Compute the portion of a particle's energy which is kinetic
	///\note Treats particles of unknown mass as massless, making their kinetic
	///      energy equal to their total energy
	double kineticEnergy(I3Particle::ParticleType type, double totalEnergy);
	
	///Figure out a particle's speed given its kinetic energy.
	///\note Treats particles of unknown mass as massless, assigning them speed c.
	double particleSpeed(I3Particle::ParticleType type, double kineticEnergy);
	
	///Guess the shape which should be associated with a given particle type.
	///Particles with long ranges (at least multiple meters) in ice are
	///considered tracks, while those with shorter ranges are labeled cascades.
	///\note Only knows about the types used by LeptonInjector.
	I3Particle::ParticleShape decideShape(I3Particle::ParticleType t);
	
	///Determine the type of neutrino which would produce the specified final state.
	///Verify that the user has provided a valid pair of particle types for
	///the final state, call log_fatal() if the settings are not valid.
	///\returns The interacting neutrino type required by the specified final state.
	I3Particle::ParticleType deduceInitialType(I3Particle::ParticleType pType1, I3Particle::ParticleType pType2);
	
	///Construct a new direction with the given relative angles with respect to
	///an existing direction.
	///\param base the existing base direction
	///\param zenith the angle of the new direction with respect to the base
	///\param azimuth the rotation of the new direction about the base
	I3Direction rotateRelative(I3Direction base, double zenith, double azimuth);
	
	///A normal I3Module can only send its output frames to the inbox associated
	///with another I3Module. This Module provides such an inbox, but instead of
	///sending output to another module's inbox it stores it in a queue of its
	///own, which it makes externally accessible.
	class OutputCollector : public I3Module{
	public:
		OutputCollector(const I3Context& ctx):I3Module(ctx){}
		void Process(){
			while(PeekFrame()){
				boost::shared_ptr<I3Frame> frame=PopFrame();
				if(!frame)
					return;
				output.push(frame);
			}
		}
		
		void DiscardOutput(){
			while(!output.empty())
				output.pop();
		}
		
		std::queue<boost::shared_ptr<I3Frame> > output;
	};
	
	void ProcessFrame(I3Module& mod, boost::shared_ptr<I3Frame> frame);
	
	class MultiLeptonInjector : public I3ConditionalModule{
	public:
		MultiLeptonInjector(const I3Context& ctx);
		///For properties appearing in both config objects, the values in rconfig will take precedence
		MultiLeptonInjector(const I3Context& ctx, RangedInjectionConfiguration rconfig, VolumeInjectionConfiguration vconfig);
		void Configure();
		void DAQ(boost::shared_ptr<I3Frame>);
	private:
		void AddParameters();
		
		I3Context innerContext;
		boost::shared_ptr<OutputCollector> collector;
		std::queue<boost::shared_ptr<I3Frame> >& results;
		std::vector<MinimalInjectionConfiguration> generatorSettings;
		std::deque<LeptonInjectorBase*> generators;
		RangedInjectionConfiguration rangedConfig;
		VolumeInjectionConfiguration volumeConfig;
	};
	
} //namespace LeptonInjector

#endif
