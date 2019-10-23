#ifndef LEPTONINJECTOR_H_INCLUDED
#define LEPTONINJECTOR_H_INCLUDED

#include <queue>

#include <earthmodel-service/EarthModelService.h>
#include <phys-services/LICrossSection.h>

#include <LeptonInjector/Particle.h>
#include <LeptonInjector/Random.h>

namespace LeptonInjector{
	
	// Generator configuration structures
	
	///Configuration parameters needed for all injection modes
	struct BasicInjectionConfiguration{
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
		ParticleType finalType1;
		///Type of second particle to be injected in the final state
		ParticleType finalType2;
		
		std::vector<char> crossSectionBlob;
		std::vector<char> totalCrossSectionBlob;
			
		void setCrossSection(const splinetable& crossSection, const splinetable& totalCrossSection);
	};
	
	///Configuration parameters for ranged injections mode
	struct RangedInjectionConfiguration : BasicInjectionConfiguration{
		RangedInjectionConfiguration();
		~RangedInjectionConfiguration();
		
		///Radius around the origin within which to target events
		double injectionRadius;
		///Length of the fixed endcaps add to the distance along which to sample interactions
		double endcapLength;
		
	};
	
	///Configuration parameters for volume injection mode
	struct VolumeInjectionConfiguration : BasicInjectionConfiguration{
		VolumeInjectionConfiguration();
		~VolumeInjectionConfiguration();
		
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
	
	
	///Construct a new direction with the given relative angles with respect to
	///an existing direction.
	///\param base the existing base direction
	///\param zenith the angle of the new direction with respect to the base
	///\param azimuth the rotation of the new direction about the base
	std::pair<double,double> rotateRelative(std::pair<double,double> base, double zenith, double azimuth);
	
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
		double seed;

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
