#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/EventProps.h>

#include <cassert>
#include <fstream>

#include <boost/math/constants/constants.hpp>
#include <boost/make_shared.hpp>

#include <icetray/open.h>
#include <icetray/I3Units.h>
#include <dataclasses/physics/I3MCTree.h>

// namespace constants = boost::math::constants;

namespace LeptonInjector{
	
		
	//--------------
	//Config objects
	
	BasicInjectionConfiguration::BasicInjectionConfiguration():
	events(1),
	energyMinimum(10*Constants::GeV),
	energyMaximum((1e9)*Constants::GeV),
	powerlawIndex(1.0),
	azimuthMinimum(0),
	azimuthMaximum(2*constants::pi<double>()),
	zenithMinimum(0),
	zenithMaximum(constants::pi<double>()),
	finalType1(ParticleType::MuMinus),
	finalType2(ParticleType::Hadrons)
	{}
	
	BasicInjectionConfiguration::~BasicInjectionConfiguration(){}
	
	RangedInjectionConfiguration::RangedInjectionConfiguration():
	injectionRadius(1200*LeptonInjector::Constants::m),
	endcapLength(1200*LeptonInjector::Constants::m)
	{}
	
	RangedInjectionConfiguration::~RangedInjectionConfiguration(){}
	
	VolumeInjectionConfiguration::VolumeInjectionConfiguration():
	cylinderRadius(1200*LeptonInjector::Constants::m),
	cylinderHeight(1200*LeptonInjector::Constants::m)
	{}
	
	VolumeInjectionConfiguration::~VolumeInjectionConfiguration(){}
	
	// TODO: update this, find out where the splinetable object is coming from 
	void BasicInjectionConfiguration::setCrossSection(const splinetable& crossSection, const splinetable& totalCrossSection){
		splinetable_buffer buf;
		buf.size=0;
		buf.mem_alloc=&malloc;
		buf.mem_realloc=&realloc;
		int result=writesplinefitstable_mem(&buf, &crossSection);
		if(result!=0){
			free(buf.data);
			log_fatal_stream("photospline error while serializing cross section: "
							 << result);
		}
		crossSectionBlob.resize(buf.size);
		std::copy((char*)buf.data,(char*)buf.data+buf.size,&crossSectionBlob[0]);
		free(buf.data);
		
		buf.size=0;
		result=writesplinefitstable_mem(&buf, &totalCrossSection);
		if(result!=0){
			free(buf.data);
			log_fatal_stream("photospline error while serializing cross section: "
							 << result);
		}
		totalCrossSectionBlob.resize(buf.size);
		std::copy((char*)buf.data,(char*)buf.data+buf.size,&totalCrossSectionBlob[0]);
		free(buf.data);
	}
	
	//--------------------

	
	//-----------
	//Module base
	
	LeptonInjectorBase::LeptonInjectorBase(BasicInjectionConfiguration& config):
	I3ConditionalModule(context),
	config(config),
	eventsGenerated(0),
	wroteConfigFrame(false),
	suspendOnCompletion(true){
		//do NOTHING with config in this constructor, as it is not yet fully constructed
		AddOutBox("OutBox");
	}
	
	LeptonInjectorBase::~LeptonInjectorBase(){
	}
	
	void LeptonInjectorBase::AddBaseParameters(){
		AddParameter("NEvents",
					 "Number of events to generate",
					 config.events);
		AddParameter("MinimumEnergy",
					 "Minimum total event energy to inject",
					 config.energyMinimum);
		AddParameter("MaximumEnergy",
					 "Maximum total event energy to inject",
					 config.energyMaximum);
		AddParameter("PowerlawIndex",
					 "Powerlaw index of the energy spectrum to inject "
					 "(should be positive)",
					 config.powerlawIndex);
		AddParameter("MinimumAzimuth",
					 "Minimum azimuth angle for injected events",
					 config.azimuthMinimum);
		AddParameter("MaximumAzimuth",
					 "Maximum azimuth angle for injected events",
					 config.azimuthMaximum);
		AddParameter("MinimumZenith",
					 "Minimum zenith angle for injected events",
					 config.zenithMinimum);
		AddParameter("MaximumZenith",
					 "Maximum zenith angle for injected events",
					 config.zenithMaximum);
		AddParameter("FinalType1",
					 "The first particle type in the final state",
					 config.finalType1);
		AddParameter("FinalType2",
					 "The seocnd particle type in the final state",
					 config.finalType2);
		AddParameter("RandomService",
					 "Name of the random service to use",
					 "I3RandomService");
		AddParameter("DoublyDifferentialCrossSectionFile",
					 "Path to the spline FITS file representing the doubly-differential cross section",
					 "");
		AddParameter("TotalCrossSectionFile",
					 "Path to the spline FITS file representing the total cross section as a function of energy"
					 " (same as DoublyDifferentialCrossSectionFile but integrated over x and y)",
					 "");
		AddParameter("SuspendOnCompletion",
					 "Suspend the tray after all events have been generated",
					 suspendOnCompletion);
	}
	
	void LeptonInjectorBase::BaseConfigure(){
		std::string randomServiceName;
		std::string dd_crossSectionFile;
		std::string total_crossSectionFile;
		
		GetParameter("NEvents",config.events);
		GetParameter("MinimumEnergy",config.energyMinimum);
		GetParameter("MaximumEnergy",config.energyMaximum);
		GetParameter("PowerlawIndex",config.powerlawIndex);
		GetParameter("MinimumAzimuth",config.azimuthMinimum);
		GetParameter("MaximumAzimuth",config.azimuthMaximum);
		GetParameter("MinimumZenith",config.zenithMinimum);
		GetParameter("MaximumZenith",config.zenithMaximum);
		GetParameter("FinalType1",config.finalType1);
		GetParameter("FinalType2",config.finalType2);
		GetParameter("RandomService",randomServiceName);
		GetParameter("DoublyDifferentialCrossSectionFile",dd_crossSectionFile);
		GetParameter("TotalCrossSectionFile",total_crossSectionFile);
		GetParameter("SuspendOnCompletion",suspendOnCompletion);
		
		if(config.events==0)
			log_fatal_stream(GetName() << ": there's no point in running this if you don't generate at least one event");
		if(config.energyMinimum<=0)
			log_fatal_stream(GetName() << ": minimum energy must be positive");
		if(config.energyMaximum<=0)
			log_fatal_stream(GetName() << ": maximum energy must be positive");
		if(config.energyMaximum<config.energyMinimum)
			log_fatal_stream(GetName() << ": maximum energy must be greater than or equal to minimum energy");
		if(config.azimuthMinimum<0.0)
			log_fatal_stream(GetName() << ": minimum azimuth angle must be greater than or equal to zero");
		if(config.azimuthMaximum>2*constants::pi<double>())
			log_fatal_stream(GetName() << ": maximum azimuth angle must be less than or equal to 2 pi");
		if(config.azimuthMinimum>config.azimuthMaximum)
			log_fatal_stream(GetName() << ": minimum azimuth angle must be less than or equal to maximum azimuth angle");
		if(config.zenithMinimum<0.0)
			log_fatal_stream(GetName() << ": minimum zenith angle must be greater than or equal to zero");
		if(config.zenithMaximum>constants::pi<double>())
			log_fatal_stream(GetName() << ": maximum zenith angle must be less than or equal to pi");
		if(config.zenithMinimum>config.zenithMaximum)
			log_fatal_stream(GetName() << ": minimum zenith angle must be less than or equal to maximum zenith angle");
		try{
			initialType=deduceInitialType(config.finalType1,config.finalType2);
		}catch(std::runtime_error& re){
			log_error_stream("While configuring " << GetName());
			throw;
		}
		random = context_.Get<boost::shared_ptr<I3RandomService> >(randomServiceName);
		if(!random)
			log_fatal_stream(GetName() << ": A random service is required");
		if(dd_crossSectionFile.empty())
			log_fatal_stream(GetName() << ": DoublyDifferentialCrossSectionFile must be specified");
		else if(total_crossSectionFile.empty())
			log_fatal_stream(GetName() << ": TotalCrossSectionFile must be specified");
		else
			crossSection.load(dd_crossSectionFile,total_crossSectionFile);
	}
	
	void LeptonInjectorBase::Finish(){
		if(eventsGenerated!=config.events)
			log_error_stream(GetName() << ": Only " << eventsGenerated <<
							 " event have been output out of a requested total of " << config.events);
	}
	
	std::array<double,3> LeptonInjectorBase::SampleFromDisk(double radius, double zenith, double azimuth){
		//choose a random point on a disk laying in the xy plane
		double t=random->Uniform(0,2*Constants::pi);
		double u=random->Uniform()+random->Uniform();
		double r=(u>1.?2.-u:u)*radius;
		std::array<double, 3>  pos = {r*cos(t) ,r*sin(t), 0.0};
		//now rotate to make the disc perpendicular to the requested normal vector
		pos = RotateY(pos, zenith );
		pos = RotateZ(pos, azimuth);
		return(pos);
	}
	
	double LeptonInjectorBase::SampleEnergy(){
		if(config.energyMinimum==config.energyMaximum)
			return(config.energyMinimum); //return the only allowed energy
			
		if(config.powerlawIndex==1.0) //sample uniformly in log space
			return(pow(10.0,random->Uniform(log10(config.energyMinimum),log10(config.energyMaximum))));
		else{
			double u=random->Uniform();
			double energyP=(1-u)*pow(config.energyMinimum,1-config.powerlawIndex) + u*pow(config.energyMaximum,1-config.powerlawIndex);
			return(pow(energyP,1/(1-config.powerlawIndex)));
		}
	}
	
    // this function returns a pair of angles
	std::pair<double,double> LeptonInjectorBase::computeFinalStateAngles(double E_total, double x, double y){
		const double M_N = crossSection.GetTargetMass();
		double theta1=0, theta2=0;
		
		//first particle is a lepton, which covers CC, NC, and leptonic GR
		if(isLepton(config.finalType1)){
			double m1=Particle(config.finalType1).GetMass();
			double E1 = (1 - y) * E_total;
			double cos_theta1, kE1;
			
			if(!isLepton(config.finalType2)){ //CC and NC have Hadrons as particle 2
				//squared kinetic energy of final state particle 1:
				double kE1sq=E1*E1 - m1*m1;
				if(kE1sq<=0){
                    throw "Negative kinetic energy. Not good";
                }
				cos_theta1=(E1 - x*y*M_N - m1*m1/(2*E_total))/sqrt(kE1sq);
				kE1=sqrt(kE1sq);
			}
			else{ //leptonic GR
				double m_e = Constants::electronMass;
				
				if(E1<=0){ throw "Bjorken Y > 1?"; }
                
				
				cos_theta1=1 - (m_e*m_e + 2*m_e*E_total - m1*m1)/(2*E_total*E1);
				kE1=E1;
			}
			
			if(cos_theta1<-1){
                // commented out until new logger is implemented
//				log_warn_stream("cos(theta) underflow (" << cos_theta1 << "); rounding up to -1"
//					"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y << ")");
				cos_theta1=-1;
			}
			else if(cos_theta1>1){
				//tell the user if the difference was large enough to plausibly not be just round-off
                
                // need new logger 
//				if((cos_theta1-1)>1e-3)
//					log_warn_stream("cos(theta) overflow (" << cos_theta1 << "); rounding down to 1"
//						"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y <<")");
				cos_theta1=1;
			}
			
			theta1=acos(cos_theta1);
			
			//longitudinal component of final state particle 2 momentum:
			double p_long=E_total-kE1*cos_theta1;
			//transverse component of final state particle 2 momentum:
			double p_trans=kE1*sin(theta1);
			theta2=atan(p_trans/p_long);
		}
		//otherwise we have hadronic GR, so both final state masses are unknown
		//and there isn't much we can do, so leave everything colinear
		// angle of particle 1 from initial dir
		return(std::make_pair(theta1,theta2));
	}
	
	// take a direction, deflect that direction by a distance /zenith/
	//		rotate the new direction around the initial direction by /azimuth/
	// So the zenith and azimuth are only what their names would suggest in the coordinate system where 
	//		/base/ is the \hat{z} axis 
	std::pair<double,double> rotateRelative(std::pair<double,double> base, double zenith, double azimuth){
		std::pair<double, double> result;
		result.first += zenith*cos(azimuth);
		result.second+= zenith*sin(azimuth);
		return(result);
	}
	
	boost::shared_ptr<I3MCTree> LeptonInjectorBase::FillTree(I3Position vertex, I3Direction dir, double energy, BasicEventProperties& properties){
		const I3CrossSection::finalStateRecord& fs=crossSection.sampleFinalState(energy,config.finalType1,random);
		
		std::pair<double,double> relativeZeniths=computeFinalStateAngles(energy,fs.x,fs.y);
		double azimuth1=random->Uniform(0,2*Constants:pi);
		double azimuth2=azimuth1+(azimuth1<Constants::pi ? 1 : -1)*Constants::pi ;
		
		//Make the first final state particle
		I3Particle p1(decideShape(config.finalType1),config.finalType1);
		p1.SetLocationType(I3Particle::InIce);
		p1.SetPos(vertex);
		p1.SetDir(rotateRelative(dir,relativeZeniths.first,azimuth1));
		p1.SetEnergy(kineticEnergy(p1.GetType(),(1-fs.y)*energy));
		p1.SetSpeed(particleSpeed(p1.GetType(),p1.GetEnergy()));
		p1.SetTime(0.0);
		
		//Make the second final state particle
		I3Particle p2(decideShape(config.finalType2),config.finalType2);
		p2.SetLocationType(I3Particle::InIce);
		p2.SetPos(vertex);
		p2.SetDir(rotateRelative(dir,relativeZeniths.second,azimuth2));
		p2.SetEnergy(kineticEnergy(p2.GetType(),fs.y*energy));
		p2.SetSpeed(particleSpeed(p2.GetType(),p2.GetEnergy()));
		p2.SetTime(0.0);
		
		boost::shared_ptr<I3MCTree> mctree(new I3MCTree);
		//Make a dummy primary
		I3Particle primary(I3Particle::Primary,initialType);
		primary.SetEnergy(kineticEnergy(primary.GetType(),energy));
		primary.SetSpeed(particleSpeed(primary.GetType(),primary.GetEnergy()));
		primary.SetDir(dir);
		primary.SetPos(vertex);
		primary.SetTime(0.0);
		primary.SetLength(0.0);
		I3MCTree::iterator primaryIt=mctree->insert(mctree->end(), primary);
		I3MCTree::iterator p1It=mctree->append_child(primaryIt, p1);
		mctree->insert_after(p1It, p2);
		
		properties.totalEnergy=energy;
		properties.zenith=dir.GetZenith();
		properties.azimuth=dir.GetAzimuth();
		properties.finalStateX=fs.x;
		properties.finalStateY=fs.y;
		properties.finalType1=p1.GetType();
		properties.finalType2=p2.GetType();
		properties.initialType=primary.GetType();
		
		return(mctree);
	}
	
	//-----------------------
	//Ranged injection module
	
	RangedLeptonInjector::RangedLeptonInjector(const I3Context& context):
	LeptonInjectorBase(context,config){
		init();
	}
	
	RangedLeptonInjector::RangedLeptonInjector(const I3Context& context, RangedInjectionConfiguration config_):
	LeptonInjectorBase(context,config),config(config_){
		init();
	}
	
	void RangedLeptonInjector::init(){
		AddBaseParameters();
		AddParameter("InjectionRadius",
					 "Radius around the origin within which to target events",
					 config.injectionRadius);
		AddParameter("EndcapLength",
					 "Length of the fixed endcaps add to each end of the distance "
					 "along which to sample interactions",
					 config.endcapLength);
		AddParameter("EarthModel",
					 "Name of the Earth model service to use",
					 "");
		AddOutBox("OutBox");
	}
	
	void RangedLeptonInjector::Configure(){
		BaseConfigure();
		GetParameter("InjectionRadius",config.injectionRadius);
		GetParameter("EndcapLength",config.endcapLength);
		std::string earthModelName;
		GetParameter("EarthModel",earthModelName);
		
		if(config.injectionRadius<0)
			log_fatal_stream(GetName() << ": InjectionRadius must be non-negative");
		if(config.endcapLength<0)
			log_fatal_stream(GetName() << ": EndcapLength must be non-negative");
		earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
		if(!earthModel)
			log_fatal_stream(GetName() << ": an Earth model service is required");
	}
	
	void RangedLeptonInjector::DAQ(boost::shared_ptr<I3Frame> frame){
		//first, make sure configuration gets written once
		if(!wroteConfigFrame){
			boost::shared_ptr<I3Frame> sframe(new I3Frame('S'));
			boost::shared_ptr<RangedInjectionConfiguration> sconfig(new RangedInjectionConfiguration(config));
			sconfig->setCrossSection(getCrossSection(),getTotalCrossSection());
			sframe->Put("LeptonInjectorProperties",sconfig);
			PushFrame(sframe);
			wroteConfigFrame=true;
		}
		if(DoneGenerating()){
			PushFrame(frame);
			return;
		}
		
		//Choose an energy
		double energy=SampleEnergy();
		
		//Pick a direction on the sphere
		I3Direction dir(acos(random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum))),
						random->Uniform(config.azimuthMinimum,config.azimuthMaximum));
		log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
		
		//decide the point of closest approach
		I3Position pca=SampleFromDisk(config.injectionRadius,dir.GetZenith(),dir.GetAzimuth());
		log_trace_stream("pca=(" << pca.GetX() << ',' << pca.GetY() << ',' << pca.GetZ() << ')');
		
		//Figure out where we want the vertex
		//Add up the column depth for the range of a muon at this energy with the
		//column depth for the fixed endcaps to ensure that the whole detector is
		//covered
		using namespace earthmodel::EarthModelCalculator;
		double totalColumnDepth=MWEtoColumnDepthCGS(GetLeptonRange(energy))
		+earthModel->GetColumnDepthInCGS(pca-config.endcapLength*dir,pca+config.endcapLength*dir);
		//See whether that much column depth actually exists along the chosen path
		{
			double maxDist=earthModel->DistanceForColumnDepthToPoint(pca+config.endcapLength*dir,dir,totalColumnDepth)-config.endcapLength;
			double actualColumnDepth=earthModel->GetColumnDepthInCGS(pca+config.endcapLength*dir,pca-maxDist*dir);
			if(actualColumnDepth<(totalColumnDepth-1)){ //if actually smaller, clip as needed, but for tiny differences we don't care
				log_debug_stream("Wanted column depth of " << totalColumnDepth << " but found only " << actualColumnDepth << " g/cm^2");
				totalColumnDepth=actualColumnDepth;
			}
		}
		//Choose how much of the total column depth this event should have to traverse
		double traversedColumnDepth=totalColumnDepth*random->Uniform();
		//endcapLength is subtracted so that dist==0 corresponds to pca
		double dist=earthModel->DistanceForColumnDepthToPoint(pca+config.endcapLength*dir,dir,totalColumnDepth-traversedColumnDepth)-config.endcapLength;
		
		{ //ensure that the point we picked is inside the atmosphere
			I3Position atmoEntry, atmoExit;
			int isect=GetIntersectionsWithSphere(earthModel->GetEarthCoordPosFromDetCoordPos(pca),
												 earthModel->GetEarthCoordDirFromDetCoordDir(dir),
												 earthModel->GetAtmoRadius(),atmoEntry,atmoExit);
			if(isect<2)
				log_fatal_stream("PCA not inside atmosphere: " << pca << " (" << earthModel->GetEarthCoordPosFromDetCoordPos(pca) << ')');
			atmoEntry=earthModel->GetDetCoordPosFromEarthCoordPos(atmoEntry);
			double atmoDist=(pca-atmoEntry).Magnitude();
			if(std::abs(dist-atmoDist)<100.0)
				dist=std::min(dist,atmoDist);
		}
		I3Position vertex=pca-dist*dir;
		
		//assemble the MCTree
		boost::shared_ptr<RangedEventProperties> properties(new RangedEventProperties);
		boost::shared_ptr<I3MCTree> mctree=FillTree(vertex,dir,energy,*properties);
		
		//set subclass properties
		properties->impactParameter=(pca-I3Position(0,0,0)).Magnitude();
		properties->totalColumnDepth=totalColumnDepth;
		
		//package up output and send it
		frame->Put(mctree);
		frame->Put("EventProperties",properties);
		PushFrame(frame);
		
		//update event count and check for completion
		eventsGenerated++;
		if(eventsGenerated==config.events && suspendOnCompletion)
			RequestSuspension();
	}
	
	I3_MODULE(RangedLeptonInjector);
	
	//-----------------------
	//Volume injection module
	
	VolumeLeptonInjector::VolumeLeptonInjector(const I3Context& context):
	LeptonInjectorBase(context,config){
		init();
	}
	
	VolumeLeptonInjector::VolumeLeptonInjector(const I3Context& context, VolumeInjectionConfiguration config_):
	LeptonInjectorBase(context,config),config(config_){
		init();
	}
	
	void VolumeLeptonInjector::init(){
		AddBaseParameters();
		AddParameter("CylinderRadius",
					 "Radius of the vertical cylinder around the origin within "
					 "which to place events",
					 config.cylinderRadius);
		AddParameter("CylinderHeight",
					 "Height of the vertical cylinder around the origin within "
					 "which to place events",
					 config.cylinderHeight);
		AddOutBox("OutBox");
	}
	
	void VolumeLeptonInjector::Configure(){
		BaseConfigure();
		GetParameter("CylinderRadius",config.cylinderRadius);
		GetParameter("CylinderHeight",config.cylinderHeight);
		if(config.cylinderRadius<0)
			log_fatal_stream(GetName() << ": CylinderRadius must be non-negative");
		if(config.cylinderHeight<0)
			log_fatal_stream(GetName() << ": CylinderHeight must be non-negative");
	}
	
	void VolumeLeptonInjector::DAQ(boost::shared_ptr<I3Frame> frame){
		//first, make sure configuration gets written once
		if(!wroteConfigFrame){
			boost::shared_ptr<I3Frame> sframe(new I3Frame('S'));
			boost::shared_ptr<VolumeInjectionConfiguration> sconfig(new VolumeInjectionConfiguration(config));
			sconfig->setCrossSection(getCrossSection(),getTotalCrossSection());
			sframe->Put("LeptonInjectorProperties",sconfig);
			PushFrame(sframe);
			wroteConfigFrame=true;
		}
		if(DoneGenerating()){
			PushFrame(frame);
			return;
		}
		
		//Choose an energy
		double energy=SampleEnergy();
		
		//Pick a direction on the sphere
		I3Direction dir(acos(random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum))),
						random->Uniform(config.azimuthMinimum,config.azimuthMaximum));
		log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
		
		//Pick a position in the xy-plane
		I3Position vertex=SampleFromDisk(config.cylinderRadius);
		//Add on the vertical component
		vertex.SetZ(random->Uniform(-config.cylinderHeight/2,config.cylinderHeight/2));
		log_trace_stream("vtx=(" << vertex.GetX() << ',' << vertex.GetY() << ',' << vertex.GetZ() << ')');
		
		//assemble the MCTree
		boost::shared_ptr<VolumeEventProperties> properties(new VolumeEventProperties);
		boost::shared_ptr<I3MCTree> mctree=FillTree(vertex,dir,energy,*properties);
		
		//set subclass properties
		properties->radius=vertex.GetRho();
		properties->z=vertex.GetZ();
	
		//package up output and send it
		frame->Put(mctree);
		frame->Put("EventProperties",properties);
		PushFrame(frame);
		
		//update event count and check for completion
		eventsGenerated++;
		if(eventsGenerated==config.events && suspendOnCompletion)
			RequestSuspension();
	}
	
	I3_MODULE(VolumeLeptonInjector);
	

	
	void MultiLeptonInjector::AddParameters(){
		AddOutBox("OutBox");
		AddParameter("MinimumEnergy",
					 "Minimum total event energy to inject",
					 rangedConfig.energyMinimum);
		AddParameter("MaximumEnergy",
					 "Maximum total event energy to inject",
					 rangedConfig.energyMaximum);
		AddParameter("PowerlawIndex",
					 "Powerlaw index of the energy spectrum to inject "
					 "(should be positive)",
					 rangedConfig.powerlawIndex);
		AddParameter("MinimumAzimuth",
					 "Minimum azimuth angle for injected events",
					 rangedConfig.azimuthMinimum);
		AddParameter("MaximumAzimuth",
					 "Maximum azimuth angle for injected events",
					 rangedConfig.azimuthMaximum);
		AddParameter("MinimumZenith",
					 "Minimum zenith angle for injected events",
					 rangedConfig.zenithMinimum);
		AddParameter("MaximumZenith",
					 "Maximum zenith angle for injected events",
					 rangedConfig.zenithMaximum);
		AddParameter("RandomService",
					 "Name of the random service to use",
					 "I3RandomService");
		
		AddParameter("InjectionRadius",
					 "Radius around the origin within which to target events",
					 rangedConfig.injectionRadius);
		AddParameter("EndcapLength",
					 "Length of the fixed endcaps add to each end of the distance "
					 "along which to sample interactions",
					 rangedConfig.endcapLength);
		AddParameter("EarthModel",
					 "Name of the Earth model service to use",
					 "");
		
		AddParameter("CylinderRadius",
					 "Radius of the vertical cylinder around the origin within "
					 "which to place events",
					 volumeConfig.cylinderRadius);
		AddParameter("CylinderHeight",
					 "Height of the vertical cylinder around the origin within "
					 "which to place events",
					 volumeConfig.cylinderHeight);
		
		AddParameter("Generators","The collection of configurations to generate",generatorSettings);
	}
	
	void MultiLeptonInjector::Configure(){
		try{
			//get the set of generators to be run
			GetParameter("Generators",generatorSettings);
			
			if(generatorSettings.empty())
				log_fatal_stream(GetName()+": There is no point in running this module without specfying at least one generator");
			
			//figure out whether there are any ranged or volume injectors
			bool hasRanged=false, hasVolume=false;
			for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=generatorSettings.begin(), end=generatorSettings.end(); genSet!=end; genSet++){
				hasRanged |= genSet->ranged;
				hasVolume |= !genSet->ranged;
			}
			
			//get the properties shared by all generators
			std::string randomServiceName;
			boost::shared_ptr<I3RandomService> random;
			//fetch each parameter directly into one configuration object,
			//and clone it into the other
			GetParameter("MinimumEnergy",rangedConfig.energyMinimum);
			volumeConfig.energyMinimum=rangedConfig.energyMinimum;
			GetParameter("MaximumEnergy",rangedConfig.energyMaximum);
			volumeConfig.energyMaximum=rangedConfig.energyMaximum;
			GetParameter("PowerlawIndex",rangedConfig.powerlawIndex);
			volumeConfig.powerlawIndex=rangedConfig.powerlawIndex;
			GetParameter("MinimumAzimuth",rangedConfig.azimuthMinimum);
			volumeConfig.azimuthMinimum=rangedConfig.azimuthMinimum;
			GetParameter("MaximumAzimuth",rangedConfig.azimuthMaximum);
			volumeConfig.azimuthMaximum=rangedConfig.azimuthMaximum;
			GetParameter("MinimumZenith",rangedConfig.zenithMinimum);
			volumeConfig.zenithMinimum=rangedConfig.zenithMinimum;
			GetParameter("MaximumZenith",rangedConfig.zenithMaximum);
			volumeConfig.zenithMaximum=rangedConfig.zenithMaximum;
			GetParameter("RandomService",randomServiceName);
			
			if(rangedConfig.energyMinimum<=0)
				log_fatal_stream(GetName() << ": minimum energy must be positive");
			if(rangedConfig.energyMaximum<=0)
				log_fatal_stream(GetName() << ": maximum energy must be positive");
			if(rangedConfig.energyMaximum<rangedConfig.energyMinimum)
				log_fatal_stream(GetName() << ": maximum energy must be greater than or equal to minimum energy");
			if(rangedConfig.azimuthMinimum<0.0)
				log_fatal_stream(GetName() << ": minimum azimuth angle must be greater than or equal to zero");
			if(rangedConfig.azimuthMaximum>2*constants::pi<double>())
				log_fatal_stream(GetName() << ": maximum azimuth angle must be less than or equal to 2 pi");
			if(rangedConfig.azimuthMinimum>rangedConfig.azimuthMaximum)
				log_fatal_stream(GetName() << ": minimum azimuth angle must be less than or equal to maximum azimuth angle");
			if(rangedConfig.zenithMinimum<0.0)
				log_fatal_stream(GetName() << ": minimum zenith angle must be greater than or equal to zero");
			if(rangedConfig.zenithMaximum>constants::pi<double>())
				log_fatal_stream(GetName() << ": maximum zenith angle must be less than or equal to pi");
			if(rangedConfig.zenithMinimum>rangedConfig.zenithMaximum)
				log_fatal_stream(GetName() << ": minimum zenith angle must be less than or equal to maximum zenith angle");
			random = context_.Get<boost::shared_ptr<I3RandomService> >(randomServiceName);
			if(!random)
				log_fatal_stream(GetName() << ": A random service is required");
			
			innerContext.Put(random,randomServiceName);
			
			//get the properties for ranged injectors
			std::string earthModelName;
			boost::shared_ptr<earthmodel::EarthModelService> earthModel;
			if(hasRanged){
				GetParameter("InjectionRadius",rangedConfig.injectionRadius);
				GetParameter("EndcapLength",rangedConfig.endcapLength);
				GetParameter("EarthModel",earthModelName);
				
				if(rangedConfig.injectionRadius<0)
					log_fatal_stream(GetName() << ": InjectionRadius must be non-negative");
				if(rangedConfig.endcapLength<0)
					log_fatal_stream(GetName() << ": EndcapLength must be non-negative");
				earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
				if(!earthModel)
					log_fatal_stream(GetName() << ": an Earth model service is required");
				
				innerContext.Put(earthModel,earthModelName);
			}
			
			//get the properties for volume injectors
			if(hasVolume){
				GetParameter("CylinderRadius",volumeConfig.cylinderRadius);
				GetParameter("CylinderHeight",volumeConfig.cylinderHeight);
				
				if(volumeConfig.cylinderRadius<0)
					log_fatal_stream(GetName() << ": CylinderRadius must be non-negative");
				if(volumeConfig.cylinderHeight<0)
					log_fatal_stream(GetName() << ": CylinderHeight must be non-negative");
			}
			
			//construct all generators
			unsigned int i=0;
			for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=generatorSettings.begin(), end=generatorSettings.end(); genSet!=end; genSet++){
				log_debug_stream("Configuring injector " << i << ":");
				LeptonInjectorBase* generator=NULL;
				try{
					if(genSet->ranged){
						log_debug_stream(" this is a ranged injector");
						generator=new RangedLeptonInjector(innerContext,rangedConfig);
						generator->GetConfiguration().Set("EarthModel",boost::python::object(earthModelName));
					}
					else{ //volume
						log_debug_stream(" this is a volume injector");
						generator=new VolumeLeptonInjector(innerContext,volumeConfig);
					}
					
					//set properties not shared with other injectors, or which are not part of the config object
					generator->GetConfiguration().Set("NEvents",boost::python::object(genSet->events));
					generator->GetConfiguration().Set("FinalType1",boost::python::object(genSet->finalType1));
					generator->GetConfiguration().Set("FinalType2",boost::python::object(genSet->finalType2));
					generator->GetConfiguration().Set("RandomService",boost::python::object(randomServiceName));
					generator->GetConfiguration().Set("DoublyDifferentialCrossSectionFile",boost::python::object(genSet->crossSectionPath));
					generator->GetConfiguration().Set("TotalCrossSectionFile",boost::python::object(genSet->totalCrossSectionPath));
					generator->GetConfiguration().Set("SuspendOnCompletion",boost::python::object(false));
					
					generator->SetName(GetName()+"_Generator_"+boost::lexical_cast<std::string>(i++));
					generator->Configure();
				}catch(...){
					delete generator;
					throw;
				}
				generators.push_back(generator);
			}
			
			//bind the first generator to the collector
			generators.front()->ConnectOutBox("OutBox",collector);
		}catch(std::runtime_error& err){
			throw std::runtime_error("While configuring "+GetName()+":\n"+err.what());
		}
	}
	
	void MultiLeptonInjector::DAQ(boost::shared_ptr<I3Frame> frame){
		if(generators.empty()){
			RequestSuspension();
			return;
		}
		while(generators.front()->DoneGenerating()){
			delete generators.front();
			generators.pop_front();
			//if there are no more generators, we are done
			if(generators.empty()){
				RequestSuspension();
				return;
			}
			//bind the next generator to the collector
			generators.front()->ConnectOutBox("OutBox",collector);
		}
		generators.front()->DAQ(frame);
		collector->Process();
		while(!collector->output.empty()){
			PushFrame(collector->output.front());
			collector->output.pop();
		}
	}
		
	I3_MODULE(MultiLeptonInjector);
	
} //namespace LeptonInjector
