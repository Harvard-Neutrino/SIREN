#include <LeptonInjector/LeptonInjector.h>

#include <fstream>
#include <tableio/I3TableRowDescription.h>

#include <boost/detail/endian.hpp>

//put all of the serialization boilerplate here in one place

namespace LeptonInjector{

	template <typename Archive>
	void BasicInjectionConfiguration::serialize(Archive &ar, unsigned version){
		ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
		
		ar & make_nvp("generatedEvents", events);
		ar & make_nvp("energyMinimum", energyMinimum);
		ar & make_nvp("energyMaximum", energyMaximum);
		ar & make_nvp("powerlawIndex", powerlawIndex);
		ar & make_nvp("azimuthMinimum", azimuthMinimum);
		ar & make_nvp("azimuthMaximum", azimuthMaximum);
		ar & make_nvp("zenithMinimum", zenithMinimum);
		ar & make_nvp("zenithMaximum", zenithMaximum);
		ar & make_nvp("finalType1", finalType1);
		ar & make_nvp("finalType2", finalType2);
		ar & make_nvp("crossSection", crossSectionBlob);
		ar & make_nvp("totalCrossSection", totalCrossSectionBlob);
	}
	
	template <typename Archive>
	void RangedInjectionConfiguration::serialize(Archive &ar, unsigned version){
		ar & make_nvp("BasicInjectionConfiguration",
					  base_object<BasicInjectionConfiguration>(*this));
		
		ar & make_nvp("injectionRadius", injectionRadius);
		ar & make_nvp("endcapLength", endcapLength);
	}
	
	template <typename Archive>
	void VolumeInjectionConfiguration::serialize(Archive &ar, unsigned version){
		ar & make_nvp("BasicInjectionConfiguration",
					  base_object<BasicInjectionConfiguration>(*this));
		
		ar & make_nvp("cylinderRadius", cylinderRadius);
		ar & make_nvp("cylinderHeight", cylinderHeight);
	}
	
	template <typename Archive>
	void BasicEventProperties::serialize(Archive &ar, unsigned version){
		ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
		
		ar & make_nvp("totalEnergy", totalEnergy);
		ar & make_nvp("zenith", zenith);
		ar & make_nvp("azimuth", azimuth);
		ar & make_nvp("finalStateX", finalStateX);
		ar & make_nvp("finalStateY", finalStateY);
		ar & make_nvp("finalType1", finalType1);
		ar & make_nvp("finalType2", finalType2);
		ar & make_nvp("initialType", initialType);
	}
	
	template <typename Archive>
	void RangedEventProperties::serialize(Archive &ar, unsigned version){
		ar & make_nvp("BasicEventProperties",
					  base_object<BasicEventProperties>(*this));
		
		ar & make_nvp("impactParameter", impactParameter);
		ar & make_nvp("totalColumnDepth", totalColumnDepth);
	}
	
	template <typename Archive>
	void VolumeEventProperties::serialize(Archive &ar, unsigned version){
		ar & make_nvp("BasicEventProperties",
					  base_object<BasicEventProperties>(*this));
		
		ar & make_nvp("radius", radius);
		ar & make_nvp("z", z);
	}
	
	template<typename T>
	struct endian_adapter{
		const T& t;
		endian_adapter(const T& t):t(t){}
	};
	
	template<typename T>
	endian_adapter<T> little_endian(const T& t){ return(endian_adapter<T>(t)); }
	
	std::ostream& endianWrite(std::ostream& os, const char* data, size_t dataSize){
#if defined(BOOST_LITTLE_ENDIAN)
		//just write bytes
		os.write(data,dataSize);
#elif defined(BOOST_BIG_ENDIAN)
		//write bytes in reverse order
		for(size_t i=1; i<=dataSize; i++)
			os.write(data+dataSize-i),1);
#elif defined(BOOST_PDP_ENDIAN)
		//complain bitterly
		#error PDP-endian systems are not supported.
#else
		#error Unable to determine machine endianness!
#endif
		return(os);
	}
	
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const endian_adapter<T>& e){
		return(endianWrite(os,(char*)&e.t,sizeof(T)));
	}
	
	//automatically write string's length before its data
	std::ostream& operator<<(std::ostream& os, const endian_adapter<std::string>& e){
		os << little_endian(e.t.size());
		return(endianWrite(os,e.t.c_str(),e.t.size()));
	}
	
	//same for vector<char>
	std::ostream& operator<<(std::ostream& os, const endian_adapter<std::vector<char> >& e){
		os << little_endian(e.t.size());
		return(endianWrite(os,&e.t[0],e.t.size()));
	}
	
	/*
	 Base:
	 uint64_t: block length
	 size_t: name length
	 char[name length]: block type name
	 uint8_t: version number
	 char[block length-17-name length]: block data
	 */
	void writeBlockHeader(std::ostream& os, uint64_t blockDataSize,
						  const std::string& blockTypeName, uint8_t version){
		size_t nameLen=blockTypeName.size();
		uint64_t totalBlockSize=8+sizeof(nameLen)+nameLen+1+blockDataSize;
		os << little_endian(totalBlockSize)
		   << little_endian(blockTypeName)
		   << little_endian(version);
		if(!os.good())
			log_fatal("Writing block header failed");
	}
	
	/*
	 EnumDef:
	 block data:
	 uint8_t: enum name length
	 char[enum name length]: enum name
	 uint32_t: number of enumerators
	 int64_t: enumerator value
	 uint8_t: enumerator name length
	 char[enumerator name length]: enumerator name
	 */
	template<typename Enum>
	void writeEnumDefBlock(std::ostream& os, const std::string& enumName,
						   const std::vector<std::pair<std::string,Enum> >& enumerators){
		//compute data size
		uint64_t dataSize=0;
		size_t nameLen=enumName.size();
		dataSize+=sizeof(nameLen);
		dataSize+=nameLen;
		if(enumerators.size()>=(1ULL<<32))
			log_fatal_stream("Number of enumerators (" << enumerators.size()
							 << ") too large");
		uint32_t numEnum=enumerators.size();
		dataSize+=sizeof(numEnum);
		for(size_t i=0; i<numEnum; i++){
			dataSize+=8+sizeof(size_t)+enumerators[i].first.size();
		}
		//write header
		writeBlockHeader(os,dataSize,"EnumDef",1);
		//write data
		os << little_endian(enumName) << little_endian(numEnum);
		for(size_t i=0; i<numEnum; i++){
			os << little_endian((int64_t)enumerators[i].second)
			   << little_endian(enumerators[i].first);
		}
		if(!os.good())
			log_fatal("Writing enum def block failed");
	}
	
	void writeConfiguration(std::ostream& os, const RangedInjectionConfiguration& config){
		//compute data size
		uint64_t dataSize=0;
		dataSize+=4; //events
		dataSize+=8; //energyMinimum
		dataSize+=8; //energyMaximum
		dataSize+=8; //powerlawIndex
		dataSize+=8; //azimuthMinimum
		dataSize+=8; //azimuthMaximum
		dataSize+=8; //zenithMinimum
		dataSize+=8; //zenithMaximum
		dataSize+=sizeof(I3Particle::ParticleType); //finalType1
		dataSize+=sizeof(I3Particle::ParticleType); //finalType2
		dataSize+=sizeof(size_t); //crossSection size
		dataSize+=config.crossSectionBlob.size(); //crossSection
        dataSize+=sizeof(size_t);
		dataSize+=config.totalCrossSectionBlob.size(); //crossSection
		dataSize+=8; //injectionRadius
		dataSize+=8; //endcapLength
		//write header
		writeBlockHeader(os,dataSize,"RangedInjectionConfiguration",1);
		//write data
		os << little_endian(config.events)
		   << little_endian(config.energyMinimum/I3Units::GeV)
		   << little_endian(config.energyMaximum/I3Units::GeV)
		   << little_endian(config.powerlawIndex)
		   << little_endian(config.azimuthMinimum/I3Units::radian)
		   << little_endian(config.azimuthMaximum/I3Units::radian)
		   << little_endian(config.zenithMinimum/I3Units::radian)
		   << little_endian(config.zenithMaximum/I3Units::radian)
		   << little_endian(config.finalType1)
		   << little_endian(config.finalType2)
		   << little_endian(config.crossSectionBlob)
		   << little_endian(config.totalCrossSectionBlob)
		   << little_endian(config.injectionRadius/I3Units::meter)
		   << little_endian(config.endcapLength/I3Units::meter);
		if(!os.good())
			log_fatal("Writing ranged injection config block failed");
	}
	
	void writeConfiguration(std::ostream& os, const VolumeInjectionConfiguration& config){
		//compute data size
		uint64_t dataSize=0;
		dataSize+=4; //events
		dataSize+=8; //energyMinimum
		dataSize+=8; //energyMaximum
		dataSize+=8; //powerlawIndex
		dataSize+=8; //azimuthMinimum
		dataSize+=8; //azimuthMaximum
		dataSize+=8; //zenithMinimum
		dataSize+=8; //zenithMaximum
		dataSize+=sizeof(I3Particle::ParticleType); //finalType1
		dataSize+=sizeof(I3Particle::ParticleType); //finalType2
		dataSize+=sizeof(size_t); //crossSection size
		dataSize+=config.crossSectionBlob.size(); //crossSection
        dataSize+=sizeof(size_t);
		dataSize+=config.totalCrossSectionBlob.size(); //crossSection
		dataSize+=8; //cylinderRadius
		dataSize+=8; //cylinderHeight
		//write header
		writeBlockHeader(os,dataSize,"VolumeInjectionConfiguration",1);
		//write data
		os << little_endian(config.events)
		   << little_endian(config.energyMinimum/I3Units::GeV)
		   << little_endian(config.energyMaximum/I3Units::GeV)
		   << little_endian(config.powerlawIndex)
		   << little_endian(config.azimuthMinimum/I3Units::radian)
		   << little_endian(config.azimuthMaximum/I3Units::radian)
		   << little_endian(config.zenithMinimum/I3Units::radian)
		   << little_endian(config.zenithMaximum/I3Units::radian)
		   << little_endian(config.finalType1)
		   << little_endian(config.finalType2)
		   << little_endian(config.crossSectionBlob)
		   << little_endian(config.totalCrossSectionBlob)
		   << little_endian(config.cylinderRadius/I3Units::meter)
		   << little_endian(config.cylinderHeight/I3Units::meter);
		if(!os.good())
			log_fatal("Writing volume injection config block failed");
	}
	
	class InjectionConfigSerializer : public I3Module{
	public:
		InjectionConfigSerializer(const I3Context& ctx):
		I3Module(ctx),wroteParticleTypes(false){
			Register(I3Frame::Stream('S'),&InjectionConfigSerializer::S);
			AddParameter("OutputPath","");
			AddOutBox("OutBox");
		}
		
		void Configure(){
			std::string outputPath;
			GetParameter("OutputPath",outputPath);
			output.open(outputPath.c_str());
			if(!output.good())
				log_fatal_stream("Failed to open " << outputPath << " for writing");
		}
		
		void S(boost::shared_ptr<I3Frame> frame){
			if(frame->Has("LeptonInjectorProperties")){
				if(!wroteParticleTypes){
					MAKE_ENUM_VECTOR(type,I3Particle,I3Particle::ParticleType,I3PARTICLE_H_I3Particle_ParticleType);
					writeEnumDefBlock(output, "I3Particle::ParticleType", type);
					wroteParticleTypes=true;
				}
				boost::shared_ptr<const I3FrameObject> config
				  =frame->Get<boost::shared_ptr<const I3FrameObject> >("LeptonInjectorProperties");
				boost::shared_ptr<const RangedInjectionConfiguration> rconfig;
				boost::shared_ptr<const VolumeInjectionConfiguration> vconfig;
				if((rconfig=boost::dynamic_pointer_cast<const RangedInjectionConfiguration>(config)))
					writeConfiguration(output, *rconfig);
				else if((vconfig=boost::dynamic_pointer_cast<const VolumeInjectionConfiguration>(config)))
					writeConfiguration(output, *vconfig);
				else
					log_fatal("Unable to write LeptonInjectorProperties with unrecognized type");
			}
			PushFrame(frame);
		}
		
	private:
		std::ofstream output;
		bool wroteParticleTypes;
	};
	
	I3_MODULE(InjectionConfigSerializer);
	
} //namespace LeptonInjector

I3_SERIALIZABLE(LeptonInjector::RangedInjectionConfiguration);
I3_SERIALIZABLE(LeptonInjector::BasicInjectionConfiguration);
I3_SERIALIZABLE(LeptonInjector::VolumeInjectionConfiguration);
I3_SERIALIZABLE(LeptonInjector::BasicEventProperties);
I3_SERIALIZABLE(LeptonInjector::RangedEventProperties);
I3_SERIALIZABLE(LeptonInjector::VolumeEventProperties);
