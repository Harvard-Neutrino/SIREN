#include <icetray/load_project.h>
#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/stream_to_string.hpp>
#include <dataclasses/ostream_overloads.hpp>
#include <tableio/converter/pybindings.h>
#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/converter/LeptonInjectionConfigurationConverter.h>

using namespace boost::python;

BOOST_PYTHON_MODULE(LeptonInjector){
	using namespace LeptonInjector;
	load_project("libLeptonInjector", false);
	load_project("libtableio", false);
	
	class_<BasicInjectionConfiguration,bases<I3FrameObject> >("BasicInjectionConfiguration")
	.def_readonly("events",&BasicInjectionConfiguration::events)
	.def_readonly("energyMinimum",&BasicInjectionConfiguration::energyMinimum)
	.def_readonly("energyMaximum",&BasicInjectionConfiguration::energyMaximum)
	.def_readonly("powerlawIndex",&BasicInjectionConfiguration::powerlawIndex)
	.def_readonly("azimuthMinimum",&BasicInjectionConfiguration::azimuthMinimum)
	.def_readonly("azimuthMaximum",&BasicInjectionConfiguration::azimuthMaximum)
	.def_readonly("zenithMinimum",&BasicInjectionConfiguration::zenithMinimum)
	.def_readonly("zenithMaximum",&BasicInjectionConfiguration::zenithMaximum)
	.def_readonly("finalType1",&BasicInjectionConfiguration::finalType1)
	.def_readonly("finalType2",&BasicInjectionConfiguration::finalType2)
	;
	
	class_<RangedInjectionConfiguration,bases<BasicInjectionConfiguration,I3FrameObject> >("RangedInjectionConfiguration")
	.def_readonly("injectionRadius",&RangedInjectionConfiguration::injectionRadius)
	.def_readonly("endcapLength",&RangedInjectionConfiguration::endcapLength)
	;
	
	class_<VolumeInjectionConfiguration,bases<BasicInjectionConfiguration,I3FrameObject> >("VolumeInjectionConfiguration")
	.def_readonly("cylinderRadius",&VolumeInjectionConfiguration::cylinderRadius)
	.def_readonly("cylinderHeight",&VolumeInjectionConfiguration::cylinderHeight)
	;
	
	class_<MinimalInjectionConfiguration, boost::shared_ptr<MinimalInjectionConfiguration> >("injector",
	  init<unsigned int,I3Particle::ParticleType,I3Particle::ParticleType,std::string,std::string,bool>(
	    (args("NEvents"),args("FinalType1"),args("FinalType2"),args("DoublyDifferentialCrossSectionFile"),args("TotalCrossSectionFile"),args("Ranged"))
	  )
	)
	.def_readwrite("events",&MinimalInjectionConfiguration::events)
	.def_readwrite("finalType1",&MinimalInjectionConfiguration::finalType1)
	.def_readwrite("finalType2",&MinimalInjectionConfiguration::finalType2)
	.def_readwrite("crossSectionPath",&MinimalInjectionConfiguration::crossSectionPath)
	.def_readwrite("totalCrossSectionPath",&MinimalInjectionConfiguration::totalCrossSectionPath)
	.def_readwrite("ranged",&MinimalInjectionConfiguration::ranged)
	;
	
	class_<std::vector<MinimalInjectionConfiguration> >("MultiConfigList")
	.def(vector_indexing_suite<std::vector<MinimalInjectionConfiguration> >())
	;
	from_python_sequence<std::vector<MinimalInjectionConfiguration>, variable_capacity_policy>();
	
	class_<BasicEventProperties,bases<I3FrameObject> >("BasicEventProperties")
	.def_readonly("totalEnergy",&BasicEventProperties::totalEnergy)
	.def_readonly("zenith",&BasicEventProperties::zenith)
	.def_readonly("azimuth",&BasicEventProperties::azimuth)
	.def_readonly("finalStateX",&BasicEventProperties::finalStateX)
	.def_readonly("finalStateY",&BasicEventProperties::finalStateY)
	.def_readonly("finalType1",&BasicEventProperties::finalType1)
	.def_readonly("finalType2",&BasicEventProperties::finalType2)
	.def_readonly("initialType",&BasicEventProperties::initialType)
	;
	
	class_<RangedEventProperties,bases<BasicEventProperties,I3FrameObject> >("RangedEventProperties")
	.def_readonly("impactParameter",&RangedEventProperties::impactParameter)
	.def_readonly("totalColumnDepth",&RangedEventProperties::totalColumnDepth)
	;
	
	class_<VolumeEventProperties,bases<BasicEventProperties,I3FrameObject> >("VolumeEventProperties")
	.def_readonly("radius",&VolumeEventProperties::radius)
	.def_readonly("z",&VolumeEventProperties::z)
	;

	{
	I3CONVERTER_NAMESPACE(LeptonInjector);
	I3CONVERTER_EXPORT_DEFAULT(EventPropertiesConverter,"Converts an EventProperties");
	}
}
