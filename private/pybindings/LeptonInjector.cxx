#include <LeptonInjector.h>
#include <Controller.h>
#include <Random.h>

// #include <converter/LeptonInjectionConfigurationConverter.h>
#include <boost/python.hpp>


using namespace boost::python;

BOOST_PYTHON_MODULE(LeptonInjector){
	using namespace LeptonInjector;

    class_<LI_random>("LI_random",init<unsigned int>(args("seed")));

    class_<Controller>("Controller", init<std::vector<MinimalInjectionConfiguration>,double,double,double,double,double,double,double,double,double,double,double>(
          (args("injectors"),args("minimum energy"),args("maximum energy"),args("spectral index"),args("minimum azimuth"),args("maximum azimuth"),args("minimum zenith"),args("maximum zenith"),args("injection radius"), args("endcap length"), args("cylinder radius"), args("cylinder height"))
          ) 
     )
     .def("Execute",&Controller::Execute)
     .def("AddInjector",&Controller::AddInjector)
     .def("Output",&Controller::NameOutfile)
     ;

    
//    class_<ParticleType>("Particle");
    

	class_<BasicInjectionConfiguration>("BasicInjectionConfiguration")
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
	
	class_<RangedInjectionConfiguration,bases<BasicInjectionConfiguration> >("RangedInjectionConfiguration")
	.def_readonly("injectionRadius",&RangedInjectionConfiguration::injectionRadius)
	.def_readonly("endcapLength",&RangedInjectionConfiguration::endcapLength)
	;
	
	class_<VolumeInjectionConfiguration,bases<BasicInjectionConfiguration> >("VolumeInjectionConfiguration")
	.def_readonly("cylinderRadius",&VolumeInjectionConfiguration::cylinderRadius)
	.def_readonly("cylinderHeight",&VolumeInjectionConfiguration::cylinderHeight)
	;
	
	class_<MinimalInjectionConfiguration, std::shared_ptr<MinimalInjectionConfiguration>>("injector",
	  init<unsigned int,ParticleType,ParticleType,std::string,std::string,bool>(
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
	

	class_<BasicEventProperties>("BasicEventProperties")
	.def_readonly("totalEnergy",&BasicEventProperties::totalEnergy)
	.def_readonly("zenith",&BasicEventProperties::zenith)
	.def_readonly("azimuth",&BasicEventProperties::azimuth)
	.def_readonly("finalStateX",&BasicEventProperties::finalStateX)
	.def_readonly("finalStateY",&BasicEventProperties::finalStateY)
	.def_readonly("finalType1",&BasicEventProperties::finalType1)
	.def_readonly("finalType2",&BasicEventProperties::finalType2)
	.def_readonly("initialType",&BasicEventProperties::initialType)
	;
	
	class_<RangedEventProperties,bases<BasicEventProperties> >("RangedEventProperties")
	.def_readonly("impactParameter",&RangedEventProperties::impactParameter)
	.def_readonly("totalColumnDepth",&RangedEventProperties::totalColumnDepth)
	;
	
	class_<VolumeEventProperties,bases<BasicEventProperties> >("VolumeEventProperties")
	.def_readonly("radius",&VolumeEventProperties::radius)
	.def_readonly("z",&VolumeEventProperties::z)
	;

}
