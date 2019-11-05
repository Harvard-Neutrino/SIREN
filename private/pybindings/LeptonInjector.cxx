#include <LeptonInjector.h>
#include <Controller.h>
#include <Random.h>
#include <Constants.h>


// #include <converter/LeptonInjectionConfigurationConverter.h>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(pylepton_injector){
	using namespace LeptonInjector;
  

    class_<LI_random>("RNG",init<unsigned int>(args("seed")=1))
        .def("Uniform",&LI_random::Uniform)
    ;

    class_<Controller>("Controller", init<MinimalInjectionConfiguration,double,double,double,double,double,double,double,double,double,double,double>(  
        (args("injectors"),args("minimum energy"),args("maximum energy"),args("spectral index"),args("minimum azimuth"),args("maximum azimuth"),args("minimum zenith"),args("maximum zenith"),args("injection radius")=1200., args("endcap length")=1200., args("cylinder radius")=1200., args("cylinder height")=1200.))
        )
        .def("Execute",&Controller::Execute)
        .def("AddInjector",&Controller::AddInjector)
        .def("Output",&Controller::NameOutfile)
     ;

    
    enum_<ParticleType>("Particle")
        .value("EPlus",ParticleType::EPlus)
        .value("EMinus",ParticleType::EMinus)
        .value("MuPlus",ParticleType::MuPlus)
        .value("MuMinus",ParticleType::MuMinus)
        .value("TauPlus",ParticleType::TauPlus)
        .value("TauMinus",ParticleType::TauMinus)
        .value("NuE",ParticleType::NuE)
        .value("NuEBar",ParticleType::NuEBar)
        .value("NuMuBar",ParticleType::NuMuBar)
        .value("NuTau",ParticleType::NuTau)
        .value("NuTauBar",ParticleType::NuTauBar)
        .value("NuMu",ParticleType::NuMu)
        .value("Hadrons",ParticleType::Hadrons)
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
}
