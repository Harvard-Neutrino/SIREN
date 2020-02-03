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
        .def("LICFile",&Controller::NameLicFile)
        .def("Seed",&Controller::setSeed)
     ;

    
    enum_<Particle::ParticleType>("Particle")
        .value("EPlus",Particle::EPlus)
        .value("EMinus",Particle::EMinus)
        .value("MuPlus",Particle::MuPlus)
        .value("MuMinus",Particle::MuMinus)
        .value("TauPlus",Particle::TauPlus)
        .value("TauMinus",Particle::TauMinus)
        .value("NuE",Particle::NuE)
        .value("NuEBar",Particle::NuEBar)
        .value("NuMuBar",Particle::NuMuBar)
        .value("NuTau",Particle::NuTau)
        .value("NuTauBar",Particle::NuTauBar)
        .value("NuMu",Particle::NuMu)
        .value("Hadrons",Particle::Hadrons)
    ;
   
    class_<MinimalInjectionConfiguration, std::shared_ptr<MinimalInjectionConfiguration>>("injector",
	  init<unsigned int,Particle::ParticleType,Particle::ParticleType,std::string,std::string,bool>(
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
