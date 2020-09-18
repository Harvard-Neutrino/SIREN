#include <string>
#include <vector>
#include <utility>

#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Controller.h>
#include <LeptonInjector/Random.h>
#include <LeptonInjector/Constants.h>


// #include <converter/LeptonInjectionConfigurationConverter.h>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "array_ref.h"
#include "array_indexing_suite.h"

using namespace boost::python;

namespace LeptonInjector{
    struct LIConstants {
    };
}

BOOST_PYTHON_MODULE(LeptonInjector){
	using namespace LeptonInjector;

    class_<array_ref<double>>( "double_array" )
        .def( array_indexing_suite<array_ref<double>>() )
    ;

    class_<std::pair<double, double> >("DoublePair")
    	.def_readwrite("first", &std::pair<double, double>::first)
    	.def_readwrite("second", &std::pair<double, double>::second)
	;

    class_<LI_random>("RNG",init<unsigned int>(args("seed")=1))
        .def("Uniform",&LI_random::Uniform)
    ;

    class_<LI_Position>("LI_Position", init<>())
        .def(init<double, double, double>())
        .def(init<const LI_Position&>())
        .def(init<std::array<double,LeptonInjector::n_dimensions>>())
        .def("at",&LI_Position::at)
        .def("Magnitude",&LI_Position::Magnitude)
        .def("GetZ",&LI_Position::GetZ)
        .def("GetY",&LI_Position::GetY)
        .def("GetX",&LI_Position::GetX)
        .def("SetZ",&LI_Position::SetZ)
        .def("SetY",&LI_Position::SetY)
        .def("SetX",&LI_Position::SetX)
        .def(self + self)
        .def(self - self)
        .def(self += self)
        .def(self == self)
        .def(self * double())
        .def(double() * self)
        .def(self * LI_Direction())
        .def(self * self)
        ;

    class_<LI_Direction>("LI_Direction", init<>())
        .def(init<double, double>())
        .def(init<std::array<double, 2>>())
        .def(init<std::pair<double, double>>())
        .def(init<const LI_Direction&>())
        .def(init<const LI_Position&>())
        .def("GetZ",&LI_Direction::GetZ)
        .def("GetY",&LI_Direction::GetY)
        .def("GetX",&LI_Direction::GetX)
        .def_readwrite("zenith", &LI_Direction::zenith)
        .def_readwrite("azimuth", &LI_Direction::azimuth)
        .def(-self)
        .def(self * double())
        .def(double() * self)
        .def(self * LI_Position())
        ;

    def("RotateX", &RotateX);
    def("RotateY", &RotateY);
    def("RotateZ", &RotateZ);
    def("rotateRelative", &rotateRelative);
    def("computeCylinderIntersections", &computeCylinderIntersections);


    class_<Controller>("Controller", init<Injector,double,double,double,double,double,double,double,double,double,double,double>(
        (args("injectors"),args("minimum energy"),args("maximum energy"),args("spectral index"),args("minimum azimuth"),args("maximum azimuth"),args("minimum zenith"),args("maximum zenith"),args("injection radius")=1200., args("endcap length")=1200., args("cylinder radius")=1200., args("cylinder height")=1200.))
        )
        .def("Execute",&Controller::Execute)
        .def("AddInjector",&Controller::AddInjector)
        .def("NameOutfile",&Controller::NameOutfile)
        .def("NameLicFile",&Controller::NameLicFile)
        .def("Overwrite",&Controller::Overwrite)
        .def("setSeed",&Controller::setSeed)
     ;


    {
    scope particle = class_<Particle>("Particle", init<>())
        .def_readwrite("type", &Particle::type)
        .def_readwrite("energy", &Particle::energy)
        .def_readwrite("direction", &Particle::direction)
        .add_property( "position",
                +[](Particle *obj) {
                    return array_ref<double>( obj->position );
                })
        .def("GetMass", &Particle::GetMass)
        .def("HasMass", &Particle::HasMass)
        .def("GetTypeString", &Particle::GetTypeString)
        ;

    enum_<Particle::ParticleType>("ParticleType")
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
    }

    def("isLepton", &isLepton);
    def("isCharged", &isCharged);
    def("particleName", &particleName);
    def("particleMass", &particleMass);
    def("kineticEnergy", &kineticEnergy);
    def("particleSpeed", &particleSpeed);
    def("decideShape", &decideShape);
    def("deduceInitialType", &deduceInitialType);
    def("getInteraction", &getInteraction);

    class_<Injector, std::shared_ptr<Injector>>("Injector",
	  init<unsigned int,Particle::ParticleType,Particle::ParticleType,std::string,std::string,bool>(
	    (args("NEvents"),args("FinalType1"),args("FinalType2"),args("DoublyDifferentialCrossSectionFile"),args("TotalCrossSectionFile"),args("Ranged"))
	  )
	)
	.def_readwrite("events",&Injector::events)
	.def_readwrite("finalType1",&Injector::finalType1)
	.def_readwrite("finalType2",&Injector::finalType2)
	.def_readwrite("crossSectionPath",&Injector::crossSectionPath)
	.def_readwrite("totalCrossSectionPath",&Injector::totalCrossSectionPath)
	.def_readwrite("ranged",&Injector::ranged)
	;

    scope constants = class_<LIConstants>("Constants");

    constants.attr("pi") = Constants::pi;
	constants.attr("tau") = Constants::tau;
	constants.attr("degrees") = Constants::degrees;
	constants.attr("deg") = Constants::deg;
	constants.attr("radian") = Constants::radian;
	constants.attr("m") = Constants::m;
	constants.attr("meter") = Constants::meter;
	constants.attr("cm") = Constants::cm;
	constants.attr("centimeter") = Constants::centimeter;
	constants.attr("second") = Constants::second;
	constants.attr("c") = Constants::c;
	constants.attr("protonMass") = Constants::protonMass;
	constants.attr("neutronMass") = Constants::neutronMass;
	constants.attr("isoscalarMass") = Constants::isoscalarMass;
	constants.attr("electronMass") = Constants::electronMass;
	constants.attr("muonMass") = Constants::muonMass;
	constants.attr("tauMass") = Constants::tauMass;
	constants.attr("s") = Constants::s;
	constants.attr("tauLifeTime") = Constants::tauLifeTime;
	constants.attr("MuonLifeTime") = Constants::MuonLifeTime;
	constants.attr("wMass") = Constants::wMass;
	constants.attr("wWidth") = Constants::wWidth;
	constants.attr("zMass") = Constants::zMass;
	constants.attr("WBranchE") = Constants::WBranchE;
	constants.attr("WBranchMuon") = Constants::WBranchMuon;
	constants.attr("WBranchTau") = Constants::WBranchTau;
	constants.attr("WBranchHadronic") = Constants::WBranchHadronic;
	constants.attr("nuEMass") = Constants::nuEMass;
	constants.attr("nuMuMass") = Constants::nuMuMass;
	constants.attr("nuTauMass") = Constants::nuTauMass;
	constants.attr("GeV") = Constants::GeV;
	constants.attr("EeV") = Constants::EeV;
	constants.attr("PeV") = Constants::PeV;
	constants.attr("TeV") = Constants::TeV;
	constants.attr("MeV") = Constants::MeV;
	constants.attr("keV") = Constants::keV;
	constants.attr("eV") = Constants::eV;
	constants.attr("Joule") = Constants::Joule;
	constants.attr("FermiConstant") = Constants::FermiConstant;
	constants.attr("avogadro") = Constants::avogadro;
	constants.attr("thetaWeinberg") = Constants::thetaWeinberg;
	constants.attr("gravConstant") = Constants::gravConstant;
	constants.attr("fineStructure") = Constants::fineStructure;
}
