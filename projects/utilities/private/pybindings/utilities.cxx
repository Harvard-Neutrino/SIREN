
#include <vector>

#include "../../public/SIREN/utilities/Random.h"
#include "../../public/SIREN/utilities/Constants.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(utilities,m) {
  using namespace siren::utilities;

  class_<SIREN_random, std::shared_ptr<SIREN_random>>(m, "SIREN_random")
    .def(init<>())
    .def(init<unsigned int>())
    .def("Uniform",&SIREN_random::Uniform)
    .def("set_seed",&SIREN_random::set_seed);

  module Constants = m.def_submodule("Constants");
  // geometry
  Constants.attr("pi") = &Constants::pi;
  Constants.attr("tau") = &Constants::tau;
  Constants.attr("degrees") = &Constants::degrees;
  // This is used since the user will enter a number in degrees, while the c trigonometric functions
  // expect angles presented in radians.
  Constants.attr("deg") = &Constants::deg;
  Constants.attr("radian") = &Constants::radian;

  // meter
  Constants.attr("m") = &Constants::m;
  Constants.attr("meter") = &Constants::meter;
  Constants.attr("cm") = &Constants::cm;
  Constants.attr("centimeter") = &Constants::centimeter;

  // second is a billion to be consistent with IceCube code
  Constants.attr("second") = &Constants::second;

  // speed of light
  Constants.attr("c") = &Constants::c;

  // masses, GeV/c^2
  Constants.attr("protonMass") = &Constants::protonMass;
  Constants.attr("neutronMass") = &Constants::neutronMass;
  Constants.attr("isoscalarMass") = &Constants::isoscalarMass;
  Constants.attr("electronMass") = &Constants::electronMass;
  Constants.attr("muonMass") = &Constants::muonMass;
  Constants.attr("tauMass") = &Constants::tauMass;
  Constants.attr("lambda0Mass") = &Constants::lambda0Mass;
  Constants.attr("Pi0Mass") = &Constants::Pi0Mass;
  Constants.attr("PiPlusMass") = &Constants::PiPlusMass;
  Constants.attr("PiMinusMass") = &Constants::PiMinusMass;
  Constants.attr("K0Mass") = &Constants::K0Mass;
  Constants.attr("KPlusMass") = &Constants::KPlusMass;
  Constants.attr("KMinusMass") = &Constants::KMinusMass;
  Constants.attr("KPrime0Mass") = &Constants::KPrime0Mass;
  Constants.attr("KPrimePlusMass") = &Constants::KPrimePlusMass;
  Constants.attr("KPrimeMinusMass") = &Constants::KPrimeMinusMass;
  Constants.attr("D0Mass") = &Constants::D0Mass;
  Constants.attr("DPlusMass") = &Constants::DPlusMass;
  Constants.attr("DMinusMass") = &Constants::DMinusMass;
  Constants.attr("DsPlusMass") = &Constants::DsPlusMass;
  Constants.attr("DsMinusMass") = &Constants::DsMinusMass;
  Constants.attr("EtaMass") = &Constants::EtaMass;
  Constants.attr("EtaPrimeMass") = &Constants::EtaPrimeMass;
  Constants.attr("Rho0Mass") = &Constants::Rho0Mass;
  Constants.attr("RhoPlusMass") = &Constants::RhoPlusMass;
  Constants.attr("RhoMinusMass") = &Constants::RhoMinusMass;
  Constants.attr("OmegaMass") = &Constants::OmegaMass;
  Constants.attr("PhiMass") = &Constants::PhiMass;

  // confusing units
  // static const double second          = 1.523e15; // [eV^-1 sec^-1]
  Constants.attr("s") = &Constants::s;
  Constants.attr("tauLifeTime") = &Constants::tauLifeTime;
  Constants.attr("MuonLifeTime") = &Constants::MuonLifeTime;

  // GeV/c^2
  Constants.attr("wMass") = &Constants::wMass;
  Constants.attr("wWidth") = &Constants::wWidth;
  Constants.attr("zMass") = &Constants::zMass;

  Constants.attr("WBranchE") = &Constants::WBranchE;
  Constants.attr("WBranchMuon") = &Constants::WBranchMuon;
  Constants.attr("WBranchTau") = &Constants::WBranchTau;
  Constants.attr("WBranchHadronic") = &Constants::WBranchHadronic;

  Constants.attr("nuEMass") = &Constants::nuEMass;
  Constants.attr("nuMuMass") = &Constants::nuMuMass;
  Constants.attr("nuTauMass") = &Constants::nuTauMass;
  /*
  W boson - http://pdg.lbl.gov/2019/listings/rpp2019-list-w-boson.pdf
  Z boson - http://pdg.lbl.gov/2018/listings/rpp2018-list-z-boson.pdf
  */

  // Unit Conversions

  Constants.attr("elementaryCharge") = &Constants::elementaryCharge;

  Constants.attr("GeV") = &Constants::GeV;
  Constants.attr("EeV") = &Constants::EeV;
  Constants.attr("PeV") = &Constants::PeV;
  Constants.attr("TeV") = &Constants::TeV;
  Constants.attr("MeV") = &Constants::MeV;
  Constants.attr("keV") = &Constants::keV;
  Constants.attr("eV") = &Constants::eV;
  Constants.attr("Joule") = &Constants::Joule;

  Constants.attr("GeV_per_amu") = &Constants::GeV_per_amu;
  Constants.attr("invGeVsq_per_cmsq") = &Constants::invGeVsq_per_cmsq;


  // may need to fix these after setting GeV to 1.0
  Constants.attr("FermiConstant") = &Constants::FermiConstant;
  Constants.attr("avogadro") = &Constants::avogadro;
  Constants.attr("thetaWeinberg") = &Constants::thetaWeinberg;
  Constants.attr("gravConstant") = &Constants::gravConstant;
  Constants.attr("fineStructure") = &Constants::fineStructure;
  Constants.attr("hbarc") = &Constants::hbarc;

  // CKM matrix elements
  // from https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf
  Constants.attr("Vud") = &Constants::Vud;
  Constants.attr("Vcd") = &Constants::Vcd;
  Constants.attr("Vtd") = &Constants::Vtd;
  Constants.attr("Vus") = &Constants::Vus;
  Constants.attr("Vcs") = &Constants::Vcs;
  Constants.attr("Vts") = &Constants::Vts;
  Constants.attr("Vub") = &Constants::Vub;
  Constants.attr("Vcb") = &Constants::Vcb;
  Constants.attr("Vtb") = &Constants::Vtb;
}
