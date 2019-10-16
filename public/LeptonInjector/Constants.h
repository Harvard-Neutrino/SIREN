#ifndef LI_CONSTANTS
#define LI_CONSTANTS

#include <math.h> // pow function

namespace LeptonInjector{ namespace Constants{

// masses, GeV/c^2
static const double protonMass      = 0.938272;
static const double neutronMass     = 0.939565;
static const double isoscalarMass   = 0.5*( protonMass + neutronMass );
static const double electronMass    = 0.000511;
static const double muonMass        = 0.105658374;
static const double tauMass         = 1.77686;

// confusing units
static const double second          = 1.523e15; // [eV^-1 sec^-1]
static const double tauLifeTime     = second*2.906e-13;
static const double MuonLifeTime    = second*2.196e-6;

// GeV/c^2
static const double wMass           = 80.379;
static const double wWidth          = 2.085; //GeV?
static const double zMass           = 91.1876;

static const double WBranchE        = 0.1071;
static const double WBranchMuon     = 0.1063;
static const double WBranchTau      = 0.1138;
static const double WBranchHadronic = 0.6741;

static const double nuEMass         = 0.;
static const double nuMuMass        = 0.;
static const double nuTauMass       = 0.;
/*
W boson - http://pdg.lbl.gov/2019/listings/rpp2019-list-w-boson.pdf
Z boson - http://pdg.lbl.gov/2018/listings/rpp2018-list-z-boson.pdf
*/

// Unit Conversions 
static const double EeV             = 1.0e18; // [eV/EeV]
static const double PeV             = 1.0e15;
static const double TeV             = 1.0e12;
static const double GeV             = 1.0e9;
static const double MeV             = 1.0e6;
static const double keV             = 1.0e3;
static const double  eV             = 1.0;
static const double Joule           = 1.0/(1.60225e-19); // eV/J

// Taken from SQUIDS, but modified because instantiating an object for constants is silly

static const double FermiConstant   = 1.16639e-23/pow(GeV,2); // [GeV^-2] 
static const double avogadro        = 6.0221415e+23; // [mol cm^-3]
static const double thetaWeinberg   = 0.2312; // dimensionless 
static const double gravConstant    = 6.6700e-11; // [m^3 kg^-1 s^-2]
static const double fineStructure   = 1.0/137.0; // dimensionless

} // namespace Constants
} // namespace nucross
#endif
