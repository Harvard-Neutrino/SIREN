#pragma once
#ifndef SIREN_Constants_H
#define SIREN_Constants_H

#include <math.h> // pow function

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

// This header defines a lot of constants - not all actually used in this code.
// This was, for the most part, just copied from previous code of mine

namespace siren {
namespace utilities {
namespace Constants {

// geometry
static const double pi         = 3.141592653589793238462643383279502884197;
static const double tau        = 2*pi; //let's not play favorites
static const double degrees         = pi/180.; // converts degrees into radians.
// This is used since the user will enter a number in degrees, while the c trigonometric functions
// expect angles presented in radians.
static const double &deg            = degrees; // adding this in too...
static const double radian          = 1.;

// meter
static const double m               = 1.;
static const double &meter          = m;
static const double cm              = 0.01*m;
static const double &centimeter     = cm;

// second is a billion to be consistent with IceCube code
static const double second          = 1.e9;

// speed of light
static const double c               = 299792458.*(m/second); // [m sec^-1]

// masses, GeV/c^2
static const double protonMass      = 0.938272;
static const double neutronMass     = 0.939565;
static const double isoscalarMass   = 0.5*( protonMass + neutronMass );
static const double electronMass    = 0.000511;
static const double muonMass        = 0.105658374;
static const double tauMass         = 1.77686;
static const double lambda0Mass     = 1.1156836; // GeV
static const double Pi0Mass         = 0.1349770; // GeV
static const double PiPlusMass      = 0.13957039; // GeV
static const double PiMinusMass     = 0.13957039; // GeV
static const double K0Mass          = 0.493614; // GeV
static const double KplusMass       = 0.493677; // GeV
static const double KminusMass      = 0.493677; // GeV
static const double EtaMass         = 0.547862; // GeV
static const double EtaPrimeMass    = 0.95778;  // GeV
static const double RhoMass         = 0.77526; // GeV
static const double OmegaMass       = 0.78266; // GeV
static const double D0Mass          = 1.86962; // GeV
static const double DPlusMass       = 1.86484; // GeV
static const double CharmMass       = 1.27; // GeV

// confusing units
// static const double second          = 1.523e15; // [eV^-1 sec^-1]
static const double &s              = second;
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

static const double elementaryCharge  = 1.602176634e-19; // C

static const double GeV             = 1.0;
static const double EeV             = (1.0e9)*GeV; // [eV/EeV]
static const double PeV             = (1.0e6)*GeV;
static const double TeV             = (1.0e3)*GeV;
static const double MeV             = (1.0e-3)*GeV;
static const double keV             = (1.0e-6)*GeV;
static const double  eV             = (1.0e-9)*GeV;
static const double Joule           = eV/(1.602176634e-19); // eV/J

static const double GeV_per_amu     = 0.9314941024171441;
static const double invGeVsq_per_cmsq = 2.568189461147501e27;


// may need to fix these after setting GeV to 1.0
static const double FermiConstant   = 1.16639e-23/pow(eV,2); // [GeV^-2]
static const double avogadro        = 6.0221415e+23; // [mol cm^-3]
static const double thetaWeinberg   = 0.2312; // dimensionless
static const double gravConstant    = 6.6700e-11; // [m^3 kg^-1 s^-2]
static const double fineStructure   = 1.0/137.0; // dimensionless
static const double hbarc           = 197.3*(1e-9)*(1e-7)*GeV*cm; // [GeV m]

//hbar
static const double hbar            = 6.58211957 * (1e-25); // GeV seconds


} // namespace Constants
} // namespace utilities
} // namespace siren

#endif // SIREN_Constants_H

