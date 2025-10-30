#include "SIREN/interactions/HNLDecay.h"

#include <cmath>

#include <gsl/gsl_integration.h>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include <CRunDec3.1/CRunDec.h>

#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/ParticleMasses.h"

#include "SIREN/math/Vector3D.h"

#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Sampling.h"

#include "SIREN/detector/MaterialModel.h"

#include "SIREN/interactions/Decay.h"

///////////////////
// HELPER FUNCTIONS
///////////////////

double lambda(double a, double b, double c){
  return a*a + b*b + c*c - 2*a*b - 2*a*c - 2*b*c;
}

// Functions to get three body decay widths
// Follows 1905.00284v2

double integrand1(double s, void * params) {
  double x = ((double*)params)[0];
  double y = ((double*)params)[1];
  double z = ((double*)params)[2];
  return 12 * (1./s) * (s - x - y) * (1 + z - s) * sqrt(lambda(s,x,y)*lambda(1,s,z));
}

double integrand2(double s, void * params) {
  double x = ((double*)params)[0];
  double y = ((double*)params)[1];
  double z = ((double*)params)[2];
  return 24 * sqrt(y*z) * (1./s) * (1 + x - s) * sqrt(lambda(s,y,z)*lambda(1,s,x));
}

double I1(double x, double y, double z) {
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &integrand1;
  double params[3] = { x, y, z };
  F.params = params;
  double result,error;
  gsl_integration_qags (&F, pow(sqrt(x)+sqrt(y),2), pow(1-sqrt(z),2), 0, 1e-5, 1000,
                        workspace, &result, &error);
  assert(error/result < 1e-3);
  gsl_integration_workspace_free (workspace);
  return result;
}

double I2(double x, double y, double z) {
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &integrand2;
  double params[3] = { x, y, z };
  F.params = params;
  double result,error;
  gsl_integration_qags (&F, pow(sqrt(y)+sqrt(z),2), pow(1-sqrt(x),2), 0, 1e-5, 1000,
                        workspace, &result, &error);
  assert(error/result < 1e-3);
  gsl_integration_workspace_free (workspace);
  return result;
}

// Functions to get quark decay widths
// All follow arXiv:1805.08567

double L(double x, double thresh=3e-3){
  double num;
  if (x<thresh) {
    num = 2*pow(x,6) + 6*pow(x,8) + 18*pow(x,10);
  }
  else {
    num = 1 - 3 * x*x - (1 - x*x) * sqrt(1 - 4 * x*x);
  }
  double denom = x*x * (1 + sqrt(1 - 4 * x*x));
  return log(num/denom);
}

// Integrand for CC decays
double integrand(double x, void * params){
  double xu = ((double*)params)[0];
  double xd = ((double*)params)[1];
  double xl = ((double*)params)[2];
  //return (x - xl*xl - xd*xd) * (1 + xu*xu - x) * sqrt(lambda(x,xl*xl,xd*xd)*lambda(1,x,xu*xu));
  return 12 * 1./x * (x - xl*xl - xd*xd) * (1 + xu*xu - x) * sqrt(lambda(x,xl*xl,xd*xd)*lambda(1,x,xu*xu));
}

double I(double xu, double xd, double xl) {
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &integrand;
  double params[3] = { xu, xd, xl };
  F.params = params;
  double result,error;
  gsl_integration_qags (&F, pow(xd+xl,2), pow(1-xu,2), 0, 1e-5, 1000,
                        workspace, &result, &error);
  assert(error/result < 1e-3);
  gsl_integration_workspace_free (workspace);
  return result;
}

double GammaQuarksCC(double U, double m_N, double mu, double md, double ml, double Nw=1) {
  double xu = mu/m_N;
  double xd = md/m_N;
  double xl = ml/m_N;
  if(xu+xd+xl>1) return 0;
  return Nw * pow(siren::utilities::Constants::FermiConstant,2) * pow(m_N,5) / (192*pow(siren::utilities::Constants::pi,3)) * U*U * I1(xu*xu,xd*xd,xl*xl);
}

double GammaQuarksNC(double U, double m_N, double mq, double Nz = 1, std::string fs = "qup") {
  double C1f, C2f;
  if (fs=="qup") {
    C1f = 1./4. * (1 - 8./3. * siren::utilities::Constants::thetaWeinberg + 32./9. * pow(siren::utilities::Constants::thetaWeinberg,2));
    C2f = 1./3. * siren::utilities::Constants::thetaWeinberg * (4./3. * siren::utilities::Constants::thetaWeinberg - 1);
  }
  else if (fs=="qdown") {
    C1f = 1./4. * (1 - 4./3. * siren::utilities::Constants::thetaWeinberg + 8./9. * pow(siren::utilities::Constants::thetaWeinberg,2));
    C2f = 1./6. * siren::utilities::Constants::thetaWeinberg * (2./3. * siren::utilities::Constants::thetaWeinberg - 1);
  }
  else return 0;

  double x = mq/m_N;
  if (2*x>=1) return 0;
  double prefactor = Nz * pow(siren::utilities::Constants::FermiConstant,2) * pow(m_N,5) / (192*pow(siren::utilities::Constants::pi,3)) * U*U;
  double factor1 = (1 - 14 * x*x - 2 * pow(x,4) - 12 * pow(x,6))*sqrt(1 - 4 * x*x) + 12 * pow(x,4) * (pow(x,4) - 1) * L(x);
  double factor2 = x*x * (2 + 10 * pow(x,2) - 12 * pow(x,4)) * sqrt(1 - 4 * x*x) + 6 * pow(x,4) * (1 - 2 * x*x + 2 * pow(x,4)) * L(x);
  return prefactor * (C1f*factor1 + 4*C2f*factor2);
}

double GammaHadronsCC(double U, double m_N, double m_l) {
  std::map<siren::dataclasses::Particle::ParticleType,double> quark_masses;
  quark_masses[siren::dataclasses::Particle::ParticleType::u] = siren::utilities::Constants::upMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::d] = siren::utilities::Constants::downMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::c] = siren::utilities::Constants::charmMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::s] = siren::utilities::Constants::strangeMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::t] = siren::utilities::Constants::topMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::b] = siren::utilities::Constants::bottomMass;
  std::map<std::pair<siren::dataclasses::Particle::ParticleType,siren::dataclasses::Particle::ParticleType>,double> V_CKM;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::u,
                  siren::dataclasses::Particle::ParticleType::d)] = siren::utilities::Constants::Vud;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::u,
                  siren::dataclasses::Particle::ParticleType::s)] = siren::utilities::Constants::Vus;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::u,
                  siren::dataclasses::Particle::ParticleType::b)] = siren::utilities::Constants::Vub;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::c,
                  siren::dataclasses::Particle::ParticleType::d)] = siren::utilities::Constants::Vcd;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::c,
                  siren::dataclasses::Particle::ParticleType::s)] = siren::utilities::Constants::Vcs;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::c,
                  siren::dataclasses::Particle::ParticleType::b)] = siren::utilities::Constants::Vcb;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::t,
                  siren::dataclasses::Particle::ParticleType::d)] = siren::utilities::Constants::Vtd;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::t,
                  siren::dataclasses::Particle::ParticleType::s)] = siren::utilities::Constants::Vts;
  V_CKM[std::make_pair(siren::dataclasses::Particle::ParticleType::t,
                  siren::dataclasses::Particle::ParticleType::b)] = siren::utilities::Constants::Vtb;
  std::vector<siren::dataclasses::Particle::ParticleType> upquarks = {siren::dataclasses::Particle::ParticleType::u,
                                                                      siren::dataclasses::Particle::ParticleType::c,
                                                                      siren::dataclasses::Particle::ParticleType::t};
  std::vector<siren::dataclasses::Particle::ParticleType> dnquarks = {siren::dataclasses::Particle::ParticleType::d,
                                                                      siren::dataclasses::Particle::ParticleType::s,
                                                                      siren::dataclasses::Particle::ParticleType::b};
  double Gamma_qq = 0;
  for(auto upquark : upquarks) {
    for(auto dnquark : dnquarks) {
      Gamma_qq += GammaQuarksCC(U, m_N, quark_masses[upquark], quark_masses[dnquark], m_l, 3*V_CKM[std::make_pair(upquark,dnquark)]);
    }
  }
  return Gamma_qq;
}

double GammaHadronsNC(double U, double m_N) {
  std::map<siren::dataclasses::Particle::ParticleType,double> quark_masses;
  quark_masses[siren::dataclasses::Particle::ParticleType::u] = siren::utilities::Constants::upMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::d] = siren::utilities::Constants::downMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::c] = siren::utilities::Constants::charmMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::s] = siren::utilities::Constants::strangeMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::t] = siren::utilities::Constants::topMass;
  quark_masses[siren::dataclasses::Particle::ParticleType::b] = siren::utilities::Constants::bottomMass;
  std::vector<siren::dataclasses::Particle::ParticleType> upquarks = {siren::dataclasses::Particle::ParticleType::u,
                                                                      siren::dataclasses::Particle::ParticleType::c,
                                                                      siren::dataclasses::Particle::ParticleType::t};
  std::vector<siren::dataclasses::Particle::ParticleType> dnquarks = {siren::dataclasses::Particle::ParticleType::d,
                                                                      siren::dataclasses::Particle::ParticleType::s,
                                                                      siren::dataclasses::Particle::ParticleType::b};

  double Gamma_qq = 0;
  for (const auto& upquark : upquarks) {
    Gamma_qq += GammaQuarksNC(U, m_N, quark_masses[upquark], 3, "qup");
  }
  for (const auto& dnquark : dnquarks) {
    double x = 1.;
    if (dnquark == siren::dataclasses::Particle::ParticleType::s) {
      x = 1 - 4 * pow(siren::utilities::Constants::KPlusMass/m_N,2);
      if (x <=0) continue;
    }
    Gamma_qq += sqrt(x) * GammaQuarksNC(U, m_N, quark_masses[dnquark], 3, "qdown");
  }
  return Gamma_qq;
}

namespace siren {
namespace interactions {

bool HNLDecay::equal(Decay const & other) const {
    const HNLDecay* x = dynamic_cast<const HNLDecay*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(
                    primary_types,
                    hnl_mass,
                    nature,
                    mixing)
            ==
            std::tie(
                    x->primary_types,
                    x->hnl_mass,
                    x->nature,
                    x->mixing);
}

double HNLDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

double HNLDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
  std::vector<dataclasses::InteractionSignature> signatures = GetPossibleSignaturesFromParent(primary);
  double gamma_tot = 0;
  dataclasses::InteractionRecord record;
  for(auto signature : signatures) {
    record.signature = signature;
    gamma_tot += TotalDecayWidthForFinalState(record);
  }
  return gamma_tot;
}

double HNLDecay::CCMesonDecayWidth(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionRecord proxyRecord = record;
  // Check lepton charge
  // negative case
  double GammaMeson = 0;
  if (proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::EMinus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::MuMinus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::TauMinus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::EPlus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::MuPlus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::TauPlus) {
    for(auto meson : PlusChargedMesons) {
      proxyRecord.signature.secondary_types[1] = meson;
      GammaMeson += TotalDecayWidthForFinalState(proxyRecord);
    }
  }
  else if (proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::EPlus ||
           proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::MuPlus ||
           proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::TauPlus) {
    for(auto meson : MinusChargedMesons) {
      proxyRecord.signature.secondary_types[1] = meson;
      GammaMeson += TotalDecayWidthForFinalState(proxyRecord);
    }
  }
  else {
    std::cerr << "Asked for CC Meson Decay width but first secondary particle is " << proxyRecord.signature.secondary_types[0] << std::endl;
  }
  return GammaMeson;
}

double HNLDecay::NCMesonDecayWidth(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionRecord proxyRecord = record;
  assert(proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuE ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuMu ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuTau ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuEBar ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuMuBar ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuTauBar);
  double GammaMeson = 0;
  for(auto meson : NeutralMesons) {
    proxyRecord.signature.secondary_types[1] = meson;
    GammaMeson += TotalDecayWidthForFinalState(proxyRecord);
  }
  return GammaMeson;
}

double HNLDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {

  // All decay widths from 2007.03701
  // and 0901.3589

  double mixing_element, width, m_alpha, m_beta;
  bool charged;

  // Let's start with 2 body decays
  if(record.signature.secondary_types.size() == 2) {

    if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuE ||
      record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuEBar) {
      mixing_element = mixing[0];
      m_alpha = 0;
      charged = false;
    }
    else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EMinus ||
            record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EPlus) {
      mixing_element = mixing[0];
      m_alpha = siren::utilities::Constants::electronMass;
      charged = true;
    }
    else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMu ||
            record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMuBar) {
      mixing_element = mixing[1];
      m_alpha = 0;
      charged = false;
    }
    else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuMinus ||
            record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuPlus) {
      mixing_element = mixing[1];
      m_alpha = siren::utilities::Constants::muonMass;
      charged = true;
    }
    else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTau ||
            record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTauBar) {
      mixing_element = mixing[2];
      m_alpha = 0;
      charged = false;
    }
    else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauMinus ||
            record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauPlus) {
      mixing_element = mixing[2];
      m_alpha = siren::utilities::Constants::tauMass;
      charged = true;
    }



    double f, m_meson, Vqq, gV;
    bool pseudoscalar;
    bool meson = true;

    // Neutral Pseudoscalar mesons: N -> P nu (P = pi0, eta, etaprime)
    if(record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Pi0) {
      assert(!charged);
      f = 0.130; // GeV
      m_meson = siren::utilities::Constants::Pi0Mass;
      pseudoscalar = true;
    }
    else if(record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Eta) {
      assert(!charged);
      f = 0.0816; // GeV
      m_meson = siren::utilities::Constants::EtaMass;
      pseudoscalar = true;
    }
    else if(record.signature.secondary_types[1]==siren::dataclasses::ParticleType::EtaPrime) {
      assert(!charged);
      f = -0.0946; // GeV
      m_meson = siren::utilities::Constants::EtaPrimeMass;
      pseudoscalar = true;
    }
    // Charged Pseudoscalar mesons: N -> P+- l-+ (P = pi, K, D, Ds)
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::PiPlus ||
             record.signature.secondary_types[1]==siren::dataclasses::ParticleType::PiMinus) {
      assert(charged);
      f = 0.130; // GeV
      m_meson = siren::utilities::Constants::PiPlusMass;
      pseudoscalar = true;
      Vqq = siren::utilities::Constants::Vud;
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::KPlus ||
             record.signature.secondary_types[1]==siren::dataclasses::ParticleType::KMinus) {
      assert(charged);
      f = 0.156; // GeV
      m_meson = siren::utilities::Constants::KPlusMass;
      pseudoscalar = true;
      Vqq = siren::utilities::Constants::Vus;
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::DPlus ||
             record.signature.secondary_types[1]==siren::dataclasses::ParticleType::DMinus) {
      assert(charged);
      f = 0.212; // GeV
      m_meson = siren::utilities::Constants::DPlusMass;
      pseudoscalar = true;
      Vqq = siren::utilities::Constants::Vcd;
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::DsPlus ||
             record.signature.secondary_types[1]==siren::dataclasses::ParticleType::DsMinus) {
      assert(charged);
      f = 0.249; // GeV
      m_meson = siren::utilities::Constants::DsPlusMass;
      pseudoscalar = true;
      Vqq = siren::utilities::Constants::Vcs;
    }
    // Neutral Vector mesons: N -> V nu (V = rho, omega, phi, K*)
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Rho0) {
      assert(!charged);
      f = 0.171; // GeV^2
      m_meson = siren::utilities::Constants::Rho0Mass;
      pseudoscalar = false;
      gV = 1 - 2 * siren::utilities::Constants::thetaWeinberg;
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Omega) {
      assert(!charged);
      f = 0.155; // GeV^2
      m_meson = siren::utilities::Constants::OmegaMass;
      pseudoscalar = false;
      gV = -2./3. * siren::utilities::Constants::thetaWeinberg;
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Phi) {
      assert(!charged);
      f = 0.232; // GeV^2
      m_meson = siren::utilities::Constants::PhiMass;
      pseudoscalar = false;
      gV = -sqrt(2)*(1./2. - 2./3.*siren::utilities::Constants::thetaWeinberg);
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::KPrime0) {
      assert(!charged);
      f = 0.178; // GeV^2
      m_meson = siren::utilities::Constants::KPrime0Mass;
      pseudoscalar = false;
      gV = 0; // Not provided in 2007.03701
    }
    // Charged Vector mesons: N -> V+- l-+ (V = rho, K*)
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::RhoPlus ||
             record.signature.secondary_types[1]==siren::dataclasses::ParticleType::RhoMinus) {
      assert(charged);
      f = 0.171; // GeV^2
      m_meson = siren::utilities::Constants::RhoPlusMass;
      pseudoscalar = false;
      Vqq = siren::utilities::Constants::Vud;
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::ParticleType::KPrimePlus ||
             record.signature.secondary_types[1]==siren::dataclasses::ParticleType::KPrimeMinus) {
      assert(charged);
      f = 0.178; // GeV^2
      m_meson = siren::utilities::Constants::KPrimePlusMass;
      pseudoscalar = false;
      Vqq = siren::utilities::Constants::Vus;
    }
    // Hadrons: N -> nu(l) + Hadrons
    // Implementation follows 1805.08567
    else if (charged && record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Hadrons) {
      meson = false;
      if (hnl_mass <= 1) return 0; // threshold for hadron decay
      // if(_GammaHadronsCC<=0) {
      //   SetGammaHadrons(GammaHadronsCC(mixing_element, hnl_mass, m_alpha),"CC");
      // }
      // width = _GammaHadronsCC;
      double cc_meson = CCMesonDecayWidth(record);
      if (nature==ChiralNature::Majorana) cc_meson *= 0.5; // make sure this is the dirac decay width for now
      width = GammaHadronsCC(mixing_element, hnl_mass, m_alpha);// - CCMesonDecayWidth(record);
      width -= cc_meson;
      // TODO: this correction should come in different a la the NC case
      // However, this treatment matches the results from 2007.03701
      width *= (1 + DeltaQCD()); // loop correction
    }
    else if (!charged && record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Hadrons) {
      meson = false;
      if (hnl_mass <= 1) return 0; // threshold for hadron decay
      // if(_GammaHadronsNC<=0) {
      //   SetGammaHadrons(GammaHadronsNC(mixing_element, hnl_mass),"NC");
      // }
      // width = _GammaHadronsNC;
      double nc_meson = NCMesonDecayWidth(record);
      if (nature==ChiralNature::Majorana) nc_meson *= 0.5; // make sure this is the dirac decay width for now
      width = (1 + DeltaQCD())*GammaHadronsNC(mixing_element, hnl_mass);// - NCMesonDecayWidth(record);
      width -= nc_meson;
      //width *= (1 + DeltaQCD()); // loop correction
    }
    // Weak Bosons
    // https://arxiv.org/abs/0901.3589v2
    else if (record.signature.secondary_types[1]==siren::dataclasses::Particle::ParticleType::Z0) {
      assert(!charged);
      meson = false;
      double xz = siren::utilities::Constants::zMass/hnl_mass;
      if (xz>=1) return 0;
      double muz = xz*xz;
      double Gamma_longitudinal = pow(siren::utilities::Constants::gweak,2) / (64 * siren::utilities::Constants::pi * pow(siren::utilities::Constants::wMass,2)) * pow(hnl_mass,3) * pow(1-muz,2);
      double cosw = cos(asin(sqrt(siren::utilities::Constants::thetaWeinberg)));
      double Gamma_transverse = pow(siren::utilities::Constants::gweak,2) / (32 * siren::utilities::Constants::pi * cosw*cosw) * hnl_mass * pow(1-muz,2);
      width = 0.5 * pow(mixing_element,2) * (Gamma_longitudinal + Gamma_transverse);
    }
    else if (record.signature.secondary_types[1]==siren::dataclasses::Particle::ParticleType::WPlus ||
             record.signature.secondary_types[1]==siren::dataclasses::Particle::ParticleType::WMinus) {
      assert(charged);
      meson = false;
      double xw = siren::utilities::Constants::wMass/hnl_mass;
      if (xw>=1) return 0;
      double muw = xw*xw;
      double Gamma_longitudinal = pow(siren::utilities::Constants::gweak,2) / (64 * siren::utilities::Constants::pi * pow(siren::utilities::Constants::wMass,2)) * pow(hnl_mass,3) * pow(1-muw,2);
      double Gamma_transverse = pow(siren::utilities::Constants::gweak,2) / (32 * siren::utilities::Constants::pi) * hnl_mass * pow(1-muw,2);
      width = pow(mixing_element,2) * (Gamma_longitudinal + Gamma_transverse);
    }
    // Signature not recognized
    else {
      std::cout << "HNL decay signature not recongized! Exiting\n";
      exit(0);
    }

    if (meson) {
      double x_meson = m_meson / hnl_mass;
      double x_alpha = m_alpha / hnl_mass;
      if(x_meson + x_alpha >= 1) return 0;
      double constant;

      if(pseudoscalar && !charged) {
        constant = pow(f,2) * pow(siren::utilities::Constants::FermiConstant,2) / (32 * siren::utilities::Constants::pi);
        width = constant * pow(hnl_mass,3) * pow(mixing_element,2) * pow(1-x_meson*x_meson,2);
      }
      else if(pseudoscalar && charged) {
        constant = pow(f,2) * pow(siren::utilities::Constants::FermiConstant,2) / (16 * siren::utilities::Constants::pi);
        width = constant * pow(hnl_mass,3) * pow(mixing_element,2) * pow(Vqq,2) * sqrt(lambda(1,pow(x_meson,2),pow(x_alpha,2))) * (1 - pow(x_meson,2) - pow(x_alpha,2)*(2 + pow(x_meson,2) - pow(x_alpha,2)));
      }
      else if(!pseudoscalar && !charged) {
        constant = pow(f,2) * pow(gV,2) * pow(siren::utilities::Constants::FermiConstant,2) / (32 * siren::utilities::Constants::pi);
        width = constant * pow(hnl_mass,3) / pow(m_meson,2) * pow(mixing_element,2) * (1 + 2 * pow(x_meson,2)) * pow(1 - pow(x_meson,2),2);
      }
      else if(!pseudoscalar && charged) {
        constant = pow(f,2) * pow(siren::utilities::Constants::FermiConstant,2) / (16 * siren::utilities::Constants::pi);
        width = constant * pow(hnl_mass,3) / pow(m_meson,2) * pow(mixing_element,2) * pow(Vqq,2) * sqrt(lambda(1,pow(x_meson,2),pow(x_alpha,2))) * ((1 - pow(x_meson,2))*(1+2*pow(x_meson,2)) + pow(x_alpha,2)*(pow(x_meson,2) + pow(x_alpha,2) - 2));
      }
      else {
        std::cout << "Could not find total decay width" << std::endl;
        exit(0);
      }
    }
  }
  else if(record.signature.secondary_types.size() == 3) {
    // three neutrino final state
    // see e.q. 3.19 of https://arxiv.org/abs/1905.00284
    if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight &&
       record.signature.secondary_types[1] == siren::dataclasses::ParticleType::NuLightBar &&
       record.signature.secondary_types[2] == siren::dataclasses::ParticleType::NuLight) {
        charged = false;
        mixing_element = pow(mixing[0],2) + pow(mixing[1],2) + pow(mixing[2],2);
        width = mixing_element * pow(siren::utilities::Constants::FermiConstant,2) * pow(hnl_mass,5) / (192*pow(siren::utilities::Constants::pi,3));
    }
    else {
      // make sure the first entry is a neutrino
      assert(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight ||
             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLightBar);
      // Find charged lepton masses
      int alpha,beta;
      if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::EMinus)
        {m_alpha = siren::utilities::Constants::electronMass; alpha = 0;}
      else if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::MuMinus)
        {m_alpha = siren::utilities::Constants::muonMass; alpha = 1;}
      else if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::TauMinus)
        {m_alpha = siren::utilities::Constants::tauMass; alpha = 2;}
      else {std::cerr << "Invalid HNL 3-body signature\n"; exit(0);}
      if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::EPlus)
        {m_beta = siren::utilities::Constants::electronMass; beta = 0;}
      else if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::MuPlus)
        {m_beta = siren::utilities::Constants::muonMass; beta = 1;}
      else if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::TauPlus)
        {m_beta = siren::utilities::Constants::tauMass; beta = 2;}
      else {std::cerr << "Invalid HNL 3-body signature\n"; exit(0);}
      double x_alpha = m_alpha/hnl_mass;
      double x_beta = m_beta/hnl_mass;
      // N -> nu l- l+ (l same flavor)
      if(int(record.signature.secondary_types[1]) == -int(record.signature.secondary_types[2])) {
        if (2*x_alpha>=1) return 0;
        charged = false; // even though there is a CC contribution, the dirac is factor of 2 larger
        double prefactor = pow(siren::utilities::Constants::FermiConstant,2) * pow(hnl_mass,5) / (192*pow(siren::utilities::Constants::pi,3));
        width = 0;
        double gL = -1./2 + siren::utilities::Constants::thetaWeinberg;
        double gR = siren::utilities::Constants::thetaWeinberg;
        for(int gamma = 0; gamma < 3; ++gamma) {
          width += pow(mixing[gamma],2) * ((gL*gR + (gamma==alpha)*gR)*I2(0,x_alpha*x_alpha,x_alpha*x_alpha) +
                                           (gL*gL + gR*gR + (gamma==alpha)*(1+2*gL))*I1(0,x_alpha*x_alpha,x_alpha*x_alpha));
        }
        width *= prefactor;
      }
      // N -> nu l- l+ (l differnt flavor)
      else {
        charged = true;
        if ((x_alpha+x_beta)>=1) return 0;
        double prefactor = pow(siren::utilities::Constants::FermiConstant,2) * pow(hnl_mass,5) / (384*pow(siren::utilities::Constants::pi,3));
        if(nature==ChiralNature::Majorana) {
          return prefactor * ((pow(mixing[alpha],2) * I1(0,x_alpha*x_alpha,x_beta*x_beta)) + (pow(mixing[beta],2) * I1(0,x_beta*x_beta,x_alpha*x_alpha)));
        }
        else if(nature==ChiralNature::Dirac) {
          if(record.signature.primary_type==siren::dataclasses::ParticleType::N4) {
            assert(record.signature.secondary_types[0]==siren::dataclasses::ParticleType::NuLight);
            return pow(mixing[alpha],2) * prefactor * I1(0,x_alpha*x_alpha,x_beta*x_beta);
          }
          else if(record.signature.primary_type==siren::dataclasses::ParticleType::N4Bar) {
            assert(record.signature.secondary_types[0]==siren::dataclasses::ParticleType::NuLightBar);
            return pow(mixing[beta],2) * prefactor * I1(0,x_beta*x_beta,x_alpha*x_alpha);
          }
        }
      }
    }
  }
  else {
    std::cout << "4+ body HNL decays not supported by this class\n";
    exit(0);
  }

  if(nature==ChiralNature::Majorana) return std::max(0.,2*width);
  else if(nature==ChiralNature::Dirac) return std::max(0.,width);
  return 0;
}

std::vector<std::string> HNLDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"}; // This is only true for two-body decays, may need some restructuring...
}


std::vector<dataclasses::InteractionSignature> HNLDecay::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> HNLDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const {

    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = siren::dataclasses::ParticleType::Decay;

    // Two body decays (from 2007.03701)
    signature.secondary_types.resize(2);
    if(primary==siren::dataclasses::ParticleType::N4) {
      // N -> nu P (P = pi0, eta, eta prime, rho0, omega, kprime0, hadrons, phi, Z)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuTau}) {
        signature.secondary_types[0] = particle;
        for (auto secondary_particle : NeutralMesons) {
          signature.secondary_types[1] = secondary_particle;
          signatures.push_back(signature);
        }
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Z0;
        signatures.push_back(signature);
      }
      // N -> l- P+ (P = pi, K , rho, K*, hadrons, D, Ds, W)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EMinus, siren::dataclasses::ParticleType::MuMinus, siren::dataclasses::ParticleType::TauMinus}) {
        signature.secondary_types[0] = particle;
        for (auto secondary_particle : PlusChargedMesons) {
          signature.secondary_types[1] = secondary_particle;
          signatures.push_back(signature);
        }
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::WPlus;
        signatures.push_back(signature);
      }
    }
    else if(primary==siren::dataclasses::ParticleType::N4Bar) {
      // N -> nubar P (P = pi0, eta, eta prime, rho0, omega, hadrons, phi, Z)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTauBar}) {
        signature.secondary_types[0] = particle;
        for (auto secondary_particle : NeutralMesons) {
          signature.secondary_types[1] = secondary_particle;
          signatures.push_back(signature);
        }
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Z0;
        signatures.push_back(signature);
      }
      // N -> l+ P- (P = pi, K , rho, K*, hadrons, D, Ds, W)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EPlus, siren::dataclasses::ParticleType::MuPlus, siren::dataclasses::ParticleType::TauPlus}) {
        signature.secondary_types[0] = particle;
        for (auto secondary_particle : MinusChargedMesons) {
          signature.secondary_types[1] = secondary_particle;
          signatures.push_back(signature);
        }
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::WMinus;
        signatures.push_back(signature);
      }
    }

    // Three body decays (including those from 2007.03701 and 1905.00284)
    // For the charged lepton decays, we follow 1905.00284 and ignore the flavor of the outgoing neutrino
    signature.secondary_types.resize(3);
    // N -> nu nu nu
    signature.secondary_types[0] = siren::dataclasses::ParticleType::NuLight;
    signature.secondary_types[1] = siren::dataclasses::ParticleType::NuLightBar;
    signature.secondary_types[2] = siren::dataclasses::ParticleType::NuLight;
    signatures.push_back(signature);

    if(primary==siren::dataclasses::ParticleType::N4)
      signature.secondary_types[0] = siren::dataclasses::ParticleType::NuLight;
    else if(primary==siren::dataclasses::ParticleType::N4Bar)
      signature.secondary_types[0] = siren::dataclasses::ParticleType::NuLightBar;

    // N -> nu l- l+ (l same flavor)
    signature.secondary_types[1] = siren::dataclasses::ParticleType::EMinus;
    signature.secondary_types[2] = siren::dataclasses::ParticleType::EPlus;
    signatures.push_back(signature);
    signature.secondary_types[1] = siren::dataclasses::ParticleType::MuMinus;
    signature.secondary_types[2] = siren::dataclasses::ParticleType::MuPlus;
    signatures.push_back(signature);
    signature.secondary_types[1] = siren::dataclasses::ParticleType::TauMinus;
    signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
    signatures.push_back(signature);
    for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EMinus, siren::dataclasses::ParticleType::MuMinus, siren::dataclasses::ParticleType::TauMinus}) {
      signature.secondary_types[1] = particle;
      // N -> nu l- l+ (l different flavor)
      signature.secondary_types[2] = siren::dataclasses::ParticleType::EPlus;
      if(particle != siren::dataclasses::ParticleType::EMinus) signatures.push_back(signature);
      signature.secondary_types[2] = siren::dataclasses::ParticleType::MuPlus;
      if(particle != siren::dataclasses::ParticleType::MuMinus) signatures.push_back(signature);
      signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
      if(particle != siren::dataclasses::ParticleType::TauMinus) signatures.push_back(signature);
    }

    return signatures;
}

double HNLDecay::GetMass(dataclasses::ParticleType const & secondary) const {
  // pseudoscalar meson, photon, hadrons
  if (secondary == dataclasses::ParticleType::Pi0) {
    return siren::utilities::Constants::Pi0Mass;
  }
  else if (secondary == dataclasses::ParticleType::Eta) {
    return siren::utilities::Constants::EtaMass;
  }
  else if (secondary == dataclasses::ParticleType::EtaPrime) {
    return siren::utilities::Constants::EtaPrimeMass;
  }
  else if (secondary == dataclasses::ParticleType::Gamma) {
    return 0;
  }
  else if (secondary == dataclasses::ParticleType::Hadrons) {
    return 0;
  }
  // Vector particles
  else if (secondary == dataclasses::ParticleType::Z0) {
    return siren::utilities::Constants::zMass;
  }
  else if (secondary == dataclasses::ParticleType::Rho0) {
    return siren::utilities::Constants::Rho0Mass;
  }
  else if (secondary == dataclasses::ParticleType::Omega) {
    return siren::utilities::Constants::OmegaMass;
  }
  else if (secondary == dataclasses::ParticleType::Phi) {
    return siren::utilities::Constants::PhiMass;
  }
  else if (secondary == dataclasses::ParticleType::KPrime0) {
    return siren::utilities::Constants::KPrime0Mass;
  }

  return 0;
}

// Follows https://arxiv.org/abs/1805.07523v1
double HNLDecay::GetAlpha(dataclasses::ParticleType const & secondary) const {
  // pseudoscalar meson, photon, hadrons
  if (secondary == dataclasses::ParticleType::Pi0 ||
      secondary == dataclasses::ParticleType::Eta ||
      secondary == dataclasses::ParticleType::EtaPrime ||
      secondary == dataclasses::ParticleType::Gamma ||
      secondary == dataclasses::ParticleType::Hadrons) {
    return 1;
  }
  // Vector particles
  else if (secondary == dataclasses::ParticleType::Z0) {
    double m = siren::utilities::Constants::zMass;
    return (hnl_mass*hnl_mass - 2*m*m) / (hnl_mass*hnl_mass + 2*m*m);
  }
  else if (secondary == dataclasses::ParticleType::Rho0) {
    double m = siren::utilities::Constants::Rho0Mass;
    return (hnl_mass*hnl_mass - 2*m*m) / (hnl_mass*hnl_mass + 2*m*m);
  }
  else if (secondary == dataclasses::ParticleType::Omega) {
    double m = siren::utilities::Constants::OmegaMass;
    return (hnl_mass*hnl_mass - 2*m*m) / (hnl_mass*hnl_mass + 2*m*m);
  }
  else if (secondary == dataclasses::ParticleType::Phi) {
    double m = siren::utilities::Constants::PhiMass;
    return (hnl_mass*hnl_mass - 2*m*m) / (hnl_mass*hnl_mass + 2*m*m);
  }
  else if (secondary == dataclasses::ParticleType::KPrime0) {
    double m = siren::utilities::Constants::KPrime0Mass;
    return (hnl_mass*hnl_mass - 2*m*m) / (hnl_mass*hnl_mass + 2*m*m);
  }
  return 0;
}

// Two body decays follow https://arxiv.org/abs/1805.07523v1 for NC decays and assume isotropy otherwise
// Three body decays follow arXiv:1905.00284v2
// TODO: update two-body decays to arXiv:1905.00284v2
double HNLDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    double DecayWidth = TotalDecayWidthForFinalState(record);
    // Check for isotropic decay
    if (nature==ChiralNature::Majorana) {
      return DecayWidth/2.; // This factor of 2 is for the cosTheta phase space, not the majorana nature :-)
    }
    if(record.secondary_momenta.size() ==2)
    {
      // Check for isotropic decay
      if (nature==ChiralNature::Majorana) return DecayWidth/2.; // This factor of 2 is for the cosTheta phase space, not the majorana nature :-)
      siren::math::Vector3D hnl_dir = siren::math::Vector3D(record.primary_momentum[0],
                                                            record.primary_momentum[1],
                                                            record.primary_momentum[2]);
      hnl_dir.normalize();
      unsigned int lep_index = (record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuE ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMu ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTau ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuEBar ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMuBar ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTauBar ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EMinus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuMinus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauMinus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EPlus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuPlus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauPlus) ? 0 : 1;
      unsigned int X_index = 1 - lep_index;
      rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
      rk::P4 pX(geom3::Vector3(record.secondary_momenta[X_index][1], record.secondary_momenta[X_index][2], record.secondary_momenta[X_index][3]), record.secondary_masses[X_index]);
      rk::Boost boost_to_HNL_rest = pHNL.restBoost();
      rk::P4 pX_HNLrest = pX.boost(boost_to_HNL_rest);

      siren::math::Vector3D X_dir = siren::math::Vector3D(pX_HNLrest.px(),
                                                          pX_HNLrest.py(),
                                                          pX_HNLrest.pz());
      X_dir.normalize();
      double CosThetaGamma = X_dir*hnl_dir; // scalar product
      double alpha = GetAlpha(record.signature.secondary_types[X_index]);
      alpha = std::copysign(alpha,record.primary_helicity); // 1 for RH, -1 for LH
      alpha = (record.signature.primary_type == siren::dataclasses::ParticleType::N4) ? -1*alpha : alpha;
      return DecayWidth/2. * (1 + alpha*CosThetaGamma);
    } // end 2 body decays
    else if (record.secondary_momenta.size()==3)
    {
      // three neutrino final state
      if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight &&
         record.signature.secondary_types[1] == siren::dataclasses::ParticleType::NuLightBar &&
         record.signature.secondary_types[2] == siren::dataclasses::ParticleType::NuLight)
      {
        return DecayWidth; // who cares about the outgoing neutrinos
      }
      else {
        // N (k1) -> nu (k2) l- (k3) l+ (k4)
        // follow eq 3.13 of arXiv:1905.00284v2
        // make sure the first entry is a neutrino
        assert(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight ||
               record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLightBar);
        // Find charged lepton masses

        int alpha,beta;
        double m_alpha,m_beta;

        if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::EMinus)
          {m_alpha = siren::utilities::Constants::electronMass; alpha = 0;}
        else if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::MuMinus)
          {m_alpha = siren::utilities::Constants::muonMass; alpha = 1;}
        else if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::TauMinus)
          {m_alpha = siren::utilities::Constants::tauMass; alpha = 2;}
        else {std::cerr << "Invalid HNL 3-body signature\n"; exit(0);}
        if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::EPlus)
          {m_beta = siren::utilities::Constants::electronMass; beta = 0;}
        else if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::MuPlus)
          {m_beta = siren::utilities::Constants::muonMass; beta = 1;}
        else if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::TauPlus)
          {m_beta = siren::utilities::Constants::tauMass; beta = 2;}
        else {std::cerr << "Invalid HNL 3-body signature\n"; exit(0);}

        rk::P4 k1(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), hnl_mass);
        rk::P4 k2(geom3::Vector3(record.secondary_momenta[0][1], record.secondary_momenta[0][2], record.secondary_momenta[0][3]), 0);
        rk::P4 k3(geom3::Vector3(record.secondary_momenta[1][1], record.secondary_momenta[1][2], record.secondary_momenta[1][3]), m_alpha);
        rk::P4 k4(geom3::Vector3(record.secondary_momenta[2][1], record.secondary_momenta[2][2], record.secondary_momenta[2][3]), m_beta);

        // let's get the phase space variables
        double s1 = (k2+k3).dot(k2+k3) / pow(hnl_mass,2);
        double s2 = (k2+k4).dot(k2+k4) / pow(hnl_mass,2);
        double CosTheta3 = k1.momentum().direction().dot(k3.momentum().direction()); // Angle between HNL and k3 in lab frame
        double CosTheta4 = k1.momentum().direction().dot(k4.momentum().direction()); // Angle between HNL and k4 in lab frame

        return ThreeBodyDifferentialDecayWidth(record, alpha, beta, m_alpha, m_beta, s1, s2, CosTheta3, CosTheta4);
      }
    }
    else
    {
      std::cerr << "Too many final state particles in HNL decay\n";
      exit(0);
    }

}

void HNLDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {

    siren::dataclasses::InteractionSignature const & signature = record.GetSignature();

    if(signature.secondary_types.size() == 2) {

      unsigned int lep_index = (record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuE ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMu ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTau ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuEBar ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMuBar ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTauBar ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EMinus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuMinus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauMinus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EPlus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuPlus ||
                                record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauPlus) ? 0 : 1;
      unsigned int X_index = 1 - lep_index;

      double CosTheta;
      double alpha = GetAlpha(record.signature.secondary_types[X_index]);
      alpha = std::copysign(alpha,record.primary_helicity); // 1 for RH, -1 for LH
      alpha = (record.signature.primary_type == siren::dataclasses::ParticleType::N4) ? -1*alpha : alpha;
      if(nature==ChiralNature::Majorana) {
        CosTheta = random->Uniform(-1,1);
      }
      else {
        double C = random->Uniform(0,1);
        CosTheta = (std::sqrt(1  - 2*alpha*(1 - alpha/2. - 2*C)) - 1)/alpha;
      }
      double SinTheta = std::sin(std::acos(CosTheta));

      rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
      rk::Boost boost_to_lab = pHNL.labBoost();

      geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
      geom3::Vector3 pHNL_mom = pHNL.momentum();
      geom3::UnitVector3 pHNL_dir = pHNL_mom.direction();
      geom3::Rotation3 x_to_pHNL_rot = geom3::rotationBetween(x_dir, pHNL_dir);

      double phi = random->Uniform(0, 2.0 * M_PI);
      geom3::Rotation3 rand_rot(pHNL_dir, phi);

      double mX = siren::dataclasses::GetParticleMass(record.signature.secondary_types[X_index]);
      double mLep = siren::dataclasses::GetParticleMass(record.signature.secondary_types[lep_index]);
      double EX_HNLrest = (hnl_mass*hnl_mass + mX*mX - mLep*mLep) / (2.*hnl_mass);
      rk::P4 pX_HNLrest(sqrt(pow(EX_HNLrest,2)-pow(mX,2))*geom3::Vector3(CosTheta,SinTheta,0),mX);
      pX_HNLrest.rotate(x_to_pHNL_rot);
      pX_HNLrest.rotate(rand_rot);

      rk::P4 pX = pX_HNLrest.boost(boost_to_lab);
      rk::P4 pLep(pHNL.momentum() - pX.momentum(),mLep); // ensures the neutrino has zero mass, avoids rounding errors

      siren::dataclasses::SecondaryParticleRecord & X = record.GetSecondaryParticleRecord(X_index);
      siren::dataclasses::SecondaryParticleRecord & lep = record.GetSecondaryParticleRecord(lep_index);

      assert(lep.type == siren::dataclasses::ParticleType::NuE ||
            lep.type == siren::dataclasses::ParticleType::NuMu ||
            lep.type == siren::dataclasses::ParticleType::NuTau ||
            lep.type == siren::dataclasses::ParticleType::NuEBar ||
            lep.type == siren::dataclasses::ParticleType::NuMuBar ||
            lep.type == siren::dataclasses::ParticleType::NuTauBar ||
            lep.type == siren::dataclasses::ParticleType::EMinus ||
            lep.type == siren::dataclasses::ParticleType::MuMinus ||
            lep.type == siren::dataclasses::ParticleType::TauMinus ||
            lep.type == siren::dataclasses::ParticleType::EPlus ||
            lep.type == siren::dataclasses::ParticleType::MuPlus ||
            lep.type == siren::dataclasses::ParticleType::TauPlus);

      X.SetFourMomentum({pX.e(), pX.px(), pX.py(), pX.pz()});
      X.SetMass(pX.m());
      X.SetHelicity(std::copysign(1.0, record.primary_helicity));

      lep.SetFourMomentum({pLep.e(), pLep.px(), pLep.py(), pLep.pz()});
      lep.SetMass(pLep.m());
      lep.SetHelicity(-1*record.primary_helicity);
    }
    else if(signature.secondary_types.size() == 3) {
      // three body decays

      if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight &&
         record.signature.secondary_types[1] == siren::dataclasses::ParticleType::NuLightBar &&
         record.signature.secondary_types[2] == siren::dataclasses::ParticleType::NuLight)
      {
        // N -> nu nubar nu
        // for now, be stupid and let all neutrinos have 1/3 of the HNL energy and be colinear
        // this never matters so it's probably ok to be stupid forever
        rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
        siren::dataclasses::SecondaryParticleRecord & nu1 = record.GetSecondaryParticleRecord(0);
        siren::dataclasses::SecondaryParticleRecord & nu2 = record.GetSecondaryParticleRecord(1);
        siren::dataclasses::SecondaryParticleRecord & nu3 = record.GetSecondaryParticleRecord(2);
        nu1.SetThreeMomentum({1./3.*pHNL.px(), 1./3.*pHNL.py(), 1./3.*pHNL.pz()});
        nu1.SetMass(0);
        nu1.SetHelicity(std::copysign(1.0, record.primary_helicity));
        nu2.SetThreeMomentum({1./3.*pHNL.px(), 1./3.*pHNL.py(), 1./3.*pHNL.pz()});
        nu2.SetMass(0);
        nu2.SetHelicity(std::copysign(1.0, record.primary_helicity));
        nu3.SetThreeMomentum({1./3.*pHNL.px(), 1./3.*pHNL.py(), 1./3.*pHNL.pz()});
        nu3.SetMass(0);
        nu3.SetHelicity(std::copysign(1.0, record.primary_helicity));
      }
      else {
        // N (k1) -> nu (k2) l- (k3) l+ (k4)

        // make sure the first entry is a neutrino
        assert(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight ||
               record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLightBar);
        // Find charged lepton masses

        int alpha,beta;
        double m_alpha,m_beta;

        if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::EMinus)
          {m_alpha = siren::utilities::Constants::electronMass; alpha = 0;}
        else if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::MuMinus)
          {m_alpha = siren::utilities::Constants::muonMass; alpha = 1;}
        else if(record.signature.secondary_types[1] == siren::dataclasses::ParticleType::TauMinus)
          {m_alpha = siren::utilities::Constants::tauMass; alpha = 2;}
        else {std::cerr << "Invalid HNL 3-body signature\n"; exit(0);}
        if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::EPlus)
          {m_beta = siren::utilities::Constants::electronMass; beta = 0;}
        else if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::MuPlus)
          {m_beta = siren::utilities::Constants::muonMass; beta = 1;}
        else if(record.signature.secondary_types[2] == siren::dataclasses::ParticleType::TauPlus)
          {m_beta = siren::utilities::Constants::tauMass; beta = 2;}
        else {std::cerr << "Invalid HNL 3-body signature\n"; exit(0);}

        siren::dataclasses::SecondaryParticleRecord & nu = record.GetSecondaryParticleRecord(0);
        siren::dataclasses::SecondaryParticleRecord & lep_alpha = record.GetSecondaryParticleRecord(1);
        siren::dataclasses::SecondaryParticleRecord & lep_beta = record.GetSecondaryParticleRecord(2);


        // Make proposal function
        auto proposal_func = [&] () {
          return ThreeBodyPhaseSpaceProposalDistribution(m_alpha,m_beta,random);
        };

        // Make likelihood funciton
        auto likelihood_func = [&] (std::vector<double> input) {
          // assumes input = {s1,s2,CosTheta3_HNLRest,Phi3_HNLRest,Phi4_HNLRest}
          return ThreeBodyDifferentialDecayWidth(record.record,alpha,beta,m_alpha,m_beta,input[0],input[1],input[2],input[3],input[4]);
        };

        std::vector<double> sampled_params = siren::utilities::MetropolisHasting_Sample(proposal_func,likelihood_func,random);

        // Go from sampled rest frame angles to lab frame
        std::pair<rk::P4,rk::P4> k3k4 = ThreeBodyPhaseSpaceConversion(record.record,m_alpha,m_beta,sampled_params[0],sampled_params[1],sampled_params[2],sampled_params[3],sampled_params[4]);
        rk::P4 k3 = k3k4.first;
        rk::P4 k4 = k3k4.second;
        rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
        rk::P4 k2 = pHNL - k3 - k4;

        nu.SetFourMomentum({k2.e(), k2.px(), k2.py(), k2.pz()});
        nu.SetMass(0); // set manually
        nu.SetHelicity(record.primary_helicity);

        lep_alpha.SetFourMomentum({k3.e(), k3.px(), k3.py(), k3.pz()});
        lep_alpha.SetMass(k3.m());
        lep_alpha.SetHelicity(record.primary_helicity); // I don't think this matters...

        lep_beta.SetFourMomentum({k4.e(), k4.px(), k4.py(), k4.pz()});
        lep_beta.SetMass(k4.m());
        lep_beta.SetHelicity(-1*record.primary_helicity); // I don't think this matters...



      }

    }

}

double HNLDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}

double HNLDecay::DeltaQCD(int nloops) const {
  // See https://arxiv.org/pdf/1703.03751
  // and https://arxiv.org/abs/2007.03701
  std::vector<double> quark_masses = {siren::utilities::Constants::upMass,
                                      siren::utilities::Constants::downMass,
                                      siren::utilities::Constants::strangeMass,
                                      siren::utilities::Constants::charmMass,
                                      siren::utilities::Constants::bottomMass,
                                      siren::utilities::Constants::topMass};
  int nflavors = 0;
  for (auto q_mass : quark_masses) {
    if (hnl_mass > q_mass) ++nflavors;
  }
  double m_ref, alpha_s_ref;
  if(nflavors <=5) {
    // Use the tau mass for reference
    m_ref = siren::utilities::Constants::tauMass;
    alpha_s_ref = 0.332;
  }
  else {
    // Use the Z mass
    m_ref = siren::utilities::Constants::zMass;
    alpha_s_ref = 0.1179;
  }
  double alpha_s = crundec->AlphasExact(alpha_s_ref, m_ref, hnl_mass, nflavors, nloops);
  return (alpha_s / siren::utilities::Constants::pi +
          5.2 * pow(alpha_s/siren::utilities::Constants::pi,2) +
          26.4 * pow(alpha_s/siren::utilities::Constants::pi,3));
}

void HNLDecay::SetGammaHadrons(double Gamma, std::string mode) {
  if (mode=="CC") {
    _GammaHadronsCC = Gamma;
  }
  else if (mode=="NC") {
    _GammaHadronsNC = Gamma;
  }
  else {
    std::cerr << "Invalid Gamma Hadron mode " << mode << std::endl;
  }
}

//////////////////////////////
// Three body decay class functions
//////////////////////////////

// follows 3.13 of 1905.00284v2 and subsequent sections
// Takes as input lab frame angles
// Process: N (k1) -> nu (k2) l_alpha- (k3) l_beta+ (k4)
double HNLDecay::ThreeBodyDifferentialDecayWidth(dataclasses::InteractionRecord const & record,
                                                 int & alpha, int & beta,
                                                 double & m_alpha, double & m_beta,
                                                 double & s1, double & s2,
                                                 double & CosTheta3, double & CosTheta4) const {
  // follow appendix b of 1905.00284v2
  double gL = -1./2 + siren::utilities::Constants::thetaWeinberg;
  double gR = siren::utilities::Constants::thetaWeinberg;

  double C1nu = 0,C2nu = 0,C3nu = 0,C4nu = 0,C5nu = 0,C6nu = 0;
  double C1nubar = 0,C2nubar = 0,C3nubar = 0,C4nubar = 0,C5nubar = 0,C6nubar = 0;

  for(int gamma = 0; gamma < 3; ++gamma) {
    C1nu += pow(mixing[gamma],2)*((alpha==beta)*pow(gL,2) + (gamma==alpha)*(1 + (alpha==beta)*gL));
    C2nu += (alpha==beta)*pow(gR,2)*pow(mixing[gamma],2);
    C3nu += (alpha==beta)*gR*pow(mixing[gamma],2)*((gamma==beta)+gL);
    C1nubar += (alpha==beta)*pow(gR,2)*pow(mixing[gamma],2);
    C2nubar += pow(mixing[gamma],2)*((alpha==beta)*pow(gL,2) + (gamma==beta)*(1 + (alpha==beta)*gL));
    C3nubar += (alpha==beta)*gR*pow(mixing[gamma],2)*((gamma==alpha)+gL);
  }
  C4nu = C1nu;
  C5nu = C2nu;
  C6nu = C3nu;
  C4nubar = -C1nubar;
  C5nubar = -C2nubar;
  C6nubar = -C3nubar;

  double C1=0,C2=0,C3=0,C4=0,C5=0,C6=0;

  if(nature==ChiralNature::Dirac) {
    if(record.signature.primary_type==siren::dataclasses::ParticleType::N4) {
      C1 = C1nu;
      C2 = C2nu;
      C3 = C3nu;
      C4 = C4nu;
      C5 = C5nu;
      C6 = C6nu;
    }
    else if(record.signature.primary_type==siren::dataclasses::ParticleType::N4Bar) {
      C1 = C1nubar;
      C2 = C2nubar;
      C3 = C3nubar;
      C4 = C4nubar;
      C5 = C5nubar;
      C6 = C6nubar;
    }
  }
  else if (nature==ChiralNature::Majorana) {
    C1 = C1nu+C1nubar;
    C2 = C2nu+C2nubar;
    C3 = C3nu+C3nubar;
    C4 = C4nu+C4nubar;
    C5 = C5nu+C5nubar;
    C6 = C6nu+C6nubar;
  }
  double x3 = m_alpha/hnl_mass;
  double x4 = m_beta/hnl_mass;
  double A0sq = C1*(s2 - x3*x3) * (1 + x4*x4 - s2) + C2*(s1 - x4*x4)*(1 + x3*x3 - s1) + 2*C3*x3*x4*(s1 + s2 - x3*x3 - x4*x4);
  double A1sq = (C4*(s2-x3*x3) - 2*C6*x3*x4)*sqrt(lambda(1,s2,x4*x4))*CosTheta4 + (C5*(s1-x4*x4) - 2*C6*x3*x4)*sqrt(lambda(1,s1,x3*x3))*CosTheta3;
  double prefactor = pow(siren::utilities::Constants::FermiConstant,2)*pow(hnl_mass,5) / (128*pow(siren::utilities::Constants::pi,5));
  double GammaPlus = prefactor*(A0sq+A1sq);
  double GammaMinus = prefactor*(A0sq-A1sq);
  if (record.primary_helicity>0) return GammaPlus;
  else if (record.primary_helicity<0) return GammaMinus;
  return 0; // no such thing is scalar HNLs!
}

// Takes as input the sampled phase space variables in the rest frame
// Returns k3 and k4 in the lab frame
// Process: N (k1) -> nu (k2) l_alpha- (k3) l_beta+ (k4)
std::pair<rk::P4,rk::P4> HNLDecay::ThreeBodyPhaseSpaceConversion(dataclasses::InteractionRecord const & record, double & m_alpha, double & m_beta, double & s1, double & s2, double & CosTheta3_HNLRest, double& Phi3_HNLRest, double & Phi4_HNLRest) const {


  // Step 1: compute CosTheta4_HNLRest
  double E3_HNLRest = (pow(hnl_mass,2) * (1 - s2) + pow(m_alpha,2)) / (2*hnl_mass);
  double E4_HNLRest = (pow(hnl_mass,2) * (1 - s1) + pow(m_beta,2)) / (2*hnl_mass);
  double p3_HNLRest = sqrt(E3_HNLRest*E3_HNLRest - m_alpha*m_alpha);
  double p4_HNLRest = sqrt(E4_HNLRest*E4_HNLRest - m_beta*m_beta);
  double s3 = 1 + (pow(m_alpha,2) + pow(m_beta,2))/pow(hnl_mass,2) - s1 - s2;
  double CosAlpha34_HNLRest = (E3_HNLRest*E4_HNLRest - 0.5 * (s3*pow(hnl_mass,2) - pow(m_alpha,2) - pow(m_beta,2))) / (p3_HNLRest*p4_HNLRest); // angle between k3 and k4 in HNL rest frame
  // solve a quadratic equation to get x = CosTheta4_HNLRest
  // A Sqrt(1 - x^2) + Bx - C = 0
  double A = sqrt(1-pow(CosTheta3_HNLRest,2)) * cos(Phi4_HNLRest-Phi3_HNLRest);
  double B = CosTheta3_HNLRest;
  double C = CosAlpha34_HNLRest;
  // from mathematica
  double term1 = B*C/(A*A + B*B);
  double term2 = sqrt(pow(A,4) + A*A*B*B - A*A*C*C) / (A*A + B*B);
  double CosTheta4_HNLRest = abs(term1+term2) <= 1 ? term1+term2 : term1-term2; // check for cos bounds
  assert(abs(CosTheta4_HNLRest)<=1); // double check
  double SinTheta3_HNLRest = sqrt(1-pow(CosTheta3_HNLRest,2));
  double SinTheta4_HNLRest = sqrt(1-pow(CosTheta4_HNLRest,2));


  // Step 2: boost to lab frame and compute lab-frame angles
  assert(abs(record.primary_mass-hnl_mass)<1e-6);
  rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
  rk::Boost boost_to_lab = pHNL.labBoost();
  rk::P4 k3_HNLRest(p3_HNLRest*geom3::Vector3(SinTheta3_HNLRest*cos(Phi3_HNLRest),SinTheta3_HNLRest*sin(Phi3_HNLRest),CosTheta3_HNLRest),m_alpha);
  rk::P4 k4_HNLRest(p4_HNLRest*geom3::Vector3(SinTheta4_HNLRest*cos(Phi4_HNLRest),SinTheta4_HNLRest*sin(Phi4_HNLRest),CosTheta4_HNLRest),m_alpha);
  rk::P4 k3 = k3_HNLRest.boost(boost_to_lab);
  rk::P4 k4 = k4_HNLRest.boost(boost_to_lab);
  return std::make_pair(k3,k4);
}

// Takes as input HNL rest frame angles
// Process: N (k1) -> nu (k2) l_alpha- (k3) l_beta+ (k4)
double HNLDecay::ThreeBodyDifferentialDecayWidth(dataclasses::InteractionRecord const & record, int & alpha, int & beta, double & m_alpha, double & m_beta, double & s1, double & s2, double & CosTheta3_HNLRest, double& Phi3_HNLRest, double & Phi4_HNLRest) const {

  rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
  std::pair<rk::P4,rk::P4> k3k4 = ThreeBodyPhaseSpaceConversion(record,m_alpha,m_beta,s1,s2,CosTheta3_HNLRest,Phi3_HNLRest,Phi4_HNLRest);
  rk::P4 k3 = k3k4.first;
  rk::P4 k4 = k3k4.second;
  double CosTheta3 = pHNL.momentum().direction().dot(k3.momentum().direction()); // Angle between HNL and k3 in lab frame
  double CosTheta4 = pHNL.momentum().direction().dot(k4.momentum().direction()); // Angle between HNL and k4 in lab frame

  // Call function based on lab frame angles
  return ThreeBodyDifferentialDecayWidth(record,alpha,beta,m_alpha,m_beta,s1,s2,CosTheta3,CosTheta4);

}

// Samples three body phase space
// s1, s2, CosTheta3, phi3, phi4 where angles are in HNL rest frame!
// s1 = (k2 + k3)^2 / mN^2
// s2 = (k2 + k4)^2 / mN^2
std::vector<double> HNLDecay::ThreeBodyPhaseSpaceProposalDistribution(double & m_alpha, double & m_beta, std::shared_ptr<siren::utilities::SIREN_random> random) const {
  double s1_min = pow(m_alpha,2)/pow(hnl_mass,2);
  double s1_max = pow(hnl_mass-m_beta,2)/pow(hnl_mass,2);

  // sample s1 uniformly first
  double s1 = random->Uniform(s1_min,s1_max);
  double m23_2 = s1*pow(hnl_mass,2);

  // Now follow section 3.2 of https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf
  // to get s2 * mN^2 = m24_2 = (k2 + k4)^2
  // considering m2 = mnu = 0
  // Namely, use eq 18 to find a quadratic eq a s2^s + b s2 + c = 0

  double a = pow(hnl_mass, 6)*s1;

  double b = ( pow(hnl_mass*hnl_mass * s1, 2)
              + m_alpha*m_alpha * (hnl_mass*hnl_mass - m_beta*m_beta)
              - hnl_mass*hnl_mass * s1 * (hnl_mass*hnl_mass + m_alpha*m_alpha + m_beta*m_beta)
            ) * (hnl_mass*hnl_mass);

  double c = (hnl_mass*hnl_mass * s1 * m_beta*m_beta * (hnl_mass*hnl_mass - m_alpha*m_alpha))
            - (m_alpha*m_alpha * m_beta*m_beta * (hnl_mass*hnl_mass - m_alpha*m_alpha - m_beta*m_beta));

  double s2_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double s2_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  double s2 = random->Uniform(s2_min,s2_max);

  // Now sample lab frame angles
  double discriminant = -1;
  size_t n_tries = 0;
  double CosTheta3,Phi3,Phi4;
  double E3,E4,p3,p4,s3,CosAlpha34,A,B,C;
  while (discriminant<0) {
    n_tries++;
    CosTheta3 = random->Uniform(-1,1);
    Phi3 = random->Uniform(0,2*siren::utilities::Constants::pi);
    Phi4 = random->Uniform(0,2*siren::utilities::Constants::pi);
    E3 = (pow(hnl_mass,2) * (1 - s2) + pow(m_alpha,2)) / (2*hnl_mass);
    E4 = (pow(hnl_mass,2) * (1 - s1) + pow(m_beta,2)) / (2*hnl_mass);
    p3 = sqrt(E3*E3 - m_alpha*m_alpha);
    p4 = sqrt(E4*E4 - m_beta*m_beta);
    s3 = 1 + (pow(m_alpha,2) + pow(m_beta,2))/pow(hnl_mass,2) - s1 - s2;
    CosAlpha34 = (E3*E4 - 0.5 * (s3*pow(hnl_mass,2) - pow(m_alpha,2) - pow(m_beta,2))) / (p3*p4); // angle between k3 and k4 in HNL rest frame
    // solve a quadratic equation to get x = CosTheta4_HNLRest
    // A Sqrt(1 - x^2) + Bx - C = 0
    A = sqrt(1-pow(CosTheta3,2)) * cos(Phi4-Phi3);
    B = CosTheta3;
    C = CosAlpha34;
    discriminant = pow(A,4) + A*A*B*B - A*A*C*C;
  }
  // std::cout << "Took " << n_tries << " to get a positive disc\n";
  return {s1,s2,CosTheta3,Phi3,Phi4};
}


} // namespace interactions
} // namespace siren

