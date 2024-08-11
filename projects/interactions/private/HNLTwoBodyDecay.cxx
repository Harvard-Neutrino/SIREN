#include "SIREN/interactions/HNLTwoBodyDecay.h"

#include <cmath>

#include <gsl/gsl_integration.h>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include <CRunDec3.1/CRunDec.h>

#include "SIREN/dataclasses/Particle.h"

#include "SIREN/math/Vector3D.h"

#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Constants.h"

#include "SIREN/detector/MaterialModel.h"

#include "SIREN/interactions/Decay.h"

// Fucntions to get quark decay widths
// All follow arXiv:1805.08567



double lambda(double a, double b, double c){
  return a*a + b*b + c*c - 2*a*b - 2*a*c - 2*b*c;
}

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
  return 1./x * (x - xl*xl - xd*xd) * (1 + xu*xu - x) * sqrt(lambda(x,xl*xl,xd*xd)*lambda(1,x,xu*xu));
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
  return 12 * result;
}

double GammaQuarksCC(double U, double m_N, double mu, double md, double ml, double Nw=1) {
  double xu = mu/m_N;
  double xd = md/m_N;
  double xl = ml/m_N;
  if(xu+xd+xl>1) return 0;
  return Nw * pow(siren::utilities::Constants::FermiConstant,2) * pow(m_N,5) / (192*pow(siren::utilities::Constants::pi,3)) * U*U * I(xu,xd,xl);
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

bool HNLTwoBodyDecay::equal(Decay const & other) const {
    const HNLTwoBodyDecay* x = dynamic_cast<const HNLTwoBodyDecay*>(&other);

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

double HNLTwoBodyDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

double HNLTwoBodyDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
  std::vector<dataclasses::InteractionSignature> signatures = GetPossibleSignaturesFromParent(primary);
  double gamma_tot = 0;
  dataclasses::InteractionRecord record;
  for(auto signature : signatures) {
    record.signature = signature;
    gamma_tot += TotalDecayWidthForFinalState(record);
  }
  return gamma_tot;
}

double HNLTwoBodyDecay::CCMesonDecayWidth(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionRecord proxyRecord = record;
  // Check lepton charge
  // negative case
  double GammaMeson = 0;
  if (proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::EMinus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::MuMinus ||
      proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::TauMinus) {
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

double HNLTwoBodyDecay::NCMesonDecayWidth(dataclasses::InteractionRecord const & record) const {
  dataclasses::InteractionRecord proxyRecord = record;
  assert(proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuE ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuMu ||
         proxyRecord.signature.secondary_types[0] == dataclasses::ParticleType::NuTau);
  double GammaMeson = 0;
  for(auto meson : NeutralMesons) {
    proxyRecord.signature.secondary_types[1] = meson;
    GammaMeson += TotalDecayWidthForFinalState(proxyRecord);
  }
  return GammaMeson;
}

double HNLTwoBodyDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {

  // All decay widths from 2007.03701
  // and 0901.3589

  double mixing_element, width, m_alpha;
  bool charged;

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

  // Let's start with 2 body decays
  if(record.signature.secondary_types.size() == 2) {

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
      width = GammaHadronsCC(mixing_element, hnl_mass, m_alpha) - CCMesonDecayWidth(record);
      width *= (1 + DeltaQCD()); // loop correction
    }
    else if (!charged && record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Hadrons) {
      meson = false;
      if (hnl_mass <= 1) return 0; // threshold for hadron decay
      // if(_GammaHadronsNC<=0) {
      //   SetGammaHadrons(GammaHadronsNC(mixing_element, hnl_mass),"NC");
      // }
      // width = _GammaHadronsNC;
      width = GammaHadronsNC(mixing_element, hnl_mass) - NCMesonDecayWidth(record);
      width *= (1 + DeltaQCD()); // loop correction
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
    if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight &&
       record.signature.secondary_types[1] == siren::dataclasses::ParticleType::NuLight &&
       record.signature.secondary_types[2] == siren::dataclasses::ParticleType::NuLight) {
        mixing_element = pow(mixing[0],2) + pow(mixing[1],2) + pow(mixing[2],2);
        width = mixing_element * pow(siren::utilities::Constants::FermiConstant,2) * pow(hnl_mass,5) / (192*pow(siren::utilities::Constants::pi,3));
    }
    else {
      std::cout << "Only 3 neutrino 3-body decay is supported\n";
      exit(0);
    }
  }
  else {
    std::cout << "4+ body HNL decays not supported by this class\n";
    exit(0);
  }

  if(nature==ChiralNature::Majorana) return 2*width;
  else if(nature==ChiralNature::Dirac) return width;
  return 0;
}

std::vector<std::string> HNLTwoBodyDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"};
}


std::vector<dataclasses::InteractionSignature> HNLTwoBodyDecay::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> HNLTwoBodyDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const {

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

    // Three body decays (only nu nu nu for now, TODO include those from 2007.03701)
    signature.secondary_types.resize(3);
    signature.secondary_types[0] = siren::dataclasses::ParticleType::NuLight;
    signature.secondary_types[1] = siren::dataclasses::ParticleType::NuLight;
    signature.secondary_types[2] = siren::dataclasses::ParticleType::NuLight;
    signatures.push_back(signature);
    // if(primary==siren::dataclasses::ParticleType::N4) {
    //   for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuTau}) {
    //     signature.secondary_types[0] = particle;
    //     // N -> nu nu nu
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::NuEBar;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::NuMuBar;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::NuTauBar;
    //     signatures.push_back(signature);
    //     // N -> nu l- l+ (l same flavor)
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::EMinus;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::EPlus;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::MuMinus;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::MuPlus;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::TauMinus;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
    //     signatures.push_back(signature);
    //   }
    //   for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EMinus, siren::dataclasses::ParticleType::MuMinus, siren::dataclasses::ParticleType::TauMinus}) {
    //     signature.secondary_types[0] = particle;
    //     // N -> nu l- l+ (l different flavor)
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::EPlus;
    //     if(particle != siren::dataclasses::ParticleType::EMinus) signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::MuPlus;
    //     if(particle != siren::dataclasses::ParticleType::MuMinus) signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
    //     if(particle != siren::dataclasses::ParticleType::TauMinus) signatures.push_back(signature);
    //   }
    // }
    // else if(primary==siren::dataclasses::ParticleType::N4Bar) {
    //   for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTauBar}) {
    //     signature.secondary_types[0] = particle;
    //     // N -> nu nu nu
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::NuEBar;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::NuMuBar;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::NuTauBar;
    //     signatures.push_back(signature);
    //     // N -> nu l- l+ (l same flavor)
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::EMinus;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::EPlus;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::MuMinus;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::MuPlus;
    //     signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::TauMinus;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
    //     signatures.push_back(signature);
    //   }
    //   for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EPlus, siren::dataclasses::ParticleType::MuPlus, siren::dataclasses::ParticleType::TauPlus}) {
    //     signature.secondary_types[0] = particle;
    //     // N -> nu l- l+ (l different flavor)
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::EMinus;
    //     if(particle != siren::dataclasses::ParticleType::EMinus) signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::MuMinus;
    //     if(particle != siren::dataclasses::ParticleType::MuMinus) signatures.push_back(signature);
    //     signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTauBar;
    //     signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
    //     if(particle != siren::dataclasses::ParticleType::TauMinus) signatures.push_back(signature);
    //   }

    return signatures;
}

double HNLTwoBodyDecay::GetMass(dataclasses::ParticleType const & secondary) const {
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
double HNLTwoBodyDecay::GetAlpha(dataclasses::ParticleType const & secondary) const {
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

// TODO: this function should follow from arXiv:1905.00284v2
// For now, follow https://arxiv.org/abs/1805.07523v1 for NC decays and assume isotropy otherwise
double HNLTwoBodyDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    double DecayWidth = TotalDecayWidthForFinalState(record);
    // Check for isotropic decay
    if (nature==ChiralNature::Majorana || record.secondary_momenta.size() >=2) {
      return DecayWidth/2.; // This factor of 2 is for the cosTheta phase space, not the majorana nature :-)
    }
    assert(record.secondary_momenta.size() ==2);
    siren::math::Vector3D hnl_dir = siren::math::Vector3D(record.primary_momentum[0],
                                                          record.primary_momentum[1],
                                                          record.primary_momentum[2]);
    hnl_dir.normalize();
    unsigned int nu_index = (record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuE ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMu ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTau ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuEBar ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMuBar ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTauBar ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight) ? 0 : 1;
    unsigned int X_index = 1 - nu_index;
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
}

void HNLTwoBodyDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {

    siren::dataclasses::InteractionSignature const & signature = record.GetSignature();

     unsigned int nu_index = (record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuE ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMu ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTau ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuEBar ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMuBar ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTauBar ||
                             record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuLight) ? 0 : 1;
    unsigned int X_index = 1 - nu_index;

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

    rk::P4 pX_HNLrest(hnl_mass/2.*geom3::Vector3(CosTheta,SinTheta,0),0);
    pX_HNLrest.rotate(x_to_pHNL_rot);
    pX_HNLrest.rotate(rand_rot);

    rk::P4 pX = pX_HNLrest.boost(boost_to_lab);
    rk::P4 pNu(pHNL.momentum() - pX.momentum(),0); // ensures the neutrino has zero mass, avoids rounding errors

    siren::dataclasses::SecondaryParticleRecord & X = record.GetSecondaryParticleRecord(X_index);
    siren::dataclasses::SecondaryParticleRecord & nu = record.GetSecondaryParticleRecord(nu_index);

    assert(nu.type == siren::dataclasses::ParticleType::NuE ||
           nu.type == siren::dataclasses::ParticleType::NuMu ||
           nu.type == siren::dataclasses::ParticleType::NuTau ||
           nu.type == siren::dataclasses::ParticleType::NuEBar ||
           nu.type == siren::dataclasses::ParticleType::NuMuBar ||
           nu.type == siren::dataclasses::ParticleType::NuTauBar);

    X.SetFourMomentum({pX.e(), pX.px(), pX.py(), pX.pz()});
    X.SetMass(pX.m());
    X.SetHelicity(std::copysign(1.0, record.primary_helicity));

    nu.SetFourMomentum({pNu.e(), pNu.px(), pNu.py(), pNu.pz()});
    nu.SetMass(pNu.m());
    nu.SetHelicity(-1*record.primary_helicity);

}

double HNLTwoBodyDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}

double HNLTwoBodyDecay::DeltaQCD(int nloops) const {
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

void HNLTwoBodyDecay::SetGammaHadrons(double Gamma, std::string mode) {
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


} // namespace interactions
} // namespace siren

