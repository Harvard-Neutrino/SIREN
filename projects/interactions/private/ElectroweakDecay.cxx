#include "SIREN/interactions/ElectroweakDecay.h"

#include <cmath>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "SIREN/dataclasses/Particle.h"

#include "SIREN/math/Vector3D.h"

#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Constants.h"

#include "SIREN/detector/MaterialModel.h"

#include "SIREN/interactions/Decay.h"

namespace siren {
namespace interactions {

double ElectroweakDecay::ZDecayWidth(double& cL, double& cR) const {
  double cV = cL + cR;
  double cA = cL - cR;
  double gZ = siren::utilities::Constants::gweak/(sqrt(1 - siren::utilities::Constants::thetaWeinberg));
  return gZ*gZ * siren::utilities::Constants::zMass / (48 * siren::utilities::Constants::pi) * (cV*cV + cA*cA);
}

void ElectroweakDecay::SetCKMMap() {
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
}

void ElectroweakDecay::SetSecondaryMassMap() {
  secondary_masses[siren::dataclasses::ParticleType::NuE] = siren::utilities::Constants::nuEMass;
  secondary_masses[siren::dataclasses::ParticleType::NuMu] = siren::utilities::Constants::nuMuMass;
  secondary_masses[siren::dataclasses::ParticleType::NuTau] = siren::utilities::Constants::nuTauMass;
  secondary_masses[siren::dataclasses::ParticleType::NuEBar] = siren::utilities::Constants::nuEMass;
  secondary_masses[siren::dataclasses::ParticleType::NuMuBar] = siren::utilities::Constants::nuMuMass;
  secondary_masses[siren::dataclasses::ParticleType::NuTauBar] = siren::utilities::Constants::nuTauMass;
  secondary_masses[siren::dataclasses::ParticleType::EMinus] = siren::utilities::Constants::electronMass;
  secondary_masses[siren::dataclasses::ParticleType::MuMinus] = siren::utilities::Constants::muonMass;
  secondary_masses[siren::dataclasses::ParticleType::TauMinus] = siren::utilities::Constants::tauMass;
  secondary_masses[siren::dataclasses::ParticleType::EPlus] = siren::utilities::Constants::electronMass;
  secondary_masses[siren::dataclasses::ParticleType::MuPlus] = siren::utilities::Constants::muonMass;
  secondary_masses[siren::dataclasses::ParticleType::TauPlus] = siren::utilities::Constants::tauMass;
  secondary_masses[siren::dataclasses::ParticleType::u] = siren::utilities::Constants::upMass;
  secondary_masses[siren::dataclasses::ParticleType::c] = siren::utilities::Constants::charmMass;
  secondary_masses[siren::dataclasses::ParticleType::t] = siren::utilities::Constants::topMass;
  secondary_masses[siren::dataclasses::ParticleType::uBar] = siren::utilities::Constants::upMass;
  secondary_masses[siren::dataclasses::ParticleType::cBar] = siren::utilities::Constants::charmMass;
  secondary_masses[siren::dataclasses::ParticleType::tBar] = siren::utilities::Constants::topMass;
  secondary_masses[siren::dataclasses::ParticleType::d] = siren::utilities::Constants::downMass;
  secondary_masses[siren::dataclasses::ParticleType::s] = siren::utilities::Constants::strangeMass;
  secondary_masses[siren::dataclasses::ParticleType::b] = siren::utilities::Constants::bottomMass;
  secondary_masses[siren::dataclasses::ParticleType::dBar] = siren::utilities::Constants::downMass;
  secondary_masses[siren::dataclasses::ParticleType::sBar] = siren::utilities::Constants::strangeMass;
  secondary_masses[siren::dataclasses::ParticleType::bBar] = siren::utilities::Constants::bottomMass;
}

bool ElectroweakDecay::equal(Decay const & other) const {
    const ElectroweakDecay* x = dynamic_cast<const ElectroweakDecay*>(&other);

    if(!x)
        return false;
    else
        return
            primary_types == x->primary_types;
}

double ElectroweakDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

double ElectroweakDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
  std::vector<dataclasses::InteractionSignature> signatures = GetPossibleSignaturesFromParent(primary);
  double gamma_tot = 0;
  dataclasses::InteractionRecord record;
  for(auto signature : signatures) {
    record.signature = signature;
    gamma_tot += TotalDecayWidthForFinalState(record);
  }
  return gamma_tot;
}

double ElectroweakDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
    double m1 = secondary_masses.at(record.signature.secondary_types[0]);
    double m2 = secondary_masses.at(record.signature.secondary_types[1]);
    if (m1+m2>=record.primary_mass) return 0;

    if (record.signature.primary_type == siren::dataclasses::ParticleType::WMinus) {
        for(auto l : Leptons) {
          for(auto nubar : AntiNus) {
            if (record.signature.secondary_types[0]==l && record.signature.secondary_types[1]==nubar) {
              return GammaW;
            }
          }
        }
        int iup = 0;
        for(auto ubar : UpAntiQuarks) {
          for(auto d : DownQuarks) {
            if (record.signature.secondary_types[0]==ubar && record.signature.secondary_types[1]==d) {
              return V_CKM.at(std::make_pair(UpQuarks[iup],d)) * GammaW;
            }
          }
          ++iup;
        }
    }
    else if (record.signature.primary_type == siren::dataclasses::ParticleType::WPlus) {
      for(auto lbar : AntiLeptons) {
          for(auto nu : Nus) {
            if (record.signature.secondary_types[0]==lbar && record.signature.secondary_types[1]==nu) {
              return GammaW;
            }
          }
        }
        for(auto u : UpQuarks) {
          int idown = 0;
          for(auto dbar : DownAntiQuarks) {
            if (record.signature.secondary_types[0]==u && record.signature.secondary_types[1]==dbar) {
              return V_CKM.at(std::make_pair(u,DownQuarks[idown])) * GammaW;
            }
            ++idown;
          }
        }
    }
    else if (record.signature.primary_type == siren::dataclasses::ParticleType::Z0) {
      double cL, cR;
      for(auto nu : Nus) {
        if (record.signature.secondary_types[0]==nu) {
          cL = 0.5;
          cR = 0;
          return ZDecayWidth(cL,cR);
        }
      }
      for(auto l : Leptons) {
        if (record.signature.secondary_types[0]==l) {
          cL = -0.27;
          cR = 0.23;
          return ZDecayWidth(cL,cR);
        }
      }
      for(auto u : UpQuarks) {
        if (record.signature.secondary_types[0]==u) {
          cL = 0.35;
          cR = -0.15;
          return ZDecayWidth(cL,cR);
        }
      }
      for(auto d : DownQuarks) {
        if (record.signature.secondary_types[0]==d) {
          cL = -0.42;
          cR = 0.08;
          return ZDecayWidth(cL,cR);
        }
      }
    }
    return 0;
}

std::vector<std::string> ElectroweakDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"};
}

double ElectroweakDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    // For now assume isotropic
    // dGamma / dCosTheta = Gamma/2
    return 1./2.*TotalDecayWidthForFinalState(record);

}

void ElectroweakDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {

    double mX = secondary_masses.at(record.signature.secondary_types[0]);
    double mY = secondary_masses.at(record.signature.secondary_types[1]);

    // For now assume isotropic
    // this is not true for polarized W
    double CosTheta = random->Uniform(-1,1);
    double SinTheta = std::sin(std::acos(CosTheta));

    rk::P4 pBoson(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::Boost boost_to_lab = pBoson.labBoost();

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 pBoson_mom = pBoson.momentum();
    geom3::UnitVector3 pBoson_dir = pBoson_mom.direction();
    geom3::Rotation3 x_to_pBoson_rot = geom3::rotationBetween(x_dir, pBoson_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(pBoson_dir, phi);

    double M = record.primary_mass;
    double EX_Bosonrest = (M*M - mY*mY + mX*mX)/(2*M);
    double X_mom_Bosonrest = sqrt(EX_Bosonrest*EX_Bosonrest - mX*mX);
    //double EY = (M*M - mX*mX + mY*mY)/(2*M);

    rk::P4 pX_Bosonrest(X_mom_Bosonrest*geom3::Vector3(CosTheta,SinTheta,0),mX);
    pX_Bosonrest.rotate(x_to_pBoson_rot);
    pX_Bosonrest.rotate(rand_rot);

    rk::P4 pX = pX_Bosonrest.boost(boost_to_lab);
    rk::P4 pY = pBoson - pX;
    assert(abs(pY.m()-mY)<1e-6);

    siren::dataclasses::SecondaryParticleRecord & X = record.GetSecondaryParticleRecord(0);
    siren::dataclasses::SecondaryParticleRecord & Y = record.GetSecondaryParticleRecord(1);

    X.SetFourMomentum({pX.e(), pX.px(), pX.py(), pX.pz()});
    X.SetMass(pX.m());
    X.SetHelicity(std::copysign(1.0, record.primary_helicity)); // TODO: treat helicity correctly

    Y.SetFourMomentum({pY.e(), pY.px(), pY.py(), pY.pz()});
    Y.SetMass(pY.m());
    Y.SetHelicity(std::copysign(1.0, record.primary_helicity)); // TODO: treat helicity correctly
}

double ElectroweakDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    double dd = DifferentialDecayWidth(record);
    double td = TotalDecayWidthForFinalState(record);
    if (dd == 0) return 0.;
    else if (td == 0) return 0.;
    else return dd/td;
}


std::vector<dataclasses::InteractionSignature> ElectroweakDecay::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> ElectroweakDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const {

    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = siren::dataclasses::ParticleType::Decay;

    signature.secondary_types.resize(2);
    if(primary==siren::dataclasses::ParticleType::WPlus) {
      // W+ -> l+ nu_l
      for (uint i = 0; i < AntiLeptons.size(); ++i) {
        signature.secondary_types[0] = AntiLeptons[i];
        signature.secondary_types[1] = Nus[i];
        signatures.push_back(signature);
      }
      // W+ -> u dbar
      for (auto u : UpQuarks) {
        for (auto d : DownAntiQuarks) {
          signature.secondary_types[0] = u;
          signature.secondary_types[1] = d;
          signatures.push_back(signature);
        }
      }
    }
    else if(primary==siren::dataclasses::ParticleType::WMinus) {
      // W- -> l- nu_l_bar
      for (uint i = 0; i < Leptons.size(); ++i) {
        signature.secondary_types[0] = Leptons[i];
        signature.secondary_types[1] = AntiNus[i];
        signatures.push_back(signature);
      }
      // W- -> ubar d
      for (auto u : UpAntiQuarks) {
        for (auto d : DownQuarks) {
          signature.secondary_types[0] = u;
          signature.secondary_types[1] = d;
          signatures.push_back(signature);
        }
      }
    }
    else if(primary==siren::dataclasses::ParticleType::Z0) {
      // Z -> nu nubar
      for (uint i = 0; i < Nus.size(); ++i) {
        signature.secondary_types[0] = Nus[i];
        signature.secondary_types[1] = AntiNus[i];
        signatures.push_back(signature);
      }
      // Z -> l- l+
      for (uint i = 0; i < Nus.size(); ++i) {
        signature.secondary_types[0] = Leptons[i];
        signature.secondary_types[1] = AntiLeptons[i];
        signatures.push_back(signature);
      }
      // Z -> u ubar
      for (uint i = 0; i < UpQuarks.size(); ++i) {
        signature.secondary_types[0] = UpQuarks[i];
        signature.secondary_types[1] = UpAntiQuarks[i];
        signatures.push_back(signature);
      }
      // Z -> d dbar
      for (uint i = 0; i < DownQuarks.size(); ++i) {
        signature.secondary_types[0] = DownQuarks[i];
        signature.secondary_types[1] = DownAntiQuarks[i];
        signatures.push_back(signature);
      }
    }

    return signatures;
}






} // namespace interactions
} // namespace siren