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

double ElectroweakDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const &) const {

}

void ElectroweakDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const {

}

double ElectroweakDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {

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