#include "SIREN/interactions/HNLDecay.h"

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
                    dipole_coupling)
            ==
            std::tie(
                    x->primary_types,
                    x->hnl_mass,
                    x->nature,
                    x->dipole_coupling);
}

double HNLDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

double HNLDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
  std::vector<dataclasses::InteractionSignature> signatures = GetPossibleSignaturesFromParent(primary); 
  double gamma_tot = 0;
  dataclasses::InteractionRecord record;
  for(auto signature : signatures) {
    record.siganture = signature;
    gamma_tot += TotalDecayWidthForFinalState(record);
  }
}

double HNLDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
  
  // All decay widths from 2007.03701
  
  double mixing_element, width;

  if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuE ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuEBar ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EMinus ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::EPlus) {
    mixing_element = mixing[0]
  }
  else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMu ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuMuBar ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuMinus ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::MuPlus) {
    mixing_element = mixing[1]
  }
  else if(record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTau ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::NuTauBar ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauMinus ||
     record.signature.secondary_types[0] == siren::dataclasses::ParticleType::TauPlus) {
    mixing_element = mixing[2]
  }
  // Let's start with 2 body decays 
  if(record.signature.size() == 2) {
    // N -> P nu (P = pi0, eta, etaprime)
    double f,x;
    if(record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Pi0) {
      f = 0.130; #GeV
      x = siren::utilities::Constants::Pi0Mass / hnl_mass;
    }
    else if(record.signature.secondary_types[1]==siren::dataclasses::ParticleType::Eta) f = 0.0816; #GeV
    else if(record.signature.secondary_types[1]==siren::dataclasses::ParticleType::EtaPrime) f = -0.0946; #GeV
    width = f * pow(siren::utilities::Constants::FermiConstant,2)*pow(f,2);

  }
}

std::vector<std::string> HNLDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"};
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
    if(primary==siren::dataclasses::ParticleType::NuF4) {
      // N -> nu P (P = pi0, eta, eta prime, rho0, omega, hadrons, phi)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuTau}) {
        signature.secondary_types[0] = particle;
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Pi0;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Eta;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Rho0;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Omega;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::EtaPrime;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Phi;
        signatures.push_back(signature);
      }
      // N -> l- P+ (P+ = pi, K , rho, K*, hadrons, D, Ds)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EMinus, siren::dataclasses::ParticleType::MuMinus, siren::dataclasses::ParticleType::TauMinus}) {
        signature.secondary_types[0] = particle;
        signature.secondary_types[1] = siren::dataclasses::ParticleType::PiPlus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::KPlus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::RhoPlus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::KStarPlus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::DPlus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::DsPlus;
        signatures.push_back(signature);
      }
    }
    else if(primary==siren::dataclasses::ParticleType::NuF4Bar) {
      // N -> nu P (P = pi0, eta, eta prime, rho0, omega, hadrons, phi)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTauBar}) {
        signature.secondary_types[0] = particle;
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Pi0;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Eta;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Rho0;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Omega;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::EtaPrime;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Phi;
        signatures.push_back(signature);
      }
      // N -> l- P+ (P+ = pi, K , rho, K*, hadrons, D, Ds)
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EPlus, siren::dataclasses::ParticleType::MuPlus, siren::dataclasses::ParticleType::TauPlus}) {
        signature.secondary_types[0] = particle;
        signature.secondary_types[1] = siren::dataclasses::ParticleType::PiMinus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::KMinus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::RhoMinus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::KStarMinus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::Hadrons;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::DMinus;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::DsMinus;
        signatures.push_back(signature);
      }
    }
    
    // Three body decays (from 2007.03701)
    signature.secondary_types.resize(3);
    if(primary==siren::dataclasses::ParticleType::NuF4) {
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuTau}) {
        signature.secondary_types[0] = particle;
        // N -> nu nu nu
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::NuEBar;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::NuMuBar;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::NuTauBar;
        signatures.push_back(signature);
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
      }
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EMinus, siren::dataclasses::ParticleType::MuMinus, siren::dataclasses::ParticleType::TauMinus}) {
        signature.secondary_types[0] = particle;
        // N -> nu l- l+ (l different flavor)
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::EPlus;
        if(particle != siren::dataclasses::ParticleType::EMinus) signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::MuPlus;
        if(particle != siren::dataclasses::ParticleType::MuMinus) signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
        if(particle != siren::dataclasses::ParticleType::TauMinus) signatures.push_back(signature);
      }
    }
    else if(primary==siren::dataclasses::ParticleType::NuF4Bar) {
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTauBar}) {
        signature.secondary_types[0] = particle;
        // N -> nu nu nu
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::NuEBar;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::NuMuBar;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::NuTauBar;
        signatures.push_back(signature);
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
      }
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::EPlus, siren::dataclasses::ParticleType::MuPlus, siren::dataclasses::ParticleType::TauPlus}) {
        signature.secondary_types[0] = particle;
        // N -> nu l- l+ (l different flavor)
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::EMinus;
        if(particle != siren::dataclasses::ParticleType::EMinus) signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::MuMinus;
        if(particle != siren::dataclasses::ParticleType::MuMinus) signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTauBar;
        signature.secondary_types[2] = siren::dataclasses::ParticleType::TauPlus;
        if(particle != siren::dataclasses::ParticleType::TauMinus) signatures.push_back(signature);
      }

    return signatures;
}
}

double HNLDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    double DecayWidth = TotalDecayWidthForFinalState(record);
    if(nature==ChiralNature::Majorana) {
      //TODO: make sure factor of 2 is correct here
      return DecayWidth/2.;
    }
    siren::math::Vector3D hnl_dir = siren::math::Vector3D(record.primary_momentum[0],
                                                    record.primary_momentum[1],
                                                    record.primary_momentum[2]);
    hnl_dir.normalize();
    unsigned int gamma_index = (record.signature.secondary_types[0] == siren::dataclasses::ParticleType::Gamma) ? 0 : 1;
    rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::P4 pGamma(geom3::Vector3(record.secondary_momenta[gamma_index][1], record.secondary_momenta[gamma_index][2], record.secondary_momenta[gamma_index][3]), record.secondary_masses[gamma_index]);
    rk::Boost boost_to_HNL_rest = pHNL.restBoost();
    rk::P4 pGamma_HNLrest = pGamma.boost(boost_to_HNL_rest);
    
    siren::math::Vector3D gamma_dir = siren::math::Vector3D(pGamma_HNLrest.px(),
                                                      pGamma_HNLrest.py(),
                                                      pGamma_HNLrest.pz());
    gamma_dir.normalize();
    double CosThetaGamma = gamma_dir*hnl_dir; // scalar product
    double alpha = std::copysign(1.0,record.primary_helicity); // 1 for RH, -1 for LH
    alpha = (record.signature.primary_type == siren::dataclasses::ParticleType::NuF4) ? -1*alpha : alpha;
    return DecayWidth/2. * (1 + alpha*CosThetaGamma);
}

void HNLDecay::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<siren::utilities::LI_random> random) const {
    
    unsigned int gamma_index = (record.signature.secondary_types[0] == siren::dataclasses::ParticleType::Gamma) ? 0 : 1;
    unsigned int nu_index = 1 - gamma_index;

    double CosTheta;
    double alpha = std::copysign(1.0,record.primary_helicity); // 1 for RH, -1 for LH
    alpha = (record.signature.primary_type == siren::dataclasses::ParticleType::NuF4) ? -1*alpha : alpha;
    if(nature==ChiralNature::Majorana) {
      CosTheta = random->Uniform(-1,1);
    }
    else {
      double X = random->Uniform(0,1);
      CosTheta = (std::sqrt(1  - 2*alpha*(1 - alpha/2. - 2*X)) - 1)/alpha;
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

    rk::P4 pGamma_HNLrest(hnl_mass/2.*geom3::Vector3(CosTheta,SinTheta,0),0);
    pGamma_HNLrest.rotate(x_to_pHNL_rot);
    pGamma_HNLrest.rotate(rand_rot);

    rk::P4 pGamma = pGamma_HNLrest.boost(boost_to_lab);
    rk::P4 pNu(pHNL.momentum() - pGamma.momentum(),0); // ensures the neutrino has zero mass, avoids rounding errors

    record.secondary_momenta.resize(2);
    record.secondary_masses.resize(2);
    record.secondary_helicities.resize(2);
    
    record.secondary_momenta[gamma_index][0] = pGamma.e(); // pGamma_energy
    record.secondary_momenta[gamma_index][1] = pGamma.px(); // pGamma_x
    record.secondary_momenta[gamma_index][2] = pGamma.py(); // pGamma_y
    record.secondary_momenta[gamma_index][3] = pGamma.pz(); // pGamma_z
    record.secondary_masses[gamma_index] = pGamma.m();
    record.secondary_helicities[gamma_index] = 0;

    record.secondary_momenta[nu_index][0] = pNu.e(); // pNu_energy
    record.secondary_momenta[nu_index][1] = pNu.px(); // pNu_x
    record.secondary_momenta[nu_index][2] = pNu.py(); // pNu_y
    record.secondary_momenta[nu_index][3] = pNu.pz(); // pNu_z
    record.secondary_masses[nu_index] = pNu.m();
    record.secondary_helicities[nu_index] = -1*record.primary_helicity;

}

double HNLDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}



} // namespace interactions
} // namespace siren

