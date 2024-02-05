#include "LeptonInjector/interactions/NeutrissimoDecay.h"

#include <array>                                              // for array
#include <cmath>                                              // for copysign
#include <tuple>                                              // for tie
#include <string>                                             // for basic_s...

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "LeptonInjector/interactions/Decay.h"               // for Decay
#include "LeptonInjector/dataclasses/InteractionRecord.h"     // for Interac...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle
#include "LeptonInjector/math/Vector3D.h"                     // for Vector3D
#include "LeptonInjector/utilities/Constants.h"               // for GeV, pi
#include "LeptonInjector/utilities/Random.h"                  // for LI_random

namespace LI {
namespace interactions {

bool NeutrissimoDecay::equal(Decay const & other) const {
    const NeutrissimoDecay* x = dynamic_cast<const NeutrissimoDecay*>(&other);

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

double NeutrissimoDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.GetPrimaryType());
}

double NeutrissimoDecay::TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const {
    double total_coupling_sq = 0;
    for(auto dc : dipole_coupling) total_coupling_sq += dc*dc;
    return total_coupling_sq * std::pow(hnl_mass,3) / (4*LI::utilities::Constants::pi) * LI::utilities::Constants::GeV;
}

double NeutrissimoDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
    LI::dataclasses::InteractionSignature const & signature = record.GetSignature();
    unsigned int gamma_index = (signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::Gamma) ? 0 : 1;
    unsigned int nu_index = 1 - gamma_index;
    double dipole_coupling_sq = 0;
    if(signature.secondary_types[nu_index]==LI::dataclasses::Particle::ParticleType::NuE ||
       signature.secondary_types[nu_index]==LI::dataclasses::Particle::ParticleType::NuEBar)
        dipole_coupling_sq = dipole_coupling[0]*dipole_coupling[0];
    else if(signature.secondary_types[nu_index]==LI::dataclasses::Particle::ParticleType::NuMu ||
            signature.secondary_types[nu_index]==LI::dataclasses::Particle::ParticleType::NuMuBar)
        dipole_coupling_sq = dipole_coupling[1]*dipole_coupling[1];
    else if(signature.secondary_types[nu_index]==LI::dataclasses::Particle::ParticleType::NuTau ||
            signature.secondary_types[nu_index]==LI::dataclasses::Particle::ParticleType::NuTauBar)
        dipole_coupling_sq = dipole_coupling[2]*dipole_coupling[2];
    return dipole_coupling_sq * std::pow(hnl_mass,3) / (4*LI::utilities::Constants::pi) * LI::utilities::Constants::GeV;
}

std::vector<std::string> NeutrissimoDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"};
}


std::vector<dataclasses::InteractionSignature> NeutrissimoDecay::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
        std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
        signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> NeutrissimoDecay::GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = LI::dataclasses::Particle::ParticleType::Decay;
    signature.secondary_types.resize(2);
    signature.secondary_types[0] = LI::dataclasses::Particle::ParticleType::Gamma;
    if(primary==LI::dataclasses::Particle::ParticleType::NuF4) {
      for(auto particle : std::vector<LI::dataclasses::Particle::ParticleType>{LI::dataclasses::Particle::ParticleType::NuE, LI::dataclasses::Particle::ParticleType::NuMu, LI::dataclasses::Particle::ParticleType::NuTau}) {
        signature.secondary_types[1] = particle;
        signatures.push_back(signature);
      }
    }
    else if(primary==LI::dataclasses::Particle::ParticleType::NuF4Bar) {
      for(auto particle : std::vector<LI::dataclasses::Particle::ParticleType>{LI::dataclasses::Particle::ParticleType::NuEBar, LI::dataclasses::Particle::ParticleType::NuMuBar, LI::dataclasses::Particle::ParticleType::NuTauBar}) {
        signature.secondary_types[1] = particle;
        signatures.push_back(signature);
      }
    }
    return signatures;
}

double NeutrissimoDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    double DecayWidth = TotalDecayWidthForFinalState(record);
    if(nature==ChiralNature::Majorana) {
      //TODO: make sure factor of 2 is correct here
      return DecayWidth/2.;
    }

    LI::dataclasses::Particle primary = record.GetPrimary();
    LI::dataclasses::InteractionSignature const & signature = record.GetSignature();
    std::vector<LI::dataclasses::Particle> secondaries = record.GetSecondaries();
    secondaries.resize(2);
    secondaries[0].type = signature.secondary_types[0];
    secondaries[1].type = signature.secondary_types[1];

    LI::math::Vector3D hnl_dir = LI::math::Vector3D(primary.momentum[0],
                                                    primary.momentum[1],
                                                    primary.momentum[2]);
    hnl_dir.normalize();
    unsigned int gamma_index = (signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::Gamma) ? 0 : 1;
    LI::dataclasses::Particle gamma = secondaries[gamma_index];
    rk::P4 pHNL(geom3::Vector3(primary.momentum[1], primary.momentum[2], primary.momentum[3]), primary.mass);
    rk::P4 pGamma(geom3::Vector3(gamma.momentum[1], gamma.momentum[2], gamma.momentum[3]), gamma.mass);
    rk::Boost boost_to_HNL_rest = pHNL.restBoost();
    rk::P4 pGamma_HNLrest = pGamma.boost(boost_to_HNL_rest);

    LI::math::Vector3D gamma_dir = LI::math::Vector3D(pGamma_HNLrest.px(),
                                                      pGamma_HNLrest.py(),
                                                      pGamma_HNLrest.pz());
    gamma_dir.normalize();
    double CosThetaGamma = gamma_dir*hnl_dir; // scalar product
    double alpha = std::copysign(1.0,record.GetPrimaryHelicity()); // 1 for RH, -1 for LH
    alpha = (signature.primary_type == LI::dataclasses::Particle::ParticleType::NuF4) ? -1*alpha : alpha;
    return DecayWidth/2. * (1 + alpha*CosThetaGamma);
}

void NeutrissimoDecay::SampleFinalState(dataclasses::InteractionRecord & record, std::shared_ptr<LI::utilities::LI_random> random) const {

    LI::dataclasses::InteractionSignature const & signature = record.GetSignature();
    LI::dataclasses::Particle primary = record.GetPrimary();
    std::vector<LI::dataclasses::Particle> secondaries = record.GetSecondaries();
    secondaries.resize(2);
    secondaries[0].type = signature.secondary_types[0];
    secondaries[1].type = signature.secondary_types[1];

    unsigned int gamma_index = (signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::Gamma) ? 0 : 1;
    unsigned int nu_index = 1 - gamma_index;

    double CosTheta;
    double alpha = std::copysign(1.0,record.GetPrimaryHelicity()); // 1 for RH, -1 for LH
    alpha = (signature.primary_type == LI::dataclasses::Particle::ParticleType::NuF4) ? -1*alpha : alpha;
    if(nature==ChiralNature::Majorana) {
      CosTheta = random->Uniform(-1,1);
    }
    else {
      double X = random->Uniform(0,1);
      CosTheta = (std::sqrt(1  - 2*alpha*(1 - alpha/2. - 2*X)) - 1)/alpha;
    }
    double SinTheta = std::sin(std::acos(CosTheta));

    rk::P4 pHNL(geom3::Vector3(primary.momentum[1], primary.momentum[2], primary.momentum[3]), primary.mass);
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

    secondaries[gamma_index].momentum[0] = pGamma.e(); // pGamma_energy
    secondaries[gamma_index].momentum[1] = pGamma.px(); // pGamma_x
    secondaries[gamma_index].momentum[2] = pGamma.py(); // pGamma_y
    secondaries[gamma_index].momentum[3] = pGamma.pz(); // pGamma_z
    secondaries[gamma_index].mass = pGamma.m();

    secondaries[nu_index].momentum[0] = pNu.e(); // pNu_energy
    secondaries[nu_index].momentum[1] = pNu.px(); // pNu_x
    secondaries[nu_index].momentum[2] = pNu.py(); // pNu_y
    secondaries[nu_index].momentum[3] = pNu.pz(); // pNu_z
    secondaries[nu_index].mass = pNu.m();
    secondaries[nu_index].helicity = -1*record.GetPrimaryHelicity();

    record.SetSecondaries(secondaries);
}

double NeutrissimoDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}



} // namespace interactions
} // namespace LI

