#include "SIREN/interactions/NeutrissimoDecay.h"

#include <array>                                              // for array
#include <cmath>                                              // for copysign
#include <tuple>                                              // for tie
#include <string>                                             // for basic_s...

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "SIREN/interactions/Decay.h"               // for Decay
#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/math/Vector3D.h"                     // for Vector3D
#include "SIREN/utilities/Constants.h"               // for GeV, pi
#include "SIREN/utilities/Random.h"                  // for LI_random

namespace siren {
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
    return TotalDecayWidth(record.signature.primary_type);
}

double NeutrissimoDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
    double total_coupling_sq = 0;
    for(auto dc : dipole_coupling) total_coupling_sq += dc*dc;
    return total_coupling_sq * std::pow(hnl_mass,3) / (4*siren::utilities::Constants::pi) * siren::utilities::Constants::GeV;
}

double NeutrissimoDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
    siren::dataclasses::InteractionSignature const & signature = record.signature;
    unsigned int gamma_index = (signature.secondary_types[0] == siren::dataclasses::ParticleType::Gamma) ? 0 : 1;
    unsigned int nu_index = 1 - gamma_index;
    double dipole_coupling_sq = 0;
    if(signature.secondary_types[nu_index]==siren::dataclasses::ParticleType::NuE ||
       signature.secondary_types[nu_index]==siren::dataclasses::ParticleType::NuEBar)
        dipole_coupling_sq = dipole_coupling[0]*dipole_coupling[0];
    else if(signature.secondary_types[nu_index]==siren::dataclasses::ParticleType::NuMu ||
            signature.secondary_types[nu_index]==siren::dataclasses::ParticleType::NuMuBar)
        dipole_coupling_sq = dipole_coupling[1]*dipole_coupling[1];
    else if(signature.secondary_types[nu_index]==siren::dataclasses::ParticleType::NuTau ||
            signature.secondary_types[nu_index]==siren::dataclasses::ParticleType::NuTauBar)
        dipole_coupling_sq = dipole_coupling[2]*dipole_coupling[2];
    return dipole_coupling_sq * std::pow(hnl_mass,3) / (4*siren::utilities::Constants::pi) * siren::utilities::Constants::GeV;
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

std::vector<dataclasses::InteractionSignature> NeutrissimoDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = siren::dataclasses::ParticleType::Decay;
    signature.secondary_types.resize(2);
    signature.secondary_types[0] = siren::dataclasses::ParticleType::Gamma;
    if(primary==siren::dataclasses::ParticleType::NuF4) {
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuTau}) {
        signature.secondary_types[1] = particle;
        signatures.push_back(signature);
      }
    }
    else if(primary==siren::dataclasses::ParticleType::NuF4Bar) {
      for(auto particle : std::vector<siren::dataclasses::ParticleType>{siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTauBar}) {
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

    siren::dataclasses::InteractionSignature const & signature = record.signature;

    siren::math::Vector3D hnl_dir = siren::math::Vector3D(record.primary_momentum[0],
                                                    record.primary_momentum[1],
                                                    record.primary_momentum[2]);
    hnl_dir.normalize();
    unsigned int gamma_index = (signature.secondary_types[0] == siren::dataclasses::ParticleType::Gamma) ? 0 : 1;
    std::array<double, 4> const & gamma_momentum = record.secondary_momenta[gamma_index];
    rk::P4 pHNL(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::P4 pGamma(geom3::Vector3(gamma_momentum[1], gamma_momentum[2], gamma_momentum[3]), record.secondary_masses[gamma_index]);
    rk::Boost boost_to_HNL_rest = pHNL.restBoost();
    rk::P4 pGamma_HNLrest = pGamma.boost(boost_to_HNL_rest);

    siren::math::Vector3D gamma_dir = siren::math::Vector3D(pGamma_HNLrest.px(),
                                                      pGamma_HNLrest.py(),
                                                      pGamma_HNLrest.pz());
    gamma_dir.normalize();
    double CosThetaGamma = gamma_dir*hnl_dir; // scalar product
    double alpha = std::copysign(1.0, record.primary_helicity); // 1 for RH, -1 for LH
    alpha = (signature.primary_type == siren::dataclasses::ParticleType::NuF4) ? -1*alpha : alpha;
    return DecayWidth/2. * (1 + alpha*CosThetaGamma);
}

void NeutrissimoDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::LI_random> random) const {

    siren::dataclasses::InteractionSignature const & signature = record.GetSignature();

    unsigned int gamma_index = (signature.secondary_types[0] == siren::dataclasses::ParticleType::Gamma) ? 0 : 1;
    unsigned int nu_index = 1 - gamma_index;

    double CosTheta;
    double alpha = std::copysign(1.0,record.GetPrimaryHelicity()); // 1 for RH, -1 for LH
    alpha = (signature.primary_type == siren::dataclasses::ParticleType::NuF4) ? -1*alpha : alpha;

    if(nature == ChiralNature::Majorana) {
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


    siren::dataclasses::SecondaryParticleRecord & gamma = record.GetSecondaryParticleRecord(gamma_index);
    siren::dataclasses::SecondaryParticleRecord & nu = record.GetSecondaryParticleRecord(nu_index);

    assert(gamma.type == siren::dataclasses::ParticleType::Gamma);
    assert(nu.type == siren::dataclasses::ParticleType::NuE ||
           nu.type == siren::dataclasses::ParticleType::NuMu ||
           nu.type == siren::dataclasses::ParticleType::NuTau ||
           nu.type == siren::dataclasses::ParticleType::NuEBar ||
           nu.type == siren::dataclasses::ParticleType::NuMuBar ||
           nu.type == siren::dataclasses::ParticleType::NuTauBar);

    gamma.SetFourMomentum({pGamma.e(), pGamma.px(), pGamma.py(), pGamma.pz()});
    gamma.SetMass(pGamma.m());
    gamma.SetHelicity(std::copysign(1.0, record.primary_helicity));

    nu.SetFourMomentum({pNu.e(), pNu.px(), pNu.py(), pNu.pz()});
    nu.SetMass(pNu.m());
    nu.SetHelicity(-1*record.primary_helicity);
}

double NeutrissimoDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}



} // namespace interactions
} // namespace siren

