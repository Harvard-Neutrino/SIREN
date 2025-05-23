#include "SIREN/interactions/ScalarDecay.h"

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
#include "SIREN/utilities/Random.h"                  // for SIREN_random

namespace siren {
namespace interactions {

bool ScalarDecay::equal(Decay const & other) const {
    const ScalarDecay* x = dynamic_cast<const ScalarDecay*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(
                    primary_types,
                    mass,
                    dipole_coupling)
            ==
            std::tie(
                    x->primary_types,
                    x->mass,
                    x->coupling);
}

double ScalarDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidth(record.signature.primary_type);
}

double ScalarDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
    return coupling*coupling * mass / (16*siren::utilities::Constants::pi) * siren::utilities::Constants::GeV;
}

double ScalarDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & record) const {
    siren::dataclasses::InteractionSignature const & signature = record.signature;
    if(signature.secondary_types.size()!=2) return 0;
    // make sure this is a dimuon decay
    if(!(signature.secondary_types[0]==siren::dataclasses::ParticleType::MuPlus ||
         signature.secondary_types[0]==siren::dataclasses::ParticleType::MuMinus)) return 0;
    if(!(signature.secondary_types[1]==siren::dataclasses::ParticleType::MuPlus ||
         signature.secondary_types[1]==siren::dataclasses::ParticleType::MuMinus)) return 0;
    // if so, return the total decay width
    // this could be updated to just the dimuon decay width in the future
    // right now, requires user to apply branching ratio on the back end
    return TotalDecayWidth(signature.primary_type);
}

std::vector<std::string> ScalarDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"};
}


std::vector<dataclasses::InteractionSignature> ScalarDecay::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
        std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
        signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end());
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> ScalarDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = siren::dataclasses::ParticleType::Decay;
    signature.secondary_types.resize(2);
    if(primary==siren::dataclasses::ParticleType::Scalar) {
        signature.secondary_types[0] = siren::dataclasses::ParticleType::MuPlus;
        signature.secondary_types[1] = siren::dataclasses::ParticleType::MuMinus;
        signatures.push_back(signature);
    }
    return signatures;
}

double ScalarDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & record) const {
    double DecayWidth = TotalDecayWidthForFinalState(record);
    return DecayWidth/2.; // for costheta between -1,1
}

void ScalarDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {

    siren::dataclasses::InteractionSignature const & signature = record.GetSignature();


    double CosTheta = random->Uniform(-1,1);
    double SinTheta = std::sin(std::acos(CosTheta));

    rk::P4 pScalar(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::Boost boost_to_lab = pScalar.labBoost();

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 pScalar_mom = pScalar.momentum();
    geom3::UnitVector3 pScalar_dir = pScalar_mom.direction();
    geom3::Rotation3 x_to_pScalar_rot = geom3::rotationBetween(x_dir, pScalar_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(pScalar_dir, phi);

    double muonMomentum = sqrt(std::pow(mass/2.,2) - std::pow(siren::utilities::Constants::muonMass,2));
    rk::P4 pMuPlus_ScalarRest(muonMomentum*geom3::Vector3(CosTheta,SinTheta,0),siren::utilities::Constants::muonMass);
    pMuPlus_ScalarRest.rotate(x_to_pScalar_rot);
    pMuPlus_ScalarRest.rotate(rand_rot);

    rk::P4 pMuPlus = pMuPlus_ScalarRest.boost(boost_to_lab);
    rk::P4 pMuMinus(pScalar.momentum() - pMuPlus.momentum(),siren::utilities::Constants::muonMass);


    siren::dataclasses::SecondaryParticleRecord & MuPlus = record.GetSecondaryParticleRecord(0);
    siren::dataclasses::SecondaryParticleRecord & MuMinus = record.GetSecondaryParticleRecord(1);

    assert(MuPlus.type == siren::dataclasses::ParticleType::MuPlus);
    assert(MuMinus.type == siren::dataclasses::ParticleType::MuMinus);


    MuPlus.SetFourMomentum({pMuPlus.e(), pMuPlus.px(), pMuPlus.py(), pMuPlus.pz()});
    MuPlus.SetMass(pMuPlus.m());
    MuPlus.SetHelicity(0.5); // right-handed

    MuMinus.SetFourMomentum({pMuMinus.e(), pMuMinus.px(), pMuMinus.py(), pMuMinus.pz()});
    MuPlus.SetMass(pMuMinus.m());
    MuPlus.SetHelicity(-0.5); // left-handed
}

double ScalarDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}

} // namespace interactions
} // namespace siren

