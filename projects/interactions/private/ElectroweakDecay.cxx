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

}

std::vector<std::string> ElectroweakDecay::DensityVariables() const {
    return std::vector<std::string>{"CosTheta"};
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
      signature.secondary_types[0] = siren::dataclasses::ParticleType::EPlus;
      signature.secondary_types[1] = siren::dataclasses::ParticleType::NuE;
      signatures.push_back(signature);
      signature.secondary_types[0] = siren::dataclasses::ParticleType::MuPlus;
      signature.secondary_types[1] = siren::dataclasses::ParticleType::NuMu;
      signatures.push_back(signature);
      signature.secondary_types[0] = siren::dataclasses::ParticleType::TauPlus;
      signature.secondary_types[1] = siren::dataclasses::ParticleType::NuTau;
      signatures.push_back(signature);

    }
    else if(primary==siren::dataclasses::ParticleType::WMinus) {

    }
    else if(primary==siren::dataclasses::ParticleType::Z0) {

    }

    return signatures;
}






} // namespace interactions
} // namespace siren