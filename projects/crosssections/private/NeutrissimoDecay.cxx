
#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/utilities/Errors.h"
#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Constants.h"

#include "LeptonInjector/detector/MaterialModel.h"

#include "LeptonInjector/crosssections/Decay.h"
#include "LeptonInjector/crosssections/NeutrissimoDecay.h"

namespace LI {
namespace crosssections {

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

double NeutrissimoDecay::TotalDecayWidth(dataclasses::InteractionRecord const & record) {
    return TotalDecayWidth(record.signature.primary_type);
}

double NeutrissimoDecay::TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary, double energy) {
    return std::pow(dipole_coupling,2) * std::pow(hnl_mass,3) / (4*LI::utilities::Constants::pi) * LI::utilities::Constants::GeV;
}


} // namespace crosssections
} // namespace LI

