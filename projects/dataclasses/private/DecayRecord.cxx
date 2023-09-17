#include "LeptonInjector/dataclasses/DecayRecord.h"

#include <tuple>

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/DecaySignature.h"

namespace LI {
namespace dataclasses {

bool DecayRecord::operator==(DecayRecord const & other) const {
    return std::tie(
        signature,
        primary_mass,
        primary_momentum,
        primary_helicity,
        decay_vertex,
        secondary_masses,
        secondary_momenta,
        secondary_helicity,
        decay_parameters)
        ==
        std::tie(
        other.signature,
        other.primary_mass,
        other.primary_momentum,
        other.primary_helicity,
        other.decay_vertex,
        other.secondary_masses,
        other.secondary_momenta,
        other.secondary_helicity,
        other.decay_parameters);
}

} // namespace dataclasses
} // namespace LI
