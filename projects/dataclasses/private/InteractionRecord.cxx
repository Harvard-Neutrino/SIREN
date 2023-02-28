#include "LeptonInjector/dataclasses/InteractionRecord.h"

namespace LI {
namespace dataclasses {

bool InteractionRecord::operator==(InteractionRecord const & other) const {
    return std::tie(
        signature,
        primary_mass,
        primary_momentum,
        primary_helicity,
        target_mass,
        target_momentum,
        target_helicity,
        interaction_vertex,
        secondary_masses,
        secondary_momenta,
        secondary_helicity,
        interaction_parameters)
        ==
        std::tie(
        other.signature,
        other.primary_mass,
        other.primary_momentum,
        other.primary_helicity,
        other.target_mass,
        other.target_momentum,
        other.target_helicity,
        other.interaction_vertex,
        other.secondary_masses,
        other.secondary_momenta,
        other.secondary_helicity,
        other.interaction_parameters);
}

} // namespace dataclasses
} // namespace LI
