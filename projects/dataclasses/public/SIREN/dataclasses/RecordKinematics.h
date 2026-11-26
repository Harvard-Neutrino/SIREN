#pragma once
#ifndef SIREN_Dataclasses_RecordKinematics_H
#define SIREN_Dataclasses_RecordKinematics_H

#include <cstddef>
#include <stdexcept>
#include <string>
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/math/FourMomentum.h"

namespace siren { namespace dataclasses {

inline math::FourMomentum PrimaryMomentum(InteractionRecord const & record) {
    return math::FourMomentum::FromComponents(record.primary_momentum);
}
inline math::FourMomentum SecondaryMomentum(InteractionRecord const & record, std::size_t index) {
    return math::FourMomentum::FromComponents(record.secondary_momenta.at(index));
}
inline void SetPrimaryMomentum(InteractionRecord & record, math::FourMomentum const & momentum) {
    record.primary_momentum = momentum.Components();
}
inline void SetSecondaryMomentum(InteractionRecord & record, std::size_t index,
                                 math::FourMomentum const & momentum) {
    record.secondary_momenta.at(index) = momentum.Components();
}
inline math::Vector3D InteractionVertex(InteractionRecord const & record) {
    return math::Vector3D(record.interaction_vertex);
}
inline bool HasSecondaryStorage(InteractionRecord const & record, std::size_t count) {
    return record.signature.secondary_types.size() == count
        && record.secondary_masses.size() == count
        && record.secondary_momenta.size() == count;
}
inline void RequireSecondaryStorage(InteractionRecord const & record, std::size_t count,
                                     char const * context) {
    if (!HasSecondaryStorage(record, count))
        throw std::runtime_error(std::string(context) + " requires exactly "
            + std::to_string(count) + " secondary types, masses, and momenta");
}

} }
#endif
