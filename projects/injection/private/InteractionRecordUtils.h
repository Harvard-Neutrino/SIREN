#pragma once
#ifndef SIREN_Injection_InteractionRecordUtils_H
#define SIREN_Injection_InteractionRecordUtils_H

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/math/Vector3D.h"

#include <cstddef>
#include <stdexcept>
#include <string>

namespace siren {
namespace injection {
namespace detail {

struct FourVector {
    double e;
    siren::math::Vector3D p;
};

template<typename Momentum>
inline FourVector ReadFourVector(Momentum const & momentum) {
    return {
        momentum[0],
        siren::math::Vector3D(momentum[1], momentum[2], momentum[3])};
}

inline FourVector ReadPrimary(
    siren::dataclasses::InteractionRecord const & record)
{
    return ReadFourVector(record.primary_momentum);
}

inline FourVector ReadSecondary(
    siren::dataclasses::InteractionRecord const & record,
    std::size_t index)
{
    return ReadFourVector(record.secondary_momenta[index]);
}

inline void WriteSecondary(
    siren::dataclasses::InteractionRecord & record,
    std::size_t index,
    FourVector const & momentum)
{
    record.secondary_momenta[index] = {
        momentum.e,
        momentum.p.GetX(),
        momentum.p.GetY(),
        momentum.p.GetZ()};
}

inline siren::math::Vector3D ReadVertex(
    siren::dataclasses::InteractionRecord const & record)
{
    return siren::math::Vector3D(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);
}

inline bool HasSecondaryStorage(
    siren::dataclasses::InteractionRecord const & record,
    std::size_t count)
{
    return record.signature.secondary_types.size() == count
        && record.secondary_masses.size() == count
        && record.secondary_momenta.size() == count;
}

inline void RequireSecondaryStorage(
    siren::dataclasses::InteractionRecord const & record,
    std::size_t count,
    char const * channel_name)
{
    if (!HasSecondaryStorage(record, count)) {
        throw std::runtime_error(
            std::string(channel_name) + " requires exactly "
            + std::to_string(count)
            + " secondary types, masses, and momenta");
    }
}

} // namespace detail
} // namespace injection
} // namespace siren

#endif // SIREN_Injection_InteractionRecordUtils_H
