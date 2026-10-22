#pragma once
#ifndef SIREN_DistributionVariable_H
#define SIREN_DistributionVariable_H

#include <set>
#include <string>
#include <stdexcept>

namespace siren {
namespace distributions {

enum class DistributionVariable {
    PrimaryMass,
    PrimaryEnergy,
    PrimaryDirection,
    PrimaryHelicity,
    PrimaryArea,
    InitialPosition,
    InteractionVertex,
    InteractionParameters,
};

inline std::string to_string(DistributionVariable var) {
    switch (var) {
        case DistributionVariable::PrimaryMass:       return "PrimaryMass";
        case DistributionVariable::PrimaryEnergy:     return "PrimaryEnergy";
        case DistributionVariable::PrimaryDirection:  return "PrimaryDirection";
        case DistributionVariable::PrimaryHelicity:   return "PrimaryHelicity";
        case DistributionVariable::PrimaryArea:       return "PrimaryArea";
        case DistributionVariable::InitialPosition:   return "InitialPosition";
        case DistributionVariable::InteractionVertex: return "InteractionVertex";
        case DistributionVariable::InteractionParameters: return "InteractionParameters";
    }
    throw std::runtime_error("Unknown DistributionVariable");
}

} // namespace distributions
} // namespace siren

#endif // SIREN_DistributionVariable_H
