#pragma once
#ifndef SIREN_Injection_OnShellDecayKinematics_H
#define SIREN_Injection_OnShellDecayKinematics_H

#include <cmath>

#include "SIREN/dataclasses/RecordKinematics.h"
#include "SIREN/utilities/Errors.h"

namespace siren {
namespace injection {
namespace detail {

constexpr double kOnShellEnergyTolerance = 2e-5;

inline bool OnShellParentValid(siren::dataclasses::InteractionRecord const & record) {
    if (!(record.primary_mass > 0.0) || !std::isfinite(record.primary_mass))
        return false;
    auto const & momentum = record.primary_momentum;
    for (double component : momentum) {
        if (!std::isfinite(component))
            return false;
    }
    double shell_energy = std::hypot(
        record.primary_mass, std::hypot(momentum[1], momentum[2], momentum[3]));
    return std::abs(momentum[0] - shell_energy) <= kOnShellEnergyTolerance * shell_energy;
}

inline siren::math::FourMomentum OnShellParent(siren::dataclasses::InteractionRecord const & record) {
    return siren::math::FourMomentum::FromMass(
        record.primary_mass, siren::dataclasses::PrimaryMomentum(record).Momentum());
}

inline void RequireOnShellParent(siren::dataclasses::InteractionRecord const & record) {
    if (!OnShellParentValid(record)) {
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::KinematicallyForbidden,
            "Decay parent needs a positive mass and an energy within a relative "
            "2e-5 of its mass shell");
    }
}

} // namespace detail
} // namespace injection
} // namespace siren

#endif // SIREN_Injection_OnShellDecayKinematics_H
