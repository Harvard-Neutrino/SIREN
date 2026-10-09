#pragma once
#ifndef SIREN_Injection_OnShellDecayKinematics_H
#define SIREN_Injection_OnShellDecayKinematics_H

#include <cmath>

#include "InteractionRecordUtils.h"
#include "SIREN/utilities/Errors.h"

namespace siren {
namespace injection {
namespace detail {

// Tabulated float32 beam rows need not satisfy E^2 - p^2 = m^2 in double
// precision, so decays boost with the declared mass and three-momentum and
// accept an input energy within this relative tolerance of that mass shell.
// The record keeps its input energy, which the upstream vertex's densities
// read.
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

// The parent four-momentum with its energy on the mass shell. Daughters
// conserve this four-momentum.
inline FourVector OnShellParent(siren::dataclasses::InteractionRecord const & record) {
    FourVector parent = ReadPrimary(record);
    parent.e = std::hypot(record.primary_mass, parent.p.magnitude());
    return parent;
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
