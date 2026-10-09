#pragma once
#include <cmath>
#include <limits>
#include "InteractionRecordUtils.h"
#include "SIREN/utilities/Errors.h"

namespace siren { namespace injection { namespace detail {
// Tabulated float32 beam rows need not satisfy E^2-p^2=m^2 in double
// precision. The declared mass and three-momentum define the on-shell boost.
// Permit only a small input-energy discrepancy. The record keeps its input
// energy: the parent belongs to the upstream vertex, whose densities read it.
// Daughters conserve the on-shell four-momentum OnShellParent(r).
inline bool OnShellParentValid(siren::dataclasses::InteractionRecord const & r) {
    auto const & p=r.primary_momentum;
    if (!(r.primary_mass>0) || !std::isfinite(r.primary_mass)
        || !std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) || !std::isfinite(p[3])) return false;
    double energy=std::hypot(r.primary_mass,std::hypot(p[1],p[2],p[3]));
    return std::abs(p[0]-energy)<=2e-5*energy;
}
inline FourVector OnShellParent(siren::dataclasses::InteractionRecord const & r) {
    auto p=ReadPrimary(r);
    p.e=std::hypot(r.primary_mass,p.p.magnitude());
    return p;
}
inline void RequireOnShellParent(siren::dataclasses::InteractionRecord const & r) {
    if (!OnShellParentValid(r)) throw siren::utilities::InjectionFailure(
        siren::utilities::FailureReason::KinematicallyForbidden,
        "Decay parent needs a positive mass and an energy within 2e-5 of its mass shell");
}
}}}
