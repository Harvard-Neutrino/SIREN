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
// Known-mass boost, avoiding an invariant mass reconstructed by subtracting
// large lab energies. Long-double intermediates limit cancellation on inversion
// where long double is wider than double. A forward boost of a backward
// rest-frame vector still cancels terms of size gamma*(E*+p*), so each lab
// component carries an absolute error of order eps*gamma*(E*+p*).
inline FourVector BoostOnShell(FourVector const & frame,double mass,FourVector const & p,bool inverse=false) {
    long double momentum=frame.p.magnitude();
    if (momentum==0) return p;
    auto axis=frame.p/static_cast<double>(momentum);
    long double parallel=static_cast<long double>(p.p.GetX())*axis.GetX()
        +static_cast<long double>(p.p.GetY())*axis.GetY()+static_cast<long double>(p.p.GetZ())*axis.GetZ();
    long double energy=std::hypot(momentum,static_cast<long double>(mass));
    long double sign=inverse ? -1 : 1;
    long double out_e=(energy*p.e+sign*momentum*parallel)/mass;
    long double out_p=(energy*parallel+sign*momentum*p.e)/mass;
    return {static_cast<double>(out_e),p.p+axis*static_cast<double>(out_p-parallel)};
}
}}}
