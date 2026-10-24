#pragma once
#ifndef SIREN_Math_FourMomentum_H
#define SIREN_Math_FourMomentum_H

#include <array>
#include <cmath>
#include <stdexcept>
#include <rk/rk.hh>

#include "SIREN/math/Vector3D.h"

namespace siren { namespace math {

// Value semantics; components are (E, px, py, pz) in GeV, metric (+,-,-,-).
class FourMomentum {
public:
    static FourMomentum FromComponents(std::array<double, 4> const & p) {
        return FourMomentum(rk::P4(p[0], geom3::Vector3(p[1], p[2], p[3])), false);
    }

    static FourMomentum FromMass(double mass, Vector3D const & p,
                                 bool negative_energy = false) {
        if (!std::isfinite(mass) || mass < 0.0
            || !std::isfinite(p.GetX()) || !std::isfinite(p.GetY())
            || !std::isfinite(p.GetZ()))
            throw std::domain_error("On-shell momentum needs finite components and a nonnegative mass");
        auto result = FourMomentum(rk::P4(
            geom3::Vector3(p.GetX(), p.GetY(), p.GetZ()), mass, negative_energy), true);
        result.RequireFinite();
        return result;
    }

    std::array<double, 4> Components() const {
        return {value_.e(), value_.px(), value_.py(), value_.pz()};
    }
    double Energy() const { return value_.e(); }
    Vector3D Momentum() const { return {value_.px(), value_.py(), value_.pz()}; }
    bool HasKnownMass() const { return known_mass_; }

    // This is the signed invariant of the stored components, not the cached mass.
    double MassSquared() const { return value_.squared(); }
    double Mass() const {
        RequireFinite();
        if (!known_mass_ && !(MassSquared() >= 0.0))
            throw std::domain_error("A spacelike four-momentum has no particle mass");
        double mass = value_.m();
        if (!std::isfinite(mass))
            throw std::domain_error("Four-momentum mass is not finite");
        return mass;
    }

    FourMomentum BoostToLab(FourMomentum const & frame) const {
        return Boost(frame, false);
    }
    FourMomentum BoostToRest(FourMomentum const & frame) const {
        return Boost(frame, true);
    }

    friend FourMomentum operator+(FourMomentum const & a, FourMomentum const & b) {
        return FourMomentum(a.value_ + b.value_, false);
    }
    friend FourMomentum operator-(FourMomentum const & a, FourMomentum const & b) {
        return FourMomentum(a.value_ - b.value_, false);
    }

    static double InvariantMass(FourMomentum const & a, FourMomentum const & b) {
        double ma = a.Mass(), mb = b.Mass();
        double squared = ma * ma + mb * mb + 2.0 * a.value_.dot(b.value_);
        if (!std::isfinite(squared) || squared < 0.0)
            throw std::domain_error("The pair has no finite timelike invariant mass");
        return rk::invMass(a.value_, b.value_);
    }

private:
    rk::P4 value_;
    bool known_mass_;

    FourMomentum(rk::P4 const & value, bool known_mass)
        : value_(value), known_mass_(known_mass) {}

    void RequireFinite() const {
        for (double component : Components())
            if (!std::isfinite(component))
                throw std::domain_error("Four-momentum components must be finite");
    }

    FourMomentum Boost(FourMomentum const & frame, bool inverse) const {
        Mass();
        if (!(frame.Mass() > 0.0) || !(frame.Energy() > 0.0))
            throw std::domain_error("A boost frame needs positive mass and energy");
        rk::P4 result = value_;
        result.boost(inverse ? frame.value_.restBoost() : frame.value_.labBoost());
        FourMomentum momentum(result, true);
        momentum.RequireFinite();
        return momentum;
    }
};

} }
#endif
