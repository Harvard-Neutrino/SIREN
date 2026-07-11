#pragma once
#ifndef SIREN_Math_Kinematics_H
#define SIREN_Math_Kinematics_H

#include <cmath>

namespace siren {
namespace math {

inline double Kallen(double a, double b, double c) {
    return a * a + b * b + c * c
         - 2.0 * (a * b + a * c + b * c);
}

// Breakup momentum for P -> A + B. The factorized form avoids cancellation
// in lambda(M^2,mA^2,mB^2) close to threshold.
inline double TwoBodyRestMomentum(double parent_mass,
                                  double mass_a,
                                  double mass_b) {
    if (!(parent_mass > 0.0)) return 0.0;
    double product =
        (parent_mass - mass_a - mass_b)
        * (parent_mass + mass_a + mass_b)
        * (parent_mass + mass_a - mass_b)
        * (parent_mass - mass_a + mass_b);
    if (!(product > 0.0)) return 0.0;
    return std::sqrt(product) / (2.0 * parent_mass);
}

} // namespace math
} // namespace siren

#endif // SIREN_Math_Kinematics_H
