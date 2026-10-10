#pragma once
#ifndef SIREN_TwoBodyKinematics_H
#define SIREN_TwoBodyKinematics_H

#include <array>
#include <cmath>
#include <utility>

#include "SIREN/math/Kinematics.h"

namespace siren {
namespace injection {

// Kallen (triangle) function: lambda(a, b, c) = a^2 + b^2 + c^2 - 2(ab + ac + bc)
inline double Kallen(double a, double b, double c) {
    return siren::math::Kallen(a, b, c);
}

// Rest-frame momentum magnitude for P -> A + B
// Returns |p*| = sqrt(lambda(M^2, mA^2, mB^2)) / (2*M)
inline double TwoBodyRestMomentum(double M, double mA, double mB) {
    return siren::math::TwoBodyRestMomentum(M, mA, mB);
}

// Rest-frame energy of daughter A in P -> A + B
// E_A* = (M^2 + mA^2 - mB^2) / (2*M)
inline double TwoBodyRestEnergy(double M, double mA, double mB) {
    return (M * M + mA * mA - mB * mB) / (2.0 * M);
}

struct TwoBodyLabSolution {
    double cos_theta_rest = 0.0;
    double p_lab = 0.0;
    double jacobian = 0.0;     // |dOmega_lab / dOmega_rest|
    bool valid = false;
};

// Given a parent moving with velocity beta (Lorentz factor gamma), a
// daughter with rest-frame momentum p_rest and energy E_rest, and a lab
// angle cos_theta_lab, solve for the rest-frame angle cos_theta_rest.
//
// For massive daughters, there can be 0, 1, or 2 solutions.
// Returns an array of up to 2 solutions.
//
// This implements the inverse of the Lorentz boost relation:
//   p_lab * cos(theta_lab) = gamma * (p_rest * cos(theta_rest) + beta * E_rest)
//   p_lab * sin(theta_lab) = p_rest * sin(theta_rest)
//
// Each solution also carries the Jacobian
//   |dOmega_lab/dOmega_rest| = gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab| / p_lab^2
//
// Parameters:
//   beta_parent:   parent velocity (|p_parent| / E_parent)
//   gamma_parent:  parent Lorentz factor (E_parent / M_parent)
//   p_rest:        daughter momentum magnitude in parent rest frame
//   E_rest:        daughter energy in parent rest frame
//   m_daughter:    daughter mass
//   cos_theta_lab: desired lab-frame angle (between daughter and parent direction)
std::array<TwoBodyLabSolution, 2> SolveLabAngle(
    double beta_parent,
    double gamma_parent,
    double p_rest,
    double E_rest,
    double m_daughter,
    double cos_theta_lab
);

// Compute the critical angle: the largest lab angle a massive daughter can
// reach. Returns cos(theta_critical). When the daughter's rest-frame velocity
// p_rest / E_rest is at least the parent velocity, as for a massless daughter,
// every lab angle is reachable and this returns -1.
double CriticalCosTheta(
    double beta_parent,
    double gamma_parent,
    double p_rest,
    double E_rest,
    double m_daughter
);

// Forward problem: given rest-frame angle, compute lab-frame angle
// and momentum. This is the standard Lorentz boost.
//
// Returns (cos_theta_lab, p_lab, E_lab).
struct LabFrameResult {
    double cos_theta_lab;
    double p_lab;
    double E_lab;
};

LabFrameResult BoostToLab(
    double beta_parent,
    double gamma_parent,
    double p_rest,
    double E_rest,
    double cos_theta_rest
);

} // namespace injection
} // namespace siren

#endif // SIREN_TwoBodyKinematics_H
