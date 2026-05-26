#pragma once
#ifndef SIREN_PhaseSpaceJacobian_H
#define SIREN_PhaseSpaceJacobian_H

#include "SIREN/injection/TwoBodyKinematics.h"

#include <cmath>
#include <limits>

namespace siren {
namespace injection {
namespace phase_space_jacobian {

inline double SafeInverse(double value) {
    if (value == 0.0 || !std::isfinite(value)) {
        return 0.0;
    }
    return 1.0 / value;
}

inline double RestFrameToLabSolidAngleJacobian(
    double beta_parent,
    double gamma_parent,
    double p_rest,
    double e_rest,
    double daughter_mass,
    double cos_theta_lab,
    int solution_index = 0)
{
    auto solutions = SolveLabAngle(
        beta_parent, gamma_parent, p_rest, e_rest,
        daughter_mass, cos_theta_lab);
    if (solution_index < 0 || solution_index >= static_cast<int>(solutions.size())) {
        return 0.0;
    }
    if (!solutions[solution_index].valid) {
        return 0.0;
    }
    return solutions[solution_index].jacobian;
}

inline double LabToRestFrameSolidAngleJacobian(
    double beta_parent,
    double gamma_parent,
    double p_rest,
    double e_rest,
    double daughter_mass,
    double cos_theta_lab,
    int solution_index = 0)
{
    return SafeInverse(RestFrameToLabSolidAngleJacobian(
        beta_parent, gamma_parent, p_rest, e_rest,
        daughter_mass, cos_theta_lab, solution_index));
}

inline double MomentumInTwoBodyRestFrame(
    double parent_mass,
    double mass_a,
    double mass_b)
{
    return TwoBodyRestMomentum(parent_mass, mass_a, mass_b);
}

// For P -> spectator + (first second), followed by
// (first second) -> first + second, the Dalitz invariant
// s_spectator_first is linear in the helicity angle of first in
// the pair rest frame.  This returns
// |d s_spectator_first / d cos(theta_first)| = 2 p_spectator p_first
// evaluated in the pair rest frame.
inline double Recursive2BodyToDalitzAbsJacobian(
    double parent_mass,
    double spectator_mass,
    double first_mass,
    double second_mass,
    double pair_mass_squared)
{
    if (pair_mass_squared <= 0.0) return 0.0;
    double pair_mass = std::sqrt(pair_mass_squared);

    double lambda_parent = Kallen(
        parent_mass * parent_mass,
        spectator_mass * spectator_mass,
        pair_mass_squared);
    double lambda_pair = Kallen(
        pair_mass_squared,
        first_mass * first_mass,
        second_mass * second_mass);
    if (lambda_parent <= 0.0 || lambda_pair <= 0.0) return 0.0;

    double p_spectator_in_pair_frame = std::sqrt(lambda_parent) / (2.0 * pair_mass);
    double p_first_in_pair_frame = std::sqrt(lambda_pair) / (2.0 * pair_mass);
    return 2.0 * p_spectator_in_pair_frame * p_first_in_pair_frame;
}

inline double Recursive2BodyDensityToDalitzDensity(
    double recursive_density,
    double parent_mass,
    double spectator_mass,
    double first_mass,
    double second_mass,
    double pair_mass_squared)
{
    double jacobian = Recursive2BodyToDalitzAbsJacobian(
        parent_mass, spectator_mass, first_mass, second_mass,
        pair_mass_squared);
    return recursive_density * SafeInverse(jacobian);
}

inline double DalitzDensityToRecursive2BodyDensity(
    double dalitz_density,
    double parent_mass,
    double spectator_mass,
    double first_mass,
    double second_mass,
    double pair_mass_squared)
{
    double jacobian = Recursive2BodyToDalitzAbsJacobian(
        parent_mass, spectator_mass, first_mass, second_mass,
        pair_mass_squared);
    return dalitz_density * jacobian;
}

// Helicity-angle and recursive two-body coordinates differ only by
// rotations of the angular reference axes when the same pair
// invariant mass is used.  Rotations preserve dOmega.
inline double Recursive2BodyToHelicityAnglesJacobian() {
    return 1.0;
}

inline double HelicityAnglesToRecursive2BodyJacobian() {
    return 1.0;
}

inline double Q2FromBjorkenXY(
    double x,
    double y,
    double target_mass,
    double incident_energy)
{
    return 2.0 * target_mass * incident_energy * x * y;
}

inline double BjorkenXFromQ2Y(
    double q2,
    double y,
    double target_mass,
    double incident_energy)
{
    double denom = 2.0 * target_mass * incident_energy * y;
    if (denom == 0.0) return std::numeric_limits<double>::quiet_NaN();
    return q2 / denom;
}

inline double MandelstamTFromQ2(double q2) {
    return -q2;
}

// Fixed incident energy relation between (x, y) and (Q2, y).
// Since t = -Q2, this is also the absolute Jacobian for replacing
// Bjorken x by Mandelstam t at fixed y.
inline double BjorkenXYToQ2YAbsJacobian(
    double y,
    double target_mass,
    double incident_energy)
{
    return std::abs(2.0 * target_mass * incident_energy * y);
}

inline double Q2YToBjorkenXYAbsJacobian(
    double y,
    double target_mass,
    double incident_energy)
{
    return SafeInverse(BjorkenXYToQ2YAbsJacobian(
        y, target_mass, incident_energy));
}

inline double BjorkenXYDensityToQ2YDensity(
    double bjorken_density,
    double y,
    double target_mass,
    double incident_energy)
{
    return bjorken_density * Q2YToBjorkenXYAbsJacobian(
        y, target_mass, incident_energy);
}

inline double Q2YDensityToBjorkenXYDensity(
    double q2_y_density,
    double y,
    double target_mass,
    double incident_energy)
{
    return q2_y_density * BjorkenXYToQ2YAbsJacobian(
        y, target_mass, incident_energy);
}

// ------------------------------------------------------------------ //
//  SolidAngleRest <-> MandelstamQ2  (Scatter2to2 topology)            //
// ------------------------------------------------------------------ //
//
// For 2->2 scattering in the CM frame:
//   Q^2 = -t = 2 p_CM^2 (1 - cos theta_CM)
//   |dQ^2 / d(cos theta_CM)| = 2 p_CM^2
//
// To convert a density from dOmega_CM to dQ^2 (integrating uniform phi):
//   rho(Q^2) = rho(Omega_CM) * 2*pi / (2 * p_CM^2)
//
// To convert a density from dQ^2 to dOmega_CM:
//   rho(Omega_CM) = rho(Q^2) * (2 * p_CM^2) / (2*pi)
//
// The factor of 2*pi accounts for the azimuthal integration: dOmega
// is 2D (cos_theta, phi) while Q^2 is 1D. If both densities include
// the azimuth explicitly, drop the 2*pi factor and use
// SolidAngleRestToMandelstamQ2AbsJacobian directly.

inline double SolidAngleRestToMandelstamQ2AbsJacobian(
    double s,
    double m_beam,
    double m_target)
{
    double p_cm_sq = Kallen(s, m_beam * m_beam, m_target * m_target)
                     / (4.0 * s);
    if (p_cm_sq <= 0.0) return 0.0;
    return 2.0 * p_cm_sq;
}

inline double MandelstamQ2ToSolidAngleRestAbsJacobian(
    double s,
    double m_beam,
    double m_target)
{
    return SafeInverse(SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target));
}

// Density conversions (azimuth-integrated: include 2*pi factor)

inline double SolidAngleRestDensityToMandelstamQ2Density(
    double solid_angle_density,
    double s,
    double m_beam,
    double m_target)
{
    double jacobian = SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target);
    if (jacobian <= 0.0) return 0.0;
    return solid_angle_density * 2.0 * M_PI / jacobian;
}

inline double MandelstamQ2DensityToSolidAngleRestDensity(
    double q2_density,
    double s,
    double m_beam,
    double m_target)
{
    double jacobian = SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target);
    return q2_density * jacobian / (2.0 * M_PI);
}

} // namespace phase_space_jacobian
} // namespace injection
} // namespace siren

#endif // SIREN_PhaseSpaceJacobian_H
