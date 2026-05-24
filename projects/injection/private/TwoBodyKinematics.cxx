#include "SIREN/injection/TwoBodyKinematics.h"

#include <cmath>
#include <algorithm>

namespace siren {
namespace injection {

std::array<TwoBodyLabSolution, 2> SolveLabAngle(
    double beta,
    double gamma,
    double p_rest,
    double E_rest,
    double m_daughter,
    double cos_theta_lab)
{
    std::array<TwoBodyLabSolution, 2> solutions;
    solutions[0].valid = false;
    solutions[1].valid = false;

    double beta2 = beta * beta;

    // Denominator: 1 - beta^2 * cos^2(theta_lab)
    double denom = 1.0 - beta2 * cos_theta_lab * cos_theta_lab;
    if (denom <= 0) return solutions;

    // Discriminant (from quadratic in p_lab, MeVPrtl Eq. 5):
    // D = E_rest^2 / gamma^2 - m^2 * (1 - beta^2 * cos^2(theta_lab))
    double m2 = m_daughter * m_daughter;
    double D = E_rest * E_rest / (gamma * gamma) - m2 * denom;
    if (D < 0) return solutions;

    double sqrtD = std::sqrt(D);
    double E_rest_beta_cos = E_rest * beta * cos_theta_lab / gamma;

    // Two solutions for lab-frame momentum (MeVPrtl Eq. 5):
    // p_pm = [E_rest * beta * cos(theta_lab) / gamma +/- sqrt(D)] / denom
    for (int sign = 0; sign < 2; ++sign) {
        double pm = (sign == 0) ? 1.0 : -1.0;
        double p_lab = (E_rest_beta_cos + pm * sqrtD) / denom;

        if (p_lab <= 0) continue;

        // Lab energy
        double E_lab = std::sqrt(p_lab * p_lab + m2);

        // Rest-frame angle from the boost relation:
        // p_rest * cos(theta_rest) = gamma * (p_lab * cos(theta_lab) - beta * E_lab)
        double cos_theta_rest = gamma * (p_lab * cos_theta_lab - beta * E_lab) / p_rest;

        if (cos_theta_rest < -1.0 - 1e-10 || cos_theta_rest > 1.0 + 1e-10) continue;
        cos_theta_rest = std::clamp(cos_theta_rest, -1.0, 1.0);

        // Jacobian |dOmega_lab / dOmega_rest|.
        //
        // MeVPrtl Eq. 6 gives dOmega_rest/dOmega_lab (note their primed
        // = rest frame convention).  We store the INVERSE since we need
        // the density transformation from rest to lab:
        //
        //   dOmega_rest/dOmega_lab = p_lab^2 / (gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab|)
        //   dOmega_lab/dOmega_rest = gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab| / p_lab^2

        double dp_factor = std::abs(p_lab - beta * E_lab * cos_theta_lab);
        double jacobian = 0.0;
        if (p_lab > 0) {
            jacobian = gamma * p_rest * dp_factor / (p_lab * p_lab);
        }

        solutions[sign].cos_theta_rest = cos_theta_rest;
        solutions[sign].p_lab = p_lab;
        solutions[sign].jacobian = jacobian;
        solutions[sign].valid = true;
    }

    return solutions;
}

double CriticalCosTheta(
    double beta,
    double gamma,
    double p_rest,
    double E_rest,
    double m_daughter)
{
    // The critical angle exists when the daughter's rest-frame velocity
    // is less than the parent's velocity (beta_daughter_rest < beta_parent).
    // beta_daughter_rest = p_rest / E_rest
    if (p_rest >= E_rest * beta) {
        // No critical angle: all lab directions are accessible
        return -1.0;
    }

    // cos(theta_critical) = sqrt(m^2 - E_rest^2/gamma^2) / (m * beta)
    double m2 = m_daughter * m_daughter;
    double arg = m2 - E_rest * E_rest / (gamma * gamma);
    if (arg < 0) return -1.0;

    return std::sqrt(arg) / (m_daughter * beta);
}

LabFrameResult BoostToLab(
    double beta,
    double gamma,
    double p_rest,
    double E_rest,
    double cos_theta_rest)
{
    LabFrameResult result;

    // Parallel and perpendicular components in rest frame
    double p_parallel = p_rest * cos_theta_rest;
    double sin_theta_rest = std::sqrt(std::max(0.0, 1.0 - cos_theta_rest * cos_theta_rest));
    double p_perp = p_rest * sin_theta_rest;

    // Boost: p_lab_parallel = gamma * (p_parallel + beta * E_rest)
    double p_lab_parallel = gamma * (p_parallel + beta * E_rest);
    double p_lab_perp = p_perp;

    result.p_lab = std::sqrt(p_lab_parallel * p_lab_parallel + p_lab_perp * p_lab_perp);
    result.E_lab = gamma * (E_rest + beta * p_parallel);

    if (result.p_lab > 0) {
        result.cos_theta_lab = p_lab_parallel / result.p_lab;
    } else {
        result.cos_theta_lab = 1.0;
    }

    return result;
}

} // namespace injection
} // namespace siren
