#include "SIREN/injection/TwoBodyKinematics.h"

#include <cmath>
#include <algorithm>
#include <limits>

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
    std::array<TwoBodyLabSolution, 2> solutions{};

    // This routine is called after several inverse-trigonometric and boost
    // operations. Reject invalid inputs explicitly: comparisons against NaN
    // are false, so the range checks below cannot safely detect NaNs on their
    // own.
    if (!std::isfinite(beta) || !std::isfinite(gamma) ||
        !std::isfinite(p_rest) || !std::isfinite(E_rest) ||
        !std::isfinite(m_daughter) || !std::isfinite(cos_theta_lab) ||
        beta < 0.0 || beta >= 1.0 || !(gamma > 0.0) ||
        !(p_rest > 0.0) || E_rest < 0.0 || m_daughter < 0.0) {
        return solutions;
    }

    // Dot products of normalized directions can overshoot by a few ulps.
    // Accept that harmless roundoff, but reject genuinely invalid angles.
    constexpr double cos_tolerance = 1e-12;
    if (cos_theta_lab < -1.0 - cos_tolerance ||
        cos_theta_lab > 1.0 + cos_tolerance) {
        return solutions;
    }
    cos_theta_lab = std::clamp(cos_theta_lab, -1.0, 1.0);

    double beta2 = beta * beta;

    // Denominator: 1 - beta^2 * cos^2(theta_lab)
    double denom = 1.0 - beta2 * cos_theta_lab * cos_theta_lab;
    if (!std::isfinite(denom) || denom <= 0.0) return solutions;

    // Discriminant (from quadratic in p_lab):
    // D = E_rest^2 / gamma^2 - m^2 * (1 - beta^2 * cos^2(theta_lab))
    double m2 = m_daughter * m_daughter;
    double energy_term = E_rest * E_rest / (gamma * gamma);
    double mass_term = m2 * denom;
    double D = energy_term - mass_term;
    if (!std::isfinite(D)) return solutions;
    if (D < 0.0) {
        double scale = std::max(std::abs(energy_term), std::abs(mass_term));
        double roundoff = 32.0 * std::numeric_limits<double>::epsilon() * scale;
        if (D < -roundoff) return solutions;
        D = 0.0;
    }

    double sqrtD = std::sqrt(D);
    double E_rest_beta_cos = E_rest * beta * cos_theta_lab / gamma;

    // Two solutions for lab-frame momentum:
    // p_pm = [E_rest * beta * cos(theta_lab) / gamma +/- sqrt(D)] / denom
    for (int sign = 0; sign < 2; ++sign) {
        double pm = (sign == 0) ? 1.0 : -1.0;
        double p_lab = (E_rest_beta_cos + pm * sqrtD) / denom;

        if (!std::isfinite(p_lab) || p_lab <= 0.0) continue;

        // Lab energy
        double E_lab = std::sqrt(p_lab * p_lab + m2);
        if (!std::isfinite(E_lab)) continue;

        // Rest-frame angle from the boost relation:
        // p_rest * cos(theta_rest) = gamma * (p_lab * cos(theta_lab) - beta * E_lab)
        double cos_theta_rest = gamma * (p_lab * cos_theta_lab - beta * E_lab) / p_rest;

        if (!std::isfinite(cos_theta_rest)) continue;
        if (cos_theta_rest < -1.0 - 1e-10 || cos_theta_rest > 1.0 + 1e-10) continue;
        cos_theta_rest = std::clamp(cos_theta_rest, -1.0, 1.0);

        // We store the inverse of |dOmega_rest / dOmega_lab| since we need the density transformation from rest to lab:
        //   dOmega_rest/dOmega_lab = p_lab^2 / (gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab|)
        //   dOmega_lab/dOmega_rest = gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab| / p_lab^2

        double dp_factor = std::abs(p_lab - beta * E_lab * cos_theta_lab);
        double jacobian = 0.0;
        if (p_lab > 0) {
            jacobian = gamma * p_rest * dp_factor / (p_lab * p_lab);
        }
        if (!std::isfinite(jacobian)) continue;

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
    if (!std::isfinite(beta) || !std::isfinite(gamma) ||
        !std::isfinite(p_rest) || !std::isfinite(E_rest) ||
        !std::isfinite(m_daughter) || !(gamma > 0.0) ||
        p_rest < 0.0 || E_rest < 0.0 || m_daughter < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

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

    double denominator = m_daughter * beta;
    if (!(denominator > 0.0) || !std::isfinite(denominator)) return -1.0;

    // The exact result is a cosine. Near threshold, p_rest is obtained from
    // the cancellation-prone Kallen function while E_rest uses a different
    // formula, so roundoff can put this ratio a few ulps above one.
    double cos_critical = std::sqrt(arg) / denominator;
    if (!std::isfinite(cos_critical)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::clamp(cos_critical, -1.0, 1.0);
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
