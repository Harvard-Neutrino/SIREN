#include "SIREN/detector/ExponentialDistribution1D.h"

#include <cmath>  // for exp

namespace siren {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Exponential-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExponentialDistribution1D::ExponentialDistribution1D()
    : sigma_(1.0), amplitude_(1.0), x0_(0.0) {}

ExponentialDistribution1D::ExponentialDistribution1D(const ExponentialDistribution1D& dist)
    : sigma_(dist.sigma_), amplitude_(dist.amplitude_), x0_(dist.x0_) {}

ExponentialDistribution1D::ExponentialDistribution1D(double sigma)
    : sigma_(sigma), amplitude_(1.0), x0_(0.0) {}

ExponentialDistribution1D::ExponentialDistribution1D(double sigma, double amplitude, double x0)
    : sigma_(sigma), amplitude_(amplitude), x0_(x0) {}

bool ExponentialDistribution1D::compare(const Distribution1D& dist) const {
    const ExponentialDistribution1D* dist_exp = dynamic_cast<const ExponentialDistribution1D*>(&dist);
    if(!dist_exp)
        return false;
    if(sigma_ != dist_exp->sigma_)
        return false;
    if(amplitude_ != dist_exp->amplitude_)
        return false;
    if(x0_ != dist_exp->x0_)
        return false;
    return true;
}

// With rho(x) = A*exp(sigma*(x-x0)), rho'(x) = rho(x)*sigma and the
// antiderivative is rho(x)/sigma, so Derivative and AntiDerivative are
// Evaluate(x)*sigma and Evaluate(x)/sigma.
double ExponentialDistribution1D::Derivative(double x) const {
    return Evaluate(x) * sigma_;
}

double ExponentialDistribution1D::AntiDerivative(double x) const {
    return Evaluate(x) / sigma_;
}

double ExponentialDistribution1D::Evaluate(double x) const {
    return amplitude_ * exp(sigma_*(x - x0_));
}

} // namespace siren
} // namespace detector

CEREAL_REGISTER_DYNAMIC_INIT(siren_ExponentialDistribution1D);

