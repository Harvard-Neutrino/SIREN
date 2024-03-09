#include "SIREN/detector/ExponentialDistribution1D.h"

#include <cmath>  // for exp

namespace siren {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Exponential-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExponentialDistribution1D::ExponentialDistribution1D()
    : sigma_(1.0) {}

ExponentialDistribution1D::ExponentialDistribution1D(const ExponentialDistribution1D& dist)
    : sigma_(dist.sigma_) {}

ExponentialDistribution1D::ExponentialDistribution1D(double sigma)
    : sigma_(sigma) {}

bool ExponentialDistribution1D::compare(const Distribution1D& dist) const {
    const ExponentialDistribution1D* dist_exp = dynamic_cast<const ExponentialDistribution1D*>(&dist);
    if(!dist_exp)
        return false;
    if(sigma_ != dist_exp->sigma_)
        return false;
    return true;
}

double ExponentialDistribution1D::Derivative(double x) const {
    return Evaluate(x) * sigma_;
}

double ExponentialDistribution1D::AntiDerivative(double x) const {
    return Evaluate(x) / sigma_;
}

double ExponentialDistribution1D::Evaluate(double x) const {
    return exp(sigma_*x);
}

} // namespace siren
} // namespace detector
