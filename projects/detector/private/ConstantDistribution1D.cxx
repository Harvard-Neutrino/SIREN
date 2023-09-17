#include "LeptonInjector/detector/ConstantDistribution1D.h"

#include <vector>

namespace LI {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Homogeneous-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConstantDistribution1D::ConstantDistribution1D()
    : val_(1e-25) {} // universe_mean_density from GEANT4

ConstantDistribution1D::ConstantDistribution1D(const ConstantDistribution1D& dist)
    : val_(dist.val_) {}

ConstantDistribution1D::ConstantDistribution1D(double val)
    : val_(val) {}

bool ConstantDistribution1D::compare(const Distribution1D& dist) const {
    const ConstantDistribution1D* dist_const = dynamic_cast<const ConstantDistribution1D*>(&dist);
    if(!dist_const)
        return false;
    if(val_ != dist_const->val_)
        return false;
    return true;
}

double ConstantDistribution1D::Derivative(double x) const {
    (void)x;
    return 0.0;
}

double ConstantDistribution1D::AntiDerivative(double x) const {
    return x*val_;
}

double ConstantDistribution1D::Evaluate(double x) const {
    (void)x;
    return val_;
}

} // namespace LI
} // namespace detector
