#include "SIREN/detector/PolynomialDistribution1D.h"

#include <vector>

#include "SIREN/math/Polynomial.h"

namespace siren {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PolynomialDistribution1D::PolynomialDistribution1D()
    : polynom_(std::vector<double>{}),
      Ipolynom_(polynom_.GetAntiderivative(0)),
      dpolynom_(polynom_.GetDerivative()) {}

PolynomialDistribution1D::PolynomialDistribution1D(const PolynomialDistribution1D& dist)
    : polynom_(dist.polynom_),
      Ipolynom_(dist.Ipolynom_),
      dpolynom_(dist.dpolynom_) {}

PolynomialDistribution1D::PolynomialDistribution1D(const math::Polynom& poly)
    : polynom_(poly),
      Ipolynom_(poly.GetAntiderivative(0)),
      dpolynom_(poly.GetDerivative()) {}

PolynomialDistribution1D::PolynomialDistribution1D(const std::vector<double>& coef)
    : polynom_(coef),
      Ipolynom_(polynom_.GetAntiderivative(0)),
      dpolynom_(polynom_.GetDerivative()) {}

bool PolynomialDistribution1D::compare(const Distribution1D& dist) const {
    const PolynomialDistribution1D* dist_poly = dynamic_cast<const PolynomialDistribution1D*>(&dist);
    if(!dist_poly)
        return false;
    if(polynom_ != dist_poly->polynom_)
        return false;
    return true;
}

double PolynomialDistribution1D::Derivative(double x) const {
    return dpolynom_.evaluate(x);
}

double PolynomialDistribution1D::AntiDerivative(double x) const {
    return Ipolynom_.evaluate(x);
}

double PolynomialDistribution1D::Evaluate(double x) const {
    return polynom_.evaluate(x);
}

} // namespace siren
} // namespace detector

CEREAL_REGISTER_DYNAMIC_INIT(siren_PolynomialDistribution1D);

