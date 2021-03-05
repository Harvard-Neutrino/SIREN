
#include <iostream>
#include <cmath>
#include <earthmodel-service/Vector3D.h>
#include <earthmodel-service/DensityDist.h>
#include <earthmodel-service/Polynomial.h>

using namespace earthmodel;

DensityDistribution::DensityDistribution() {}

DensityDistribution::DensityDistribution(const DensityDistribution& density_distr) {}

bool DensityDistribution::operator==(const DensityDistribution& dens_distr) const
{
    if (!this->compare(dens_distr) )
        return false;
    return true;
}


bool DensityDistribution::operator!=(const DensityDistribution& dens_distr) const {
    return !(*this == dens_distr);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%        Axis        %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Axis1D::Axis1D() {}

Axis1D::Axis1D(const Vector3D& fAxis, const Vector3D& fp0) : fAxis_(fAxis), fp0_(fp0) {}

Axis1D::Axis1D(const Axis1D& axis) : fAxis_(axis.fAxis_), fp0_(axis.fp0_) {}

bool Axis1D::operator==(const Axis1D& axis) const {
    if(fAxis_ != axis.fAxis_)
        return false;
    if(fp0_ != axis.fp0_)
        return false;
    return true;
}

bool Axis1D::operator!=(const Axis1D& axis) const {
    return !(*this == axis);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Radial       %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialAxis1D::RadialAxis1D() : Axis1D() {
    fp0_.SetCartesianCoordinates(0, 0, 0);
    fAxis_.SetSphericalCoordinates(1, 0, 0);
}

RadialAxis1D::RadialAxis1D(const Vector3D& fAxis, const Vector3D& fp0) : Axis1D(fAxis, fp0) {}

double RadialAxis1D::GetX(const Vector3D& xi) const {
    return (xi - fp0_).magnitude();
}

double RadialAxis1D::GetdX(const Vector3D& xi, const Vector3D& direction) const {
    Vector3D aux{xi - fp0_};
    aux.normalize();

    return aux * direction;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Cartesian     %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CartesianAxis1D::CartesianAxis1D() : Axis1D() {
    fAxis_.SetCartesianCoordinates(1, 0, 0);
    fp0_.SetCartesianCoordinates(0, 0, 0);
}

CartesianAxis1D::CartesianAxis1D(const Vector3D& fAxis, const Vector3D& fp0) : Axis1D(fAxis, fp0) {}

double CartesianAxis1D::GetX(const Vector3D& xi) const {
    return fAxis_ * (xi - fp0_);
}

double CartesianAxis1D::GetdX(const Vector3D& xi, const Vector3D& direction) const {
    (void)xi;

    return fAxis_ * direction;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Homogeneous-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

Density_homogeneous::Density_homogeneous()
    : DensityDistribution(), correction_factor_(1.0) {}

Density_homogeneous::Density_homogeneous(double correction_factor)
    : DensityDistribution(), correction_factor_(correction_factor) {}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& dens_distr)
    : DensityDistribution(dens_distr),
      correction_factor_(dens_distr.correction_factor_) {}


bool Density_homogeneous::compare(const DensityDistribution& dens_distr) const {
    const Density_homogeneous* dens_homogen = dynamic_cast<const Density_homogeneous*>(&dens_distr);
    if(!dens_homogen)
        return false;
    if(correction_factor_ != dens_homogen->correction_factor_ )
        return false;
    return true;
}

double Density_homogeneous::InverseIntegral(const Vector3D& xi,
                                            const Vector3D& direction,
                                            double res,
                                            double distance_to_border) const {
    (void)xi;
    (void)direction;
    (void)distance_to_border;

    return res / correction_factor_;
}

double Density_homogeneous::AntiDerivative(const Vector3D& xi,
                                           const Vector3D& direction) const {
    (void)xi;
    (void)direction;

    return correction_factor_;
}

double Density_homogeneous::Integral(const Vector3D& xi,
                                     const Vector3D& direction,
                                     double distance) const {
    return AntiDerivative(xi+direction*distance, direction) - AntiDerivative(xi, direction);
}

double Density_homogeneous::Evaluate(const Vector3D& xi) const {
    (void)xi;

    return correction_factor_;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PolynomialDistribution1D::PolynomialDistribution1D(const PolynomialDistribution1D& dist)
    : polynom_(dist.polynom_),
      Ipolynom_(dist.Ipolynom_),
      dpolynom_(dist.dpolynom_) {}

PolynomialDistribution1D::PolynomialDistribution1D(const Polynom& poly)
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

Density_polynomial::Density_polynomial(const Axis1D& axis, const Polynom& polynom)
    : polynom_(polynom),
      Polynom_(polynom_.GetAntiderivative(0)),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

Density_polynomial::Density_polynomial(const Density_polynomial& dens)
    : polynom_(dens.polynom_),
      Polynom_(dens.Polynom_),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

Density_polynomial::~Density_polynomial() {}

bool Density_polynomial::compare(const DensityDistribution& dens_distr) const {
    const Density_polynomial* dens_poly = dynamic_cast<const Density_polynomial*>(&dens_distr);
    if(!dens_poly)
        return false;
    if( polynom_ != dens_poly->polynom_ )
        return false;
    return true;
}

double Density_polynomial::Helper_function(const Vector3D& xi,
                                           const Vector3D& direction,
                                           double res,
                                           double l) const {
    return AntiDerivative(xi, direction) - AntiDerivative(xi + direction*l, direction) + res;
}

double Density_polynomial::helper_function(const Vector3D& xi,
                                           const Vector3D& direction,
                                           double res,
                                           double l) const {
    (void)res;

    return Evaluate(xi) - Evaluate(xi + l * direction);
}

double Density_polynomial::InverseIntegral(const Vector3D& xi,
                                           const Vector3D& direction,
                                           double res,
                                           double distance_to_border) const {
    std::function<double(double)> F =
        std::bind(&Density_polynomial::Helper_function, this, xi, direction,
                  res, std::placeholders::_1);

    std::function<double(double)> dF =
        std::bind(&Density_polynomial::helper_function, this, xi, direction,
                  res, std::placeholders::_1);

    // check if direction * axis larger or less than zero
    // direction * fAxis_

    try {
        res =
            NewtonRaphson(F, dF, 0, distance_to_border, distance_to_border / 2);
    } catch (MathException& e) {
        throw DensityException("Next interaction point lies in infinite.");
    }

    return res;
}

double Density_polynomial::AntiDerivative(const Vector3D& xi,
                                          const Vector3D& direction) const {
    return 0.0;
}

double Density_polynomial::Integral(const Vector3D& xi,
                                    const Vector3D& direction,
                                    double distance) const {
    return AntiDerivative(xi+direction*distance, direction) - AntiDerivative(xi, direction);
}

double Density_polynomial::Evaluate(const Vector3D& xi) const {
    return density_distribution(axis_->GetX(xi));
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Exponential-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

Density_exponential::Density_exponential(const Axis1D& axis, double sigma)
    : sigma_(sigma) {}

double Density_exponential::GetDepth(const Vector3D& xi) const {
    return axis_->GetX(xi) / sigma_;
}

bool Density_exponential::compare(const DensityDistribution& dens_distr) const {
    const Density_exponential* dens_exp = dynamic_cast<const Density_exponential*>(&dens_distr);
    if(!dens_exp)
        return false;
    if(sigma_ != dens_exp->sigma_ )
        return false;
    return true;
}

double Density_exponential::GetEffectiveDistance(const Vector3D& xi,
                                                 const Vector3D& direction) const {
    return axis_->GetdX(xi, direction) / sigma_;
}

double Density_exponential::InverseIntegral(const Vector3D& xi,
                                            const Vector3D& direction,
                                            double res,
                                            double distance_to_border) const {
    (void)distance_to_border;

    double phi = GetDepth(xi);
    double delta = GetEffectiveDistance(xi, direction);

    double aux = 1. / delta * std::log(1 + std::exp(-phi) * res * delta);

    if (std::isnan(aux))
        throw DensityException("Next interaction point lies in infinite.");

    return aux;
}

double Density_exponential::AntiDerivative(const Vector3D& xi,
                                           const Vector3D& direction) const {
    return 0.0;
}

double Density_exponential::Integral(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double distance) const {
    return AntiDerivative(xi+direction*distance, direction) - AntiDerivative(xi, direction);
}

double Density_exponential::Evaluate(const Vector3D& xi) const {
    return std::exp(GetDepth(xi));
}
