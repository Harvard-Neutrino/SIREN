
#include <iostream>
#include <cmath>
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/DensityDist.h"
#include "earthmodel-service/Polynomial.h"

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
    if (!this->compare(axis) )
        return false;
    return true;
}

bool Axis1D::operator!=(const Axis1D& axis) const {
    return !(*this == axis);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%    Distribution    %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool Distribution1D::operator==(const Distribution1D& dist) const {
    if (!this->compare(dist) )
        return false;
    return true;
}

bool Distribution1D::operator!=(const Distribution1D& dist) const {
    return !(*this == dist);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Radial       %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialAxis1D::RadialAxis1D() : Axis1D() {
    fp0_.SetCartesianCoordinates(0, 0, 0);
    fAxis_.SetSphericalCoordinates(1, 0, 0);
}

bool RadialAxis1D::compare(const Axis1D& ax) const {
    const RadialAxis1D* r_ax = dynamic_cast<const RadialAxis1D*>(&ax);
    if(!r_ax)
        return false;
    if(fp0_ != r_ax->fp0_)
        return false;
    return true;
}

RadialAxis1D::RadialAxis1D(const Vector3D& fAxis, const Vector3D& fp0) : Axis1D(fAxis, fp0) {}

RadialAxis1D::RadialAxis1D(const Vector3D& fp0) : Axis1D(Vector3D(), fp0) {}

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

bool CartesianAxis1D::compare(const Axis1D& ax) const {
    const CartesianAxis1D* c_ax = dynamic_cast<const CartesianAxis1D*>(&ax);
    if(!c_ax)
        return false;
    if(fp0_ != c_ax->fp0_)
        return false;
    if(fAxis_ != c_ax->fAxis_)
        return false;
    return true;
}

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

ConstantDistribution1D::ConstantDistribution1D()
    : val_(1e-25) {}

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

