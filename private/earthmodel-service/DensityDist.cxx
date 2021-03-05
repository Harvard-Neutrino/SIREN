
#include <iostream>
#include <cmath>
#include <earthmodel-service/Vector3D.h>
#include <earthmodel-service/DensityDist.h>
#include <earthmodel-service/Polynomial.h>

using namespace earthmodel;

DensityDistribution::DensityDistribution() : axis_(CartesianAxis().clone()) {}

DensityDistribution::DensityDistribution(const DensityDistribution& density_distr)
    : axis_(density_distr.axis_->clone()) {}

DensityDistribution::DensityDistribution(const Axis& axis) : axis_(axis.clone()) {}

bool DensityDistribution::operator==(const DensityDistribution& dens_distr) const
{
    if (*axis_ != *dens_distr.axis_)
        return false;
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

Axis::Axis() {}

Axis::Axis(const Vector3D& fAxis, const Vector3D& fp0) : fAxis_(fAxis), fp0_(fp0) {}

Axis::Axis(const Axis& axis) : fAxis_(axis.fAxis_), fp0_(axis.fp0_) {}

bool Axis::operator==(const Axis& axis) const {
    if(fAxis_ != axis.fAxis_)
        return false;
    if(fp0_ != axis.fp0_)
        return false;
    return true;
}

bool Axis::operator!=(const Axis& axis) const {
    return !(*this == axis);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Radial       %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialAxis::RadialAxis() : Axis() {
    fp0_.SetCartesianCoordinates(0, 0, 0);
    fAxis_.SetSphericalCoordinates(1, 0, 0);
}

RadialAxis::RadialAxis(const Vector3D& fAxis, const Vector3D& fp0) : Axis(fAxis, fp0) {}

double RadialAxis::GetDepth(const Vector3D& xi) const {
    return (xi - fp0_).magnitude();
}

double RadialAxis::GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const {
    Vector3D aux{xi - fp0_};
    aux.normalise();

    return -aux * direction;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Cartesian     %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CartesianAxis::CartesianAxis() : Axis() {
    fAxis_.SetCartesianCoordinates(1, 0, 0);
    fp0_.SetCartesianCoordinates(0, 0, 0);
}

CartesianAxis::CartesianAxis(const Vector3D& fAxis, const Vector3D& fp0) : Axis(fAxis, fp0) {}

double CartesianAxis::GetDepth(const Vector3D& xi) const {
    return fAxis_ * (xi - fp0_);
}

double CartesianAxis::GetEffectiveDistance(const Vector3D& xi,
                                           const Vector3D& direction) const {
    (void)xi;

    return fAxis_ * direction;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Homogeneous-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

double Density_homogeneous::Correct(const Vector3D& xi,
                                    const Vector3D& direction,
                                    double res,
                                    double distance_to_border) const {
    (void)xi;
    (void)direction;
    (void)distance_to_border;

    return res / correction_factor_;
}

double Density_homogeneous::Integrate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double l) const {
    (void)xi;
    (void)direction;

    return correction_factor_ * l;
}

double Density_homogeneous::Calculate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_homogeneous::Evaluate(const Vector3D& xi) const {
    (void)xi;

    return correction_factor_;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_polynomial::Density_polynomial(const Axis& axis, const Polynom& polynom)
    : DensityDistribution(axis),
      polynom_(polynom),
      Polynom_(polynom_.GetAntiderivative(0)),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

Density_polynomial::Density_polynomial(const Density_polynomial& dens)
    : DensityDistribution(dens),
      polynom_(dens.polynom_),
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
    return Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res;
}

double Density_polynomial::helper_function(const Vector3D& xi,
                                           const Vector3D& direction,
                                           double res,
                                           double l) const {
    (void)res;

    return Evaluate(xi) - Evaluate(xi + l * direction);
}

double Density_polynomial::Correct(const Vector3D& xi,
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

double Density_polynomial::Integrate(const Vector3D& xi,
                                     const Vector3D& direction,
                                     double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return antiderived_density_distribution(axis_->GetDepth(xi) + l * delta) /
           (delta * delta);
}

double Density_polynomial::Calculate(const Vector3D& xi,
                                     const Vector3D& direction,
                                     double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_polynomial::Evaluate(const Vector3D& xi) const {
    return density_distribution(axis_->GetDepth(xi));
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%% Exponential-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_exponential::Density_exponential(const Axis& axis, double sigma)
    : DensityDistribution(axis), sigma_(sigma) {}

double Density_exponential::GetDepth(const Vector3D& xi) const {
    return axis_->GetDepth(xi) / sigma_;
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
    return axis_->GetEffectiveDistance(xi, direction) / sigma_;
}

double Density_exponential::Correct(const Vector3D& xi,
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

double Density_exponential::Integrate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double l) const {
    double delta = GetEffectiveDistance(xi, direction);

    return std::exp(GetDepth(xi) + l * delta) / delta;
}

double Density_exponential::Calculate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_exponential::Evaluate(const Vector3D& xi) const {
    return std::exp(GetDepth(xi));
}
