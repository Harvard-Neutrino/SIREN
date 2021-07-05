
#include <cmath>
#include <tuple>
#include <iostream>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Quaternion.h"
#include "earthmodel-service/Conversions.h"
#include "earthmodel-service/EulerAngles.h"

using namespace earthmodel;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Quaternion::Quaternion() :
    x_(0.0),
    y_(0.0),
    z_(0.0),
    w_(1.0)
{
}

Quaternion::Quaternion(
    double x,
    double y,
    double z,
    double w
) :
    x_(x),
    y_(y),
    z_(z),
    w_(w)
{
}

// copy constructor
Quaternion::Quaternion(const Quaternion& quaternion) :
    x_(quaternion.x_),
    y_(quaternion.y_),
    z_(quaternion.z_),
    w_(quaternion.w_)
{
}

Quaternion::Quaternion(Quaternion&& other) :
    x_(std::move(other.x_)),
    y_(std::move(other.y_)),
    z_(std::move(other.z_)),
    w_(std::move(other.w_))
{
}

// destructor
Quaternion::~Quaternion() {}

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//

Quaternion& Quaternion::operator=(Quaternion const & quaternion) {
    if (this != &quaternion)
    {
        Quaternion tmp(quaternion);
        swap(tmp);
    }
    return *this;
}

Quaternion& Quaternion::operator=(Quaternion && other) {
    x_ = std::move(other.x_);
    y_ = std::move(other.y_);
    z_ = std::move(other.z_);
    w_ = std::move(other.w_);
    return *this;
}

Quaternion& Quaternion::operator=(Quaternion const && other) {
    x_ = other.x_;
    y_ = other.y_;
    z_ = other.z_;
    w_ = other.w_;
    return *this;
}

bool Quaternion::operator==(const Quaternion& quaternion) const
{
    return (this == &quaternion) or (
        x_ == quaternion.x_ and
        y_ == quaternion.y_ and
        z_ == quaternion.z_ and
        w_ == quaternion.w_);
}

bool Quaternion::operator!=(const Quaternion& quaternion) const
{
    return !(*this == quaternion);
}

void Quaternion::swap(Quaternion& quaternion)
{
    using std::swap;

    swap(x_, quaternion.x_);
    swap(y_, quaternion.y_);
    swap(z_, quaternion.z_);
    swap(w_, quaternion.w_);
}

namespace earthmodel {
std::ostream& operator<<(std::ostream& os, Quaternion const& quaternion)
{
    std::stringstream ss;
    ss << " Quaternion (" << &quaternion << ") ";
    os << ss.str() << '\n';
    return os;
}
} // namespace earthmodel

Quaternion Quaternion::operator*(Quaternion const & other) const
{
    Quaternion product;
    product.x_ = (other.w_ * x_) + (other.x_ * w_) + (other.y_ * z_) - (other.z_ * y_);
    product.y_ = (other.w_ * y_) + (other.y_ * w_) + (other.z_ * x_) - (other.x_ * z_);
    product.z_ = (other.w_ * z_) + (other.z_ * w_) + (other.x_ * y_) - (other.y_ * x_);
    product.w_ = (other.w_ * w_) - (other.x_ * x_) - (other.y_ * y_) - (other.z_ * z_);
    return product;
}

Quaternion & Quaternion::operator*=(Quaternion const & other)
{
    return (*this = other * (*this));
}

Quaternion Quaternion::operator*(double factor) const
{
    return Quaternion(factor * x_, factor * y_, factor * z_, factor * w_);
}

Quaternion & Quaternion::operator*=(double factor)
{
    x_ *= factor;
    y_ *= factor;
    z_ *= factor;
    w_ *= factor;
}

Quaternion Quaternion::operator+(Quaternion const & other) const
{
    Quaternion sum;
    sum.x_ = other.x_ + x_;
    sum.y_ = other.y_ + y_;
    sum.z_ = other.z_ + z_;
    sum.w_ = other.w_ + w_;
    return sum;
}

Quaternion & Quaternion::operator+=(Quaternion const & other)
{
    return (*this = other * (*this));
}

Quaternion Quaternion::operator+(double factor) const
{
    return Quaternion(factor * x_, factor * y_, factor * z_, factor * w_);
}

Quaternion & Quaternion::operator+=(double factor)
{
    x_ *= factor;
    y_ *= factor;
    z_ *= factor;
    w_ *= factor;
}

void Quaternion::GetMatrix(Matrix3D & dest) const
{
    dest.SetXX(1 - 2 * y_ * y_ - 2 * z_ * z_);
    dest.SetXY(2 * x_ * y_ + 2 * z_ * w_);
    dest.SetXZ(2 * x_ * z_ - 2 * y_ * w_);
    dest.SetYX(2 * x_ * y_ - 2 * z_ * w_);
    dest.SetYY(1 - 2 * x_ * x_ - 2 * z_ * z_);
    dest.SetYZ(2 * z_ * y_ + 2 * x_ * w_);
    dest.SetZX(2 * x_ * z_ + 2 * y_ * w_);
    dest.SetZY(2 * z_ * y_ - 2 * x_ * w_);
    dest.SetZZ(1 - 2 * x_ * x_ - 2 * y_ * y_);
}

Matrix3D Quaternion::GetMatrix() const
{
    Matrix3D mat;
    GetMatrix(mat);
    return mat;
}

Quaternion & Quaternion::invert() {
    x_ = -x_;
    y_ = -y_;
    z_ = -z_;
    return *this;
}

void Quaternion::SetPosition(Vector3D const & vec)
{
    w_ = 0.0;
    x_ = vec.GetX();
    y_ = vec.GetY();
    z_ = vec.GetZ();
}

Quaternion & Quaternion::normalize()
{
    double norm = x_ * x_ + y_ * y_ + z_ * z_ + w_ * w_;

    if(norm == 1)
        return *this;

    return (*this *= 1.0 / sqrt(norm));
}

double Quaternion::DotProduct(Quaternion const & qu) const
{
    return x_ * qu.x_ + y_ * qu.y_ + z_ * qu.z_ + w_ * qu.w_;
}

Quaternion Quaternion::lerp(Quaternion const & q1, Quaternion const & q2, double t)
{
    const double s = 1.0 - t;
    return (q1 * s) + (q2 * t);
}

Quaternion Quaternion::slerp(Quaternion const & q1, Quaternion const & q2, double t)
{
    double alpha = q1.DotProduct(q2);
    unsigned int sign = 1;

    if(alpha < 0) {
        //q1 *= -1;
        sign = -1;
        alpha *= -1;
    }

    double theta = acos(alpha);
    double istheta = 1.0 / sin(theta);
    double S = sin(theta * (1.0 - t)) * istheta * sign;
    double T = sin(theta * t) * istheta;
    return (q1 * S) + (q2 * T);
}

void Quaternion::SetAxisAngle(Vector3D const & axis, double angle)
{
    double x = angle / 2.0;
    double s = sin(x);

    x_ = s * axis.GetX();
    y_ = s * axis.GetY();
    z_ = s * axis.GetZ();
    w_ = cos(x);
}

void Quaternion::GetAxisAngle(Vector3D & axis, double & angle) const
{
    double scale = sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    if(scale == 0 or w_ > 1.0 or w_ < -1.0) {
        angle = 0;
        axis.SetCartesianCoordinates(0, 0, 1);
    } else {
        axis.SetCartesianCoordinates(x_ / scale, y_ / scale, z_ / scale);
    }
}

std::tuple<Vector3D, double> Quaternion::GetAxisAngle() const
{
    std::tuple<Vector3D, double> result;
    GetAxisAngle(std::get<0>(result), std::get<1>(result));
    return result;
}

void Quaternion::SetAngles(EulerAngles const & euler)
{
    (*this) = earthmodel::QuaternionFromEulerAngles(euler);
}

void Quaternion::GetAnglesEulerZXZ(double & alpha, double & beta, double & gamma) const
{
    double s2 = x_ * x_ + y_ * y_; // sin(beta)^2
    double c2 = w_ * w_ + z_ * z_; // cos(beta)^2
    double s = atan(z_ / w_); // (gamma+alpha)/2
    double d = atan2(y_, x_); // (gamma-alpha)/2
    alpha = s - d;
    gamma = s + d;

    if (c2 != 0.0)
        beta = 2.0 * atan(sqrt(s2 / c2));
    else
        beta = (0.5 > s2) ? 0 : M_PI;
}

void Quaternion::SetAnglesEulerZXZ(double alpha, double beta, double gamma)
{
    double theta = beta/2;

    const double sb = sin(theta);
    const double cb = cos(theta);

    theta = (gamma + alpha)/2;
    const double ss = sin(theta);
    const double cs = cos(theta);

    theta = (gamma - alpha)/2;
    const double sd = sin(theta);
    const double cd = cos(theta);

    x_ = sb * cd;
    y_ = sb * sd;
    z_ = cb * ss;
    w_ = cb * cs;

    normalize();
}

void Quaternion::GetAnglesTaitBryanZYX(double & yaw, double & pitch, double & roll) const
{
    double t0 = x_ * x_ - z_ * z_;
    double t1 = w_ * w_ - y_ * y_;
    double xx = 0.5 * (t0 + t1);
    double xy = x_ * y_ + w_ * z_;
    double xz = w_ * y_ - x_ * z_;
    double yz = 2.0 * (y_ * z_ + w_ * x_);
    double t  = xx * xx + xy * xy;

    yaw = atan2(xy, xx);
    pitch = atan(xz / sqrt(t));

    if (t != 0) {
        roll = atan2(yz, t1 - t0);
    } else {
        roll = 2.0 * atan2(x_, w_) - std::copysign(1.0, xz) * yaw;
    }
}

void Quaternion::SetAnglesTaitBryanZYX(double yaw, double pitch, double roll)
{
    
}

void Quaternion::SetAnglesTaitBryanZXY(double yaw, double pitch, double roll)
{
    double theta = yaw/2;
    const double sa = sin(theta);
    const double ca = cos(theta);
    theta = pitch/2;
    const double sb = sin(theta);
    const double cb = cos(theta);
    theta = roll/2;
    const double sg = sin(theta);
    const double cg = cos(theta);

    x_ = ca * sb * cg - sa * cb * sg;
    y_ = sa * cb * cg + ca * sb * sg;
    z_ = sa * sb * cg + ca * cb * sg;
    w_ = ca * cb * cg - sa * sb * sg;

    normalize();
}
