
#include <cmath>
#include <tuple>
#include <iostream>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Quaternion.h"
#include "earthmodel-service/Conversions.h"
#include "earthmodel-service/EulerAngles.h"
#include "earthmodel-service/EulerQuaternionConversions.h"

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

Quaternion Quaternion::compose(Quaternion const & p) const
{
    double w,x,y,z;

    double const & w0 = w_;
    double const & x0 = x_;
    double const & y0 = y_;
    double const & z0 = z_;

    double w1 = p.GetW();
    double x1 = p.GetX();
    double y1 = p.GetY();
    double z1 = p.GetZ();

    double wx = w_ * x1, wy = w_ * y1, wz = w_ * z1;
    double xx = x_ * x1, xy = x_ * y1, xz = x_ * z1;
    double yx = y_ * x1, yy = y_ * y1, yz = y_ * z1;
    double zx = z_ * x1, zy = z_ * y1, zz = z_ * z1;

    double w2 = w0 * w0;
    double x2 = x0 * x0;
    double y2 = y0 * y0;
    double z2 = z0 * z0;

    w = w1 * (w2 + x2 + y2 + z2);
    x = x1 * (w2 + x2 - y2 - z2) + 2 * (x0 * (yy + zz) + w0 * (yz - zy));
    y = y1 * (w2 - x2 + y2 - z2) + 2 * (y0 * (xx + zz) + w0 * (zx - xz));
    z = z1 * (w2 - x2 - y2 + x2) + 2 * (z0 * (xx + yy) + w0 * (xy - yx));
    return Quaternion(x, y, z, w);
}

Vector3D Quaternion::compose(Vector3D const & p) const
{
    double w,x,y,z;

    double const & w0 = w_;
    double const & x0 = x_;
    double const & y0 = y_;
    double const & z0 = z_;

    double x1 = p.GetX();
    double y1 = p.GetY();
    double z1 = p.GetZ();

    double wx = w_ * x1, wy = w_ * y1, wz = w_ * z1;
    double xx = x_ * x1, xy = x_ * y1, xz = x_ * z1;
    double yx = y_ * x1, yy = y_ * y1, yz = y_ * z1;
    double zx = z_ * x1, zy = z_ * y1, zz = z_ * z1;

    double w2 = w0 * w0;
    double x2 = x0 * x0;
    double y2 = y0 * y0;
    double z2 = z0 * z0;

    x = x1 * (w2 + x2 - y2 - z2) + 2 * (x0 * (yy + zz) + w0 * (yz - zy));
    y = y1 * (w2 - x2 + y2 - z2) + 2 * (y0 * (xx + zz) + w0 * (zx - xz));
    z = z1 * (w2 - x2 - y2 + x2) + 2 * (z0 * (xx + yy) + w0 * (xy - yx));
    return Vector3D(x, y, z);
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

double Quaternion::magnitude() const
{
    return w_ * w_ + x_ * x_ + y_ * y_ + z_ * z_;
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

void Quaternion::SetEulerAngles(EulerAngles const & euler)
{
    (*this) = earthmodel::QuaternionFromEulerAngles(euler);
}

void Quaternion::GetEulerAngles(EulerAngles & euler, EulerOrder order) const
{
    euler = earthmodel::EulerAnglesFromQuaternion(*this, order);
}

void Quaternion::GetEulerAnglesZXZs(double & alpha, double & beta, double & gamma) const
{
    EulerAngles euler = earthmodel::ZXZsFromQ(*this);
    alpha = euler.GetAlpha();
    beta = euler.GetBeta();
    gamma = euler.GetGamma();
}

void Quaternion::SetEulerAnglesZXZs(double alpha, double beta, double gamma)
{
    (*this) = earthmodel::QFromZXZs(alpha, beta, gamma);
}

void Quaternion::GetEulerAnglesXYZs(double & alpha, double & beta, double & gamma) const
{
    EulerAngles euler = earthmodel::XYZsFromQ(*this);
    alpha = euler.GetAlpha();
    beta = euler.GetBeta();
    gamma = euler.GetGamma();
}

void Quaternion::SetEulerAnglesXYZs(double alpha, double beta, double gamma)
{
    (*this) = earthmodel::QFromXYZs(alpha, beta, gamma);
}
