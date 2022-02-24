
#include <cmath>
#include <iostream>

#include "earthmodel-service/Vector3D.h"

using namespace earthmodel;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Vector3D::Vector3D()
    : cartesian_(0,0,0)
    , spherical_(0,0,0)
{
}

Vector3D::Vector3D(const double x, const double y, const double z)
    : cartesian_(x, y, z)
    , spherical_(0,0,0)
{
}

// copy constructor
Vector3D::Vector3D(const Vector3D& vector_3d)
    : cartesian_(vector_3d.cartesian_)
    , spherical_(vector_3d.spherical_)
{
}

Vector3D::Vector3D(Vector3D&& other)
    : cartesian_(std::move(other.cartesian_))
    , spherical_(std::move(other.spherical_))
{
}

Vector3D::Vector3D(std::array<double, 3> const & vec)
    : cartesian_(vec[0], vec[1], vec[2])
    , spherical_(0,0,0)
{
}

Vector3D::Vector3D(geom3::UnitVector3 const & vec)
    : cartesian_(vec.x(), vec.y(), vec.z())
    , spherical_(0,0,0)
{
}

Vector3D::Vector3D(geom3::Vector3 const & vec)
    : cartesian_(vec.x(), vec.y(), vec.z())
    , spherical_(0,0,0)
{
}

Vector3D::Vector3D(geom3::Point3 const & vec)
    : cartesian_(vec.x(), vec.y(), vec.z())
    , spherical_(0,0,0)
{
}

// destructor
Vector3D::~Vector3D() {}

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//

Vector3D& Vector3D::operator=(Vector3D const & vector_3d) {
    if (this != &vector_3d)
    {
        Vector3D tmp(vector_3d);
        swap(tmp);
    }
    return *this;
}

Vector3D& Vector3D::operator=(Vector3D && other) {
    cartesian_ = std::move(other.cartesian_);
    spherical_ = std::move(other.spherical_);
    return *this;
}

Vector3D& Vector3D::operator=(Vector3D const && other) {
    cartesian_ = other.cartesian_;
    spherical_ = other.spherical_;
    return *this;
}

Vector3D::operator std::array<double, 3>() const {
    return std::array<double, 3>{cartesian_.x_, cartesian_.y_, cartesian_.z_};
}

Vector3D::operator geom3::UnitVector3() const {
    return geom3::UnitVector3(cartesian_.x_, cartesian_.y_, cartesian_.z_, true);
}

Vector3D::operator geom3::Vector3() const {
    return geom3::Vector3(cartesian_.x_, cartesian_.y_, cartesian_.z_);
}

Vector3D::operator geom3::Point3() const {
    return geom3::Point3(cartesian_.x_, cartesian_.y_, cartesian_.z_);
}

bool Vector3D::operator==(const Vector3D& vector_3d) const
{
    if (cartesian_.x_ != vector_3d.cartesian_.x_)
        return false;
    else if (cartesian_.y_ != vector_3d.cartesian_.y_)
        return false;
    else if (cartesian_.z_ != vector_3d.cartesian_.z_)
        return false;
    else if (spherical_.radius_ != vector_3d.spherical_.radius_)
        return false;
    else if (spherical_.azimuth_ != vector_3d.spherical_.azimuth_)
        return false;
    else if (spherical_.zenith_ != vector_3d.spherical_.zenith_)
        return false;

    return true;
}

bool Vector3D::operator<(const Vector3D& vector_3d) const {
    return (this != &vector_3d) and
        std::tie(cartesian_.x_, cartesian_.y_, cartesian_.z_, spherical_.radius_, spherical_.azimuth_, spherical_.zenith_)
        <
        std::tie(vector_3d.cartesian_.x_, vector_3d.cartesian_.y_, vector_3d.cartesian_.z_, vector_3d.spherical_.radius_, vector_3d.spherical_.azimuth_, vector_3d.spherical_.zenith_);
}

bool Vector3D::operator!=(const Vector3D& vector_3d) const
{
    return !(*this == vector_3d);
}

void Vector3D::swap(Vector3D& vector_3d)
{
    using std::swap;

    swap(cartesian_.x_, vector_3d.cartesian_.x_);
    swap(cartesian_.y_, vector_3d.cartesian_.y_);
    swap(cartesian_.z_, vector_3d.cartesian_.z_);
    swap(spherical_.radius_, vector_3d.spherical_.radius_);
    swap(spherical_.azimuth_, vector_3d.spherical_.azimuth_);
    swap(spherical_.zenith_, vector_3d.spherical_.zenith_);
}

namespace earthmodel {
std::ostream& operator<<(std::ostream& os, Vector3D const& vector_3d)
{
    std::stringstream ss;
    ss << "Vector3D (" << &vector_3d << ") ";
    os << ss.str() << '\n';

    os << "Cartesian Coordinates (x[cm],y[cm],z[cm]):\n"
       << vector_3d.cartesian_.x_ << "\t" << vector_3d.cartesian_.y_ << "\t" << vector_3d.cartesian_.z_ << std::endl;
    os << "Spherical Coordinates (radius[cm],azimuth[rad],zenith[rad]):\n"
       << vector_3d.spherical_.radius_ << "\t" << vector_3d.spherical_.azimuth_ << "\t" << vector_3d.spherical_.zenith_
       << std::endl;

    return os;
}
} // namespace earthmodel

//----------------------------------------------------------------------//
//-----------------------operator basic arithmetic ---------------------//
//----------------------------------------------------------------------//

Vector3D Vector3D::operator-() const
{
    Vector3D vector_3d;
    vector_3d.cartesian_.x_ = -cartesian_.x_;
    vector_3d.cartesian_.y_ = -cartesian_.y_;
    vector_3d.cartesian_.z_ = -cartesian_.z_;
    return vector_3d;
}

namespace earthmodel {

Vector3D operator+(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D vector_sum;
    vector_sum.cartesian_.x_ = vec1.cartesian_.x_ + vec2.cartesian_.x_;
    vector_sum.cartesian_.y_ = vec1.cartesian_.y_ + vec2.cartesian_.y_;
    vector_sum.cartesian_.z_ = vec1.cartesian_.z_ + vec2.cartesian_.z_;
    return vector_sum;
}

Vector3D& operator+=(Vector3D& vec1, const Vector3D& vec2)
{
    vec1.cartesian_.x_ += vec2.cartesian_.x_;
    vec1.cartesian_.y_ += vec2.cartesian_.y_;
    vec1.cartesian_.z_ += vec2.cartesian_.z_;
    return vec1;
}

Vector3D operator-(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D vector_diff;
    vector_diff.cartesian_.x_ = vec1.cartesian_.x_ - vec2.cartesian_.x_;
    vector_diff.cartesian_.y_ = vec1.cartesian_.y_ - vec2.cartesian_.y_;
    vector_diff.cartesian_.z_ = vec1.cartesian_.z_ - vec2.cartesian_.z_;
    return vector_diff;
}

Vector3D& operator-=(Vector3D& vec1, const Vector3D& vec2)
{
    vec1.cartesian_.x_ -= vec2.cartesian_.x_;
    vec1.cartesian_.y_ -= vec2.cartesian_.y_;
    vec1.cartesian_.z_ -= vec2.cartesian_.z_;
    return vec1;
}

double operator*(const Vector3D& vec1, const Vector3D& vec2)
{
    return scalar_product(vec1, vec2);
}

Vector3D operator*(const Vector3D& vec1, const double factor1)
{
    Vector3D product;
    product.cartesian_.x_ = factor1 * vec1.cartesian_.x_;
    product.cartesian_.y_ = factor1 * vec1.cartesian_.y_;
    product.cartesian_.z_ = factor1 * vec1.cartesian_.z_;
    return product;
}

Vector3D& operator*=(Vector3D& vec1, const double factor1)
{
    vec1.cartesian_.x_ *= factor1;
    vec1.cartesian_.y_ *= factor1;
    vec1.cartesian_.z_ *= factor1;
    return vec1;
}

Vector3D operator*(const double factor1, const Vector3D& vec1)
{
    Vector3D product;
    product.cartesian_.x_ = factor1 * vec1.cartesian_.x_;
    product.cartesian_.y_ = factor1 * vec1.cartesian_.y_;
    product.cartesian_.z_ = factor1 * vec1.cartesian_.z_;
    return product;
}

Vector3D operator/(const Vector3D& vec1, const double factor1)
{
    Vector3D product;
    product.cartesian_.x_ = vec1.cartesian_.x_ / factor1;
    product.cartesian_.y_ = vec1.cartesian_.y_ / factor1;
    product.cartesian_.z_ = vec1.cartesian_.z_ / factor1;
    return product;
}

Vector3D& operator/=(Vector3D& vec1, const double factor1)
{
    vec1.cartesian_.x_ /= factor1;
    vec1.cartesian_.y_ /= factor1;
    vec1.cartesian_.z_ /= factor1;
    return vec1;
}

double scalar_product(const Vector3D& vec1, const Vector3D& vec2)
{
    return vec1.cartesian_.x_ * vec2.cartesian_.x_ + vec1.cartesian_.y_ * vec2.cartesian_.y_ + vec1.cartesian_.z_ * vec2.cartesian_.z_;
}

Vector3D vector_product(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D product;
    product.cartesian_.x_ = vec1.cartesian_.y_ * vec2.cartesian_.z_ - vec1.cartesian_.z_ * vec2.cartesian_.y_;
    product.cartesian_.y_ = vec1.cartesian_.z_ * vec2.cartesian_.x_ - vec1.cartesian_.x_ * vec2.cartesian_.z_;
    product.cartesian_.z_ = vec1.cartesian_.x_ * vec2.cartesian_.y_ - vec1.cartesian_.y_ * vec2.cartesian_.x_;
    return product;
}

Vector3D cross_product(const Vector3D& vec1, const Vector3D& vec2)
{
    return vector_product(vec1, vec2);
}

} // namespace earthmodel

double Vector3D::magnitude() const
{
    return std::sqrt(cartesian_.x_ * cartesian_.x_ + cartesian_.y_ * cartesian_.y_ + cartesian_.z_ * cartesian_.z_);
}

void Vector3D::normalize()
{
    double length   = std::sqrt(cartesian_.x_ * cartesian_.x_ + cartesian_.y_ * cartesian_.y_ + cartesian_.z_ * cartesian_.z_);
    cartesian_.x_              = cartesian_.x_ / length;
    cartesian_.y_              = cartesian_.y_ / length;
    cartesian_.z_              = cartesian_.z_ / length;
    spherical_.radius_ = 1;
}

Vector3D Vector3D::normalized() const {
    Vector3D res(*this);
    res.normalize();
    return res;
}

void Vector3D::deflect(const double cosphi_deflect, const double theta_deflect)
{
    if(cosphi_deflect != 1 || theta_deflect != 0)
    {
        CalculateSphericalCoordinates();

        double sinphi_deflect = std::sqrt( std::max(0., (1. - cosphi_deflect) * (1. + cosphi_deflect) ));
        double tx = sinphi_deflect * std::cos(theta_deflect);
        double ty = sinphi_deflect * std::sin(theta_deflect);
        double tz = std::sqrt(std::max(1. - tx * tx - ty * ty, 0.));
        if(cosphi_deflect < 0. ){
            // Backward deflection
            tz = -tz;
        }

        double sinth, costh, sinph, cosph;
        sinth = std::sin(spherical_.zenith_);
        costh = std::cos(spherical_.zenith_);
        sinph = std::sin(spherical_.azimuth_);
        cosph = std::cos(spherical_.azimuth_);

        const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
        const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

        // Rotation towards all tree axes
        Vector3D new_direction( tz * *this + tx * rotate_vector_x + ty * rotate_vector_y );

        *this = new_direction;
    }
}

void Vector3D::invert() {
    cartesian_.x_ = -cartesian_.x_;
    cartesian_.y_ = -cartesian_.y_;
    cartesian_.z_ = -cartesian_.z_;
}

Vector3D Vector3D::inverted() const {
    Vector3D res(*this);
    res.invert();
    return res;
}

//----------------------------------------------------------------------//
//---------------Spherical and cylindrical coordinates------------------//
//----------------------------------------------------------------------//

void Vector3D::CalculateCartesianFromSpherical()
{
    cartesian_.x_ = spherical_.radius_ * std::cos(spherical_.azimuth_) * std::sin(spherical_.zenith_);
    cartesian_.y_ = spherical_.radius_ * std::sin(spherical_.azimuth_) * std::sin(spherical_.zenith_);
    cartesian_.z_ = spherical_.radius_ * std::cos(spherical_.zenith_);
}

void Vector3D::CalculateSphericalCoordinates()
{
    spherical_.radius_  = std::sqrt(cartesian_.x_ * cartesian_.x_ + cartesian_.y_ * cartesian_.y_ + cartesian_.z_ * cartesian_.z_);
    spherical_.azimuth_ = std::atan2(cartesian_.y_, cartesian_.x_);
    if (spherical_.radius_ > 0.)
    {
        spherical_.zenith_ = std::acos(cartesian_.z_ / spherical_.radius_);
    } else if (spherical_.radius_ == 0.)
    {
        // log_warn("If the radius is zero, the zenith is not defined! Zero is returned!");
        spherical_.zenith_ = 0.;
    } else
    {
        // log_fatal("The radius is negativ, which is not possible!");
    }
}

// void Vector3D::CalculateCartesianFromCylindrical()
// {
//     x_ = cylindric_radius_*std::cos(cylindric_azimuth_);
//     y_ = cylindric_radius_*std::sin(cylindric_azimuth_);
//     z_ = cylindric_height_;
// }

// void Vector3D::CalculateZylindricalCoordinates()
// {
//     cylindric_radius_  = std::sqrt(x_*x_ + y_*y_);
//     cylindric_azimuth_ = std::atan2(y_, x_);
//     cylindric_height_  = z_;
// }

Vector3D::CartesianCoordinates::CartesianCoordinates(double x, double y, double z):
    x_(x), y_(y), z_(z)
{}

Vector3D::CartesianCoordinates::CartesianCoordinates(const CartesianCoordinates& coordinates):
    x_(coordinates.x_), y_(coordinates.y_), z_(coordinates.z_)
{}

Vector3D::SphericalCoordinates::SphericalCoordinates(double radius, double azimuth, double zenith):
    radius_(radius), azimuth_(azimuth), zenith_(zenith)
{}

Vector3D::SphericalCoordinates::SphericalCoordinates(const SphericalCoordinates& coordinates):
    radius_(coordinates.radius_), azimuth_(coordinates.azimuth_), zenith_(coordinates.zenith_)
{}

//----------------------------------------------------------------------//
//---------------------private member function--------------------------//
//----------------------------------------------------------------------//
