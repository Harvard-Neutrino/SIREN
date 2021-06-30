
#include <cmath>
#include <iostream>

#include "earthmodel-service/Quaternion.h"
#include "earthmodel-service/Vector3D.h"

using namespace earthmodel;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Quaternion::Quaternion() :
    x_(0.0),
    y_(0.0),
    z_(0.0),
    w_(1.0),
{
}

Quaternion::Quaternion(
    double x,
    double y,
    double z,
    double w,
) :
    x_(x),
    y_(y),
    z_(z),
    w_(w),
{
}

// copy constructor
Quaternion::Quaternion(const Quaternion& quaternion) :
    x_(quaternion.x_),
    y_(quaternion.y_),
    z_(quaternion.z_),
    w_(quaternion.w_),
{
}

Quaternion::Quaternion(Quaternion&& other) :
    x_(std::move(other.x_)),
    y_(std::move(other.y_)),
    z_(std::move(other.z_)),
    w_(std::move(other.w_)),
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

