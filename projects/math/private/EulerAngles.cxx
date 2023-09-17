#include "LeptonInjector/math/EulerAngles.h"

#include <utility>
#include <iostream>

#include "LeptonInjector/math/Vector3D.h"

using namespace LI::math;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
EulerAngles::EulerAngles() :
    order_(EulerOrder::ZXZr),
    alpha_(0.0),
    beta_(0.0),
    gamma_(0.0)
{
}

EulerAngles::EulerAngles(
    const EulerOrder order,
    const double alpha,
    const double beta,
    const double gamma
) :
    order_(order),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma)
{
}

// copy constructor
EulerAngles::EulerAngles(const EulerAngles& euler) :
    order_(euler.order_),
    alpha_(euler.alpha_),
    beta_(euler.beta_),
    gamma_(euler.gamma_)
{
}

EulerAngles::EulerAngles(EulerAngles&& other) :
    order_(std::move(other.order_)),
    alpha_(std::move(other.alpha_)),
    beta_(std::move(other.beta_)),
    gamma_(std::move(other.gamma_))
{
}

// destructor
EulerAngles::~EulerAngles() {}

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//

EulerAngles& EulerAngles::operator=(EulerAngles const & euler) {
    if (this != &euler)
    {
        EulerAngles tmp(euler);
        swap(tmp);
    }
    return *this;
}

EulerAngles& EulerAngles::operator=(EulerAngles && other) {
    order_ = std::move(other.order_);
    alpha_ = std::move(other.alpha_);
    beta_ = std::move(other.beta_);
    gamma_ = std::move(other.gamma_);
    return *this;
}

EulerAngles& EulerAngles::operator=(EulerAngles const && other) {
    order_ = other.order_;
    alpha_ = other.alpha_;
    beta_ = other.beta_;
    gamma_ = other.gamma_;
    return *this;
}

bool EulerAngles::operator==(const EulerAngles& euler) const
{
    return (this == &euler) or (
        order_ == euler.order_ and
        alpha_ == euler.alpha_ and
        beta_ == euler.beta_ and
        gamma_ == euler.gamma_);
}

bool EulerAngles::operator!=(const EulerAngles& euler) const
{
    return !(*this == euler);
}

void EulerAngles::swap(EulerAngles& euler)
{
    using std::swap;

    swap(order_, euler.order_);
    swap(alpha_, euler.alpha_);
    swap(beta_, euler.beta_);
    swap(gamma_, euler.gamma_);
}

namespace LI {
namespace math {
std::ostream& operator<<(std::ostream& os, EulerAngles const& euler)
{
    std::stringstream ss;
    ss << " EulerAngles (" << &euler << ") ";
    os << ss.str() << '\n';
    return os;
}

} // namespace math
} // namespace LI

