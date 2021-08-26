#include <cmath>
#include <tuple>
#include <iostream>

#include "earthmodel-service/Placement.h"
#include "earthmodel-service/EulerQuaternionConversions.h"

using namespace earthmodel;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

Placement::Placement() :
	position_(Vector3D(0,0,0)),
	quaternion_(QFromZXZs(0,0,0))
{
}

Placement::Placement(Vector3D position, Quaternion quaternion) :
	position_(position),
	quaternion_(quaternion)
{
}

Placement::Placement(const Placement& placement) :
	position_(placement.position_),
	quaternion(placement.quaternion_)
{
}

Placement::Placement(Placement&& other) :
    position_(std::move(other.position_)),
    quaternion_(std::move(other.quaternion_))
{
}

//destructor
Placement::~Placement() {};

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//


Placement& Placement::operator=(Placement const & placement) {
		if (this != &placement)
		{
			Placement tmp(placement);
			swap(tmp);
		}
		return *this;
}

Placement& Placement::operator=(Placement && other) {
    position_ = std::move(other.position_);
    quaternion_ = std::move(other.quaternion_);
    return *this;
}

Placement& Placement::operator=(Placement const && other) {
    position_ = std::move(other.position_);
    quaternion_ = std::move(other.quaternion_);
    return *this;
}


bool Placement::operator==(const Placement& placement) const
{
    return (this == &placement) or (
        position_ = placement.position_ and
        quaternion_ = placement.quaternion_);
}

bool Placement::operator!=(const Placement& placement) const
{
    return !(*this == placement);
}

void Placement::swap(Placement& placement)
{
    using std::swap;

    swap(position_, placement.position_);
    swap(quaternion_, placement.quaternion_);
}

friend std::ostream& operator<<(std::ostream& os, Placement const& placement);

std::shared_ptr<const Placement> create() const { return std::shared_ptr<const Placement>( new Placement(*this) ); }

//-------------------------------------//
// getter and setter functions
Vector3D Placement::GetPosition() { return position_;};
Quaternion Placement::GetQuaternion() { return quaternion_;};

Vector3D Placement::compose(Vector3D const & p0, bool inv) const
{
	return quaternion_.compose(p0,inv);
}

private:
    Vector3D position_;
    Quaternion quaternion_;
};

