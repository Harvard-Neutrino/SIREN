#include "SIREN/geometry/Geometry.h"

#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <typeinfo>
#include <typeindex>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Placement.h"

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/


namespace SI {
namespace geometry {

std::ostream& operator<<(std::ostream& os, SI::geometry::Geometry const& geometry)
{
    os << "Geometry(" << &geometry << ")" << std::endl;
    os << geometry.placement_ << std::endl;
    geometry.print(os);
    return os;
}

/******************************************************************************
 *                                  Geometry                                   *
 ******************************************************************************/

Geometry::Geometry()
{
}

Geometry::Geometry(const std::string name)
    : name_(name)
{
}

Geometry::Geometry(const std::string name, Placement const & placement)
    : name_(name)
      , placement_(placement)
{
}

Geometry::Geometry(Placement const & placement)
    : name_()
      , placement_(placement)
{
}

Geometry::Geometry(const Geometry& geometry)
    : name_(geometry.name_)
      , placement_(geometry.placement_)
{
}

// ------------------------------------------------------------------------- //
void Geometry::swap(Geometry& geometry)
{
    name_.swap(geometry.name_);
    placement_.swap(geometry.placement_);
}



// ------------------------------------------------------------------------- //
Geometry& Geometry::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        name_     = geometry.name_;
        placement_     = geometry.placement_;
    }

    return *this;
}

// ------------------------------------------------------------------------- //
bool Geometry::operator==(const Geometry& geometry) const
{
    if (name_.compare(geometry.name_) != 0 or (placement_ != geometry.placement_))
        return false;
    else
        return this->equal(geometry);
}

// ------------------------------------------------------------------------- //
bool Geometry::operator<(const Geometry& geometry) const
{
    if(typeid(this) == typeid(&geometry)) {
        if(name_ != geometry.name_)
            return name_ < geometry.name_;
        else if(placement_ != geometry.placement_)
            return placement_ < geometry.placement_;
        else
            return this->less(geometry);
    } else
        return std::type_index(typeid(this)) < std::type_index(typeid(&geometry));
}

// ------------------------------------------------------------------------- //
bool Geometry::operator!=(const Geometry& geometry) const
{
    return !(*this == geometry);
}

// ------------------------------------------------------------------------- //
// Member functions
// ------------------------------------------------------------------------- //

bool Geometry::IsInside(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const
{
    bool is_inside = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first > 0 && dist.second < 0)
    {
        is_inside = true;
    }
    return is_inside;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsInfront(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const
{
    bool is_infront = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first > 0 && dist.second > 0)
    {
        is_infront = true;
    }
    return is_infront;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsBehind(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const
{
    bool is_behind = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first < 0 && dist.second < 0)
    {
        is_behind = true;
    }
    return is_behind;
}

Geometry::ParticleLocation::Enum Geometry::GetLocation(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const {
    if(IsInfront(position, direction))
        return Geometry::ParticleLocation::InfrontGeometry;
    if(IsInside(position, direction))
        return Geometry::ParticleLocation::InsideGeometry;
    else
        return Geometry::ParticleLocation::BehindGeometry;
}

// ------------------------------------------------------------------------- //
double Geometry::DistanceToClosestApproach(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const
{
    SI::math::Vector3D pos = GlobalToLocalPosition(position);
    SI::math::Vector3D dir = GlobalToLocalPosition(direction);
    return scalar_product(-pos, dir);
}

SI::math::Vector3D Geometry::LocalToGlobalPosition(SI::math::Vector3D const & p0) const
{
    return placement_.LocalToGlobalPosition(p0);
}

SI::math::Vector3D Geometry::LocalToGlobalDirection(SI::math::Vector3D const & p0) const
{
    return placement_.LocalToGlobalDirection(p0);
}

SI::math::Vector3D Geometry::GlobalToLocalPosition(SI::math::Vector3D const & p0) const
{
    return placement_.GlobalToLocalPosition(p0);
}

SI::math::Vector3D Geometry::GlobalToLocalDirection(SI::math::Vector3D const & p0) const
{
    return placement_.GlobalToLocalDirection(p0);
}

std::pair<double, double> Geometry::DistanceToBorder(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const {
    SI::math::Vector3D local_position = GlobalToLocalPosition(position);
    SI::math::Vector3D local_direction = GlobalToLocalDirection(direction);
    return ComputeDistanceToBorder(position, direction);
}
std::vector<Geometry::Intersection> Geometry::Intersections(SI::math::Vector3D const & position, SI::math::Vector3D const & direction) const {
    SI::math::Vector3D local_position = GlobalToLocalPosition(position);
    SI::math::Vector3D local_direction = GlobalToLocalDirection(direction);
    std::vector<Geometry::Intersection> intersections = ComputeIntersections(local_position, local_direction);
    for(auto & intersection : intersections) {
        intersection.position = LocalToGlobalPosition(intersection.position);
    }
    return intersections;
}

} // namespace geometry
} // namespace SI
