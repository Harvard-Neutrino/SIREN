#include "SIREN/geometry/Geometry.h"

#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <typeinfo>
#include <typeindex>
#include <algorithm>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Placement.h"

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/


namespace siren {
namespace geometry {

std::ostream& operator<<(std::ostream& os, siren::geometry::Geometry const& geometry)
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

bool Geometry::IsInside(const siren::math::Vector3D& position, const siren::math::Vector3D& direction) const
{
    siren::math::Vector3D local_pos = GlobalToLocalPosition(position);
    siren::math::Vector3D local_dir = GlobalToLocalDirection(direction);
    std::vector<Intersection> hits = ComputeIntersections(local_pos, local_dir);

    // The first forward intersection determines containment:
    // if it is an exit (entering == false), the point is inside.
    for(auto const & hit : hits) {
        if(hit.distance > GEOMETRY_PRECISION) {
            return !hit.entering;
        }
    }
    return false;
}

bool Geometry::IsInside(const siren::math::Vector3D& position) const
{
    return IsInside(position, siren::math::Vector3D(0, 0, 1));
}

// ------------------------------------------------------------------------- //
double Geometry::DistanceToClosestApproach(const siren::math::Vector3D& position, const siren::math::Vector3D& direction) const
{
    siren::math::Vector3D pos = GlobalToLocalPosition(position);
    siren::math::Vector3D dir = GlobalToLocalPosition(direction);
    return scalar_product(-pos, dir);
}

siren::math::Vector3D Geometry::LocalToGlobalPosition(siren::math::Vector3D const & p0) const
{
    return placement_.LocalToGlobalPosition(p0);
}

siren::math::Vector3D Geometry::LocalToGlobalDirection(siren::math::Vector3D const & p0) const
{
    return placement_.LocalToGlobalDirection(p0);
}

siren::math::Vector3D Geometry::GlobalToLocalPosition(siren::math::Vector3D const & p0) const
{
    return placement_.GlobalToLocalPosition(p0);
}

siren::math::Vector3D Geometry::GlobalToLocalDirection(siren::math::Vector3D const & p0) const
{
    return placement_.GlobalToLocalDirection(p0);
}

AABB Geometry::GetWorldBoundingBox() const {
    AABB local_box = GetBoundingBox();
    // Generate the 8 corners of the local AABB
    double x0 = local_box.min_corner.GetX();
    double y0 = local_box.min_corner.GetY();
    double z0 = local_box.min_corner.GetZ();
    double x1 = local_box.max_corner.GetX();
    double y1 = local_box.max_corner.GetY();
    double z1 = local_box.max_corner.GetZ();

    AABB world_box;
    // Transform each corner to global coordinates and expand the world AABB
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x0, y0, z0)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x0, y0, z1)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x0, y1, z0)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x0, y1, z1)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x1, y0, z0)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x1, y0, z1)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x1, y1, z0)));
    world_box.ExpandToInclude(LocalToGlobalPosition(siren::math::Vector3D(x1, y1, z1)));

    return world_box;
}

std::vector<Geometry::Intersection> Geometry::Intersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    siren::math::Vector3D local_position = GlobalToLocalPosition(position);
    siren::math::Vector3D local_direction = GlobalToLocalDirection(direction);
    std::vector<Geometry::Intersection> intersections = ComputeIntersections(local_position, local_direction);
    for(auto & intersection : intersections) {
        intersection.position = LocalToGlobalPosition(intersection.position);
    }
    return intersections;
}

std::vector<Geometry::Intersection> Geometry::Intersections(LocalPosition const & position, LocalDirection const & direction) const {
    std::vector<Geometry::Intersection> intersections = ComputeIntersections(position.get(), direction.get());
    for(auto & intersection : intersections) {
        intersection.position = LocalToGlobalPosition(intersection.position);
    }
    return intersections;
}

std::pair<LocalPosition, LocalDirection> Geometry::ToLocal(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    return {LocalPosition(GlobalToLocalPosition(position)), LocalDirection(GlobalToLocalDirection(direction))};
}

} // namespace geometry
} // namespace siren

CEREAL_REGISTER_DYNAMIC_INIT(siren_Geometry);

