#include "SIREN/geometry/Box.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <math.h>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <functional>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"
#include "GeometryMacros.h"

namespace siren {
namespace geometry {

Box::Box() : Geometry("Box"), full_width_x_(0), full_width_y_(0), full_width_z_(0) { RecomputeWorldAABB(); }
Box::Box(double x, double y, double z) : Geometry("Box"), full_width_x_(x), full_width_y_(y), full_width_z_(z) { RecomputeWorldAABB(); }
Box::Box(Placement const & p) : Geometry("Box", p), full_width_x_(0), full_width_y_(0), full_width_z_(0) { RecomputeWorldAABB(); }
Box::Box(Placement const & p, double full_width_x, double full_width_y, double full_width_z) : Geometry("Box", p), full_width_x_(full_width_x), full_width_y_(full_width_y), full_width_z_(full_width_z) { RecomputeWorldAABB(); }
Box::Box(const Box& o) : Geometry(o), full_width_x_(o.full_width_x_), full_width_y_(o.full_width_y_), full_width_z_(o.full_width_z_) { RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(Box, full_width_x_, full_width_y_, full_width_z_)
SIREN_GEOMETRY_ASSIGN(Box)
SIREN_GEOMETRY_EQUAL(Box, full_width_x_, full_width_y_, full_width_z_)
SIREN_GEOMETRY_LESS(Box, full_width_x_, full_width_y_, full_width_z_)

void Box::print(std::ostream& os) const {
    os << "Box(" << full_width_x_ << ", " << full_width_y_ << ", " << full_width_z_ << ")\n";
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Box::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    double px = position.GetX(), py = position.GetY(), pz = position.GetZ();
    double dx = direction.GetX(), dy = direction.GetY(), dz = direction.GetZ();

    double hx = 0.5 * full_width_x_, hy = 0.5 * full_width_y_, hz = 0.5 * full_width_z_;

    // Use the slab method: find the ray parameter ranges where the ray is
    // inside each pair of parallel planes, then intersect the three ranges.
    double t_enter, t_exit;

    // Initialize with the full line
    t_enter = -std::numeric_limits<double>::infinity();
    t_exit  =  std::numeric_limits<double>::infinity();

    // X slab
    if(dx != 0) {
        double inv = 1.0 / dx;
        double t1 = (-hx - px) * inv; // -x face
        double t2 = ( hx - px) * inv; // +x face
        if(t1 > t2) { std::swap(t1, t2); }
        if(t1 > t_enter) { t_enter = t1; }
        if(t2 < t_exit)  { t_exit  = t2; }
    } else {
        if(px < -hx || px > hx) return {};
    }

    // Y slab
    if(dy != 0) {
        double inv = 1.0 / dy;
        double t1 = (-hy - py) * inv;
        double t2 = ( hy - py) * inv;
        if(t1 > t2) { std::swap(t1, t2); }
        if(t1 > t_enter) { t_enter = t1; }
        if(t2 < t_exit)  { t_exit  = t2; }
    } else {
        if(py < -hy || py > hy) return {};
    }

    // Z slab
    if(dz != 0) {
        double inv = 1.0 / dz;
        double t1 = (-hz - pz) * inv;
        double t2 = ( hz - pz) * inv;
        if(t1 > t2) { std::swap(t1, t2); }
        if(t1 > t_enter) { t_enter = t1; }
        if(t2 < t_exit)  { t_exit  = t2; }
    } else {
        if(pz < -hz || pz > hz) return {};
    }

    // No intersection if the entry point is past the exit point
    if(t_enter > t_exit) return {};

    // A particle on the border moving inside has one intersection (exit),
    // a particle on the border moving outside has no intersection.
    if(t_enter > 0 && t_enter < GEOMETRY_PRECISION) t_enter = 0;
    if(t_exit > 0 && t_exit < GEOMETRY_PRECISION) t_exit = 0;

    Intersection hit_enter, hit_exit;
    hit_enter.distance = t_enter;
    hit_enter.hierarchy = 0;
    hit_enter.entering = true;
    hit_enter.position = siren::math::Vector3D(px + t_enter * dx, py + t_enter * dy, pz + t_enter * dz);

    hit_exit.distance = t_exit;
    hit_exit.hierarchy = 0;
    hit_exit.entering = false;
    hit_exit.position = siren::math::Vector3D(px + t_exit * dx, py + t_exit * dy, pz + t_exit * dz);

    return {hit_enter, hit_exit};
}

// ------------------------------------------------------------------------- //
AABB Box::GetBoundingBox() const {
    double hx = full_width_x_ * 0.5;
    double hy = full_width_y_ * 0.5;
    double hz = full_width_z_ * 0.5;
    return AABB(
        math::Vector3D(-hx, -hy, -hz),
        math::Vector3D( hx,  hy,  hz)
    );
}

} // namespace geometry
} // namespace siren
