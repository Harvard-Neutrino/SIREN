#include "SIREN/geometry/Cylinder.h"

#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

Cylinder::Cylinder()
    : Geometry((std::string)("Cylinder"))
    , radius_(0.0)
    , inner_radius_(0.0)
      , z_(0.0)
{
    // Do nothing here
}

Cylinder::Cylinder(double radius, double inner_radius, double z)
    : Geometry((std::string)("Cylinder"))
    , radius_(radius)
    , inner_radius_(inner_radius)
      , z_(z)
{
    if (inner_radius_ > radius_)
    {
        //log_error("Inner radius %f is greater then radius %f (will be swaped)", inner_radius_, radius_);
        std::swap(inner_radius_, radius_);
    }
    if (inner_radius_ == radius_)
    {
        //log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
    }
}

Cylinder::Cylinder(Placement const & placement)
    : Geometry((std::string)("Cylinder"), placement)
    , radius_(0.0)
    , inner_radius_(0.0)
      , z_(0.0)
{
    // Do nothing here
}

Cylinder::Cylinder(Placement const & placement, double radius, double inner_radius, double z)
    : Geometry((std::string)("Cylinder"), placement)
    , radius_(radius)
    , inner_radius_(inner_radius)
      , z_(z)
{
    if (inner_radius_ > radius_)
    {
        //log_error("Inner radius %f is greater then radius %f (will be swaped)", inner_radius_, radius_);
        std::swap(inner_radius_, radius_);
    }
    if (inner_radius_ == radius_)
    {
        //log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
    }
}

Cylinder::Cylinder(const Cylinder& cylinder)
    : Geometry(cylinder)
    , radius_(cylinder.radius_)
    , inner_radius_(cylinder.inner_radius_)
      , z_(cylinder.z_)
{
    // Nothing to do here
}

/*Cylinder::Cylinder(const nlohmann::json& config)
  : Geometry(config)
  {
  assert(config["outer_radius"].is_number());
  assert(config["height"].is_number());


  radius_ = config["outer_radius"].get<double>(); // m
  inner_radius_ = config.value("inner_radius", 0); // m
  z_ = config["height"].get<double>(); // m

  assert(inner_radius_>=0);
  assert(radius_>inner_radius_);
  assert(z_>0);
  }*/

// ------------------------------------------------------------------------- //
void Cylinder::swap(Geometry& geometry)
{
    Cylinder* cylinder = dynamic_cast<Cylinder*>(&geometry);
    if (!cylinder)
    {
        //log_warn("Cannot swap Cylinder!");
        return;
    }

    Geometry::swap(*cylinder);

    std::swap(inner_radius_, cylinder->inner_radius_);
    std::swap(radius_, cylinder->radius_);
    std::swap(z_, cylinder->z_);
}

//------------------------------------------------------------------------- //
Cylinder& Cylinder::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);
        if (!cylinder)
        {
            //log_warn("Cannot assign Cylinder!");
            return *this;
        }
        Cylinder tmp(*cylinder);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Cylinder::equal(const Geometry& geometry) const
{
    const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);

    if (!cylinder)
        return false;
    else if (inner_radius_ != cylinder->inner_radius_)
        return false;
    else if (radius_ != cylinder->radius_)
        return false;
    else if (z_ != cylinder->z_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Cylinder::less(const Geometry& geometry) const
{
    const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);
    if(!cylinder) return false;

    return
        std::tie(inner_radius_, radius_, z_)
        <
        std::tie(cylinder->inner_radius_, cylinder->radius_, cylinder->z_);
}

void Cylinder::print(std::ostream& os) const
{
    os << "Radius: " << radius_ << "\tInnner radius: " << inner_radius_ << " Height: " << z_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Cylinder::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Calculate intersection of particle trajectory and the cylinder
    // cylinder barrel (x1 + x0)^2 + (x2 + y0)^2  = radius^2 [ z0_-0.5*z_ <
    // particle->z <z0_ - 0.5*z_ ]
    // top/bottom surface:
    // E1: x3   =   z0_ + 0.5*z
    // E2: x3   =   z0_ - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    // (-1/-1) cylinder is behind particle or particle is on border but moving
    // outside
    // ( dist_1 / dist_2 ) cylinder is infront of the particle
    // ( dist_1 / -1 ) particle is inside the cylinder or on border and moving
    // inside

    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();
    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();

    double hz = 0.5 * z_;
    double r2_outer = radius_ * radius_;
    double r2_inner = inner_radius_ * inner_radius_;

    // At most 8 intersections (2 barrel + 2 caps + 2 inner barrel + 2 inner caps)
    // but realistically at most 4 for a single hollow cylinder
    Intersection hits[8];
    int n_hits = 0;

    auto add_hit = [&](double t, double ix, double iy, double iz, bool entering) {
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        hits[n_hits].distance = t;
        hits[n_hits].hierarchy = 0;
        hits[n_hits].entering = entering;
        hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
        n_hits++;
    };

    // Test barrel surface for a given radius squared; invert_entering flips
    // the radial entering logic for inner surfaces
    auto test_barrel = [&](double r2, bool invert_entering) {
        if(dx == 0 && dy == 0) return;
        double C = dx*dx + dy*dy;
        double B_half = px*dx + py*dy;
        double A = px*px + py*py - r2;
        double det = B_half*B_half - C*A;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double inv_C = 1.0 / C;
        double t1 = (-B_half - sq) * inv_C;
        double t2 = (-B_half + sq) * inv_C;
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? t1 : t2;
            double iz = pz + t * dz;
            if(iz > -hz && iz < hz) {
                double ix = px + t * dx;
                double iy = py + t * dy;
                bool radial_entering = (ix * dx + iy * dy) < 0;
                add_hit(t, ix, iy, iz, invert_entering ? !radial_entering : radial_entering);
            }
        }
    };

    // Test endcap at z_cap for a given annular ring [r2_lo, r2_hi]
    auto test_cap = [&](double z_cap, double r2_lo, double r2_hi, bool enter) {
        if(dz == 0) return;
        double t = (z_cap - pz) / dz;
        double ix = px + t * dx;
        double iy = py + t * dy;
        double r2_hit = ix*ix + iy*iy;
        if(r2_hit <= r2_hi && r2_hit >= r2_lo) {
            add_hit(t, ix, iy, pz + t * dz, enter);
        }
    };

    // Outer barrel
    test_barrel(r2_outer, false);
    // Endcaps (annular disk between inner and outer radius)
    test_cap( hz, r2_inner, r2_outer, dz < 0);
    test_cap(-hz, r2_inner, r2_outer, dz > 0);
    // Inner barrel (if hollow)
    if(inner_radius_ > 0) {
        test_barrel(r2_inner, true);
    }

    if(n_hits == 0) return {};
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    return {hits, hits + n_hits};
}

// ------------------------------------------------------------------------- //
AABB Cylinder::GetBoundingBox() const {
    double hz = z_ * 0.5;
    return AABB(
        math::Vector3D(-radius_, -radius_, -hz),
        math::Vector3D( radius_,  radius_,  hz)
    );
}

} // namespace geometry
} // namespace siren
