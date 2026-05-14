#include "SIREN/geometry/Sphere.h"

#include <cmath>
#include <tuple>
#include <math.h>
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

Sphere::Sphere()
    : Geometry((std::string)("Sphere"))
    , radius_(0.0)
      , inner_radius_(0.0)
{
    // Do nothing here
}

Sphere::Sphere(double radius, double inner_radius)
    : Geometry("Sphere")
    , radius_(radius)
      , inner_radius_(inner_radius)
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

Sphere::Sphere(Placement const & placement)
    : Geometry((std::string)("Sphere"), placement)
    , radius_(0.0)
      , inner_radius_(0.0)
{
    // Do nothing here
}

Sphere::Sphere(Placement const & placement, double radius, double inner_radius)
    : Geometry((std::string)("Sphere"), placement)
    , radius_(radius)
      , inner_radius_(inner_radius)
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

Sphere::Sphere(const Sphere& sphere)
    : Geometry(sphere)
    , radius_(sphere.radius_)
      , inner_radius_(sphere.inner_radius_)
{
    // Nothing to do here
}



/*Sphere::Sphere(const nlohmann::json& config)
  : Geometry(config)
  {
  if(not config.is_object()) throw std::invalid_argument("No json object found.");
  if(not config.at("outer_radius").is_number())
  throw std::invalid_argument("Outer radius is not a number.");

  config["outer_radius"].get_to(radius_);
  inner_radius_ = config.value("inner_radius", 0);

  if(inner_radius_ < 0) throw std::logic_error("inner radius must be >= 0");
  if(radius_ < inner_radius_)
  throw std::logic_error("radius must be larger than inner radius");
  }*/

// ------------------------------------------------------------------------- //
void Sphere::swap(Geometry& geometry)
{
    Sphere* sphere = dynamic_cast<Sphere*>(&geometry);
    if (!sphere)
    {
        //log_warn("Cannot swap Sphere!");
        return;
    }

    Geometry::swap(*sphere);

    std::swap(inner_radius_, sphere->inner_radius_);
    std::swap(radius_, sphere->radius_);
}

//------------------------------------------------------------------------- //
Sphere& Sphere::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);
        if (!sphere)
        {
            //log_warn("Cannot assign Sphere!");
            return *this;
        }

        Sphere tmp(*sphere);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Sphere::equal(const Geometry& geometry) const
{
    const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);

    if (!sphere)
        return false;
    else if (inner_radius_ != sphere->inner_radius_)
        return false;
    else if (radius_ != sphere->radius_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Sphere::less(const Geometry& geometry) const
{
    const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);

    return
        std::tie(inner_radius_, radius_)
        <
        std::tie(sphere->inner_radius_, sphere->radius_);
}

// ------------------------------------------------------------------------- //
void Sphere::print(std::ostream& os) const
{
    os << "Radius: " << radius_ << "\tInner radius: " << inner_radius_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Sphere::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Ray-sphere: C*t^2 + 2*B*t + A = 0  (C=1 for unit direction)
    double difference_length_squared = position.GetX()*position.GetX()
                                     + position.GetY()*position.GetY()
                                     + position.GetZ()*position.GetZ();
    double B = scalar_product(position, direction);

    // At most 4 intersections (2 outer + 2 inner for hollow sphere)
    Intersection hits[4];
    int n_hits = 0;

    auto add_hit = [&](double t, bool entering) {
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        hits[n_hits].distance = t;
        hits[n_hits].hierarchy = 0;
        hits[n_hits].entering = entering;
        hits[n_hits].position = position + t * direction;
        n_hits++;
    };

    // Outer sphere
    double A = difference_length_squared - radius_ * radius_;
    double determinant = B * B - A;
    if(determinant > 0) {
        double sq = std::sqrt(determinant);
        double t1 = -B - sq;
        double t2 = -B + sq;
        add_hit(t1, true);
        add_hit(t2, false);

        // Inner sphere (hollow shell)
        if(inner_radius_ > 0) {
            A = difference_length_squared - inner_radius_ * inner_radius_;
            determinant = B * B - A;
            if(determinant > 0) {
                sq = std::sqrt(determinant);
                t1 = -B - sq;
                t2 = -B + sq;
                add_hit(t1, false);
                add_hit(t2, true);
            }
        }
    }

    if(n_hits == 0) return {};
    // Sort by distance
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    return {hits, hits + n_hits};
}

// ------------------------------------------------------------------------- //
AABB Sphere::GetBoundingBox() const {
    return AABB(
        math::Vector3D(-radius_, -radius_, -radius_),
        math::Vector3D( radius_,  radius_,  radius_)
    );
}

} // namespace geometry
} // namespace siren

CEREAL_REGISTER_DYNAMIC_INIT(siren_Sphere);

