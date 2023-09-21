#include "LeptonInjector/geometry/Sphere.h"

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

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/Placement.h"

namespace LI {
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
std::vector<Geometry::Intersection> Sphere::ComputeIntersections(LI::math::Vector3D const & position, LI::math::Vector3D const & direction) const {
    // Calculate intersection of particle trajectory and the sphere
    // sphere (x1 + x0)^2 + (x2 + y0)^2 + (x3 + z0)^2 = radius^2
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // length of direction vector =1 => C = 1
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    double A, B, t1, t2, difference_length_squared;

    double determinant;

    std::vector<Intersection> dist;

    LI::math::Vector3D intersection;

    std::function<void(double, bool)> save = [&](double t, bool entering){
        Intersection i;
        i.position = intersection;
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };

    difference_length_squared = std::pow((position).magnitude(), 2);
    A                         = difference_length_squared - radius_ * radius_;

    B = scalar_product(position, direction);

    determinant = B * B - A;

    if (determinant > 0) // determinant == 0 (boundery point) is ignored
    {
        t1 = -1 * B + std::sqrt(determinant);
        t2 = -1 * B - std::sqrt(determinant);

        // Computer precision controll
        if (t1 > 0 && t1 < GEOMETRY_PRECISION)
            t1 = 0;
        if (t2 > 0 && t2 < GEOMETRY_PRECISION)
            t2 = 0;

        if (t2 < t1)
        {
            std::swap(t1, t2);
        }

        intersection = position + t1*direction;
        save(t1, true);
        intersection = position + t2*direction;
        save(t2, false);

        if (inner_radius_ > 0)
        {
            A = difference_length_squared - inner_radius_ * inner_radius_;

            determinant = B * B - A;

            if (determinant > 0) // determinant == 0 (boundery point) is ignored
            {
                t1 = -1 * B + std::sqrt(determinant);
                t2 = -1 * B - std::sqrt(determinant);

                // Computer precision controll
                if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1 = 0;
                if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                    t2 = 0;

                if (t2 < t1)
                {
                    std::swap(t1, t2);
                }

                intersection = position + t1*direction;
                save(t1, false);
                intersection = position + t2*direction;
                save(t2, true);
            }
        }
    }

    std::function<bool(Intersection const &, Intersection const &)> comp = [](Intersection const & a, Intersection const & b){
        return a.distance < b.distance;
    };

    std::sort(dist.begin(), dist.end(), comp);
    return dist;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Sphere::ComputeDistanceToBorder(const LI::math::Vector3D& position, const LI::math::Vector3D& direction) const
{
    // Compute the surface intersections
    std::vector<Intersection> intersections = Intersections(position, direction);
    std::vector<double> dist;
    bool first = true;
    for(unsigned int i=0; i<intersections.size(); ++i) {
        Intersection const & obj = intersections[i];
        if(obj.distance > 0) {
            if(first) {
                first = false;
                dist.push_back(obj.distance);
                if(not obj.entering) {
                    break;
                }
            }
            else {
                if(not obj.entering) {
                    dist.push_back(obj.distance);
                    break;
                }
                else {
                    throw(std::runtime_error("There should never be two \"entering\" intersections in a row!"));
                }
            }
        }
    }

    std::pair<double, double> distance;

    // No intersection with the outer cylinder
    if (dist.size() < 1)
    {
        distance.first  = -1;
        distance.second = -1;
        //    return distance;
    } else if (dist.size() == 1) // particle is inside the cylinder
    {
        distance.first  = dist.at(0);
        distance.second = -1;

    } else if (dist.size() == 2) // cylinder is infront of the particle
    {
        distance.first  = dist.at(0);
        distance.second = dist.at(1);

        if (distance.second < distance.first)
        {
            std::swap(distance.first, distance.second);
        }

    } else
    {
        //log_error("This point should never be reached");
    }
    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside

    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

} // namespace geometry
} // namespace LI
