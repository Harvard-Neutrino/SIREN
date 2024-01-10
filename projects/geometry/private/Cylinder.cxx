#include "LeptonInjector/geometry/Cylinder.h"

#include <cmath>
#include <tuple>
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
            //log_warn("Cannot assign Sphere!");
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
std::vector<Geometry::Intersection> Cylinder::ComputeIntersections(LI::math::Vector3D const & position, LI::math::Vector3D const & direction) const {
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

    double A, B, C, t1, t2, t;
    double dir_vec_x = direction.GetX();
    double dir_vec_y = direction.GetY();
    double dir_vec_z = direction.GetZ();

    double determinant;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    std::vector<Intersection> dist;

    std::function<void(double, bool)> save = [&](double t, bool entering){
        Intersection i;
        i.position = LI::math::Vector3D(intersection_x,intersection_y,intersection_z);
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };

    std::function<bool()> entering_radial = [&]() {
        return LI::math::Vector3D(intersection_x, intersection_y, 0) * direction < 0;
    };

    double z_calc_pos = 0.5 * z_;
    double z_calc_neg = -0.5 * z_;

    if (!(dir_vec_x == 0 && dir_vec_y == 0)) // Otherwise the particle
        // trajectory is parallel to
        // cylinder barrel
    {

        A = std::pow((position.GetX()), 2) +
            std::pow((position.GetY()), 2) -
            radius_*radius_;

        B = 2 * ((position.GetX()) * dir_vec_x + (position.GetY()) * dir_vec_y);

        C = dir_vec_x * dir_vec_x + dir_vec_y * dir_vec_y;

        B /= C;
        A /= C;

        determinant = 0.25 * B*B - A;

        if (determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1 = -1 * B / 2 + std::sqrt(determinant);
            t2 = -1 * B / 2 - std::sqrt(determinant);

            // Computer precision controll
            if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                t1 = 0;
            if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                t2 = 0;

            intersection_z = position.GetZ() + t1 * dir_vec_z;
            // is inside the borders
            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
            {
                intersection_x = position.GetX() + t1 * dir_vec_x;
                intersection_y = position.GetY() + t1 * dir_vec_y;
                save(t1, entering_radial());
            }

            intersection_z = position.GetZ() + t2 * dir_vec_z;
            // is inside the borders
            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
            {
                intersection_x = position.GetX() + t2 * dir_vec_x;
                intersection_y = position.GetY() + t2 * dir_vec_y;
                save(t2, entering_radial());
            }
        }
    }

    // intersection with E1
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel
        // to E1 (Should not happen)
    {
        t = (z_calc_pos - position.GetZ()) / dir_vec_z;
        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_y = position.GetY() + t * dir_vec_y;

        if (std::sqrt(std::pow((intersection_x), 2) +
                    std::pow((intersection_y), 2)) <=
                radius_ &&
                std::sqrt(std::pow((intersection_x), 2) +
                    std::pow((intersection_y), 2)) >=
                inner_radius_)
        {
            intersection_z = position.GetZ() + t * dir_vec_z;
            save(t, direction.GetZ() < 0);
        }
    }

    // intersection with E2
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel
        // to E2 (Should not happen)
    {
        t = (z_calc_neg - position.GetZ()) / dir_vec_z;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_y = position.GetY() + t * dir_vec_y;

        if (std::sqrt(std::pow((intersection_x), 2) +
                    std::pow((intersection_y), 2)) <=
                radius_ &&
                std::sqrt(std::pow((intersection_x), 2) +
                    std::pow((intersection_y), 2)) >=
                inner_radius_)
        {
            intersection_z = position.GetZ() + t * dir_vec_z;
            save(t, direction.GetZ() > 0);
        }
    }

    // This cylinder might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner cylinder.

    if (inner_radius_ > 0)
    {
        if (!(dir_vec_x == 0 && dir_vec_y == 0))
        {

            A = std::pow((position.GetX()), 2) +
                std::pow((position.GetY()), 2) -
                inner_radius_*inner_radius_;

            B = 2 *
                ((position.GetX()) * dir_vec_x + (position.GetY()) * dir_vec_y);

            C = dir_vec_x * dir_vec_x + dir_vec_y * dir_vec_y;

            B /= C;
            A /= C;

            determinant = 0.25 * B*B - A;

            if (determinant > 0) // determinant == 0 (boundery point) is ignored
            {
                t1 = -1 * B / 2 + std::sqrt(determinant);
                t2 = -1 * B / 2 - std::sqrt(determinant);

                // Computer precision controll
                if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1 = 0;
                if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                    t2 = 0;

                // Ok we have a intersection with the inner cylinder

                intersection_z = position.GetZ() + t1 * dir_vec_z;
                // is inside the borders
                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                {
                    intersection_x = position.GetX() + t1 * dir_vec_x;
                    intersection_y = position.GetY() + t1 * dir_vec_y;
                    save(t1, not entering_radial());
                }

                intersection_z = position.GetZ() + t2 * dir_vec_z;
                // is inside the borders
                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                {
                    intersection_x = position.GetX() + t2 * dir_vec_x;
                    intersection_y = position.GetY() + t2 * dir_vec_y;
                    save(t2, not entering_radial());
                }
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
std::pair<double, double> Cylinder::ComputeDistanceToBorder(const LI::math::Vector3D& position, const LI::math::Vector3D& direction) const
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
