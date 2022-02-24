/*
 * Geometry.cxx
 *
 *  Created on: 05.06.2013
 *      Author: koehne
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/Placement.h"

using namespace earthmodel;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace earthmodel {

std::ostream& operator<<(std::ostream& os, Geometry const& geometry)
{
    os << "Geometry(" << &geometry << ")" << std::endl;
    os << geometry.placement_ << std::endl;
    geometry.print(os);
    return os;
}

} // namespace earthmodel

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

bool Geometry::IsInside(const Vector3D& position, const Vector3D& direction) const
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
bool Geometry::IsInfront(const Vector3D& position, const Vector3D& direction) const
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
bool Geometry::IsBehind(const Vector3D& position, const Vector3D& direction) const
{
    bool is_behind = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first < 0 && dist.second < 0)
    {
        is_behind = true;
    }
    return is_behind;
}

Geometry::ParticleLocation::Enum Geometry::GetLocation(const Vector3D& position, const Vector3D& direction) const {
    if(IsInfront(position, direction))
        return Geometry::ParticleLocation::InfrontGeometry;
    if(IsInside(position, direction))
        return Geometry::ParticleLocation::InsideGeometry;
    else
        return Geometry::ParticleLocation::BehindGeometry;
}

// ------------------------------------------------------------------------- //
double Geometry::DistanceToClosestApproach(const Vector3D& position, const Vector3D& direction) const
{
    Vector3D pos = GlobalToLocalPosition(position);
    Vector3D dir = GlobalToLocalPosition(direction);
    return scalar_product(-pos, dir);
}

Vector3D Geometry::LocalToGlobalPosition(Vector3D const & p0) const
{
    return placement_.LocalToGlobalPosition(p0);
}

Vector3D Geometry::LocalToGlobalDirection(Vector3D const & p0) const
{
    return placement_.LocalToGlobalDirection(p0);
}

Vector3D Geometry::GlobalToLocalPosition(Vector3D const & p0) const
{
    return placement_.GlobalToLocalPosition(p0);
}

Vector3D Geometry::GlobalToLocalDirection(Vector3D const & p0) const
{
    return placement_.GlobalToLocalDirection(p0);
}

std::pair<double, double> Geometry::DistanceToBorder(const Vector3D& position, const Vector3D& direction) const {
    Vector3D local_position = GlobalToLocalPosition(position);
    Vector3D local_direction = GlobalToLocalDirection(direction);
    return ComputeDistanceToBorder(position, direction);
}
std::vector<Geometry::Intersection> Geometry::Intersections(Vector3D const & position, Vector3D const & direction) const {
    Vector3D local_position = GlobalToLocalPosition(position);
    Vector3D local_direction = GlobalToLocalDirection(direction);
    std::vector<Geometry::Intersection> intersections = ComputeIntersections(local_position, local_direction);
    for(auto & intersection : intersections) {
        intersection.position = LocalToGlobalPosition(intersection.position);
    }
    return intersections;
}

Box::Box()
    : Geometry((std::string)("Box"))
    , x_(0.0)
    , y_(0.0)
      , z_(0.0)
{
    // Do nothing here
}

Box::Box(double x, double y, double z)
    : Geometry("Box")
    , x_(x)
    , y_(y)
      , z_(z)
{
    // Do nothing here
}

Box::Box(Placement const & placement)
    : Geometry((std::string)("Box"), placement)
    , x_(0.0)
    , y_(0.0)
      , z_(0.0)
{
    // Do nothing here
}

Box::Box(Placement const & placement, double x, double y, double z)
    : Geometry((std::string)("Box"), placement)
    , x_(x)
    , y_(y)
      , z_(z)
{
    // Do nothing here
}

Box::Box(const Box& box)
    : Geometry(box)
    , x_(box.x_)
    , y_(box.y_)
      , z_(box.z_)
{
    // Nothing to do here
}

/*Box::Box(const nlohmann::json& config)
  : Geometry(config)
  {
  if(not config.at("length").is_number())
  throw std::invalid_argument("Length is not a number.");
  if(not config.at("width").is_number())
  throw std::invalid_argument("Width is not a number.");
  if(not config.at("height").is_number())
  throw std::invalid_argument("Height is not a number.");

  x_ = config["length"].get<double>();
  y_ = config["width"].get<double>();
  z_ = config["height"].get<double>();

  if(x_ < 0) throw std::logic_error("lenght must be > 0");
  if(y_ < 0) throw std::logic_error("width must be > 0");
  if(z_ < 0) throw std::logic_error("height must be > 0");
  }*/

// ------------------------------------------------------------------------- //
void Box::swap(Geometry& geometry)
{
    Box* box = dynamic_cast<Box*>(&geometry);
    if (!box)
    {
        //log_warn("Cannot swap Box!");
        return;
    }

    Geometry::swap(*box);

    std::swap(x_, box->x_);
    std::swap(y_, box->y_);
    std::swap(z_, box->z_);
}

//------------------------------------------------------------------------- //
Box& Box::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Box* box = dynamic_cast<const Box*>(&geometry);
        if (!box)
        {
            //log_warn("Cannot assign Sphere!");
            return *this;
        }

        Box tmp(*box);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Box::equal(const Geometry& geometry) const
{
    const Box* box = dynamic_cast<const Box*>(&geometry);

    if (!box)
        return false;
    else if (x_ != box->x_)
        return false;
    else if (y_ != box->y_)
        return false;
    else if (z_ != box->z_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Box::less(const Geometry& geometry) const
{
    const Box* box = dynamic_cast<const Box*>(&geometry);

    return
        std::tie(x_, y_, z_)
        <
        std::tie(box->x_, box->y_, box->z_);
}

// ------------------------------------------------------------------------- //
void Box::print(std::ostream& os) const
{
    os << "Width_x: " << x_ << "\tWidth_y " << y_ << "\tHeight: " << z_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Box::ComputeIntersections(Vector3D const & position, Vector3D const & direction) const {
    // Calculate intersection of particle trajectory and the box
    // Surface of the box is defined by six planes:
    // E1: x1   =   position.GetX() + 0.5*x
    // E2: x1   =   position.GetX() - 0.5*x
    // E3: x2   =   position.GetY() + 0.5*y
    // E4: x2   =   position.GetY() - 0.5*y
    // E5: x3   =   position.GetZ() + 0.5*z
    // E6: x3   =   position.GetZ() - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    double dir_vec_x = direction.GetX();
    double dir_vec_y = direction.GetY();
    double dir_vec_z = direction.GetZ();

    double t;
    double intersection_x;
    double intersection_y;
    double intersection_z;
    bool entering;

    std::vector<Intersection> dist;

    std::function<void()> save = [&](){
        Intersection i;
        i.position = Vector3D(intersection_x,intersection_y,intersection_z);
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };

    double x_calc_pos =   0.5 * x_;
    double x_calc_neg = - 0.5 * x_;
    double y_calc_pos =   0.5 * y_;
    double y_calc_neg = - 0.5 * y_;
    double z_calc_pos =   0.5 * z_;
    double z_calc_neg = - 0.5 * z_;

    // intersection with E1
    if (dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E1
    {
        t = (x_calc_pos - position.GetX()) / dir_vec_x;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_y = position.GetY() + t * dir_vec_y;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_y >= y_calc_neg && intersection_y <= y_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_x = position.GetX() + t * dir_vec_x;
            entering = direction.GetX() < 0;
            save();
        }
    }

    // intersection with E2
    if (dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E2
    {
        t = (x_calc_neg - position.GetX()) / dir_vec_x;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_y = position.GetY() + t * dir_vec_y;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_y >= y_calc_neg && intersection_y <= y_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_x = position.GetX() + t * dir_vec_x;
            entering = direction.GetX() > 0;
            save();
        }
    }

    // intersection with E3
    if (dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E3
    {
        t = (y_calc_pos - position.GetY()) / dir_vec_y;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_y = position.GetY() + t * dir_vec_y;
            entering = direction.GetY() < 0;
            save();
        }
    }

    // intersection with E4
    if (dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E4
    {
        t = (y_calc_neg - position.GetY()) / dir_vec_y;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_y = position.GetY() + t * dir_vec_y;
            entering = direction.GetY() > 0;
            save();
        }
    }

    // intersection with E5
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E5
    {
        t = (z_calc_pos - position.GetZ()) / dir_vec_z;

        // Computer precision controll
        if (std::fabs(t) < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_y = position.GetY() + t * dir_vec_y;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_y >= y_calc_neg &&
                intersection_y <= y_calc_pos)
        {
            intersection_z = position.GetZ() + t * dir_vec_z;
            entering = direction.GetZ() < 0;
            save();
        }
    }

    // intersection with E6
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E6
    {
        t = (z_calc_neg - position.GetZ()) / dir_vec_z;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_y = position.GetY() + t * dir_vec_y;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_y >= y_calc_neg &&
                intersection_y <= y_calc_pos)
        {
            intersection_z = position.GetZ() + t * dir_vec_z;
            entering = direction.GetZ() > 0;
            save();
        }
    }

    std::function<bool(Intersection const &, Intersection const &)> comp = [](Intersection const & a, Intersection const & b){
        return a.distance < b.distance;
    };

    std::sort(dist.begin(), dist.end(), comp);
    return dist;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Box::ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const
{
    // Compute the surface intersections
    std::vector<Intersection> intersections = Intersections(position, direction);
    std::vector<double> dist;
    for(unsigned int i=0; i<intersections.size(); ++i) {
        if(intersections[i].distance > 0) {
            dist.push_back(intersections[i].distance);
        }
    }

    std::pair<double, double> distance;

    if (dist.size() < 1) // No intersection with the box
    {
        distance.first  = -1;
        distance.second = -1;
    } else if (dist.size() == 1) // Particle is inside the box and we have one
        // intersection in direction of the particle
        // trajectory
    {
        distance.first  = dist.at(0);
        distance.second = -1;
    } else if (dist.size() == 2) // Particle is outside and the box is infront
        // of the particle trajectory ( two
        // intersections).
    {
        distance.first  = dist.at(0);
        distance.second = dist.at(1);
        if (distance.second < distance.first)
        {
            std::swap(distance.first, distance.second);
        }

    } else
    {
        //log_error("This point should nerver be reached... (-1/-1) is returned");

        distance.first  = -1;
        distance.second = -1;
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
std::vector<Geometry::Intersection> Cylinder::ComputeIntersections(Vector3D const & position, Vector3D const & direction) const {
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
    bool entering;


    std::vector<Intersection> dist;

    std::function<void(double, bool)> save = [&](double t, bool entering){
        Intersection i;
        i.position = Vector3D(intersection_x,intersection_y,intersection_z);
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };

    std::function<bool()> entering_radial = [&]() {
        return Vector3D(intersection_x, intersection_y, 0) * direction < 0;
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
std::pair<double, double> Cylinder::ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const
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
std::vector<Geometry::Intersection> Sphere::ComputeIntersections(Vector3D const & position, Vector3D const & direction) const {
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

    Vector3D intersection;

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
std::pair<double, double> Sphere::ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const
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

ExtrPoly::ExtrPoly()
    : Geometry((std::string)("ExtrPoly"))
      , polygon_({})
, zsections_({})
{
    // Do nothing here
}


ExtrPoly::ExtrPoly(const std::vector<std::vector<double>>& polygon,
        const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry("ExtrPoly")
    , polygon_(polygon)
      , zsections_(zsections)
{
    if (polygon.size() < 3)
    {
        std::cout << "Need 3 polygon vertices at least!! Give it another shot";
        return;
    }
    ComputeLateralPlanes();
}

ExtrPoly::ExtrPoly(Placement const & placement)
    : Geometry((std::string)("ExtrPoly"), placement)
      , polygon_({})
, zsections_({})
{
    ComputeLateralPlanes();
}


ExtrPoly::ExtrPoly(Placement const & placement, const std::vector<std::vector<double>>& polygon,
        const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry((std::string)("ExtrPoly"), placement)
    , polygon_(polygon)
      , zsections_(zsections)
{
    if (polygon.size() < 3)
    {
        std::cout << "Need 3 polygon vertices at least!! Give it another shot";
        return;
    }
    ComputeLateralPlanes();
}

ExtrPoly::ExtrPoly(const ExtrPoly& extr)
    : Geometry(extr)
    , polygon_(extr.polygon_)
      , zsections_(extr.zsections_)
{
    ComputeLateralPlanes();
}


// ------------------------------------------------------------------------- //
void ExtrPoly::swap(Geometry& geometry)
{
    ExtrPoly* extr = dynamic_cast<ExtrPoly*>(&geometry);
    if (!extr)
    {
        //log_warn("Cannot swap ExtrPoly!");
        return;
    }

    Geometry::swap(*extr);

    std::swap(polygon_, extr->polygon_);
    std::swap(zsections_, extr->zsections_);
}

//------------------------------------------------------------------------- //
ExtrPoly& ExtrPoly::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);
        if (!extr)
        {
            //log_warn("Cannot assign ExtrPoly!");
            return *this;
        }

        ExtrPoly tmp(*extr);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::equal(const Geometry& geometry) const
{
    const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);

    if (!extr)
        return false;
    else if (polygon_ != extr->polygon_)
        return false;
    else if (!(zsections_ == extr->zsections_))
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::less(const Geometry& geometry) const
{
    const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);

    return
        std::tie(polygon_, zsections_)
        <
        std::tie(extr->polygon_, extr->zsections_);
}

// ------------------------------------------------------------------------- //
void ExtrPoly::print(std::ostream& os) const
{
    os << "Not implemented yet\n";
    //os << "Polygon: " << polygon_ << "\tZsections: " << zsections_ << '\n';
}

// ------------------------------------------------------------------------- //
void ExtrPoly::ComputeLateralPlanes()
{
    int Nv = polygon_.size();
    planes_.resize(Nv);
    for (int i=0, k=Nv-1; i<Nv; k = i++)
    {
        std::vector<double> dir = {polygon_[i][0] - polygon_[k][0],polygon_[i][1] - polygon_[k][1]};
        double norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
        dir[0]/=norm; dir[1]/=norm;
        planes_[i].a = -dir[1];
        planes_[i].b = dir[0];
        planes_[i].c = 0;
        planes_[i].d = dir[1]*polygon_[i][0] - dir[0]*polygon_[i][1];
    }
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> ExtrPoly::ComputeIntersections(Vector3D const & position, Vector3D const & direction) const {
    // Calculate intersection of particle trajectory and the extr poly
    // Implementation follows that of Geant4, see here:
    //
    // https://gitlab.cern.ch/geant4/geant4/-/blob/master/source/geometry/solids/specific/src/G4ExtrudedSolid.cc
    //
    // NOTE: Only works for convex right prisms at the moment

    std::vector<Geometry::Intersection> dist;

    Vector3D intersection;

    std::function<void(double, bool)> save = [&](double t, bool entering){
        Intersection i;
        i.position = Vector3D(position.GetX() + direction.GetX()*t,
                position.GetY() + direction.GetY()*t,
                position.GetZ() + direction.GetZ()*t);
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };

    
    int Nz = zsections_.size();
    double z0 = zsections_[0].zpos;
    double z1 = zsections_[Nz-1].zpos;

    if ((position.GetZ() <= z0 + GEOMETRY_PRECISION) && direction.GetZ() <= 0) return dist;
    if ((position.GetZ() >= z1 - GEOMETRY_PRECISION) && direction.GetZ() >= 0) return dist;

    // Intersection with Z planes
    double dz = (z1 - z0)*0.5;
    double pz = position.GetZ() - dz - z0;

    double invz = (direction.GetZ() == 0) ? DBL_MAX : -1./direction.GetZ();
    double ddz = (invz < 0) ? dz : -dz;
    double tzmin = (pz + ddz)*invz;
    double tzmax = (pz - ddz)*invz;

    // Intersection with lateral planes
    int np = planes_.size();
    double txmin = tzmin, txmax = tzmax;
    for (int i=0; i<np; ++i)
    {
        double cosa = planes_[i].a*direction.GetX()+planes_[i].b*direction.GetY();
        double distnce = planes_[i].a*position.GetX()+planes_[i].b*position.GetY()+planes_[i].d;
        // case 1: particle is outside of outward-facing normal vector of plane in XY projection
        if (distnce >= -GEOMETRY_PRECISION)
        {
            if (cosa >= 0) { return dist; } // If particle is currently moving away from any face, it will never intersect
            double tmp  = -distnce/cosa;
            if (txmin < tmp)  { txmin = tmp; }

        }
        // case 2: particle is inside of outward-facing normal vector of plane
        else if (cosa > 0)
        {
            double tmp  = -distnce/cosa;
            if (txmax > tmp)  { txmax = tmp; }
        }
    }
    double tmin = txmin, tmax = txmax;
		if (tmax <= tmin + GEOMETRY_PRECISION)   // touch or no hit
      {
        return dist;
      }


    save(tmin,true);
    save(tmax,false);

    std::function<bool(Intersection const &, Intersection const &)> comp = [](Intersection const & a, Intersection const & b){
        return a.distance < b.distance;
    };

    std::sort(dist.begin(), dist.end(), comp);
    return dist;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> ExtrPoly::ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const
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
