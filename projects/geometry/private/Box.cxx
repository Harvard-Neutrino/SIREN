#include "SIREN/geometry/Box.h"

#include <cmath>
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

namespace siren {
namespace geometry {

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
std::vector<Geometry::Intersection> Box::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    double px = position.GetX(), py = position.GetY(), pz = position.GetZ();
    double dx = direction.GetX(), dy = direction.GetY(), dz = direction.GetZ();

    double hx = 0.5 * x_, hy = 0.5 * y_, hz = 0.5 * z_;

    // At most 2 intersections for a convex shape; use a fixed-size local array
    Intersection hits[2];
    int n_hits = 0;

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

    // Precision control: Note on boundary (t near 0), the original code
    // treats |t| < GEOMETRY_PRECISION as "on the border" and sets to 0.
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
std::pair<double, double> Box::ComputeDistanceToBorder(const siren::math::Vector3D& position, const siren::math::Vector3D& direction) const
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

// ------------------------------------------------------------------------- //
AABB Box::GetBoundingBox() const {
    double hx = x_ * 0.5;
    double hy = y_ * 0.5;
    double hz = z_ * 0.5;
    return AABB(
        math::Vector3D(-hx, -hy, -hz),
        math::Vector3D( hx,  hy,  hz)
    );
}

} // namespace geometry
} // namespace siren
