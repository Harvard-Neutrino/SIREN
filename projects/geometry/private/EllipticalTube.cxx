#include "SIREN/geometry/EllipticalTube.h"

#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

EllipticalTube::EllipticalTube()
    : Geometry((std::string)("EllipticalTube"))
    , dx_(0.0)
    , dy_(0.0)
    , dz_(0.0) {
}

EllipticalTube::EllipticalTube(double dx, double dy, double dz)
    : Geometry((std::string)("EllipticalTube"))
    , dx_(dx)
    , dy_(dy)
    , dz_(dz) {
    if(dx_ <= 0 || dy_ <= 0 || dz_ <= 0) {
        throw std::invalid_argument("EllipticalTube: dx, dy, dz must all be positive!");
    }
}

EllipticalTube::EllipticalTube(Placement const & placement)
    : Geometry((std::string)("EllipticalTube"), placement)
    , dx_(0.0)
    , dy_(0.0)
    , dz_(0.0) {
}

EllipticalTube::EllipticalTube(Placement const & placement, double dx, double dy, double dz)
    : Geometry((std::string)("EllipticalTube"), placement)
    , dx_(dx)
    , dy_(dy)
    , dz_(dz) {
    if(dx_ <= 0 || dy_ <= 0 || dz_ <= 0) {
        throw std::invalid_argument("EllipticalTube: dx, dy, dz must all be positive!");
    }
}

EllipticalTube::EllipticalTube(const EllipticalTube& other)
    : Geometry(other)
    , dx_(other.dx_)
    , dy_(other.dy_)
    , dz_(other.dz_) {
}

// ------------------------------------------------------------------------- //
void EllipticalTube::swap(Geometry& geometry) {
    EllipticalTube* et = dynamic_cast<EllipticalTube*>(&geometry);
    if(!et) return;

    Geometry::swap(*et);

    std::swap(dx_, et->dx_);
    std::swap(dy_, et->dy_);
    std::swap(dz_, et->dz_);
}

// ------------------------------------------------------------------------- //
EllipticalTube& EllipticalTube::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const EllipticalTube* et = dynamic_cast<const EllipticalTube*>(&geometry);
        if(!et) return *this;

        EllipticalTube tmp(*et);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool EllipticalTube::equal(const Geometry& geometry) const {
    const EllipticalTube* et = dynamic_cast<const EllipticalTube*>(&geometry);
    if(!et) return false;
    if(dx_ != et->dx_) return false;
    if(dy_ != et->dy_) return false;
    if(dz_ != et->dz_) return false;
    return true;
}

// ------------------------------------------------------------------------- //
bool EllipticalTube::less(const Geometry& geometry) const {
    const EllipticalTube* et = dynamic_cast<const EllipticalTube*>(&geometry);
    if(!et) return false;
    return std::tie(dx_, dy_, dz_) < std::tie(et->dx_, et->dy_, et->dz_);
}

// ------------------------------------------------------------------------- //
void EllipticalTube::print(std::ostream& os) const {
    os << "dx: " << dx_ << "\tdy: " << dy_ << "\tdz: " << dz_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> EllipticalTube::ComputeIntersections(
        siren::math::Vector3D const & position,
        siren::math::Vector3D const & direction) const {

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dirx = direction.GetX();
    double diry = direction.GetY();
    double dirz = direction.GetZ();

    Intersection hits[4];
    int nhits = 0;

    double idx2 = 1.0 / (dx_ * dx_);
    double idy2 = 1.0 / (dy_ * dy_);

    // Origin shift: move to closest approach point to reduce magnitudes and
    // prevent catastrophic cancellation in the discriminant at far-field.
    double t_shift = -(px*dirx + py*diry + pz*dirz);
    double qx = px + t_shift * dirx;
    double qy = py + t_shift * diry;

    // Elliptical barrel: (x/dx)^2 + (y/dy)^2 = 1
    // Substituting x = qx + t*dirx, y = qy + t*diry (shifted origin):
    // A*t^2 + 2*B*t + C = 0
    double A = dirx * dirx * idx2 + diry * diry * idy2;
    double B = qx * dirx * idx2 + qy * diry * idy2;
    double C = qx * qx * idx2 + qy * qy * idy2 - 1.0;

    double disc = B * B - A * C;

    if(disc > 0 && A > GEOMETRY_PRECISION * GEOMETRY_PRECISION) {
        double sqrt_disc = std::sqrt(disc);
        double inv_A = 1.0 / A;

        double t1 = (-B - sqrt_disc) * inv_A + t_shift;
        double t2 = (-B + sqrt_disc) * inv_A + t_shift;

        // Check t1
        {
            double iz = pz + t1 * dirz;
            if(std::fabs(iz) <= dz_) {
                double ix = px + t1 * dirx;
                double iy = py + t1 * diry;
                // Outward normal on ellipse: (x/dx^2, y/dy^2, 0)
                double ndot = ix * dirx * idx2 + iy * diry * idy2;
                bool entering = (ndot < 0);
                hits[nhits].distance = t1;
                hits[nhits].hierarchy = 0;
                hits[nhits].entering = entering;
                hits[nhits].position = siren::math::Vector3D(ix, iy, iz);
                ++nhits;
            }
        }

        // Check t2
        {
            double iz = pz + t2 * dirz;
            if(std::fabs(iz) <= dz_) {
                double ix = px + t2 * dirx;
                double iy = py + t2 * diry;
                double ndot = ix * dirx * idx2 + iy * diry * idy2;
                bool entering = (ndot < 0);
                hits[nhits].distance = t2;
                hits[nhits].hierarchy = 0;
                hits[nhits].entering = entering;
                hits[nhits].position = siren::math::Vector3D(ix, iy, iz);
                ++nhits;
            }
        }
    }

    // Endcaps: z = +dz and z = -dz
    if(std::fabs(dirz) > GEOMETRY_PRECISION) {
        // Top cap: z = +dz
        {
            double t = (dz_ - pz) / dirz;
            double ix = px + t * dirx;
            double iy = py + t * diry;
            double ellipse_val = ix * ix * idx2 + iy * iy * idy2;
            if(ellipse_val <= 1.0 + GEOMETRY_PRECISION) {
                bool entering = (dirz < 0);
                hits[nhits].distance = t;
                hits[nhits].hierarchy = 0;
                hits[nhits].entering = entering;
                hits[nhits].position = siren::math::Vector3D(ix, iy, dz_);
                ++nhits;
            }
        }
        // Bottom cap: z = -dz
        {
            double t = (-dz_ - pz) / dirz;
            double ix = px + t * dirx;
            double iy = py + t * diry;
            double ellipse_val = ix * ix * idx2 + iy * iy * idy2;
            if(ellipse_val <= 1.0 + GEOMETRY_PRECISION) {
                bool entering = (dirz > 0);
                hits[nhits].distance = t;
                hits[nhits].hierarchy = 0;
                hits[nhits].entering = entering;
                hits[nhits].position = siren::math::Vector3D(ix, iy, -dz_);
                ++nhits;
            }
        }
    }

    if(nhits == 0) return {};

    std::sort(hits, hits + nhits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });

    return std::vector<Intersection>(hits, hits + nhits);
}

// ------------------------------------------------------------------------- //
AABB EllipticalTube::GetBoundingBox() const {
    return AABB(
        siren::math::Vector3D(-dx_, -dy_, -dz_),
        siren::math::Vector3D( dx_,  dy_,  dz_)
    );
}

} // namespace geometry
} // namespace siren
