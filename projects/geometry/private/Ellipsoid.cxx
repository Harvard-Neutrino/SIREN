#include "SIREN/geometry/Ellipsoid.h"

#include <cmath>
#include <tuple>
#include <math.h>
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

Ellipsoid::Ellipsoid()
    : Geometry((std::string)("Ellipsoid"))
    , ax_(0.0)
    , by_(0.0)
    , cz_(0.0)
    , zcut1_(0.0)
      , zcut2_(0.0) {
    // Do nothing here
}

Ellipsoid::Ellipsoid(double ax, double by, double cz)
    : Geometry((std::string)("Ellipsoid"))
    , ax_(ax)
    , by_(by)
    , cz_(cz)
    , zcut1_(-cz)
      , zcut2_(cz) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) {
        throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    }
}

Ellipsoid::Ellipsoid(double ax, double by, double cz, double zcut1, double zcut2)
    : Geometry((std::string)("Ellipsoid"))
    , ax_(ax)
    , by_(by)
    , cz_(cz)
    , zcut1_(zcut1)
      , zcut2_(zcut2) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) {
        throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    }
    // Clamp z-cuts to valid range
    if(zcut1_ < -cz_) zcut1_ = -cz_;
    if(zcut1_ >  cz_) zcut1_ =  cz_;
    if(zcut2_ < -cz_) zcut2_ = -cz_;
    if(zcut2_ >  cz_) zcut2_ =  cz_;
    if(zcut1_ >= zcut2_) {
        throw std::invalid_argument("Ellipsoid zcut1 must be less than zcut2!");
    }
}

Ellipsoid::Ellipsoid(Placement const & placement)
    : Geometry((std::string)("Ellipsoid"), placement)
    , ax_(0.0)
    , by_(0.0)
    , cz_(0.0)
    , zcut1_(0.0)
      , zcut2_(0.0) {
    // Do nothing here
}

Ellipsoid::Ellipsoid(Placement const & placement, double ax, double by, double cz)
    : Geometry((std::string)("Ellipsoid"), placement)
    , ax_(ax)
    , by_(by)
    , cz_(cz)
    , zcut1_(-cz)
      , zcut2_(cz) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) {
        throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    }
}

Ellipsoid::Ellipsoid(Placement const & placement, double ax, double by, double cz, double zcut1, double zcut2)
    : Geometry((std::string)("Ellipsoid"), placement)
    , ax_(ax)
    , by_(by)
    , cz_(cz)
    , zcut1_(zcut1)
      , zcut2_(zcut2) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) {
        throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    }
    if(zcut1_ < -cz_) zcut1_ = -cz_;
    if(zcut1_ >  cz_) zcut1_ =  cz_;
    if(zcut2_ < -cz_) zcut2_ = -cz_;
    if(zcut2_ >  cz_) zcut2_ =  cz_;
    if(zcut1_ >= zcut2_) {
        throw std::invalid_argument("Ellipsoid zcut1 must be less than zcut2!");
    }
}

Ellipsoid::Ellipsoid(const Ellipsoid& other)
    : Geometry(other)
    , ax_(other.ax_)
    , by_(other.by_)
    , cz_(other.cz_)
    , zcut1_(other.zcut1_)
      , zcut2_(other.zcut2_) {
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Ellipsoid::swap(Geometry& geometry) {
    Ellipsoid* other = dynamic_cast<Ellipsoid*>(&geometry);
    if(!other) {
        //log_warn("Cannot swap Ellipsoid!");
        return;
    }

    Geometry::swap(*other);

    std::swap(ax_, other->ax_);
    std::swap(by_, other->by_);
    std::swap(cz_, other->cz_);
    std::swap(zcut1_, other->zcut1_);
    std::swap(zcut2_, other->zcut2_);
}

//------------------------------------------------------------------------- //
Ellipsoid& Ellipsoid::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Ellipsoid* other = dynamic_cast<const Ellipsoid*>(&geometry);
        if(!other) {
            //log_warn("Cannot assign Ellipsoid!");
            return *this;
        }

        Ellipsoid tmp(*other);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Ellipsoid::equal(const Geometry& geometry) const
{
    const Ellipsoid* other = dynamic_cast<const Ellipsoid*>(&geometry);

    if(!other)
        return false;
    else if(ax_ != other->ax_)
        return false;
    else if(by_ != other->by_)
        return false;
    else if(cz_ != other->cz_)
        return false;
    else if(zcut1_ != other->zcut1_)
        return false;
    else if(zcut2_ != other->zcut2_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Ellipsoid::less(const Geometry& geometry) const
{
    const Ellipsoid* other = dynamic_cast<const Ellipsoid*>(&geometry);
    if(!other) return false;

    return
        std::tie(ax_, by_, cz_, zcut1_, zcut2_)
        <
        std::tie(other->ax_, other->by_, other->cz_, other->zcut1_, other->zcut2_);
}

// ------------------------------------------------------------------------- //
void Ellipsoid::print(std::ostream& os) const
{
    os << "ax: " << ax_ << "\tby: " << by_ << "\tcz: " << cz_
       << "\tzcut1: " << zcut1_ << "\tzcut2: " << zcut2_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Ellipsoid::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Ray-ellipsoid intersection.
    // Surface: (x/ax)^2 + (y/by)^2 + (z/cz)^2 = 1
    // Ray: P + t*D
    // Substituting gives A*t^2 + 2*B*t + C = 0

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    double inv_ax2 = 1.0 / (ax_ * ax_);
    double inv_by2 = 1.0 / (by_ * by_);
    double inv_cz2 = 1.0 / (cz_ * cz_);

    double A = dx * dx * inv_ax2 + dy * dy * inv_by2 + dz * dz * inv_cz2;
    double B = px * dx * inv_ax2 + py * dy * inv_by2 + pz * dz * inv_cz2;
    double C = px * px * inv_ax2 + py * py * inv_by2 + pz * pz * inv_cz2 - 1.0;

    double disc = B * B - A * C;

    Intersection hits[4];
    int n_hits = 0;

    // Ellipsoid surface intersections
    if(disc >= 0.0) {
        double sqrt_disc = std::sqrt(disc);
        double inv_A = 1.0 / A;

        double t1 = (-B - sqrt_disc) * inv_A;
        double t2 = (-B + sqrt_disc) * inv_A;

        // Process both roots
        for(int i = 0; i < 2; ++i) {
            double t = (i == 0) ? t1 : t2;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            double ix = px + t * dx;
            double iy = py + t * dy;
            double iz = pz + t * dz;

            // Check z-cut bounds
            if(iz >= zcut1_ - GEOMETRY_PRECISION && iz <= zcut2_ + GEOMETRY_PRECISION) {
                // Outward normal: gradient of (x/ax)^2 + (y/by)^2 + (z/cz)^2
                // = (2x/ax^2, 2y/by^2, 2z/cz^2), factor of 2 irrelevant for sign
                double nx = ix * inv_ax2;
                double ny = iy * inv_by2;
                double nz = iz * inv_cz2;
                bool entering = (nx * dx + ny * dy + nz * dz) < 0;

                hits[n_hits].distance = t;
                hits[n_hits].hierarchy = 0;
                hits[n_hits].entering = entering;
                hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
                n_hits++;
            }
        }
    }

    // Z-cut plane intersections
    // Bottom cut: z = zcut1_
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
        double t = (zcut1_ - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        double ix = px + t * dx;
        double iy = py + t * dy;
        double iz = pz + t * dz;

        // Check if point is inside or on the ellipsoid at this z
        double val = ix * ix * inv_ax2 + iy * iy * inv_by2 + iz * iz * inv_cz2;
        if(val <= 1.0 + GEOMETRY_PRECISION) {
            // Entering bottom cut when moving upward (dz > 0)
            bool entering = dz > 0;
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = entering;
            hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
            n_hits++;
        }
    }

    // Top cut: z = zcut2_
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
        double t = (zcut2_ - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        double ix = px + t * dx;
        double iy = py + t * dy;
        double iz = pz + t * dz;

        double val = ix * ix * inv_ax2 + iy * iy * inv_by2 + iz * iz * inv_cz2;
        if(val <= 1.0 + GEOMETRY_PRECISION) {
            // Entering top cut when moving downward (dz < 0)
            bool entering = dz < 0;
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = entering;
            hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
            n_hits++;
        }
    }

    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    return {hits, hits + n_hits};
}

// ------------------------------------------------------------------------- //
AABB Ellipsoid::GetBoundingBox() const {
    return AABB(
        math::Vector3D(-ax_, -by_, zcut1_),
        math::Vector3D( ax_,  by_, zcut2_)
    );
}

} // namespace geometry
} // namespace siren
