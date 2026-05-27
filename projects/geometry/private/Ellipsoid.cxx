#include "SIREN/geometry/Ellipsoid.h"

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
#include "GeometryMacros.h"

namespace siren {
namespace geometry {

Ellipsoid::Ellipsoid() : Geometry("Ellipsoid"), ax_(0), by_(0), cz_(0), zcut1_(0), zcut2_(0) { RecomputeWorldAABB(); }
Ellipsoid::Ellipsoid(double ax, double by, double cz) : Geometry("Ellipsoid"), ax_(ax), by_(by), cz_(cz), zcut1_(-cz), zcut2_(cz) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    RecomputeWorldAABB();
}
Ellipsoid::Ellipsoid(double ax, double by, double cz, double zcut1, double zcut2)
    : Geometry("Ellipsoid"), ax_(ax), by_(by), cz_(cz), zcut1_(zcut1), zcut2_(zcut2) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    if(zcut1_ < -cz_) zcut1_ = -cz_;
    if(zcut1_ >  cz_) zcut1_ =  cz_;
    if(zcut2_ < -cz_) zcut2_ = -cz_;
    if(zcut2_ >  cz_) zcut2_ =  cz_;
    if(zcut1_ >= zcut2_) throw std::invalid_argument("Ellipsoid zcut1 must be less than zcut2!");
    RecomputeWorldAABB();
}
Ellipsoid::Ellipsoid(Placement const & p) : Geometry("Ellipsoid", p), ax_(0), by_(0), cz_(0), zcut1_(0), zcut2_(0) { RecomputeWorldAABB(); }
Ellipsoid::Ellipsoid(Placement const & p, double ax, double by, double cz) : Geometry("Ellipsoid", p), ax_(ax), by_(by), cz_(cz), zcut1_(-cz), zcut2_(cz) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    RecomputeWorldAABB();
}
Ellipsoid::Ellipsoid(Placement const & p, double ax, double by, double cz, double zcut1, double zcut2)
    : Geometry("Ellipsoid", p), ax_(ax), by_(by), cz_(cz), zcut1_(zcut1), zcut2_(zcut2) {
    if(ax_ <= 0 || by_ <= 0 || cz_ <= 0) throw std::invalid_argument("Ellipsoid semi-axes must be positive!");
    if(zcut1_ < -cz_) zcut1_ = -cz_;
    if(zcut1_ >  cz_) zcut1_ =  cz_;
    if(zcut2_ < -cz_) zcut2_ = -cz_;
    if(zcut2_ >  cz_) zcut2_ =  cz_;
    if(zcut1_ >= zcut2_) throw std::invalid_argument("Ellipsoid zcut1 must be less than zcut2!");
    RecomputeWorldAABB();
}
Ellipsoid::Ellipsoid(const Ellipsoid& o) : Geometry(o), ax_(o.ax_), by_(o.by_), cz_(o.cz_), zcut1_(o.zcut1_), zcut2_(o.zcut2_) { RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(Ellipsoid, ax_, by_, cz_, zcut1_, zcut2_)
SIREN_GEOMETRY_ASSIGN(Ellipsoid)
SIREN_GEOMETRY_EQUAL(Ellipsoid, ax_, by_, cz_, zcut1_, zcut2_)
SIREN_GEOMETRY_LESS(Ellipsoid, ax_, by_, cz_, zcut1_, zcut2_)

void Ellipsoid::print(std::ostream& os) const {
    os << "Ellipsoid(" << ax_ << ", " << by_ << ", " << cz_ << ", " << zcut1_ << ", " << zcut2_ << ")\n";
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Ellipsoid::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Ray-ellipsoid intersection.
    // Surface: (x/ax)^2 + (y/by)^2 + (z/cz)^2 = 1
    // Ray: P + t*D
    // Substituting gives A*t^2 + 2*B*t + C = 0
    //
    // For far-field rays (|P| >> ellipsoid size), B^2 ~ A*C ~ |P|^2 and the
    // discriminant B^2 - A*C suffers catastrophic cancellation. We shift the
    // ray origin to the closest approach point along the ray to the coordinate
    // origin, keeping the quadratic coefficients small.

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    double inv_ax2 = 1.0 / (ax_ * ax_);
    double inv_by2 = 1.0 / (by_ * by_);
    double inv_cz2 = 1.0 / (cz_ * cz_);

    // Origin shift: move ray origin to closest approach to the coordinate origin.
    // This minimizes |q|^2 and keeps B and C well-conditioned for far-field rays.
    double t_shift = -(px * dx + py * dy + pz * dz);
    double qx = px + t_shift * dx;
    double qy = py + t_shift * dy;
    double qz = pz + t_shift * dz;

    double A = dx * dx * inv_ax2 + dy * dy * inv_by2 + dz * dz * inv_cz2;
    double B = qx * dx * inv_ax2 + qy * dy * inv_by2 + qz * dz * inv_cz2;
    double C = qx * qx * inv_ax2 + qy * qy * inv_by2 + qz * qz * inv_cz2 - 1.0;

    double disc = B * B - A * C;

    Intersection hits[4];
    int n_hits = 0;

    // Ellipsoid surface intersections
    if(disc >= 0.0) {
        double sqrt_disc = std::sqrt(disc);
        double inv_A = 1.0 / A;

        double t1 = (-B - sqrt_disc) * inv_A + t_shift;
        double t2 = (-B + sqrt_disc) * inv_A + t_shift;

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
    // Deduplicate hits at the same distance (can occur where the ellipsoid
    // surface meets a z-cut plane — both branches report the same point).
    int n_unique = 0;
    for(int i = 0; i < n_hits; ++i) {
        if(n_unique > 0 && std::fabs(hits[i].distance - hits[n_unique - 1].distance) < GEOMETRY_PRECISION) {
            continue;
        }
        hits[n_unique++] = hits[i];
    }
    return {hits, hits + n_unique};
}

// ------------------------------------------------------------------------- //
AABB Ellipsoid::GetBoundingBox() const {
    // Tighten x/y extents when z-cuts are active.
    // Max x extent at height z is ax * sqrt(1 - z^2/cz^2).
    // The maximum over [zcut1, zcut2] is at z=0 if 0 is in range,
    // otherwise at whichever cut is closer to z=0.
    double z_for_max = std::clamp(0.0, zcut1_, zcut2_);
    double frac_sq = (z_for_max * z_for_max) / (cz_ * cz_);
    double scale = std::sqrt(std::max(0.0, 1.0 - frac_sq));
    double max_x = ax_ * scale;
    double max_y = by_ * scale;
    return AABB(
        math::Vector3D(-max_x, -max_y, zcut1_),
        math::Vector3D( max_x,  max_y, zcut2_)
    );
}

} // namespace geometry
} // namespace siren
