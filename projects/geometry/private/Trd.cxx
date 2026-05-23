#include "SIREN/geometry/Trd.h"

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

Trd::Trd() : Geometry("Trd"), dx1_(0), dx2_(0), dy1_(0), dy2_(0), dz_(0) { RecomputeWorldAABB(); }
Trd::Trd(double dx1, double dx2, double dy1, double dy2, double dz) : Geometry("Trd"), dx1_(dx1), dx2_(dx2), dy1_(dy1), dy2_(dy2), dz_(dz) {
    if(dz_ <= 0) throw std::invalid_argument("Trd half-height dz must be positive!");
    if(dx1_ < 0 || dx2_ < 0 || dy1_ < 0 || dy2_ < 0) throw std::invalid_argument("Trd half-widths (dx1, dx2, dy1, dy2) must be non-negative!");
    RecomputeWorldAABB();
}
Trd::Trd(Placement const & p) : Geometry("Trd", p), dx1_(0), dx2_(0), dy1_(0), dy2_(0), dz_(0) { RecomputeWorldAABB(); }
Trd::Trd(Placement const & p, double dx1, double dx2, double dy1, double dy2, double dz) : Geometry("Trd", p), dx1_(dx1), dx2_(dx2), dy1_(dy1), dy2_(dy2), dz_(dz) {
    if(dz_ <= 0) throw std::invalid_argument("Trd half-height dz must be positive!");
    if(dx1_ < 0 || dx2_ < 0 || dy1_ < 0 || dy2_ < 0) throw std::invalid_argument("Trd half-widths (dx1, dx2, dy1, dy2) must be non-negative!");
    RecomputeWorldAABB();
}
Trd::Trd(const Trd& o) : Geometry(o), dx1_(o.dx1_), dx2_(o.dx2_), dy1_(o.dy1_), dy2_(o.dy2_), dz_(o.dz_) { RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(Trd, dx1_, dx2_, dy1_, dy2_, dz_)
SIREN_GEOMETRY_ASSIGN(Trd)
SIREN_GEOMETRY_EQUAL(Trd, dx1_, dx2_, dy1_, dy2_, dz_)
SIREN_GEOMETRY_LESS(Trd, dx1_, dx2_, dy1_, dy2_, dz_)

void Trd::print(std::ostream& os) const {
    os << "Trd(" << dx1_ << ", " << dx2_ << ", " << dy1_ << ", " << dy2_ << ", " << dz_ << ")\n";
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Trd::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Calculate intersection of ray with the trapezoid (frustum of rectangular pyramid).
    // The Trd has 6 faces:
    //   Top face:    z = +dz, rectangle [-dx2, dx2] x [-dy2, dy2]
    //   Bottom face: z = -dz, rectangle [-dx1, dx1] x [-dy1, dy1]
    //   4 side faces: planar trapezoids
    //
    // At any z in [-dz, +dz], the cross-section is a rectangle with half-widths:
    //   x_half(z) = dx1 + (dx2 - dx1) * (z + dz) / (2*dz)
    //   y_half(z) = dy1 + (dy2 - dy1) * (z + dz) / (2*dz)

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    double t;
    double ix, iy, iz; // intersection point
    bool entering;

    // Max 6 intersections: 6 faces
    Intersection hits[6];
    int n_hits = 0;

    // Helper: at a given z-level, compute the half-widths of the Trd cross-section
    // x_half(z) = dx1_ + (dx2_ - dx1_) * (z + dz_) / (2*dz_)
    // y_half(z) = dy1_ + (dy2_ - dy1_) * (z + dz_) / (2*dz_)
    double inv2dz = (dz_ > 0.0) ? 1.0 / (2.0 * dz_) : 0.0;

    // --- Top face: z = +dz_ ---
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
        t = (dz_ - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        ix = px + t * dx;
        iy = py + t * dy;
        if(ix >= -dx2_ && ix <= dx2_ && iy >= -dy2_ && iy <= dy2_) {
            iz = pz + t * dz;
            entering = direction.GetZ() < 0;
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = entering;
            hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
            n_hits++;
        }
    }

    // --- Bottom face: z = -dz_ ---
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
        t = (-dz_ - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        ix = px + t * dx;
        iy = py + t * dy;
        if(ix >= -dx1_ && ix <= dx1_ && iy >= -dy1_ && iy <= dy1_) {
            iz = pz + t * dz;
            entering = direction.GetZ() > 0;
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = entering;
            hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
            n_hits++;
        }
    }

    // --- Side faces ---
    // The 4 side faces are planar. Each side face connects a top edge to a
    // bottom edge. We derive the plane equation for each face and intersect
    // the ray with that plane, then check bounds.
    //
    // +x face: passes through (dx1_, y, -dz_) and (dx2_, y, +dz_) for all valid y.
    //   Two points on this face:
    //     P1 = (dx1_, 0, -dz_), P2 = (dx2_, 0, +dz_)
    //   The face also extends in the y-direction, so a third point is:
    //     P3 = (dx1_, 1, -dz_)
    //   Edge vectors: v1 = P2-P1 = (dx2_-dx1_, 0, 2*dz_), v2 = P3-P1 = (0, 1, 0)
    //   Normal = v2 x v1 = (1*2*dz_ - 0, 0 - 0*(dx2_-dx1_), 0 - 1*0)
    //          = (2*dz_, 0, -(dx2_-dx1_))
    //   Unnormalized outward normal for +x face: n = (2*dz_, 0, dx1_-dx2_)
    //   Plane equation: n . (r - P1) = 0
    //     2*dz_*(x - dx1_) + (dx1_ - dx2_)*(z + dz_) = 0

    // +x face
    {
        double nx = 2.0 * dz_;
        double nz = dx1_ - dx2_;
        double denom = nx * dx + nz * dz; // n . direction
        if(std::fabs(denom) > GEOMETRY_PRECISION) {
            // n . (P1 - position) = nx*(dx1_ - px) + nz*(-dz_ - pz)
            double num = nx * (dx1_ - px) + nz * (-dz_ - pz);
            t = num / denom;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            iz = pz + t * dz;
            if(iz >= -dz_ && iz <= dz_) {
                ix = px + t * dx;
                // Check y bounds at this z
                double yh = dy1_ + (dy2_ - dy1_) * (iz + dz_) * inv2dz;
                iy = py + t * dy;
                if(iy >= -yh && iy <= yh) {
                    entering = denom < 0;
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = entering;
                    hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
                    n_hits++;
                }
            }
        }
    }

    // -x face: outward normal is (-2*dz_, 0, dx1_-dx2_), point on face (-dx1_, 0, -dz_)
    {
        double nx = -2.0 * dz_;
        double nz = dx1_ - dx2_;
        double denom = nx * dx + nz * dz;
        if(std::fabs(denom) > GEOMETRY_PRECISION) {
            double num = nx * (-dx1_ - px) + nz * (-dz_ - pz);
            t = num / denom;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            iz = pz + t * dz;
            if(iz >= -dz_ && iz <= dz_) {
                ix = px + t * dx;
                double yh = dy1_ + (dy2_ - dy1_) * (iz + dz_) * inv2dz;
                iy = py + t * dy;
                if(iy >= -yh && iy <= yh) {
                    entering = denom < 0;
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = entering;
                    hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
                    n_hits++;
                }
            }
        }
    }

    // +y face: passes through (x, dy1_, -dz_) and (x, dy2_, +dz_).
    //   P1 = (0, dy1_, -dz_), P2 = (0, dy2_, +dz_), P3 = (1, dy1_, -dz_)
    //   v1 = P2-P1 = (0, dy2_-dy1_, 2*dz_), v2 = P3-P1 = (1, 0, 0)
    //   Normal = v1 x v2 = (0*0 - 2*dz_*0, 2*dz_*1 - 0*0, 0*0 - (dy2_-dy1_)*1)
    //          = (0, 2*dz_, dy1_-dy2_)
    //   Outward normal for +y face: n = (0, 2*dz_, dy1_-dy2_)
    {
        double ny = 2.0 * dz_;
        double nz = dy1_ - dy2_;
        double denom = ny * dy + nz * dz;
        if(std::fabs(denom) > GEOMETRY_PRECISION) {
            double num = ny * (dy1_ - py) + nz * (-dz_ - pz);
            t = num / denom;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            iz = pz + t * dz;
            if(iz >= -dz_ && iz <= dz_) {
                iy = py + t * dy;
                double xh = dx1_ + (dx2_ - dx1_) * (iz + dz_) * inv2dz;
                ix = px + t * dx;
                if(ix >= -xh && ix <= xh) {
                    entering = denom < 0;
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = entering;
                    hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
                    n_hits++;
                }
            }
        }
    }

    // -y face: outward normal is (0, -2*dz_, dy1_-dy2_), point on face (0, -dy1_, -dz_)
    {
        double ny = -2.0 * dz_;
        double nz = dy1_ - dy2_;
        double denom = ny * dy + nz * dz;
        if(std::fabs(denom) > GEOMETRY_PRECISION) {
            double num = ny * (-dy1_ - py) + nz * (-dz_ - pz);
            t = num / denom;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            iz = pz + t * dz;
            if(iz >= -dz_ && iz <= dz_) {
                iy = py + t * dy;
                double xh = dx1_ + (dx2_ - dx1_) * (iz + dz_) * inv2dz;
                ix = px + t * dx;
                if(ix >= -xh && ix <= xh) {
                    entering = denom < 0;
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = entering;
                    hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
                    n_hits++;
                }
            }
        }
    }

    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });

    // Net-parity dedup: collapse coincident hits from shared face boundaries.
    std::vector<Intersection> result;
    int i = 0;
    while(i < n_hits) {
        int j = i + 1;
        while(j < n_hits && std::fabs(hits[j].distance - hits[i].distance) <= GEOMETRY_PRECISION)
            ++j;
        int net = 0;
        for(int k = i; k < j; ++k)
            net += hits[k].entering ? 1 : -1;
        if(net > 0) {
            hits[i].entering = true;
            result.push_back(hits[i]);
        } else if(net < 0) {
            hits[i].entering = false;
            result.push_back(hits[i]);
        }
        i = j;
    }
    return result;
}

// ------------------------------------------------------------------------- //
AABB Trd::GetBoundingBox() const {
    double max_dx = std::max(dx1_, dx2_);
    double max_dy = std::max(dy1_, dy2_);
    return AABB(
        math::Vector3D(-max_dx, -max_dy, -dz_),
        math::Vector3D( max_dx,  max_dy,  dz_)
    );
}

} // namespace geometry
} // namespace siren
