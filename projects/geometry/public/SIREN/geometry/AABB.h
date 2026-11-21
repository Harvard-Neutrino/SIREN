#pragma once
#ifndef SIREN_AABB_H
#define SIREN_AABB_H

#include <algorithm>
#include <limits>
#include <cmath>
#include <array>

#include "SIREN/math/Vector3D.h"

namespace siren {
namespace geometry {

// Axis-Aligned Bounding Box
struct AABB {
    math::Vector3D min_corner;
    math::Vector3D max_corner;

    AABB()
        : min_corner(std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max())
        , max_corner(-std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()) {}

    AABB(math::Vector3D const & min_c, math::Vector3D const & max_c)
        : min_corner(min_c), max_corner(max_c) {}

    // Expand this AABB to include a point
    void ExpandToInclude(math::Vector3D const & p) {
        min_corner.SetCartesianCoordinates(
            std::min(min_corner.GetX(), p.GetX()),
            std::min(min_corner.GetY(), p.GetY()),
            std::min(min_corner.GetZ(), p.GetZ())
        );
        max_corner.SetCartesianCoordinates(
            std::max(max_corner.GetX(), p.GetX()),
            std::max(max_corner.GetY(), p.GetY()),
            std::max(max_corner.GetZ(), p.GetZ())
        );
    }

    // Expand this AABB to include another AABB
    void ExpandToInclude(AABB const & other) {
        ExpandToInclude(other.min_corner);
        ExpandToInclude(other.max_corner);
    }

    // Compute the centroid of the AABB
    math::Vector3D Centroid() const {
        return math::Vector3D(
            0.5 * (min_corner.GetX() + max_corner.GetX()),
            0.5 * (min_corner.GetY() + max_corner.GetY()),
            0.5 * (min_corner.GetZ() + max_corner.GetZ())
        );
    }

    // Compute surface area (used for SAH if needed)
    double SurfaceArea() const {
        double dx = max_corner.GetX() - min_corner.GetX();
        double dy = max_corner.GetY() - min_corner.GetY();
        double dz = max_corner.GetZ() - min_corner.GetZ();
        return 2.0 * (dx * dy + dy * dz + dz * dx);
    }

    // Return the index of the largest axis (0=x, 1=y, 2=z)
    int LargestAxis() const {
        double dx = max_corner.GetX() - min_corner.GetX();
        double dy = max_corner.GetY() - min_corner.GetY();
        double dz = max_corner.GetZ() - min_corner.GetZ();
        if(dx >= dy && dx >= dz) return 0;
        if(dy >= dz) return 1;
        return 2;
    }

    // Get extent along a given axis (0=x, 1=y, 2=z)
    double GetExtent(int axis) const {
        switch(axis) {
            case 0: return max_corner.GetX() - min_corner.GetX();
            case 1: return max_corner.GetY() - min_corner.GetY();
            case 2: return max_corner.GetZ() - min_corner.GetZ();
            default: return 0.0;
        }
    }

    // Get centroid coordinate along a given axis
    double GetCentroidAxis(int axis) const {
        switch(axis) {
            case 0: return 0.5 * (min_corner.GetX() + max_corner.GetX());
            case 1: return 0.5 * (min_corner.GetY() + max_corner.GetY());
            case 2: return 0.5 * (min_corner.GetZ() + max_corner.GetZ());
            default: return 0.0;
        }
    }

    // Check if the AABB is valid (min <= max on all axes)
    bool IsValid() const {
        return min_corner.GetX() <= max_corner.GetX()
            && min_corner.GetY() <= max_corner.GetY()
            && min_corner.GetZ() <= max_corner.GetZ();
    }

    // Check if a point is inside this AABB
    bool Contains(math::Vector3D const & p) const {
        return p.GetX() >= min_corner.GetX() && p.GetX() <= max_corner.GetX()
            && p.GetY() >= min_corner.GetY() && p.GetY() <= max_corner.GetY()
            && p.GetZ() >= min_corner.GetZ() && p.GetZ() <= max_corner.GetZ();
    }
};

// Slab-method ray-AABB test. inv_direction = 1.0/direction (precomputed).
// tmin/tmax are set to entry/exit distances along the ray.
// When a direction component is zero, inv_direction is +/-inf. If the ray
// origin sits exactly on a slab boundary, 0*inf = NaN, which collapses
// the t-interval via fmin/fmax. We guard each axis with an isinf check:
// a parallel ray inside the slab imposes no constraint; outside, it misses.
inline bool RayAABBIntersect(
    math::Vector3D const & origin,
    math::Vector3D const & inv_direction,
    AABB const & box,
    double & tmin,
    double & tmax) {
    tmin = -std::numeric_limits<double>::infinity();
    tmax =  std::numeric_limits<double>::infinity();

    double lo, hi, t1, t2;

    lo = box.min_corner.GetX() - origin.GetX();
    hi = box.max_corner.GetX() - origin.GetX();
    if(std::isinf(inv_direction.GetX())) {
        if(lo > 0 || hi < 0) return false;
    } else {
        t1 = lo * inv_direction.GetX();
        t2 = hi * inv_direction.GetX();
        tmin = std::fmax(tmin, std::fmin(t1, t2));
        tmax = std::fmin(tmax, std::fmax(t1, t2));
    }

    lo = box.min_corner.GetY() - origin.GetY();
    hi = box.max_corner.GetY() - origin.GetY();
    if(std::isinf(inv_direction.GetY())) {
        if(lo > 0 || hi < 0) return false;
    } else {
        t1 = lo * inv_direction.GetY();
        t2 = hi * inv_direction.GetY();
        tmin = std::fmax(tmin, std::fmin(t1, t2));
        tmax = std::fmin(tmax, std::fmax(t1, t2));
    }

    lo = box.min_corner.GetZ() - origin.GetZ();
    hi = box.max_corner.GetZ() - origin.GetZ();
    if(std::isinf(inv_direction.GetZ())) {
        if(lo > 0 || hi < 0) return false;
    } else {
        t1 = lo * inv_direction.GetZ();
        t2 = hi * inv_direction.GetZ();
        tmin = std::fmax(tmin, std::fmin(t1, t2));
        tmax = std::fmin(tmax, std::fmax(t1, t2));
    }

    // Do NOT check tmax >= 0: SIREN uses the full line (both directions)
    // for intersection queries, including negative-distance hits needed
    // by SectorLoop for containment resolution.
    return tmax >= tmin;
}

// Simplified overload: returns true/false without tmin/tmax output
inline bool RayAABBIntersect(
    math::Vector3D const & origin,
    math::Vector3D const & inv_direction,
    AABB const & box) {
    double tmin, tmax;
    return RayAABBIntersect(origin, inv_direction, box, tmin, tmax);
}

} // namespace geometry
} // namespace siren

#endif // SIREN_AABB_H
