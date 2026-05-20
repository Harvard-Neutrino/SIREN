// Shared phi-range utilities for shapes with azimuthal angular extent.
// Included by Sphere, Cylinder, Cone, CutTube, Polycone, GenericPolycone,
// Polyhedra, Torus.

#pragma once
#ifndef SIREN_PhiUtils_H
#define SIREN_PhiUtils_H

#include <cmath>

namespace siren {
namespace geometry {
namespace phi_utils {

inline constexpr double TWO_PI = 2.0 * M_PI;

// Normalize angle to [0, 2*pi)
inline double NormalizePhi(double phi) {
    phi = std::fmod(phi, TWO_PI);
    if(phi < 0) phi += TWO_PI;
    return phi;
}

// Check if the azimuthal angle of point (x,y) falls within
// [start_phi, start_phi + delta_phi]. Handles wraparound.
inline bool PhiInRange(double x, double y, double start_phi, double delta_phi) {
    double phi = NormalizePhi(std::atan2(y, x));
    double sp = NormalizePhi(start_phi);
    double ep = sp + delta_phi;
    if(ep <= TWO_PI + 1e-9) {
        return phi >= sp - 1e-9 && phi <= ep + 1e-9;
    } else {
        return phi >= sp - 1e-9 || phi <= NormalizePhi(ep) + 1e-9;
    }
}

// Determine the initial phi-wedge state at t=-infinity for a ray.
// At t=-infinity the ray position is far from the z-axis (when the ray
// has any transverse component), so phi is well-defined there:
//   phi(-inf) = atan2(-dy, -dx)
// This avoids the atan2 singularity at (0,0) that plagues rays
// passing through or near the z-axis.
inline bool InitialPhiState(double px, double py, double dx, double dy,
                            double start_phi, double delta_phi) {
    double rxy2_dir = dx*dx + dy*dy;
    if(rxy2_dir > 1e-20) {
        return PhiInRange(-dx, -dy, start_phi, delta_phi);
    }
    double rxy2_orig = px*px + py*py;
    if(rxy2_orig > 1e-20) {
        return PhiInRange(px, py, start_phi, delta_phi);
    }
    return delta_phi >= M_PI - 1e-9;
}

// Correct entering flag for a wedge hit at the z-axis.
// Both phi-cut half-planes pass through the z-axis, so two wedge
// hits coincide there with potentially wrong entering flags.
// After crossing the z-axis the ray moves in direction (dx,dy),
// so the post-crossing phi is atan2(dy, dx).
inline bool ZAxisWedgeEntering(double dx, double dy,
                               double start_phi, double delta_phi) {
    if(dx*dx + dy*dy > 1e-20) {
        return PhiInRange(dx, dy, start_phi, delta_phi);
    }
    return delta_phi >= M_PI - 1e-9;
}

} // namespace phi_utils
} // namespace geometry
} // namespace siren

#endif // SIREN_PhiUtils_H
