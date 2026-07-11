#pragma once
#ifndef SIREN_Injection_GeometryVolume_H
#define SIREN_Injection_GeometryVolume_H

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

// Exact volume for geometry primitives supported by directed channels.
// Returns NaN when no exact analytic implementation is available.
double ExactGeometryVolume(siren::geometry::Geometry const & geometry);

// Resolve the normalization volume shared by every detector-directed channel.
// A positive caller-supplied value is authoritative and bypasses the AABB fill
// fraction guard. Volume-mode channels otherwise require a supported exact
// primitive whose fill fraction is large enough for rejection sampling.
double ResolveDetectorDirectedVolume(
    siren::geometry::Geometry const & geometry,
    bool volume_mode,
    double supplied_volume = -1.0);

} // namespace injection
} // namespace siren

#endif // SIREN_Injection_GeometryVolume_H
