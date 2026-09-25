#pragma once
#ifndef SIREN_Injection_GeometryVolume_H
#define SIREN_Injection_GeometryVolume_H

#include <string>

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

// Exact volume for geometry primitives supported by directed channels: box,
// cylinder, sphere, ellipsoid, elliptical tube, cone, torus, Trd,
// parallelepiped, polycone and generic polycone with a simple R-Z profile
// (checked with rounding-aware orientations; the moment is taken with z
// measured from the profile), including their cuts and azimuthal segments,
// and a closed triangle mesh of one connected piece whose edge-connected
// shells are all outward (divergence theorem about its box centre, from edge
// vectors, compensated). Returns NaN when no exact implementation is available.
double ExactGeometryVolume(siren::geometry::Geometry const & geometry);

// Independent estimate of any solid's volume by chord integration: along each
// of 48 fixed directions spread over a hemisphere (skew to the axes of the
// solid's own frame, where its box is tight), a jittered 48 x 48 grid of rays
// over the smallest rectangle around the box's projection,
// each contributing the length of its path inside the solid (from the solid's
// Intersections). The volume is the mean over directions and the standard
// error their scatter, so a direction parallel to a thin wall widens the error
// rather than biasing the result. The estimate is inconclusive (resolved =
// false, with the reason) when fewer than 20 rays meet the solid along more
// than 4 of the directions (a speck, an empty or an extremely thin solid seen
// edge-on from many directions); when a ray's crossings are not
// finite, or more than a few do not alternate entering/exiting (an open,
// doubled, nested or self-intersecting surface); when the solid's IsInside,
// which the channels sample with, disagrees with its crossings; when more than
// 1% of chords are within 5 GEOMETRY_PRECISION, where crossings are merged;
// when parts of the solid may have been missed by more than 1% of the result
// (a small core inside a thin shell, which the rays can miss entirely): each
// primitive with an exact volume inside the solid's box is estimated on its
// own and its unexplained shortfall counts, and any other part (a connected
// mesh piece, an intersection or subtraction node) met by fewer than 10 rays
// counts with its bounding box. Material that a Boolean combination of broad
// operands selects (by intersection or subtraction) inside a part the rays do
// meet, such as a small boss on a thin plate, can still go unseen: for Boolean
// solids the check is a safeguard, not a guarantee; or
// when the result is not finite, exceeds the bounding box, or has a standard
// error above 5%. Otherwise [lower, upper] is the accepted range, max(2%, 5
// standard errors) around the result. The ray sequence is fixed, so the
// estimate is reproducible.
struct GeometryVolumeEstimate {
    double volume = 0.0;
    double standard_error = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    long long points = 0;
    bool resolved = false;
    bool malformed = false;  // the reason is a surface that is not a solid's
    std::string reason;
};
GeometryVolumeEstimate EstimateGeometryVolume(siren::geometry::Geometry const & geometry);

// Resolve the normalization volume shared by every detector-directed channel.
// In volume mode a positive caller-supplied value must agree with the exact
// volume to 1e-4 where one exists (the exact value is then used; an exact
// volume above the bounding box is refused, and a mesh's only counts when the
// rays resolve the mesh and their estimate agrees with it); otherwise it
// must not exceed the bounding box (beyond the box's rounding) and must lie in
// EstimateGeometryVolume's accepted range (about 2% for well-resolved solids,
// wider for thin curved shells). A supplied volume for a solid the estimate
// cannot resolve is refused with the reason: its normalization cannot be
// checked. A supplied value bypasses the AABB fill fraction guard. Any volume
// is refused when the layer IsInside samples as outside (the last
// GEOMETRY_PRECISION before each exit along +z) exceeds 0.5% of it: the
// channel's samples would not fill the solid (a solid a few nm thick).
// Volume-mode channels without one require a supported exact primitive whose
// fill fraction is large enough for rejection sampling.
double ResolveDetectorDirectedVolume(
    siren::geometry::Geometry const & geometry,
    bool volume_mode,
    double supplied_volume = -1.0);

// Check a restored (archived) normalization volume by the same rules, without
// changing it: an archive whose volume the constructor would reject is refused
// instead of being reread as a different density. Only volume mode uses it.
void ValidateArchivedDetectorDirectedVolume(
    siren::geometry::Geometry const * geometry,
    bool volume_mode,
    double archived_volume);

} // namespace injection
} // namespace siren

#endif // SIREN_Injection_GeometryVolume_H
