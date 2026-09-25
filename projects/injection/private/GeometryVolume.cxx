#include "SIREN/injection/GeometryVolume.h"

#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/math/Vector3D.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace siren {
namespace injection {

namespace {

constexpr double kTwoPi = 2.0 * M_PI;
constexpr double kMinimumViableFillFraction = 1e-4;

double AABBVolume(siren::geometry::Geometry const & geometry) {
    auto aabb = geometry.GetWorldBoundingBox();
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    return extent.GetX() * extent.GetY() * extent.GetZ();
}

} // namespace

double ExactGeometryVolume(siren::geometry::Geometry const & geometry) {
    if (auto const * cylinder =
            dynamic_cast<siren::geometry::Cylinder const *>(&geometry)) {
        double outer = cylinder->GetRadius();
        double inner = cylinder->GetInnerRadius();
        double delta_phi = std::clamp(
            cylinder->GetDeltaPhi(), 0.0, kTwoPi);
        return 0.5 * delta_phi
            * (outer * outer - inner * inner) * cylinder->GetZ();
    }
    if (auto const * sphere =
            dynamic_cast<siren::geometry::Sphere const *>(&geometry)) {
        double outer = sphere->GetRadius();
        double inner = sphere->GetInnerRadius();
        double delta_phi = std::clamp(sphere->GetDeltaPhi(), 0.0, kTwoPi);
        double theta_start = std::clamp(
            sphere->GetStartTheta(), 0.0, M_PI);
        double theta_end = std::clamp(
            theta_start + sphere->GetDeltaTheta(), theta_start, M_PI);
        double angular_integral = delta_phi
            * (std::cos(theta_start) - std::cos(theta_end));
        return (outer * outer * outer - inner * inner * inner)
            * angular_integral / 3.0;
    }
    if (auto const * box =
            dynamic_cast<siren::geometry::Box const *>(&geometry)) {
        return box->GetX() * box->GetY() * box->GetZ();
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double ResolveDetectorDirectedVolume(
    siren::geometry::Geometry const & geometry,
    bool volume_mode,
    double supplied_volume)
{
    double aabb_volume = AABBVolume(geometry);
    bool has_supplied_volume = supplied_volume > 0.0;
    double target_volume = has_supplied_volume
        ? supplied_volume
        : ExactGeometryVolume(geometry);

    if (!volume_mode) return target_volume;
    if (!(aabb_volume > 0.0) || !std::isfinite(aabb_volume)) {
        throw std::runtime_error("Target bounding box has zero volume");
    }
    if (!(target_volume > 0.0) || !std::isfinite(target_volume)) {
        throw std::runtime_error(
            "Volume mode requires an exact caller-supplied volume for this geometry");
    }
    if (!has_supplied_volume &&
        target_volume / aabb_volume <= kMinimumViableFillFraction) {
        throw std::runtime_error(
            "Target volume is too small relative to its bounding box for sampling to be viable");
    }
    return target_volume;
}

} // namespace injection
} // namespace siren
