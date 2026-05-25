#pragma once
#ifndef SIREN_DetectorDirectedChannelUtils_H
#define SIREN_DetectorDirectedChannelUtils_H

#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <memory>
#include <utility>
#include <algorithm>

namespace siren {
namespace injection {
namespace detail {

static const double kTwoPi = 2.0 * M_PI;
static const double kFourPi = 4.0 * M_PI;

inline double AABBVolume(siren::geometry::Geometry const & geo) {
    auto aabb = geo.GetWorldBoundingBox();
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    return extent.GetX() * extent.GetY() * extent.GetZ();
}

inline double ComputeGeometryVolume(siren::geometry::Geometry const & geo, double aabb_volume) {
    if (auto const * cyl = dynamic_cast<siren::geometry::Cylinder const *>(&geo)) {
        double R = cyl->GetRadius();
        double r = cyl->GetInnerRadius();
        double H = cyl->GetZ();
        return M_PI * (R * R - r * r) * H;
    }
    if (auto const * sph = dynamic_cast<siren::geometry::Sphere const *>(&geo)) {
        double R = sph->GetRadius();
        double r = sph->GetInnerRadius();
        return (4.0 / 3.0) * M_PI * (R * R * R - r * r * r);
    }
    if (auto const * box = dynamic_cast<siren::geometry::Box const *>(&geo)) {
        return box->GetX() * box->GetY() * box->GetZ();
    }

    auto aabb = geo.GetWorldBoundingBox();
    int N = 10000;
    int n_inside = 0;
    for (int i = 0; i < N; ++i) {
        int ix = i;
        double fx = 0.0, fy = 0.0, fz = 0.0;
        double bx = 1.0 / 2.0, by = 1.0 / 3.0, bz = 1.0 / 5.0;
        for (int j = 0; j < 15; ++j) {
            fx += (ix % 2) * bx;
            bx /= 2.0;
            ix /= 2;
        }
        ix = i;
        for (int j = 0; j < 10; ++j) {
            fy += (ix % 3) * by;
            by /= 3.0;
            ix /= 3;
        }
        ix = i;
        for (int j = 0; j < 7; ++j) {
            fz += (ix % 5) * bz;
            bz /= 5.0;
            ix /= 5;
        }
        double x = aabb.min_corner.GetX() + fx * (aabb.max_corner.GetX() - aabb.min_corner.GetX());
        double y = aabb.min_corner.GetY() + fy * (aabb.max_corner.GetY() - aabb.min_corner.GetY());
        double z = aabb.min_corner.GetZ() + fz * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ());
        if (geo.IsInside(siren::math::Vector3D(x, y, z))) {
            ++n_inside;
        }
    }
    if (n_inside == 0) return aabb_volume;
    return aabb_volume * static_cast<double>(n_inside) / N;
}

inline double ConeSolidAngle(
    siren::geometry::Geometry const & target,
    siren::math::Vector3D const & position)
{
    auto aabb = target.GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D offset = center - position;
    double dist = offset.magnitude();
    if (dist < bounding_radius) return kFourPi;

    double sin_theta = bounding_radius / dist;
    double cos_theta = std::sqrt(std::max(0.0, 1.0 - sin_theta * sin_theta));
    return kTwoPi * (1.0 - cos_theta);
}

inline std::pair<siren::math::Vector3D, double> SampleConeDirection(
    siren::geometry::Geometry const & target,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position)
{
    auto aabb = target.GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D to_center = center - position;
    double dist = to_center.magnitude();

    if (dist < bounding_radius) {
        double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
        double phi = kTwoPi * random->Uniform(0, 1);
        double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        return {
            siren::math::Vector3D(
                sin_theta * std::cos(phi),
                sin_theta * std::sin(phi),
                cos_theta),
            kFourPi
        };
    }

    double sin_cone = bounding_radius / dist;
    double cos_cone = std::sqrt(std::max(0.0, 1.0 - sin_cone * sin_cone));
    double solid_angle = kTwoPi * (1.0 - cos_cone);

    double cos_theta = cos_cone + (1.0 - cos_cone) * random->Uniform(0, 1);
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double phi = kTwoPi * random->Uniform(0, 1);

    siren::math::Vector3D axis = to_center;
    axis.normalize();

    siren::math::Vector3D perp1, perp2;
    if (std::abs(axis.GetX()) < 0.9)
        perp1 = siren::math::Vector3D(1, 0, 0);
    else
        perp1 = siren::math::Vector3D(0, 1, 0);
    perp2 = siren::math::vector_product(axis, perp1);
    perp2.normalize();
    perp1 = siren::math::vector_product(perp2, axis);
    perp1.normalize();

    siren::math::Vector3D direction = axis * cos_theta
        + perp1 * (sin_theta * std::cos(phi))
        + perp2 * (sin_theta * std::sin(phi));
    direction.normalize();

    return {direction, solid_angle};
}

inline bool DirectionHitsTarget(
    siren::geometry::Geometry const & target,
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction)
{
    auto intersections = target.Intersections(position, direction);
    for (auto const & isec : intersections) {
        if (isec.distance > 0 && isec.entering) return true;
    }
    return false;
}

inline siren::math::Vector3D SampleVolumePoint(
    siren::geometry::Geometry const & target,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    auto aabb = target.GetWorldBoundingBox();
    for (int attempt = 0; attempt < 10000; ++attempt) {
        double x = aabb.min_corner.GetX()
            + random->Uniform(0, 1) * (aabb.max_corner.GetX() - aabb.min_corner.GetX());
        double y = aabb.min_corner.GetY()
            + random->Uniform(0, 1) * (aabb.max_corner.GetY() - aabb.min_corner.GetY());
        double z = aabb.min_corner.GetZ()
            + random->Uniform(0, 1) * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ());
        siren::math::Vector3D point(x, y, z);
        if (target.IsInside(point)) return point;
    }
    throw std::runtime_error("Failed to sample a point inside target volume after 10000 attempts");
}

inline double SolidAngleDensity(
    siren::geometry::Geometry const & target,
    double target_volume,
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction)
{
    auto intersections = target.Intersections(position, direction);
    std::vector<siren::geometry::Geometry::Intersection> pos_intersections;
    for (auto const & isec : intersections) {
        if (isec.distance > 0) {
            pos_intersections.push_back(isec);
        }
    }
    std::sort(pos_intersections.begin(), pos_intersections.end(),
              [](auto const & a, auto const & b) { return a.distance < b.distance; });

    double r2_integral = 0.0;
    if (!pos_intersections.empty()) {
        size_t start_idx = 0;
        if (!pos_intersections[0].entering) {
            // Origin is inside the volume.
            // First segment is from 0 to the first exit.
            double r_exit = pos_intersections[0].distance;
            r2_integral += (r_exit * r_exit * r_exit) / 3.0;
            start_idx = 1;
        }

        for (size_t i = start_idx; i < pos_intersections.size(); ++i) {
            if (pos_intersections[i].entering) {
                double r_enter = pos_intersections[i].distance;
                double r_exit = r_enter;
                for (size_t j = i + 1; j < pos_intersections.size(); ++j) {
                    if (!pos_intersections[j].entering) {
                        r_exit = pos_intersections[j].distance;
                        break;
                    }
                }
                if (r_exit > r_enter) {
                    r2_integral += (r_exit * r_exit * r_exit
                                  - r_enter * r_enter * r_enter) / 3.0;
                }
            }
        }
    }

    if (r2_integral <= 0.0 || target_volume <= 0.0) return 0.0;
    return r2_integral / target_volume;
}

inline siren::math::Vector3D SampleDirectedDirection(
    siren::geometry::Geometry const & target,
    DetectorDirected2BodyChannel::Mode mode,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position)
{
    if (mode == DetectorDirected2BodyChannel::Mode::Cone) {
        return SampleConeDirection(target, random, position).first;
    }
    siren::math::Vector3D target_point;
    siren::math::Vector3D direction;
    double diff_mag = 0.0;
    for (int attempt = 0; attempt < 100; ++attempt) {
        target_point = SampleVolumePoint(target, random);
        direction = target_point - position;
        diff_mag = direction.magnitude();
        if (diff_mag > 1e-6) {
            break;
        }
    }
    if (diff_mag <= 1e-6) {
        throw std::runtime_error("Failed to sample a target point sufficiently separated from decay vertex after 100 attempts");
    }
    direction.normalize();
    return direction;
}

inline double AngularDensity(
    siren::geometry::Geometry const & target,
    double target_volume,
    DetectorDirected2BodyChannel::Mode mode,
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction)
{
    if (mode == DetectorDirected2BodyChannel::Mode::Cone) {
        if (!DirectionHitsTarget(target, position, direction)) return 0.0;
        return 1.0 / ConeSolidAngle(target, position);
    }
    return SolidAngleDensity(target, target_volume, position, direction);
}

} // namespace detail
} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirectedChannelUtils_H
