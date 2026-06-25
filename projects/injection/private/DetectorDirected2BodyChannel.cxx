#include "SIREN/injection/DetectorDirected2BodyChannel.h"

#include "DetectorDirectedChannelUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace siren {
namespace injection {

static const double TWO_PI = 2.0 * M_PI;
static const double FOUR_PI = 4.0 * M_PI;

DetectorDirected2BodyChannel::DetectorDirected2BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int daughter_index,
    Mode mode,
    double volume)
    : target_(std::move(target))
    , daughter_index_(daughter_index)
    , mode_(mode)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected2BodyChannel requires a target geometry");
    }
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    aabb_volume_ = extent.GetX() * extent.GetY() * extent.GetZ();

    if (volume > 0.0) {
        // Caller-supplied volume: trust it, skip the MC estimate and the
        // viability guard (the caller asserts the tile is samplable).
        target_volume_ = volume;
        return;
    }

    target_volume_ = detail::ComputeGeometryVolume(*target_, aabb_volume_);

    if (mode_ == Mode::Volume) {
        if (aabb_volume_ <= 0.0) {
            throw std::runtime_error("Target bounding box has zero volume");
        }
        if (target_volume_ / aabb_volume_ < 1e-4) {
            throw std::runtime_error("Target volume is too small relative to its bounding box for sampling to be viable");
        }
    }
}

void DetectorDirected2BodyChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

bool DetectorDirected2BodyChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const
{
    if (record.signature.secondary_types.size() != 2) return false;
    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);
    auto geo = detail::ClassifyDirectedRegime(
        record.primary_momentum[0], record.primary_momentum[1],
        record.primary_momentum[2], record.primary_momentum[3],
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_);
    return detail::IsDirectedStepActive(geo.regime, geo.inside_geometry);
}

// ================================================================ //
//  Cone mode helpers                                                 //
// ================================================================ //

double DetectorDirected2BodyChannel::ConeSolidAngle(
    siren::math::Vector3D const & position) const
{
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D offset = center - position;
    double dist = offset.magnitude();

    if (dist < bounding_radius) return FOUR_PI;

    double sin_theta = bounding_radius / dist;
    double cos_theta = std::sqrt(1.0 - sin_theta * sin_theta);
    return TWO_PI * (1.0 - cos_theta);
}

std::pair<siren::math::Vector3D, double>
DetectorDirected2BodyChannel::SampleConeDirection(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position) const
{
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D to_center = center - position;
    double dist = to_center.magnitude();

    if (dist < bounding_radius) {
        double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
        double phi = TWO_PI * random->Uniform(0, 1);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        return {
            siren::math::Vector3D(
                sin_theta * std::cos(phi),
                sin_theta * std::sin(phi),
                cos_theta),
            FOUR_PI
        };
    }

    double sin_cone = bounding_radius / dist;
    double cos_cone = std::sqrt(1.0 - sin_cone * sin_cone);
    double solid_angle = TWO_PI * (1.0 - cos_cone);

    double cos_theta = cos_cone + (1.0 - cos_cone) * random->Uniform(0, 1);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = TWO_PI * random->Uniform(0, 1);

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

bool DetectorDirected2BodyChannel::DirectionHitsTarget(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const
{
    auto intersections = target_->Intersections(position, direction);
    for (auto const & isec : intersections) {
        if (isec.distance > 0 && isec.entering) return true;
    }
    return false;
}

// ================================================================ //
//  Volume mode helpers                                               //
// ================================================================ //

siren::math::Vector3D DetectorDirected2BodyChannel::SampleVolumePoint(
    std::shared_ptr<siren::utilities::SIREN_random> random) const
{
    auto aabb = target_->GetWorldBoundingBox();

    for (int attempt = 0; attempt < 10000; ++attempt) {
        double x = aabb.min_corner.GetX()
            + random->Uniform(0, 1) * (aabb.max_corner.GetX() - aabb.min_corner.GetX());
        double y = aabb.min_corner.GetY()
            + random->Uniform(0, 1) * (aabb.max_corner.GetY() - aabb.min_corner.GetY());
        double z = aabb.min_corner.GetZ()
            + random->Uniform(0, 1) * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ());
        siren::math::Vector3D point(x, y, z);
        if (target_->IsInside(point)) return point;
    }
    throw std::runtime_error("Failed to sample a point inside target volume after 10000 attempts");
}

double DetectorDirected2BodyChannel::SolidAngleDensity(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const
{
    auto intersections = target_->Intersections(position, direction);

    std::vector<siren::geometry::Geometry::Intersection> pos_intersections;
    for (auto const & isec : intersections) {
        if (isec.distance > 0) {
            pos_intersections.push_back(isec);
        }
    }
    std::sort(pos_intersections.begin(), pos_intersections.end(),
              [](auto const & a, auto const & b) { return a.distance < b.distance; });

    // g(omega) = (1/V) * integral of r^2 dr over intersection segments
    //          = (1/V) * sum_segments (r_exit^3 - r_enter^3) / 3
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

    if (r2_integral <= 0) return 0.0;
    return r2_integral / target_volume_;
}

// Regime classification, sampling, and density evaluation are in
// DetectorDirectedChannelUtils.h (detail::ClassifyDirectedRegime,
// detail::SampleDirectedStep, detail::DensityDirectedStep).

// ================================================================ //
//  Sample and Density — thin adapters over shared step functions    //
// ================================================================ //

void DetectorDirected2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 2) {
        throw std::runtime_error(
            "DetectorDirected2BodyChannel requires exactly 2 secondaries");
    }

    double m_A = record.secondary_masses[daughter_index_];
    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];

    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    auto result = detail::SampleDirectedStep(
        E_parent, px_parent, py_parent, pz_parent,
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_, target_volume_, mode_, random);

    record.secondary_momenta[daughter_index_] = {
        result.E_lab,
        result.p_lab * result.lab_dir.GetX(),
        result.p_lab * result.lab_dir.GetY(),
        result.p_lab * result.lab_dir.GetZ()
    };
    record.secondary_momenta[1 - daughter_index_] = {
        E_parent - result.E_lab,
        px_parent - result.p_lab * result.lab_dir.GetX(),
        py_parent - result.p_lab * result.lab_dir.GetY(),
        pz_parent - result.p_lab * result.lab_dir.GetZ()
    };
}

double DetectorDirected2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (record.signature.secondary_types.size() != 2) return 0.0;

    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    return detail::DensityDirectedStep(
        record.primary_momentum[0],
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3],
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        record.secondary_momenta[daughter_index_][0],
        record.secondary_momenta[daughter_index_][1],
        record.secondary_momenta[daughter_index_][2],
        record.secondary_momenta[daughter_index_][3],
        decay_pos, *target_, target_volume_, mode_);
}

} // namespace injection
} // namespace siren
