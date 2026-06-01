#pragma once
#ifndef SIREN_DetectorDirected2BodyChannel_H
#define SIREN_DetectorDirected2BodyChannel_H

#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/math/Vector3D.h"

#include <memory>
#include <string>
#include <utility>

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

// Biased 2-body decay channel that forces one daughter toward a
// target geometry (typically the detector fiducial volume).
//
// Two sampling modes:
//
// CONE mode (simple):
//   Samples uniformly on a bounding cone around the target.
//   Density = 1/Omega_cone for directions inside the cone, 0 outside.
//   Directions that miss the actual geometry get density 0; the
//   multi-channel isotropic term handles those events.
//
// VOLUME mode (accurate):
//   Samples a point uniformly in the target volume. The direction
//   from the decay position to that point defines the lab direction.
//   Density = (1/V) * integral(r^2 dr) over intersection segments,
//   which is the exact angular density of the volume-uniform
//   distribution.  Every sampled direction is guaranteed to hit the
//   target by construction.
//
// In both modes, after choosing a lab direction, the 2-body
// kinematics are solved to find the rest-frame angle(s) and the
// Jacobian |dOmega_lab/dOmega_rest|.  If the direction is
// kinematically forbidden (beyond the critical angle for massive
// daughters), the event falls back to isotropic sampling and the
// Density at that point returns 0 from the biased channel.
class DetectorDirected2BodyChannel : public PhaseSpaceChannel {
public:
    enum class Mode { Cone, Volume };

    // If `volume > 0`, it is used directly as the target volume for the
    // Volume-mode chord-depth density, skipping the internal Monte-Carlo
    // volume estimate and the viability guard.  This lets a caller supply an
    // accurate volume for a composite/thin tile (e.g. a BooleanGeometry
    // subtraction shell) whose AABB-rejection fill fraction is small.
    DetectorDirected2BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int daughter_index = 0,
        Mode mode = Mode::Volume,
        double volume = -1.0
    );

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override { return "DetectorDirected2Body"; }
    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Decay2Body;
    }
    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::SolidAngleRest();
    }

    // Set the true volume of the target geometry (for Volume mode).
    // If not called, the AABB volume is used as an approximation.
    void SetVolume(double volume);

    // True if this channel actually directs (a non-isotropic proposal) at the
    // phase-space point in `record`, vs falling back to isotropic 1/4pi
    // (Disjoint / KinInBound / parent-at-rest-inside).  Lets a diagnostic
    // attribute the channel's variance into directing vs the shared fallback.
    bool DirectingActive(
        siren::dataclasses::InteractionRecord const & record) const;

private:
    std::shared_ptr<siren::geometry::Geometry const> target_;
    int daughter_index_;
    Mode mode_;
    double aabb_volume_;
    double target_volume_;

    // ---- Cone mode helpers ----

    // Sample a direction uniformly on a bounding cone.
    // Returns (direction, cone_solid_angle).
    std::pair<siren::math::Vector3D, double> SampleConeDirection(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        siren::math::Vector3D const & position
    ) const;

    double ConeSolidAngle(siren::math::Vector3D const & position) const;

    bool DirectionHitsTarget(
        siren::math::Vector3D const & position,
        siren::math::Vector3D const & direction
    ) const;

    // ---- Volume mode helpers ----

    // Sample a point uniformly inside the target volume (AABB rejection).
    siren::math::Vector3D SampleVolumePoint(
        std::shared_ptr<siren::utilities::SIREN_random> random
    ) const;

    // Exact solid angle density at a direction: (1/V) * sum(r_exit^3 - r_enter^3)/3
    double SolidAngleDensity(
        siren::math::Vector3D const & position,
        siren::math::Vector3D const & direction
    ) const;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirected2BodyChannel_H
