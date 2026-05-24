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
// Sampling strategy:
//   1. Sample a lab-frame direction uniformly over the solid angle
//      subtended by the target geometry as seen from the decay point.
//   2. Solve the 2-body kinematics to find the rest-frame angle(s)
//      that produce this lab direction (0, 1, or 2 solutions for
//      massive daughters).
//   3. When 2 solutions exist, pick one with probability proportional
//      to the kinematic weight (Jacobian).
//   4. Construct the full final-state 4-momenta.
//
// Density at any phase-space point:
//   g(Omega_lab) * |dOmega_lab / dOmega_rest|
//
// where g is the lab-frame angular density (1/solid_angle for
// uniform sampling over the target) and the Jacobian converts
// between the rest-frame and lab-frame solid angle measures.
//
// For phase-space points where the lab direction does NOT intersect
// the target geometry, the density is 0.
class DetectorDirected2BodyChannel : public PhaseSpaceChannel {
public:
    // target: the geometry to direct the daughter toward (e.g.,
    //         the detector fiducial volume)
    // daughter_index: which secondary (0 or 1) to direct
    DetectorDirected2BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int daughter_index = 0
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

private:
    std::shared_ptr<siren::geometry::Geometry const> target_;
    int daughter_index_;

    // Sample a direction uniformly over the solid angle subtended
    // by the target geometry as seen from the given position.
    // Returns (direction, solid_angle).
    std::pair<siren::math::Vector3D, double> SampleTargetDirection(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        siren::math::Vector3D const & position
    ) const;

    // Compute the solid angle subtended by the target geometry
    // as seen from the given position.
    double TargetSolidAngle(
        siren::math::Vector3D const & position
    ) const;

    // Check whether a ray from position in the given direction
    // intersects the target geometry.
    bool DirectionHitsTarget(
        siren::math::Vector3D const & position,
        siren::math::Vector3D const & direction
    ) const;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirected2BodyChannel_H
