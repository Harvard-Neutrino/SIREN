#pragma once
#ifndef SIREN_DetectorDirected2BodyChannel_H
#define SIREN_DetectorDirected2BodyChannel_H

#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>

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
//   The bounding cone includes some directions that can miss the exact
//   geometry; those directions retain the cone density because they are part
//   of the proposal support.
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
// Jacobian |dOmega_lab/dOmega_rest|. Regimes where directing provides no
// usable support fall back deterministically to isotropic sampling. If an
// active directed sampler exhausts its kinematic solutions, the injection
// attempt fails rather than emitting a sample absent from Density.
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
    // Prefer the constructor argument for composite geometries (whose volume
    // cannot be derived analytically) and very thin geometries so the
    // viability check can use it before construction completes.
    void SetVolume(double volume);

    // True if this channel actually directs (a non-isotropic proposal) at the
    // phase-space point in `record`, vs falling back to isotropic 1/4pi
    // (Disjoint / KinInBound / parent-at-rest-inside).  Lets a diagnostic
    // attribute the channel's variance into directing vs the shared fallback.
    bool DirectingActive(
        siren::dataclasses::InteractionRecord const & record) const override;

private:
    std::shared_ptr<siren::geometry::Geometry const> target_;
    int daughter_index_;
    Mode mode_;
    double target_volume_;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirected2BodyChannel_H
