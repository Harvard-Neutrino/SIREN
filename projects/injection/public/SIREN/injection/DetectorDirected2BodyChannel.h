#pragma once
#ifndef SIREN_DetectorDirected2BodyChannel_H
#define SIREN_DetectorDirected2BodyChannel_H

#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

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
class OnShellCascadeChannel;

class DetectorDirected2BodyChannel : public PhaseSpaceChannel {
public:
    enum class Mode { Cone, Volume };

    // If `volume > 0`, it is the target volume for the Volume-mode chord-depth
    // density, in place of the analytic volume and the viability guard. It is
    // checked first (GeometryVolume.h): against the analytic volume where one
    // exists, otherwise against a chord-integration estimate of the solid. This lets
    // a caller supply the volume of a composite/thin tile (e.g. a
    // BooleanGeometry subtraction shell) whose AABB fill fraction is small.
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

    // Set the true volume of the target geometry (for Volume mode). It is
    // validated like the constructor argument (see GeometryVolume.h).
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
    friend class cereal::access;
    friend class OnShellCascadeChannel;

    DetectorDirected2BodyChannel() = default;

    // A volume-mode component whose owner has already checked the volume, so
    // it is neither checked again nor replaced.
    struct CheckedVolume { double value; };
    DetectorDirected2BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int daughter_index,
        CheckedVolume volume);

    std::shared_ptr<siren::geometry::Geometry const> target_;
    int daughter_index_;
    Mode mode_;
    double target_volume_;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 1) {
            int mode = static_cast<int>(mode_);
            archive(::cereal::make_nvp("Target", target_));
            archive(::cereal::make_nvp("DaughterIndex", daughter_index_));
            archive(::cereal::make_nvp("Mode", mode));
            archive(::cereal::make_nvp("TargetVolume", target_volume_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "DetectorDirected2BodyChannel: legacy directed-decay density archive rejected; regenerate events with the current sampler");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 1) {
            int mode;
            archive(::cereal::make_nvp("Target", target_));
            archive(::cereal::make_nvp("DaughterIndex", daughter_index_));
            archive(::cereal::make_nvp("Mode", mode));
            archive(::cereal::make_nvp("TargetVolume", target_volume_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
            if(mode != static_cast<int>(Mode::Cone)
               && mode != static_cast<int>(Mode::Volume)) {
                throw std::runtime_error(
                    "DetectorDirected2BodyChannel: invalid Mode value "
                    + std::to_string(mode) + " in archive");
            }
            mode_ = static_cast<Mode>(mode);
            ValidateArchivedDetectorDirectedVolume(
                target_.get(), mode_ == Mode::Volume, target_volume_);
        } else {
            throw std::runtime_error(
                "DetectorDirected2BodyChannel: legacy directed-decay density archive rejected; regenerate events with the current sampler");
        }
    }
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::DetectorDirected2BodyChannel, 1);
CEREAL_REGISTER_TYPE(siren::injection::DetectorDirected2BodyChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::DetectorDirected2BodyChannel);

#endif // SIREN_DetectorDirected2BodyChannel_H
