#pragma once
#ifndef SIREN_DetectorDirected3BodyChannel_H
#define SIREN_DetectorDirected3BodyChannel_H

#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

// Recursive detector-directed 3-body decay channel.
//
// The channel factorizes P -> spectator + X and X -> a + b.
// It samples the intermediate invariant mass s_X, then points X at
// the target geometry, and finally points one daughter of X at the
// same target geometry. This is the MadGraph-style recursive
// two-body decomposition used by the biasing design.
class DetectorDirected3BodyChannel : public PhaseSpaceChannel {
public:
    enum class InvariantMassMode { Uniform, BreitWigner, PowerLaw };

    DetectorDirected3BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int spectator_index = 0,
        int pair_first_index = 1,
        int pair_second_index = 2,
        int directed_pair_index = 1,
        InvariantMassMode mass_mode = InvariantMassMode::Uniform,
        double resonance_mass = 0.0,
        double resonance_width = 0.0,
        double power_law_nu = 0.8,
        double power_law_offset = 0.0,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume
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

    std::string Name() const override { return "DetectorDirected3Body"; }
    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Decay3Body;
    }
    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::Recursive2Body(
            spectator_index_, pair_first_index_, pair_second_index_);
    }

    void SetVolume(double volume);

private:
    std::shared_ptr<siren::geometry::Geometry const> target_;
    int spectator_index_;
    int pair_first_index_;
    int pair_second_index_;
    int directed_pair_index_;
    InvariantMassMode mass_mode_;
    double resonance_mass_;
    double resonance_width_;
    double power_law_nu_;
    double power_law_offset_;
    DetectorDirected2BodyChannel::Mode mode_;
    double target_volume_;

    double SampleInvariantMassSquared(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        double s_min,
        double s_max
    ) const;

    double InvariantMassDensity(double s, double s_min, double s_max) const;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirected3BodyChannel_H
