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

// Detector-directed 3-body phase-space channel.
//
// Two factorization modes:
//
//   Direct (default):  P -> directed + X,  X -> other_a + other_b
//     Samples the directed daughter's lab direction toward the target
//     using the parent's boost.  The complementary system X decays
//     isotropically.  The invariant mass variable is s_X = M^2_{a+b}.
//     Use this for highly-boosted parents (e.g., pion 3-body decay).
//
//   Recursive:  P -> spectator + pair,  pair -> directed + other
//     Samples the pair direction isotropically, then directs one
//     daughter within the pair's frame.  The invariant mass variable
//     is s_pair = M^2_{directed+other}.  Use this when the pair has
//     a resonance structure (e.g., off-shell chi' with BW sampling)
//     or when the parent boost is small.
class DetectorDirected3BodyChannel : public PhaseSpaceChannel {
public:
    enum class InvariantMassMode { Uniform, BreitWigner, PowerLaw };
    enum class Factorization { Direct, Recursive };

    // Direct mode constructor: just specify which daughter to bias.
    // The other two indices are inferred (ascending order from {0,1,2}).
    DetectorDirected3BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int directed_index,
        InvariantMassMode mass_mode = InvariantMassMode::Uniform,
        double resonance_mass = 0.0,
        double resonance_width = 0.0,
        double power_law_nu = 0.8,
        double power_law_offset = 0.0,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay3Body
    );

    // Recursive mode constructor (backward compatible).
    DetectorDirected3BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int spectator_index,
        int pair_first_index,
        int pair_second_index,
        int directed_pair_index,
        InvariantMassMode mass_mode = InvariantMassMode::Uniform,
        double resonance_mass = 0.0,
        double resonance_width = 0.0,
        double power_law_nu = 0.8,
        double power_law_offset = 0.0,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay3Body
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
        return topology_;
    }
    PhaseSpaceMeasure Measure() const override {
        if (factorization_ == Factorization::Direct) {
            return PhaseSpaceMeasure::Recursive2Body(
                other_a_index_, other_b_index_, directed_index_);
        }
        return PhaseSpaceMeasure::Recursive2Body(
            spectator_index_, pair_first_index_, pair_second_index_);
    }

    void SetVolume(double volume);

private:
    Factorization factorization_;
    std::shared_ptr<siren::geometry::Geometry const> target_;

    // Direct mode indices
    int directed_index_;
    int other_a_index_;
    int other_b_index_;

    // Recursive mode indices (also used for backward compat)
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
    PhaseSpaceTopology topology_;
    double target_volume_;

    void SampleDirect(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record) const;

    double DensityDirect(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record) const;

    void SampleRecursive(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record) const;

    double DensityRecursive(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record) const;

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
