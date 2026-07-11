#pragma once
#ifndef SIREN_DetectorDirected3BodyChannel_H
#define SIREN_DetectorDirected3BodyChannel_H

#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>
#include <vector>

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

struct TabulatedMappingTable;

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
    enum class InvariantMassMode { Uniform, BreitWigner, PowerLaw, Tabulated };
    enum class Factorization { Direct, Recursive };

    // Direct mode constructor: just specify which daughter to bias.
    // The other two indices are inferred (ascending order from {0,1,2}).
    // A positive `volume` supplies the exact target volume for composite or
    // very thin geometries, as in DetectorDirected2BodyChannel.
    DetectorDirected3BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int directed_index,
        InvariantMassMode mass_mode = InvariantMassMode::Uniform,
        double resonance_mass = 0.0,
        double resonance_width = 0.0,
        double power_law_nu = 0.8,
        double power_law_offset = 0.0,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay3Body,
        std::vector<double> mass_cdf_nodes = {},
        std::vector<double> mass_cdf_values = {},
        double volume = -1.0
    );

    // Recursive mode constructor: explicit spectator/pair/directed indices. A
    // positive `volume` supplies the exact target volume for thin/composite
    // geometries.
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
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay3Body,
        std::vector<double> mass_cdf_nodes = {},
        std::vector<double> mass_cdf_values = {},
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

    std::string Name() const override { return "DetectorDirected3Body"; }
    PhaseSpaceTopology Topology() const override {
        return topology_;
    }
    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::Recursive2Body(
            roles_.outer, roles_.pair_first, roles_.pair_second);
    }

    void SetVolume(double volume);

    // True if the directed sub-step actually directs (a non-isotropic proposal)
    // at the point in `record`, vs the isotropic 1/4pi fallback.  Lets a
    // diagnostic attribute the channel's variance into directing vs fallback.
    bool DirectingActive(
        siren::dataclasses::InteractionRecord const & record) const override;

private:
    struct RoleIndices {
        int outer;
        int pair_first;
        int pair_second;
        int directed;

        int PairOther() const {
            return directed == pair_first ? pair_second : pair_first;
        }
    };

    std::shared_ptr<siren::geometry::Geometry const> target_;
    RoleIndices roles_;

    using SampleImplementation = void (DetectorDirected3BodyChannel::*)(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        siren::dataclasses::InteractionRecord &) const;
    using DensityImplementation = double (DetectorDirected3BodyChannel::*)(
        std::shared_ptr<siren::detector::DetectorModel const>,
        siren::dataclasses::InteractionRecord const &) const;
    using ActiveImplementation = bool (DetectorDirected3BodyChannel::*)(
        siren::dataclasses::InteractionRecord const &) const;

    SampleImplementation sample_implementation_;
    DensityImplementation density_implementation_;
    ActiveImplementation active_implementation_;

    InvariantMassMode mass_mode_;
    double resonance_mass_;
    double resonance_width_;
    double power_law_nu_;
    double power_law_offset_;
    std::shared_ptr<TabulatedMappingTable const> mass_cdf_table_;
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

    bool DirectingActiveDirect(
        siren::dataclasses::InteractionRecord const & record) const;
    bool DirectingActiveRecursive(
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
