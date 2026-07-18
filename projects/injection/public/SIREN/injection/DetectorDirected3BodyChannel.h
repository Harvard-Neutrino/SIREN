#pragma once
#ifndef SIREN_DetectorDirected3BodyChannel_H
#define SIREN_DetectorDirected3BodyChannel_H

#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

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
    friend cereal::access;

    DetectorDirected3BodyChannel() = default;

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

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (version != 0) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel only supports version <= 0!");
        }

        // Archive the installed factorization, never method pointers.
        Factorization factorization;
        if (sample_implementation_ == &DetectorDirected3BodyChannel::SampleDirect &&
            density_implementation_ == &DetectorDirected3BodyChannel::DensityDirect &&
            active_implementation_ == &DetectorDirected3BodyChannel::DirectingActiveDirect) {
            factorization = Factorization::Direct;
        } else if (
            sample_implementation_ == &DetectorDirected3BodyChannel::SampleRecursive &&
            density_implementation_ == &DetectorDirected3BodyChannel::DensityRecursive &&
            active_implementation_ == &DetectorDirected3BodyChannel::DirectingActiveRecursive) {
            factorization = Factorization::Recursive;
        } else {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel has inconsistent implementations");
        }

        archive(::cereal::make_nvp("Target", target_));
        archive(::cereal::make_nvp("OuterIndex", roles_.outer));
        archive(::cereal::make_nvp("PairFirstIndex", roles_.pair_first));
        archive(::cereal::make_nvp("PairSecondIndex", roles_.pair_second));
        archive(::cereal::make_nvp("DirectedIndex", roles_.directed));
        archive(::cereal::make_nvp("Factorization", static_cast<int>(factorization)));
        archive(::cereal::make_nvp("InvariantMassMode", static_cast<int>(mass_mode_)));
        archive(::cereal::make_nvp("ResonanceMass", resonance_mass_));
        archive(::cereal::make_nvp("ResonanceWidth", resonance_width_));
        archive(::cereal::make_nvp("PowerLawNu", power_law_nu_));
        archive(::cereal::make_nvp("PowerLawOffset", power_law_offset_));
        archive(::cereal::make_nvp("MassCdfTable", mass_cdf_table_));
        archive(::cereal::make_nvp("Mode", static_cast<int>(mode_)));
        archive(::cereal::make_nvp("Topology", static_cast<int>(topology_)));
        archive(::cereal::make_nvp("TargetVolume", target_volume_));
        archive(cereal::virtual_base_class<PhaseSpaceChannel>(this));
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if (version != 0) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel only supports version <= 0!");
        }

        int factorization_int;
        int mass_mode_int;
        int mode_int;
        int topology_int;
        archive(::cereal::make_nvp("Target", target_));
        archive(::cereal::make_nvp("OuterIndex", roles_.outer));
        archive(::cereal::make_nvp("PairFirstIndex", roles_.pair_first));
        archive(::cereal::make_nvp("PairSecondIndex", roles_.pair_second));
        archive(::cereal::make_nvp("DirectedIndex", roles_.directed));
        archive(::cereal::make_nvp("Factorization", factorization_int));
        archive(::cereal::make_nvp("InvariantMassMode", mass_mode_int));
        archive(::cereal::make_nvp("ResonanceMass", resonance_mass_));
        archive(::cereal::make_nvp("ResonanceWidth", resonance_width_));
        archive(::cereal::make_nvp("PowerLawNu", power_law_nu_));
        archive(::cereal::make_nvp("PowerLawOffset", power_law_offset_));
        archive(::cereal::make_nvp("MassCdfTable", mass_cdf_table_));
        archive(::cereal::make_nvp("Mode", mode_int));
        archive(::cereal::make_nvp("Topology", topology_int));
        archive(::cereal::make_nvp("TargetVolume", target_volume_));
        archive(cereal::virtual_base_class<PhaseSpaceChannel>(this));

        if (factorization_int != static_cast<int>(Factorization::Direct) &&
            factorization_int != static_cast<int>(Factorization::Recursive)) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel: invalid Factorization value " +
                std::to_string(factorization_int) + " in archive");
        }
        if (mass_mode_int != static_cast<int>(InvariantMassMode::Uniform) &&
            mass_mode_int != static_cast<int>(InvariantMassMode::BreitWigner) &&
            mass_mode_int != static_cast<int>(InvariantMassMode::PowerLaw) &&
            mass_mode_int != static_cast<int>(InvariantMassMode::Tabulated)) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel: invalid InvariantMassMode value " +
                std::to_string(mass_mode_int) + " in archive");
        }
        if (mode_int != static_cast<int>(DetectorDirected2BodyChannel::Mode::Cone) &&
            mode_int != static_cast<int>(DetectorDirected2BodyChannel::Mode::Volume)) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel: invalid Mode value " +
                std::to_string(mode_int) + " in archive");
        }
        if (topology_int != static_cast<int>(PhaseSpaceTopology::Decay2Body) &&
            topology_int != static_cast<int>(PhaseSpaceTopology::Decay3Body) &&
            topology_int != static_cast<int>(PhaseSpaceTopology::DecayNBody) &&
            topology_int != static_cast<int>(PhaseSpaceTopology::Scatter2to2) &&
            topology_int != static_cast<int>(PhaseSpaceTopology::Scatter2to3) &&
            topology_int != static_cast<int>(PhaseSpaceTopology::Unspecified)) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel: invalid PhaseSpaceTopology value " +
                std::to_string(topology_int) + " in archive");
        }

        mass_mode_ = static_cast<InvariantMassMode>(mass_mode_int);
        mode_ = static_cast<DetectorDirected2BodyChannel::Mode>(mode_int);
        topology_ = static_cast<PhaseSpaceTopology>(topology_int);
        // Tabulated mode owns exactly one immutable table.
        if ((mass_mode_ == InvariantMassMode::Tabulated) !=
            static_cast<bool>(mass_cdf_table_)) {
            throw std::runtime_error(
                "DetectorDirected3BodyChannel: invalid tabulated table state");
        }

        Factorization factorization = static_cast<Factorization>(factorization_int);
        if (factorization == Factorization::Direct) {
            sample_implementation_ = &DetectorDirected3BodyChannel::SampleDirect;
            density_implementation_ = &DetectorDirected3BodyChannel::DensityDirect;
            active_implementation_ = &DetectorDirected3BodyChannel::DirectingActiveDirect;
        } else {
            sample_implementation_ = &DetectorDirected3BodyChannel::SampleRecursive;
            density_implementation_ = &DetectorDirected3BodyChannel::DensityRecursive;
            active_implementation_ = &DetectorDirected3BodyChannel::DirectingActiveRecursive;
        }
    }

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

CEREAL_CLASS_VERSION(siren::injection::DetectorDirected3BodyChannel, 0);
CEREAL_REGISTER_TYPE(siren::injection::DetectorDirected3BodyChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::DetectorDirected3BodyChannel);

#endif // SIREN_DetectorDirected3BodyChannel_H
