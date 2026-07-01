#pragma once
#ifndef SIREN_PhysicalChannelAdapters_H
#define SIREN_PhysicalChannelAdapters_H

#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>

namespace siren { namespace interactions { class Decay; } }
namespace siren { namespace interactions { class CrossSection; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }

namespace siren {
namespace injection {

// Wraps an existing Decay as a PhaseSpaceChannel.
//
// This is the "physical channel" in the multi-channel framework:
// it samples from the physical angular distribution (including
// polarization, matrix element structure, resonances) and
// evaluates the physical density via FinalStateProbability.
//
// Unlike the standalone Isotropic2BodyChannel, this channel
// captures the actual physics of the decay model it wraps.
//
// Usage in multi-channel:
//   auto physical = make_shared<PhysicalDecayChannel>(my_hnl_decay);
//   auto biased = make_shared<DetectorDirected2BodyChannel>(target);
//   mc.channels = {physical, biased};
//   mc.weights = {0.01, 0.99};
class PhysicalDecayChannel : public PhaseSpaceChannel {
public:
    explicit PhysicalDecayChannel(
        std::shared_ptr<siren::interactions::Decay> decay);
    PhysicalDecayChannel(
        std::shared_ptr<siren::interactions::Decay> decay,
        siren::dataclasses::InteractionSignature const & signature);
    PhysicalDecayChannel(
        std::shared_ptr<siren::interactions::Decay> decay,
        PhaseSpaceConvention convention);

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override;
    PhaseSpaceTopology Topology() const override;
    PhaseSpaceMeasure Measure() const override;
    PhaseSpaceConvention Convention() const override;

    std::shared_ptr<siren::interactions::Decay> GetDecay() const;

private:
    std::shared_ptr<siren::interactions::Decay> decay_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;
};

// Wraps an existing CrossSection as a PhaseSpaceChannel.
//
// Same adapter pattern as PhysicalDecayChannel but for
// scattering processes.
class PhysicalCrossSectionChannel : public PhaseSpaceChannel {
public:
    explicit PhysicalCrossSectionChannel(
        std::shared_ptr<siren::interactions::CrossSection> cross_section);
    PhysicalCrossSectionChannel(
        std::shared_ptr<siren::interactions::CrossSection> cross_section,
        siren::dataclasses::InteractionSignature const & signature);
    PhysicalCrossSectionChannel(
        std::shared_ptr<siren::interactions::CrossSection> cross_section,
        PhaseSpaceConvention convention);

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override;
    PhaseSpaceTopology Topology() const override;
    PhaseSpaceMeasure Measure() const override;
    PhaseSpaceConvention Convention() const override;

    std::shared_ptr<siren::interactions::CrossSection> GetCrossSection() const;

private:
    std::shared_ptr<siren::interactions::CrossSection> cross_section_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;
};

} // namespace injection
} // namespace siren

#endif // SIREN_PhysicalChannelAdapters_H
