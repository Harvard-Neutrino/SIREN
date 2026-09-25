#pragma once
#ifndef SIREN_PhysicalChannelAdapters_H
#define SIREN_PhysicalChannelAdapters_H

#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"

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
    // Explicit convention override for models whose legacy metadata cannot
    // describe the density precisely. The topology must agree with signature;
    // the measure is trusted as the caller's declaration of the returned
    // FinalStateProbability density.
    PhysicalDecayChannel(
        std::shared_ptr<siren::interactions::Decay> decay,
        siren::dataclasses::InteractionSignature const & signature,
        PhaseSpaceConvention const & convention);

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

    std::shared_ptr<siren::interactions::Decay> GetDecay() const;

private:
    friend class cereal::access;

    PhysicalDecayChannel() = default;

    std::shared_ptr<siren::interactions::Decay> decay_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            int topology = static_cast<int>(topology_);
            archive(::cereal::make_nvp("Decay", decay_));
            archive(::cereal::make_nvp("Topology", topology));
            archive(::cereal::make_nvp("Measure", measure_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "PhysicalDecayChannel only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            int topology;
            archive(::cereal::make_nvp("Decay", decay_));
            archive(::cereal::make_nvp("Topology", topology));
            archive(::cereal::make_nvp("Measure", measure_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
            if(topology != static_cast<int>(PhaseSpaceTopology::Decay2Body)
               && topology != static_cast<int>(PhaseSpaceTopology::Decay3Body)
               && topology != static_cast<int>(PhaseSpaceTopology::DecayNBody)
               && topology != static_cast<int>(PhaseSpaceTopology::Scatter2to2)
               && topology != static_cast<int>(PhaseSpaceTopology::Scatter2to3)
               && topology != static_cast<int>(PhaseSpaceTopology::Unspecified)) {
                throw std::runtime_error(
                    "PhysicalDecayChannel: invalid PhaseSpaceTopology value "
                    + std::to_string(topology) + " in archive");
            }
            topology_ = static_cast<PhaseSpaceTopology>(topology);
        } else {
            throw std::runtime_error(
                "PhysicalDecayChannel only supports version <= 0!");
        }
    }
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
    // Explicit convention override; see PhysicalDecayChannel above.
    PhysicalCrossSectionChannel(
        std::shared_ptr<siren::interactions::CrossSection> cross_section,
        siren::dataclasses::InteractionSignature const & signature,
        PhaseSpaceConvention const & convention);

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

    std::shared_ptr<siren::interactions::CrossSection> GetCrossSection() const;

private:
    friend class cereal::access;

    PhysicalCrossSectionChannel() = default;

    std::shared_ptr<siren::interactions::CrossSection> cross_section_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            int topology = static_cast<int>(topology_);
            archive(::cereal::make_nvp("CrossSection", cross_section_));
            archive(::cereal::make_nvp("Topology", topology));
            archive(::cereal::make_nvp("Measure", measure_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "PhysicalCrossSectionChannel only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            int topology;
            archive(::cereal::make_nvp("CrossSection", cross_section_));
            archive(::cereal::make_nvp("Topology", topology));
            archive(::cereal::make_nvp("Measure", measure_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
            if(topology != static_cast<int>(PhaseSpaceTopology::Decay2Body)
               && topology != static_cast<int>(PhaseSpaceTopology::Decay3Body)
               && topology != static_cast<int>(PhaseSpaceTopology::DecayNBody)
               && topology != static_cast<int>(PhaseSpaceTopology::Scatter2to2)
               && topology != static_cast<int>(PhaseSpaceTopology::Scatter2to3)
               && topology != static_cast<int>(PhaseSpaceTopology::Unspecified)) {
                throw std::runtime_error(
                    "PhysicalCrossSectionChannel: invalid PhaseSpaceTopology value "
                    + std::to_string(topology) + " in archive");
            }
            topology_ = static_cast<PhaseSpaceTopology>(topology);
        } else {
            throw std::runtime_error(
                "PhysicalCrossSectionChannel only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::PhysicalDecayChannel, 0);
CEREAL_REGISTER_TYPE(siren::injection::PhysicalDecayChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::PhysicalDecayChannel);

CEREAL_CLASS_VERSION(siren::injection::PhysicalCrossSectionChannel, 0);
CEREAL_REGISTER_TYPE(siren::injection::PhysicalCrossSectionChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::PhysicalCrossSectionChannel);

#endif // SIREN_PhysicalChannelAdapters_H
