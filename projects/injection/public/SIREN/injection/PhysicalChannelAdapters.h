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

namespace detail {

void CheckPhysicalChannelArchiveVersion(std::uint32_t version, char const * channel_name);
PhaseSpaceTopology ReadPhysicalChannelTopology(int topology, char const * channel_name);

template<class Archive, class Channel, class Interaction, class Measure>
void ArchivePhysicalChannel(
        Archive & archive, Channel const * channel, char const * interaction_name,
        Interaction & interaction, int & topology, Measure & measure) {
    archive(::cereal::make_nvp(interaction_name, interaction));
    archive(::cereal::make_nvp("Topology", topology));
    archive(::cereal::make_nvp("Measure", measure));
    archive(::cereal::virtual_base_class<PhaseSpaceChannel>(channel));
}

} // namespace detail

// Adapts a Decay sampler and FinalStateProbability to a proposal channel.
class PhysicalDecayChannel : public PhaseSpaceChannel {
public:
    explicit PhysicalDecayChannel(
        std::shared_ptr<siren::interactions::Decay> decay);
    PhysicalDecayChannel(
        std::shared_ptr<siren::interactions::Decay> decay,
        siren::dataclasses::InteractionSignature const & signature);
    // Declares the convention of a model that does not declare its own. The
    // topology must match the signature; the measure is taken as the caller's
    // statement of what FinalStateProbability returns.
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
        detail::CheckPhysicalChannelArchiveVersion(version, "PhysicalDecayChannel");
        int topology = static_cast<int>(topology_);
        detail::ArchivePhysicalChannel(
            archive, this, "Decay", decay_, topology, measure_);
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        detail::CheckPhysicalChannelArchiveVersion(version, "PhysicalDecayChannel");
        int topology;
        detail::ArchivePhysicalChannel(
            archive, this, "Decay", decay_, topology, measure_);
        topology_ = detail::ReadPhysicalChannelTopology(topology, "PhysicalDecayChannel");
    }
};

class PhysicalCrossSectionChannel : public PhaseSpaceChannel {
public:
    explicit PhysicalCrossSectionChannel(
        std::shared_ptr<siren::interactions::CrossSection> cross_section);
    PhysicalCrossSectionChannel(
        std::shared_ptr<siren::interactions::CrossSection> cross_section,
        siren::dataclasses::InteractionSignature const & signature);
    // Declares the convention; see PhysicalDecayChannel.
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
        detail::CheckPhysicalChannelArchiveVersion(version, "PhysicalCrossSectionChannel");
        int topology = static_cast<int>(topology_);
        detail::ArchivePhysicalChannel(
            archive, this, "CrossSection", cross_section_, topology, measure_);
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        detail::CheckPhysicalChannelArchiveVersion(version, "PhysicalCrossSectionChannel");
        int topology;
        detail::ArchivePhysicalChannel(
            archive, this, "CrossSection", cross_section_, topology, measure_);
        topology_ = detail::ReadPhysicalChannelTopology(topology, "PhysicalCrossSectionChannel");
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
