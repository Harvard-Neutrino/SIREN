#pragma once
#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

namespace siren { namespace injection {

// Uniform rest-solid-angle proposal on a conservative preimage of the target
// bounding cone. Empty or full-sphere bounds use normalized isotropic sampling.
class RestFrameEnvelope2BodyChannel : public PhaseSpaceChannel {
public:
    RestFrameEnvelope2BodyChannel(
        std::shared_ptr<siren::geometry::Geometry const> target, int daughter_index=0);
    void Sample(std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        siren::dataclasses::InteractionRecord &) const override;
    double Density(std::shared_ptr<siren::detector::DetectorModel const>,
        siren::dataclasses::InteractionRecord const &) const override;
    bool DirectingActive(siren::dataclasses::InteractionRecord const &) const override;
    std::string Name() const override { return "RestFrameEnvelope2Body"; }
    PhaseSpaceTopology Topology() const override { return PhaseSpaceTopology::Decay2Body; }
    PhaseSpaceMeasure Measure() const override { return PhaseSpaceMeasure::SolidAngleRest(); }

    template<class Archive> void save(Archive & ar, std::uint32_t version) const {
        if (version!=1) throw std::runtime_error("RestFrameEnvelope2BodyChannel: legacy boost convention rejected; regenerate the archive");
        ar(cereal::make_nvp("Target", target_), cereal::make_nvp("DaughterIndex", daughter_index_),
           cereal::virtual_base_class<PhaseSpaceChannel>(this));
    }
    template<class Archive> void load(Archive & ar, std::uint32_t version) {
        if (version!=1) throw std::runtime_error("RestFrameEnvelope2BodyChannel: legacy boost convention rejected; regenerate the archive");
        ar(cereal::make_nvp("Target", target_), cereal::make_nvp("DaughterIndex", daughter_index_),
           cereal::virtual_base_class<PhaseSpaceChannel>(this));
        if (!target_ || daughter_index_<0 || daughter_index_>1)
            throw std::runtime_error("RestFrameEnvelope2BodyChannel: invalid archived target or daughter");
    }
private:
    friend class cereal::access;
    RestFrameEnvelope2BodyChannel() = default;
    std::shared_ptr<siren::geometry::Geometry const> target_;
    int daughter_index_=0;
};
}}
CEREAL_CLASS_VERSION(siren::injection::RestFrameEnvelope2BodyChannel, 1);
CEREAL_REGISTER_TYPE(siren::injection::RestFrameEnvelope2BodyChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::injection::PhaseSpaceChannel,siren::injection::RestFrameEnvelope2BodyChannel);
