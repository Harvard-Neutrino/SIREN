#pragma once
#ifndef SIREN_Isotropic2BodyChannel_H
#define SIREN_Isotropic2BodyChannel_H

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
#include <cereal/types/polymorphic.hpp>

namespace siren {
namespace injection {

// Uniform proposal with density 1/(4*pi) in rest-frame solid angle.
class Isotropic2BodyChannel : public PhaseSpaceChannel {
public:
    // daughter_index: which secondary (0 or 1) in the record is
    // the "daughter A" whose direction we track.
    explicit Isotropic2BodyChannel(int daughter_index = 0);

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override { return "Isotropic2Body"; }
    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Decay2Body;
    }
    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::SolidAngleRest();
    }

private:
    friend class cereal::access;

    int daughter_index_;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("DaughterIndex", daughter_index_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "Isotropic2BodyChannel only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DaughterIndex", daughter_index_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "Isotropic2BodyChannel only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::Isotropic2BodyChannel, 0);
CEREAL_REGISTER_TYPE(siren::injection::Isotropic2BodyChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::Isotropic2BodyChannel);

#endif // SIREN_Isotropic2BodyChannel_H
