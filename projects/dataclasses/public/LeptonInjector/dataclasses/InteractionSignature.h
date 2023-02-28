#pragma once
#ifndef LI_InteractionSignature_H
#define LI_InteractionSignature_H

#include <vector>

#include "LeptonInjector/dataclasses/Particle.h"

namespace LI {
namespace dataclasses {

struct InteractionSignature {
    LI::dataclasses::Particle::ParticleType primary_type;
    LI::dataclasses::Particle::ParticleType target_type;
    std::vector<LI::dataclasses::Particle::ParticleType> secondary_types;
    bool operator==(InteractionSignature const & other) const;
    bool operator<(InteractionSignature const & other) const;
    friend std::ostream& operator<<(std::ostream& os, InteractionSignature const& signature);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("TargetType", target_type));
            archive(cereal::make_nvp("SecondaryTypes", secondary_types));
        } else {
            throw std::runtime_error("InteractionSignature only supports version <= 0!");
        }
    }
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::InteractionSignature, 0);

#endif // LI_InteractionSignature_H
