#ifndef LI_DecaySignature_H
#define LI_DecaySignature_H

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include <vector>
#include <ostream>
#include <stdint.h>
#include <stdexcept>

#include "LeptonInjector/dataclasses/Particle.h"

namespace LI {
namespace dataclasses {

struct DecaySignature {
    LI::dataclasses::Particle::ParticleType primary_type;
    std::vector<LI::dataclasses::Particle::ParticleType> secondary_types;
    bool operator==(DecaySignature const & other) const;
    friend std::ostream& operator<<(std::ostream& os, DecaySignature const& signature);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("SecondaryTypes", secondary_types));
        } else {
            throw std::runtime_error("DecaySignature only supports version <= 0!");
        }
    }
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::DecaySignature, 0);

#endif // LI_DecaySignature_H
