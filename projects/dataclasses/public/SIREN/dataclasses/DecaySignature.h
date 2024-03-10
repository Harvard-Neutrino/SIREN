#ifndef SIREN_DecaySignature_H
#define SIREN_DecaySignature_H

#include <cstdint>                                // for uint32_t
#include <vector>                                 // for vector
#include <ostream>                                // for ostream
#include <stdexcept>                              // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"

namespace siren {
namespace dataclasses {

struct DecaySignature {
    siren::dataclasses::ParticleType primary_type;
    std::vector<siren::dataclasses::ParticleType> secondary_types;

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
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::DecaySignature, 0);

#endif // SIREN_DecaySignature_H
