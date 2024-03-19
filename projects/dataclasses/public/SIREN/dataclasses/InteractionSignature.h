#pragma once
#ifndef SIREN_InteractionSignature_H
#define SIREN_InteractionSignature_H

#include <iosfwd>                                 // for ostream
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t
#include <stdexcept>                              // for runtime_error
#include <cereal/cereal.hpp>                      // for make_nvp, CEREAL_CL...

#include "SIREN/dataclasses/Particle.h"  // for Particle
                                                  //

namespace siren { namespace dataclasses { struct InteractionSignature; } }

std::ostream& operator<<(std::ostream& os, siren::dataclasses::InteractionSignature const& signature);

namespace siren {
namespace dataclasses {

struct InteractionSignature {
    siren::dataclasses::ParticleType primary_type = siren::dataclasses::ParticleType::unknown;
    siren::dataclasses::ParticleType target_type = siren::dataclasses::ParticleType::unknown;
    std::vector<siren::dataclasses::ParticleType> secondary_types;

    bool operator==(InteractionSignature const & other) const;
    bool operator<(InteractionSignature const & other) const;
    friend std::ostream& ::operator<<(std::ostream& os, InteractionSignature const& signature);
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
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::InteractionSignature, 0);

#endif // SIREN_InteractionSignature_H
