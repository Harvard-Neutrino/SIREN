#ifndef LI_ParticleID_H
#define LI_ParticleID_H

#include <cstdint>                                      // for uint32_t
#include <ostream>                                      // for ostream
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace LI {
namespace dataclasses {

struct ParticleID {
    uint64_t major_id;
    int64_t minor_id;

    ParticleID() : major_id(0), minor_id(0){};
    ParticleID(uint64_t major, int32_t minor) : 
        major_id(major), minor_id(minor) {};

    static ParticleID GenerateID();

    bool operator<(ParticleID const & other) const;

    bool operator==(ParticleID const & other) const;
    friend std::ostream& operator<<(std::ostream& os, ParticleID const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("MajorID", major_id));
            archive(::cereal::make_nvp("MinorID", minor_id));
        } else {
            throw std::runtime_error("ParticleID only supports version <= 0!");
        }
    };
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::ParticleID, 0);

#endif // LI_ParticleID_H
