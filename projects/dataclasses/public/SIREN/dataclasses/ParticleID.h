#ifndef SIREN_ParticleID_H
#define SIREN_ParticleID_H

#include <cstdint>                                      // for uint32_t
#include <ostream>                                      // for ostream
#include <stdexcept>                                    // for runtime_error
#include <sstream>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace siren { namespace dataclasses { class ParticleID; } }

std::ostream& operator<<(std::ostream& os, siren::dataclasses::ParticleID const& record);

namespace siren {
namespace dataclasses {

class ParticleID {
    bool id_set = false;
    uint64_t major_id = 0;
    int64_t minor_id = 0;
public:
    static ParticleID GenerateID();

    ParticleID();
    ParticleID(uint64_t major, int32_t minor);
    ParticleID & operator=(ParticleID const & other) = default;
    ParticleID(ParticleID const & other) = default;
    ParticleID & operator=(ParticleID && other) = default;
    ParticleID(ParticleID && other) = default;

    bool IsSet() const;
    operator bool() const;

    uint64_t GetMajorID() const;
    int32_t GetMinorID() const;

    void SetID(uint64_t major, int32_t minor);

    bool operator<(ParticleID const & other) const;

    bool operator==(ParticleID const & other) const;
    friend std::ostream& ::operator<<(std::ostream& os, ParticleID const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("IDSet", id_set));
            archive(::cereal::make_nvp("MajorID", major_id));
            archive(::cereal::make_nvp("MinorID", minor_id));
        } else {
            throw std::runtime_error("ParticleID only supports version <= 0!");
        }
    };
};

} // namespace dataclasses
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::ParticleID, 0);

#endif // SIREN_ParticleID_H
