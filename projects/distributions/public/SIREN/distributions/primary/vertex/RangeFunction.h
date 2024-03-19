#pragma once
#ifndef SIREN_RangeFunction_H
#define SIREN_RangeFunction_H

#include <cstdint>            // for uint32_t
#include <stdexcept>          // for runtime_error

#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace siren { namespace dataclasses { enum class ParticleType : int32_t; } }

namespace siren {
namespace distributions {

class RangeFunction {
friend cereal::access;
public:
    virtual ~RangeFunction() {};
public:
    RangeFunction();
    virtual double operator()(siren::dataclasses::ParticleType const & primary_type, double energy) const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            throw std::runtime_error("RangeFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("RangeFunction only supports version <= 0!");
        }
    }
    bool operator==(RangeFunction const & distribution) const;
    bool operator<(RangeFunction const & distribution) const;
protected:
    virtual bool equal(RangeFunction const & distribution) const = 0;
    virtual bool less(RangeFunction const & distribution) const = 0;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::RangeFunction, 0);

#endif // SIREN_RangeFunction_H
