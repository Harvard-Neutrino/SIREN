#pragma once
#ifndef SIREN_DepthFunction_H
#define SIREN_DepthFunction_H

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

class DepthFunction {
friend cereal::access;
public:
    virtual ~DepthFunction() {};
public:
    DepthFunction();
    virtual double operator()(siren::dataclasses::ParticleType const & primary_type, double energy) const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            throw std::runtime_error("DepthFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("DepthFunction only supports version <= 0!");
        }
    }
    bool operator==(DepthFunction const & distribution) const;
    bool operator<(DepthFunction const & distribution) const;
protected:
    virtual bool equal(DepthFunction const & distribution) const = 0;
    virtual bool less(DepthFunction const & distribution) const = 0;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::DepthFunction, 0);
CEREAL_REGISTER_TYPE(siren::distributions::DepthFunction);

#endif // SIREN_DepthFunction_H
