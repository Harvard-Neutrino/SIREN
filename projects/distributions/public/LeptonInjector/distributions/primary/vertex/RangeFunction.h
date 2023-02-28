#pragma once
#ifndef LI_RangeFunction_H
#define LI_RangeFunction_H

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace LI {
namespace dataclasses {
struct InteractionSignature;
}
namespace distributions {

class RangeFunction {
friend cereal::access;
public:
    virtual ~RangeFunction() {};
public:
    RangeFunction();
    virtual double operator()(LI::dataclasses::InteractionSignature const & signature, double energy) const;
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
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::RangeFunction, 0);

#endif // LI_RangeFunction_H
