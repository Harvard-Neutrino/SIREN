#pragma once
#ifndef LI_DecayRangeFunction_H
#define LI_DecayRangeFunction_H

#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/primary/vertex/RangeFunction.h"

namespace LI { namespace dataclasses { enum class ParticleType : int32_t; } }

namespace LI {
namespace distributions {

class DecayRangeFunction : virtual public RangeFunction {
friend cereal::access;
protected:
    DecayRangeFunction() {};
private:
    double particle_mass; // GeV
    double decay_width; // GeV
    double multiplier;
    double max_distance;
public:
    DecayRangeFunction(double particle_mass, double decay_width, double multiplier, double max_distance);
    static double DecayLength(double mass, double width, double energy);
    double operator()(LI::dataclasses::ParticleType const & primary_type, double energy) const override;
    double DecayLength(LI::dataclasses::ParticleType const & primary_type, double energy) const;
    double Range(LI::dataclasses::ParticleType const & primary_type, double energy) const;
    double Multiplier() const;
    double ParticleMass() const;
    double DecayWidth() const;
    double MaxDistance() const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("ParticleMass", particle_mass));
            archive(::cereal::make_nvp("DecayWidth", decay_width));
            archive(::cereal::make_nvp("Multiplier", multiplier));
            archive(::cereal::make_nvp("MaxDistance", max_distance));
            archive(cereal::virtual_base_class<RangeFunction>(this));
        } else {
            throw std::runtime_error("DecayRangeFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<DecayRangeFunction> & construct, std::uint32_t const version) {
        if(version == 0) {
            double mass;
            double width;
            double multiplier;
            double max_distance;
            archive(::cereal::make_nvp("ParticleMass", mass));
            archive(::cereal::make_nvp("DecayWidth", width));
            archive(::cereal::make_nvp("Multiplier", multiplier));
            archive(::cereal::make_nvp("MaxDistance", max_distance));
            construct(mass, width, multiplier, max_distance);
            archive(cereal::virtual_base_class<RangeFunction>(construct.ptr()));
        } else {
            throw std::runtime_error("DecayRangeFunction only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(RangeFunction const & distribution) const override;
    virtual bool less(RangeFunction const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::DecayRangeFunction, 0);
CEREAL_REGISTER_TYPE(LI::distributions::DecayRangeFunction);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::RangeFunction, LI::distributions::DecayRangeFunction);

#endif // LI_DecayRangeFunction_H
