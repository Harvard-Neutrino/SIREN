#pragma once
#ifndef SIREN_DecayRangeFunction_H
#define SIREN_DecayRangeFunction_H

#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/primary/vertex/RangeFunction.h"

namespace siren { namespace dataclasses { enum class ParticleType : int32_t; } }

namespace siren {
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
    double operator()(siren::dataclasses::ParticleType const & primary_type, double energy) const override;
    double DecayLength(siren::dataclasses::ParticleType const & primary_type, double energy) const;
    double Range(siren::dataclasses::ParticleType const & primary_type, double energy) const;
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
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::DecayRangeFunction, 0);
CEREAL_REGISTER_TYPE(siren::distributions::DecayRangeFunction);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::RangeFunction, siren::distributions::DecayRangeFunction);

#endif // SIREN_DecayRangeFunction_H
