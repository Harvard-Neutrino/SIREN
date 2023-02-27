#pragma once
#ifndef LI_DecayRangePositionDistribution_H
#define LI_DecayRangePositionDistribution_H

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/utilities/Particle.h"
#include "LeptonInjector/utilities/Interpolator.h"

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/VertexPositionDistribution.h"

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace crosssections {
struct InteractionRecord;
struct InteractionSignature;
class CrossSectionCollection;
} // namespace crosssections
} // namespace LeptonInjector

namespace LI {
namespace distributions {

class RangeFunction {
friend cereal::access;
public:
    RangeFunction();
    virtual double operator()(LI::crosssections::InteractionSignature const & signature, double energy) const;
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
    double operator()(LI::crosssections::InteractionSignature const & signature, double energy) const override;
    double DecayLength(LI::crosssections::InteractionSignature const & signature, double energy) const;
    double Range(LI::crosssections::InteractionSignature const & signature, double energy) const;
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

class RangePositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    double radius;
    double endcap_length;
    std::shared_ptr<RangeFunction> range_function;
    std::set<LI::utilities::Particle::ParticleType> target_types;

    LI::math::Vector3D SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const;

    LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    RangePositionDistribution();
    RangePositionDistribution(const RangePositionDistribution &) = default;
    RangePositionDistribution(double radius, double endcap_length, std::shared_ptr<RangeFunction> range_function, std::set<LI::utilities::Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("RangeFunction", range_function));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("RangePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<RangePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::set<LI::utilities::Particle::ParticleType> t;
            std::shared_ptr<RangeFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("RangeFunction", f));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, f, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("RangePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class DecayRangePositionDistribution : virtual public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DecayRangeFunction> range_function;
    std::set<LI::utilities::Particle::ParticleType> target_types;

    LI::math::Vector3D SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const;

    LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    DecayRangePositionDistribution();
    DecayRangePositionDistribution(const DecayRangePositionDistribution &) = default;
    DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function, std::set<LI::utilities::Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    virtual bool AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> second_cross_sections) const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("DecayRangeFunction", range_function));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("DecayRangePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<DecayRangePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::set<LI::utilities::Particle::ParticleType> t;
            std::shared_ptr<DecayRangeFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("DecayRangeFunction", f));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, f, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("DecayRangePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::RangeFunction, 0);

CEREAL_CLASS_VERSION(LI::distributions::DecayRangeFunction, 0);
CEREAL_REGISTER_TYPE(LI::distributions::DecayRangeFunction);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::RangeFunction, LI::distributions::DecayRangeFunction);

CEREAL_CLASS_VERSION(LI::distributions::RangePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::RangePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::VertexPositionDistribution, LI::distributions::RangePositionDistribution);

CEREAL_CLASS_VERSION(LI::distributions::DecayRangePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::DecayRangePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::VertexPositionDistribution, LI::distributions::DecayRangePositionDistribution);

#endif // LI_DecayRangePositionDistribution_H
