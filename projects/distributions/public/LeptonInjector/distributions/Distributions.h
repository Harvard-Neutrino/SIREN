#pragma once
#ifndef LI_Distributions_H
#define LI_Distributions_H

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

#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/utilities/Interpolator.h"

#include "LeptonInjector/math/Vector3D.h"

namespace LI {
namespace utilities {
class LI_random;
}

namespace detector {
class EarthModel;
} // namespace detector

namespace dataclasses {
struct InteractionRecord;
struct InteractionSignature;
}
namespace crosssections {
class CrossSectionCollection;
}

} // namespace LeptonInjector

namespace LI {
namespace distributions {

class PhysicallyNormalizedDistribution {
friend cereal::access;
protected:
    double normalization = 1.0;
public:
    PhysicallyNormalizedDistribution();
    PhysicallyNormalizedDistribution(double norm);
    virtual void SetNormalization(double norm);
    virtual double GetNormalization() const;
    virtual bool IsNormalizationSet() const;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Normalization", normalization));
        } else {
            throw std::runtime_error("PhysicallyNormalizedDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Normalization", normalization));
        } else {
            throw std::runtime_error("PhysicallyNormalizedDistribution only supports version <= 0!");
        }
    }
};

class WeightableDistribution {
friend cereal::access;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const;
    virtual std::string Name() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            throw std::runtime_error("WeightableDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("WeightableDistribution only supports version <= 0!");
        }
    }
    bool operator==(WeightableDistribution const & distribution) const;
    bool operator<(WeightableDistribution const & distribution) const;
    virtual bool AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> second_cross_sections) const;
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class NormalizationConstant : virtual public WeightableDistribution, virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
protected:
    NormalizationConstant();
public:
    NormalizationConstant(double norm);
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const;
    virtual std::string Name() const;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("NormalizationConstant only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("NormalizationConstant only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const;
    virtual bool less(WeightableDistribution const & distribution) const;
};

class InjectionDistribution : virtual public WeightableDistribution {
friend cereal::access;
private:
public:
    virtual void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("InjectionDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("InjectionDistribution only supports version <= 0!");
        }
    }
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PhysicallyNormalizedDistribution, 0);

CEREAL_CLASS_VERSION(LI::distributions::WeightableDistribution, 0);

CEREAL_CLASS_VERSION(LI::distributions::NormalizationConstant, 0);
CEREAL_REGISTER_TYPE(LI::distributions::NormalizationConstant);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::WeightableDistribution, LI::distributions::NormalizationConstant);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PhysicallyNormalizedDistribution, LI::distributions::NormalizationConstant);

CEREAL_CLASS_VERSION(LI::distributions::InjectionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::InjectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::WeightableDistribution, LI::distributions::InjectionDistribution);

#endif // LI_Distributions_H


