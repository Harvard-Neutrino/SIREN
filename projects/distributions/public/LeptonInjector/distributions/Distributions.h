#pragma once
#ifndef LI_Distributions_H
#define LI_Distributions_H

#include <string>                                        // for string
#include <memory>                                        // for shared_ptr
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error
#include <vector>                                        // for vector

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/InteractionTree.h"  // for InteractionT...

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace utilities {
class LI_random;
}

namespace detector {
class DetectorModel;
} // namespace detector

namespace dataclasses {
class InteractionRecord;
struct InteractionSignature;
}
namespace interactions {
class InteractionCollection;
}

} // namespace LeptonInjector

namespace LI {
namespace distributions {

class PhysicallyNormalizedDistribution {
friend cereal::access;
protected:
    bool normalization_set = false;
    double normalization = 1.0;
public:
    PhysicallyNormalizedDistribution();
    virtual ~PhysicallyNormalizedDistribution() {};
    PhysicallyNormalizedDistribution(double norm);
    virtual void SetNormalization(double norm);
    virtual double GetNormalization() const;
    virtual bool IsNormalizationSet() const;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("NormalizationSet", normalization_set));
            archive(::cereal::make_nvp("Normalization", normalization));
        } else {
            throw std::runtime_error("PhysicallyNormalizedDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("NormalizationSet", normalization_set));
            archive(::cereal::make_nvp("Normalization", normalization));
        } else {
            throw std::runtime_error("PhysicallyNormalizedDistribution only supports version <= 0!");
        }
    }
};

class WeightableDistribution {
friend cereal::access;
public:
    virtual ~WeightableDistribution() {};
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const = 0;
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
    virtual bool AreEquivalent(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::DetectorModel const> second_detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> second_interactions) const;
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class NormalizationConstant : virtual public WeightableDistribution, virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
protected:
    NormalizationConstant();
public:
    virtual ~NormalizationConstant() {};
    NormalizationConstant(double norm);
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override;
    virtual std::string Name() const override;
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
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class InjectionDistribution : virtual public WeightableDistribution {
friend cereal::access;
private:
public:
    virtual ~InjectionDistribution() {};
    virtual void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const;
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


