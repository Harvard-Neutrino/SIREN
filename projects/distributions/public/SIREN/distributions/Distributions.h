#pragma once
#ifndef SIREN_Distributions_H
#define SIREN_Distributions_H

#include <string>                                        // for string
#include <memory>                                        // for shared_ptr
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error
#include <vector>                                        // for vector

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/InteractionTree.h"  // for InteractionT...

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace utilities {
class SIREN_random;
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

} // namespace sirenREN

namespace siren {
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
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const = 0;
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
    virtual bool AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const;
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
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
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

class PrimaryInjectionDistribution : virtual public WeightableDistribution {
friend cereal::access;
private:
public:
    virtual ~PrimaryInjectionDistribution() {};
    virtual void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const = 0;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryInjectionDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryInjectionDistribution only supports version <= 0!");
        }
    }
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

class SecondaryInjectionDistribution : virtual public WeightableDistribution {
friend cereal::access;
private:
public:
    virtual ~SecondaryInjectionDistribution() {};
    virtual void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const = 0;
    virtual std::shared_ptr<SecondaryInjectionDistribution> clone() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("SecondaryInjectionDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("SecondaryInjectionDistribution only supports version <= 0!");
        }
    }
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PhysicallyNormalizedDistribution, 0);

CEREAL_CLASS_VERSION(siren::distributions::WeightableDistribution, 0);

CEREAL_CLASS_VERSION(siren::distributions::NormalizationConstant, 0);
CEREAL_REGISTER_TYPE(siren::distributions::NormalizationConstant);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::WeightableDistribution, siren::distributions::NormalizationConstant);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PhysicallyNormalizedDistribution, siren::distributions::NormalizationConstant);

CEREAL_CLASS_VERSION(siren::distributions::PrimaryInjectionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::PrimaryInjectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::WeightableDistribution, siren::distributions::PrimaryInjectionDistribution);

CEREAL_CLASS_VERSION(siren::distributions::SecondaryInjectionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::SecondaryInjectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::WeightableDistribution, siren::distributions::SecondaryInjectionDistribution);

#endif // SIREN_Distributions_H


