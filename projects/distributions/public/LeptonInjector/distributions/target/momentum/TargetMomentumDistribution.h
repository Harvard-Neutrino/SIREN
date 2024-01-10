#pragma once
#ifndef LI_TargetMomentumDistribution_H
#define LI_TargetMomentumDistribution_H

#include <array>                                         // for array
#include <memory>
#include <string>
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/Distributions.h"  // for InjectionDis...

namespace LI { namespace crosssections { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class EarthModel; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class TargetMomentumDistribution : virtual public InjectionDistribution {
friend cereal::access;
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("TargetMomentumDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("TargetMomentumDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

class TargetAtRest : virtual public TargetMomentumDistribution {
friend cereal::access;
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<TargetMomentumDistribution>(this));
        } else {
            throw std::runtime_error("TargetAtRest only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<TargetMomentumDistribution>(this));
        } else {
            throw std::runtime_error("TargetAtRest only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::TargetMomentumDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::TargetMomentumDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::TargetMomentumDistribution);

CEREAL_CLASS_VERSION(LI::distributions::TargetAtRest, 0);
CEREAL_REGISTER_TYPE(LI::distributions::TargetAtRest);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::TargetMomentumDistribution, LI::distributions::TargetAtRest);

#endif // LI_TargetMomentumDistribution_H
