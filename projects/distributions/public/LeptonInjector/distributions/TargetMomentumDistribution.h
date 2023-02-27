#pragma once
#ifndef LI_TargetMomentumDistribution_H
#define LI_TargetMomentumDistribution_H

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

class TargetMomentumDistribution : virtual public InjectionDistribution {
friend cereal::access;
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
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
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class TargetAtRest : virtual public TargetMomentumDistribution {
friend cereal::access;
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
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
