#pragma once
#ifndef LI_PrimaryEnergyDistribution_H
#define LI_PrimaryEnergyDistribution_H

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

#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace detector {
class EarthModel;
} // namespace detector

namespace crosssections {
struct InteractionRecord;
struct InteractionSignature;
class CrossSectionCollection;
} // namespace crosssections
} // namespace LI

namespace LI {
namespace distributions {

class PrimaryEnergyDistribution : virtual public InjectionDistribution, virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
private:
    virtual double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PrimaryEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PhysicallyNormalizedDistribution, LI::distributions::PrimaryEnergyDistribution);

#endif // LI_PrimaryEnergyDistribution_H
