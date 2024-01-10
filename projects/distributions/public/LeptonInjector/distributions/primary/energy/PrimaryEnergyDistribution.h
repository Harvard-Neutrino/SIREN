#pragma once
#ifndef LI_PrimaryEnergyDistribution_H
#define LI_PrimaryEnergyDistribution_H

#include <memory>                                        // for shared_ptr
#include <string>                                        // for string
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/Distributions.h"  // for WeightableDi...

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace utilities { class LI_random; } }
namespace cereal { class access; }

namespace LI {
namespace distributions {

class PrimaryEnergyDistribution : virtual public InjectionDistribution, virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
public:
    virtual ~PrimaryEnergyDistribution() {};
private:
    virtual double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const override = 0;
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
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PrimaryEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PhysicallyNormalizedDistribution, LI::distributions::PrimaryEnergyDistribution);

#endif // LI_PrimaryEnergyDistribution_H
