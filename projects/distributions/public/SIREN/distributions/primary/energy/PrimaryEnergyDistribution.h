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

#include "SIREN/distributions/Distributions.h"  // for WeightableDi...

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace utilities { class LI_random; } }
namespace cereal { class access; }

namespace SI {
namespace distributions {

class PrimaryEnergyDistribution : virtual public PrimaryInjectionDistribution, virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
public:
    virtual ~PrimaryEnergyDistribution() {};
private:
    virtual double SampleEnergy(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const = 0;
public:
    void Sample(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override = 0;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
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
} // namespace SI

CEREAL_CLASS_VERSION(SI::distributions::PrimaryEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(SI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::PrimaryInjectionDistribution, SI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::PhysicallyNormalizedDistribution, SI::distributions::PrimaryEnergyDistribution);

#endif // LI_PrimaryEnergyDistribution_H
