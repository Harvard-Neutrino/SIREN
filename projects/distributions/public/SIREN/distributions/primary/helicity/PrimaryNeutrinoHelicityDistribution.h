#pragma once
#ifndef LI_PrimaryNeutrinoHelicityDistribution_H
#define LI_PrimaryNeutrinoHelicityDistribution_H

#include <memory>                                        // for shared_ptr
#include <string>                                        // for string
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/Distributions.h"  // for InjectionDis...

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace distributions {

class PrimaryNeutrinoHelicityDistribution : virtual public PrimaryInjectionDistribution {
friend cereal::access;
public:
    virtual void Sample(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override;
	PrimaryNeutrinoHelicityDistribution();
	PrimaryNeutrinoHelicityDistribution(const PrimaryNeutrinoHelicityDistribution &) = default;
    virtual std::vector<std::string> DensityVariables() const override;
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryNeutrinoHelicityDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryNeutrinoHelicityDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace SI

CEREAL_CLASS_VERSION(SI::distributions::PrimaryNeutrinoHelicityDistribution, 0);
CEREAL_REGISTER_TYPE(SI::distributions::PrimaryNeutrinoHelicityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::PrimaryInjectionDistribution, SI::distributions::PrimaryNeutrinoHelicityDistribution);

#endif // LI_PrimaryNeutrinoHelicityDistribution_H
