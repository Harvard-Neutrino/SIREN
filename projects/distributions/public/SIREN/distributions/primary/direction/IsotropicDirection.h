#pragma once
#ifndef LI_IsotropicDirection_H
#define LI_IsotropicDirection_H

#include <memory>
#include <string>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace SI { namespace distributions { class WeightableDistribution; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace math {
class Vector3D;
}
namespace distributions {

class IsotropicDirection : virtual public PrimaryDirectionDistribution {
friend cereal::access;
public:
    SI::math::Vector3D SampleDirection(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("IsotropicDirection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("IsotropicDirection only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace SI


CEREAL_CLASS_VERSION(SI::distributions::IsotropicDirection, 0);
CEREAL_REGISTER_TYPE(SI::distributions::IsotropicDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::PrimaryDirectionDistribution, SI::distributions::IsotropicDirection);

#endif // LI_IsotropicDirection_H
