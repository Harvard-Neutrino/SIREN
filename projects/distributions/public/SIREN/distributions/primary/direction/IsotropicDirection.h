#pragma once
#ifndef SIREN_IsotropicDirection_H
#define SIREN_IsotropicDirection_H

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

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace math {
class Vector3D;
}
namespace distributions {

class IsotropicDirection : virtual public PrimaryDirectionDistribution {
friend cereal::access;
public:
    siren::math::Vector3D SampleDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
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
} // namespace siren


CEREAL_CLASS_VERSION(siren::distributions::IsotropicDirection, 0);
CEREAL_REGISTER_TYPE(siren::distributions::IsotropicDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryDirectionDistribution, siren::distributions::IsotropicDirection);

#endif // SIREN_IsotropicDirection_H
