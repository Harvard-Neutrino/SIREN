#pragma once
#ifndef LI_VertexPositionDistribution_H
#define LI_VertexPositionDistribution_H

#include <tuple>                                         // for tuple
#include <memory>                                        // for shared_ptr
#include <string>                                        // for string
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/InteractionTree.h"  // for InteractionT...
#include "SIREN/distributions/Distributions.h"  // for WeightableDi...
#include "SIREN/math/Vector3D.h"                // for Vector3D

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace distributions {

class VertexPositionDistribution : virtual public PrimaryInjectionDistribution {
friend cereal::access;
public:
    virtual ~VertexPositionDistribution() {};
private:
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> SamplePosition(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const = 0;
public:
    virtual void Sample(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override = 0;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override = 0;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> InjectionBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & interaction) const = 0;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> InjectionBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionTreeDatum const & datum) const;
    virtual bool AreEquivalent(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<SI::detector::DetectorModel const> second_detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> second_interactions) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("VertexPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("VertexPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace SI

CEREAL_CLASS_VERSION(SI::distributions::VertexPositionDistribution, 0);
CEREAL_REGISTER_TYPE(SI::distributions::VertexPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::PrimaryInjectionDistribution, SI::distributions::VertexPositionDistribution);

#endif // LI_VertexPositionDistribution_H
