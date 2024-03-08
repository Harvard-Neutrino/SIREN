#pragma once
#ifndef LI_OrientedPositionDistribution_H
#define LI_OrientedPositionDistribution_H

#include <tuple>
#include <memory>
#include <string>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace SI { namespace distributions { class WeightableDistribution; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace distributions {

class OrientedCylinderPositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    double radius;

    SI::math::Vector3D SampleFromDisk(std::shared_ptr<SI::utilities::LI_random> rand, SI::math::Vector3D const & dir) const;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> GetBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record, SI::math::Vector3D & point_of_closest_approach) const = 0;
    virtual void SampleWithinBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record, std::tuple<SI::math::Vector3D, SI::math::Vector3D> boundaries) const = 0;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> SamplePosition(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double PositionProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record, std::tuple<SI::math::Vector3D, SI::math::Vector3D> boundaries) const = 0;
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override;
    virtual std::string Name() const override = 0;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override = 0;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> InjectionBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & interaction) const override;
    virtual bool AreEquivalent(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<SI::detector::DetectorModel const> second_detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> second_interactions) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("OrientedCylinderPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("OrientedCylinderPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace SI

CEREAL_CLASS_VERSION(SI::distributions::OrientedCylinderPositionDistribution, 0);
CEREAL_REGISTER_TYPE(SI::distributions::OrientedCylinderPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::VertexPositionDistribution, SI::distributions::OrientedCylinderPositionDistribution);

#endif // LI_OrientedPositionDistribution_H
