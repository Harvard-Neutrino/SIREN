#pragma once
#ifndef LI_CylinderVolumePositionDistribution_H
#define LI_CylinderVolumePositionDistribution_H

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
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace SI { namespace distributions { class WeightableDistribution; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace distributions {

class CylinderVolumePositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
protected:
    CylinderVolumePositionDistribution() {};
private:
    SI::geometry::Cylinder cylinder;
    std::tuple<SI::math::Vector3D, SI::math::Vector3D> SamplePosition(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override;
    CylinderVolumePositionDistribution(SI::geometry::Cylinder);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> InjectionBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & interaction) const override;
    virtual bool AreEquivalent(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<SI::detector::DetectorModel const> second_detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> second_interactions) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Cylinder", cylinder));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("CylinderVolumePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<CylinderVolumePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            SI::geometry::Cylinder c;
            archive(::cereal::make_nvp("Cylinder", c));
            construct(c);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("CylinderVolumePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace SI

CEREAL_CLASS_VERSION(SI::distributions::CylinderVolumePositionDistribution, 0);
CEREAL_REGISTER_TYPE(SI::distributions::CylinderVolumePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::VertexPositionDistribution, SI::distributions::CylinderVolumePositionDistribution);

#endif // LI_CylinderVolumePositionDistribution_H
