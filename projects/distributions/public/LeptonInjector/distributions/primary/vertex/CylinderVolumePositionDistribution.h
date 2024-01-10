#pragma once
#ifndef LI_CylinderVolumePositionDistribution_H
#define LI_CylinderVolumePositionDistribution_H

#include <memory>
#include <string>
#include <cstdint>
#include <utility>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"
#include "LeptonInjector/geometry/Cylinder.h"
#include "LeptonInjector/math/Vector3D.h"

namespace LI { namespace crosssections { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class EarthModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class CylinderVolumePositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
protected:
    CylinderVolumePositionDistribution() {};
private:
    LI::geometry::Cylinder cylinder;
    LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override;
    CylinderVolumePositionDistribution(LI::geometry::Cylinder);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & interaction) const override;
    virtual bool AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> second_cross_sections) const override;
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
            LI::geometry::Cylinder c;
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
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::CylinderVolumePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::CylinderVolumePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::VertexPositionDistribution, LI::distributions::CylinderVolumePositionDistribution);

#endif // LI_CylinderVolumePositionDistribution_H
