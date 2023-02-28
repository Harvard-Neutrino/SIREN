#pragma once
#ifndef LI_OrientedPositionDistribution_H
#define LI_OrientedPositionDistribution_H

#include <memory>
#include <string>
#include <utility>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace detector {
class EarthModel;
} // namespace detector

namespace dataclasses {
struct InteractionRecord;
struct InteractionSignature;
}
namespace crosssections {
class CrossSectionCollection;
} // namespace crosssections
} // namespace LeptonInjector

namespace LI {
namespace distributions {

class OrientedCylinderPositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    double radius;

    LI::math::Vector3D SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> GetBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record, LI::math::Vector3D & point_of_closest_approach) const = 0;
    virtual void SampleWithinBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record, std::pair<LI::math::Vector3D, LI::math::Vector3D> boundaries) const = 0;
    virtual LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const;
public:
    virtual double PositionProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record, std::pair<LI::math::Vector3D, LI::math::Vector3D> boundaries) const = 0;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const;
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & interaction) const;
    virtual bool AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> second_cross_sections) const;
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
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::OrientedCylinderPositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::OrientedCylinderPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::VertexPositionDistribution, LI::distributions::OrientedCylinderPositionDistribution);

#endif // LI_OrientedPositionDistribution_H
