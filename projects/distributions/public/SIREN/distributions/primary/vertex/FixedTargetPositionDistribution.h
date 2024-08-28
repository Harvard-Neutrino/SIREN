#pragma once
#ifndef SIREN_FixedTargetPositionDistribution_H
#define SIREN_FixedTargetPositionDistribution_H

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

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class FixedTargetPositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    siren::geometry::Cylinder cylinder;
    std::shared_ptr<siren::geometry::Geometry> fiducial_volume = nullptr;
    double max_length = std::numeric_limits<double>::infinity();
    std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    FixedTargetPositionDistribution() {};
    FixedTargetPositionDistribution(siren::geometry::Cylinder cylinder, std::shared_ptr<siren::geometry::Geometry> fiducial_volume);
    FixedTargetPositionDistribution(siren::geometry::Cylinder cylinder, std::shared_ptr<siren::geometry::Geometry> fiducial_volume, double max_length);
    FixedTargetPositionDistribution(siren::geometry::Cylinder cylinder, double max_length);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    virtual bool AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Cylinder", cylinder));
            archive(::cereal::make_nvp("MaxLength", max_length));
            archive(::cereal::make_nvp("FidVol", fiducial_volume));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("FixedTargetPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<FixedTargetPositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            siren::geometry::Cylinder c;
            double max_length;
            std::shared_ptr<siren::geometry::Geometry> fiducial_volume;
            archive(::cereal::make_nvp("Cylinder", c));
            archive(::cereal::make_nvp("MaxLength", max_length));
            archive(::cereal::make_nvp("FidVol", fiducial_volume));
            construct(c,fiducial_volume,max_length);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("FixedTargetPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::FixedTargetPositionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::FixedTargetPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::FixedTargetPositionDistribution);

#endif // SIREN_FixedTargetPositionDistribution_H
