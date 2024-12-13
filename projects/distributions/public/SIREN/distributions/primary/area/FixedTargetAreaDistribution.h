#pragma once
#ifndef SIREN_FixedTargetAreaDistribution_H
#define SIREN_FixedTargetAreaDistribution_H

#include <tuple>
#include <memory>
#include <string>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/primary/area/PrimaryAreaDistribution.h"
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

class FixedTargetAreaDistribution : virtual public PrimaryAreaDistribution {
friend cereal::access;
private:
    siren::geometry::Cylinder cylinder;
    siren::math::Vector3D SamplePointOfClosestApproach(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    FixedTargetAreaDistribution() {};
    FixedTargetAreaDistribution(siren::geometry::Cylinder cylinder);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual bool AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Cylinder", cylinder));
            archive(cereal::virtual_base_class<PrimaryAreaDistribution>(this));
        } else {
            throw std::runtime_error("FixedTargetAreaDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<FixedTargetAreaDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            siren::geometry::Cylinder c;
            double max_length;
            std::shared_ptr<siren::geometry::Geometry> fiducial_volume;
            archive(::cereal::make_nvp("Cylinder", c));
            construct(c);
            archive(cereal::virtual_base_class<PrimaryAreaDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("FixedTargetAreaDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::FixedTargetAreaDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::FixedTargetAreaDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryAreaDistribution, siren::distributions::FixedTargetAreaDistribution);

#endif // SIREN_FixedTargetAreaDistribution_H
