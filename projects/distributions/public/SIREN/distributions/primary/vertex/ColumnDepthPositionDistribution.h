#pragma once
#ifndef LI_ColumnDepthPositionDistribution_H
#define LI_ColumnDepthPositionDistribution_H

#include <set>
#include <tuple>
#include <memory>
#include <string>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/distributions/primary/vertex/DepthFunction.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace SI { namespace distributions { class WeightableDistribution; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace distributions {

class ColumnDepthPositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
protected:
    ColumnDepthPositionDistribution() {};
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DepthFunction> depth_function;
    std::set<SI::dataclasses::ParticleType> target_types;

    SI::math::Vector3D SampleFromDisk(std::shared_ptr<SI::utilities::LI_random> rand, SI::math::Vector3D const & dir) const;

    std::tuple<SI::math::Vector3D, SI::math::Vector3D> SamplePosition(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    std::tuple<SI::math::Vector3D, SI::math::Vector3D> GetSamplePosition(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record); 
    virtual double GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const override;
    ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::set<SI::dataclasses::ParticleType> target_types);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> InjectionBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("DepthFunction", depth_function));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("ColumnDepthPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<ColumnDepthPositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::set<SI::dataclasses::ParticleType> t;
            std::shared_ptr<DepthFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("DepthFunction", f));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, f, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("ColumnDepthPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace SI

CEREAL_CLASS_VERSION(SI::distributions::ColumnDepthPositionDistribution, 0);
CEREAL_REGISTER_TYPE(SI::distributions::ColumnDepthPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::distributions::VertexPositionDistribution, SI::distributions::ColumnDepthPositionDistribution);

#endif // LI_ColumnDepthPositionDistribution_H
