#pragma once
#ifndef SIREN_PointSourcePositionDistribution_H
#define SIREN_PointSourcePositionDistribution_H

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
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class PointSourcePositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    siren::math::Vector3D origin;
    double max_distance;
    std::set<siren::dataclasses::ParticleType> target_types;

    siren::math::Vector3D SampleFromDisk(std::shared_ptr<siren::utilities::SIREN_random> rand, siren::math::Vector3D const & dir) const;

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    PointSourcePositionDistribution();
    PointSourcePositionDistribution(const PointSourcePositionDistribution &) = default;
    PointSourcePositionDistribution(siren::math::Vector3D origin, double max_distance, std::set<siren::dataclasses::ParticleType> target_types);
    std::string Name() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Origin", origin));
            archive(::cereal::make_nvp("MaxDistance", max_distance));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("PointSourcePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PointSourcePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            siren::math::Vector3D r;
            double l;
            std::set<siren::dataclasses::ParticleType> t;
            archive(::cereal::make_nvp("Origin", r));
            archive(::cereal::make_nvp("MaxDistance", l));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PointSourcePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PointSourcePositionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::PointSourcePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::PointSourcePositionDistribution);

#endif // SIREN_PointSourcePositionDistribution_H
