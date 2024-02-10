#pragma once
#ifndef LI_PointSourcePositionDistribution_H
#define LI_PointSourcePositionDistribution_H

#include <set>
#include <memory>
#include <string>
#include <utility>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"
#include "LeptonInjector/math/Vector3D.h"

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class PointSourcePositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    LI::math::Vector3D origin;
    double max_distance;
    std::set<LI::dataclasses::Particle::ParticleType> target_types;

    LI::math::Vector3D SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const;

    LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override;
    PointSourcePositionDistribution();
    PointSourcePositionDistribution(const PointSourcePositionDistribution &) = default;
    PointSourcePositionDistribution(LI::math::Vector3D origin, double max_distance, std::set<LI::dataclasses::Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
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
            LI::math::Vector3D r;
            double l;
            std::set<LI::dataclasses::Particle::ParticleType> t;
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
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PointSourcePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PointSourcePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::VertexPositionDistribution, LI::distributions::PointSourcePositionDistribution);

#endif // LI_PointSourcePositionDistribution_H
