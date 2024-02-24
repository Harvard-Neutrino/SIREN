#pragma once
#ifndef LI_SecondaryPhysicalVertexDistribution_H
#define LI_SecondaryPhysicalVertexDistribution_H

#include <tuple>
#include <limits>
#include <memory>
#include <string>
#include <cstdint>
#include <stddef.h>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/InteractionTree.h"
#include "LeptonInjector/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"
#include "LeptonInjector/math/Vector3D.h"

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace geometry { class Geometry; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class SecondaryPhysicalVertexDistribution : virtual public SecondaryVertexPositionDistribution {
friend cereal::access;
public:

    SecondaryPhysicalVertexDistribution();
    SecondaryPhysicalVertexDistribution(const SecondaryPhysicalVertexDistribution &) = default;

    virtual void SampleVertex(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::SecondaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override;

    std::string Name() const override;
    virtual std::shared_ptr<SecondaryInjectionDistribution> clone() const override;
    virtual std::tuple<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & interaction) const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<SecondaryVertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("SecondaryPhysicalVertexDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<SecondaryPhysicalVertexDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            construct();
            archive(cereal::virtual_base_class<SecondaryVertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("SecondaryPhysicalVertexDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::SecondaryPhysicalVertexDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::SecondaryPhysicalVertexDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::SecondaryVertexPositionDistribution, LI::distributions::SecondaryPhysicalVertexDistribution);

#endif // LI_SecondaryPhysicalVertexDistribution_H
