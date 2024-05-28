#pragma once
#ifndef SIREN_SecondaryPhysicalVertexDistribution_H
#define SIREN_SecondaryPhysicalVertexDistribution_H

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
#include <cereal/types/memory.hpp>

#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class InjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace geometry { class Geometry; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class SecondaryPhysicalVertexDistribution : virtual public SecondaryVertexPositionDistribution {
friend cereal::access;
public:

    SecondaryPhysicalVertexDistribution();
    SecondaryPhysicalVertexDistribution(const SecondaryPhysicalVertexDistribution &) = default;

    virtual void SampleVertex(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;

    std::string Name() const override;
    virtual std::shared_ptr<SecondaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;

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
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::SecondaryPhysicalVertexDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::SecondaryPhysicalVertexDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::SecondaryVertexPositionDistribution, siren::distributions::SecondaryPhysicalVertexDistribution);

#endif // SIREN_SecondaryPhysicalVertexDistribution_H
