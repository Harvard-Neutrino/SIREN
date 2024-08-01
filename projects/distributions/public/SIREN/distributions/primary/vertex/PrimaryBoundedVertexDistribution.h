#pragma once
#ifndef SIREN_PrimaryBoundedVertexDistribution_H
#define SIREN_PrimaryBoundedVertexDistribution_H

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
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace geometry { class Geometry; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class PrimaryBoundedVertexDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
private:
    std::shared_ptr<siren::geometry::Geometry> fiducial_volume = nullptr;
    double max_length = std::numeric_limits<double>::infinity();
public:

    PrimaryBoundedVertexDistribution();
    PrimaryBoundedVertexDistribution(const PrimaryBoundedVertexDistribution &) = default;
    PrimaryBoundedVertexDistribution(double max_length);
    PrimaryBoundedVertexDistribution(std::shared_ptr<siren::geometry::Geometry> fiducial_volume);
    PrimaryBoundedVertexDistribution(std::shared_ptr<siren::geometry::Geometry> fiducial_volume, double max_length);

    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;

    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("MaxLength", max_length));
            archive(::cereal::make_nvp("FidVol", fiducial_volume));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryBoundedVertexDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PrimaryBoundedVertexDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double max_length;
            std::shared_ptr<siren::geometry::Geometry> fiducial_volume;
            archive(::cereal::make_nvp("MaxLength", max_length));
            archive(::cereal::make_nvp("FidVol", fiducial_volume));
            construct(fiducial_volume,max_length);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PrimaryBoundedVertexDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PrimaryBoundedVertexDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::PrimaryBoundedVertexDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::PrimaryBoundedVertexDistribution);

#endif // SIREN_PrimaryBoundedVertexDistribution_H
