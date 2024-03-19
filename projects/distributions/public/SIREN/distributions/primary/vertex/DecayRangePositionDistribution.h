#pragma once
#ifndef SIREN_DecayRangePositionDistribution_H
#define SIREN_DecayRangePositionDistribution_H

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
#include "SIREN/distributions/primary/vertex/DecayRangeFunction.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class DecayRangePositionDistribution : virtual public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DecayRangeFunction> range_function;

    siren::math::Vector3D SampleFromDisk(std::shared_ptr<siren::utilities::SIREN_random> rand, siren::math::Vector3D const & dir) const;

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    DecayRangePositionDistribution();
    DecayRangePositionDistribution(const DecayRangePositionDistribution &) = default;
    DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function);
    std::string Name() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual bool AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("DecayRangeFunction", range_function));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("DecayRangePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<DecayRangePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::shared_ptr<DecayRangeFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("DecayRangeFunction", f));
            construct(r, l, f);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("DecayRangePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::DecayRangePositionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::DecayRangePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::DecayRangePositionDistribution);

#endif // SIREN_DecayRangePositionDistribution_H
