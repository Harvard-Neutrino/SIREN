#pragma once
#ifndef LI_VertexPositionDistribution_H
#define LI_VertexPositionDistribution_H

#include <memory>                                        // for shared_ptr
#include <string>                                        // for string
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <utility>                                       // for pair
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/InteractionTree.h"  // for InteractionT...
#include "LeptonInjector/distributions/Distributions.h"  // for WeightableDi...
#include "LeptonInjector/math/Vector3D.h"                // for Vector3D

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class VertexPositionDistribution : virtual public InjectionDistribution {
friend cereal::access;
public:
    virtual ~VertexPositionDistribution() {};
private:
    virtual LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const = 0;
    virtual LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum & datum) const;
public:
    virtual void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const override;
    virtual bool IsPositionDistribution() const override {return true;}
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const override = 0;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & interaction) const = 0;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum const & datum) const;
    virtual bool AreEquivalent(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::DetectorModel const> second_earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> second_cross_sections) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("VertexPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("VertexPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::VertexPositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::VertexPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::VertexPositionDistribution);

#endif // LI_VertexPositionDistribution_H
