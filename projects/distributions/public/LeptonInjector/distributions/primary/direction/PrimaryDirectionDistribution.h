#pragma once
#ifndef LI_PrimaryDirectionDistribution_H
#define LI_PrimaryDirectionDistribution_H

#include <string>                                        // for string
#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include "LeptonInjector/distributions/Distributions.h"  // for WeightableDi...

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace math { class Vector3D; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace math {
class Vector3D;
} // namespace math

namespace detector {
class DetectorModel;
} // namespace detector

namespace dataclasses {
struct InteractionRecord;
struct InteractionSignature;
}
namespace interactions {
class InteractionCollection;
} // namespace interactions
} // namespace LeptonInjector

namespace LI {
namespace distributions {

class PrimaryDirectionDistribution : virtual public InjectionDistribution {
friend cereal::access;
public:
    virtual ~PrimaryDirectionDistribution() {};
protected:
    PrimaryDirectionDistribution() {};
private:
    virtual LI::math::Vector3D SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PrimaryDirectionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryDirectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::PrimaryDirectionDistribution);

#endif // LI_PrimaryDirectionDistribution_H
