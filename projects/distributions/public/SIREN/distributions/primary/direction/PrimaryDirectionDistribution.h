#pragma once
#ifndef SIREN_PrimaryDirectionDistribution_H
#define SIREN_PrimaryDirectionDistribution_H

#include <string>                                        // for string
#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include "SIREN/distributions/Distributions.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace math { class Vector3D; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace utilities {
class SIREN_random;
} // namespace utilities

namespace math {
class Vector3D;
} // namespace math

namespace detector {
class DetectorModel;
} // namespace detector

namespace dataclasses {
class InteractionRecord;
struct InteractionSignature;
}
namespace interactions {
class InteractionCollection;
} // namespace interactions
} // namespace sirenREN

namespace siren {
namespace distributions {

class PrimaryDirectionDistribution : virtual public PrimaryInjectionDistribution {
friend cereal::access;
public:
    virtual ~PrimaryDirectionDistribution() {};
protected:
    PrimaryDirectionDistribution() {};
private:
    virtual siren::math::Vector3D SampleDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const = 0;
public:
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const override = 0;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PrimaryDirectionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::PrimaryDirectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryInjectionDistribution, siren::distributions::PrimaryDirectionDistribution);

#endif // SIREN_PrimaryDirectionDistribution_H
