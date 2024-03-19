#ifndef SIREN_FixedDirection_H
#define SIREN_FixedDirection_H

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class FixedDirection : virtual public PrimaryDirectionDistribution {
friend cereal::access;
protected:
    FixedDirection() {};
private:
    siren::math::Vector3D dir;
public:
    FixedDirection(siren::math::Vector3D dir) : dir(dir) {};
private:
    siren::math::Vector3D SampleDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Direction", dir));
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("FixedDirection only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<FixedDirection> & construct, std::uint32_t const version) {
        if(version == 0) {
            siren::math::Vector3D d;
            archive(::cereal::make_nvp("Direction", d));
            construct(d);
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("FixedDirection only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::FixedDirection, 0);
CEREAL_REGISTER_TYPE(siren::distributions::FixedDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryDirectionDistribution, siren::distributions::FixedDirection);

#endif // SIREN_FixedDirection_H
