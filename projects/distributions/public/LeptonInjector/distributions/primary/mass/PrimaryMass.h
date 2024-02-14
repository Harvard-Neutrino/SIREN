#pragma once
#ifndef LI_PrimaryMass_H
#define LI_PrimaryMass_H

#include <string>                                        // for string
#include <memory>                                        // for shared_ptr
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error
#include <vector>                                        // for vector

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/Particle.h"         // for Particle
#include "LeptonInjector/distributions/Distributions.h"  // for InjectionDis...

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class PrimaryMass : virtual public PrimaryInjectionDistribution {
friend cereal::access;
protected:
    PrimaryMass() {};
private:
    double primary_mass;
public:
    PrimaryMass(double primary_mass = 0);
    double GetPrimaryMass() const;
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryMass only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PrimaryMass> & construct, std::uint32_t const version) {
        if(version == 0) {
            double m;
            archive(::cereal::make_nvp("PrimaryMass", m));
            construct(m);
            archive(cereal::virtual_base_class<PrimaryInjectionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PrimaryMass only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PrimaryMass, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryMass);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryInjectionDistribution, LI::distributions::PrimaryMass);

#endif // LI_PrimaryMass_H
