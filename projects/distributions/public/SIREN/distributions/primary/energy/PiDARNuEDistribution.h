#pragma once
#ifndef SIREN_PiDARNuEDistribution_H
#define SIREN_PiDARNuEDistribution_H

#include <memory>
#include <string>
#include <cstdint>
#include <stddef.h>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/utilities/Constants.h"
#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class PiDARNuEDistribution : virtual public PrimaryEnergyDistribution {
friend cereal::access;
private:
    static double const & mu_mass;
public:
    PiDARNuEDistribution() = default;
    double pdf(double energy) const;
    double cdf(double energy) const;
    double inv_cdf(double u) const;
    double SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("PiDARNuEDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PiDARNuEDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            construct();
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PiDARNuEDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PiDARNuEDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::PiDARNuEDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryEnergyDistribution, siren::distributions::PiDARNuEDistribution);

#endif // SIREN_PiDARNuEDistribution_H
