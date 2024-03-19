#pragma once
#ifndef SIREN_ModifiedMoyalPlusExponentialEnergyDistribution_H
#define SIREN_ModifiedMoyalPlusExponentialEnergyDistribution_H

#include <memory>
#include <string>
#include <cstdint>
#include <stddef.h>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class ModifiedMoyalPlusExponentialEnergyDistribution : virtual public PrimaryEnergyDistribution {
friend cereal::access;
protected:
    ModifiedMoyalPlusExponentialEnergyDistribution() {};
private:
    double energyMin;
    double energyMax;
    double mu;
    double sigma;
    double A;
    double l;
    double B;
    double integral;
    const size_t burnin = 40;
    double unnormed_pdf(double energy) const;
    double pdf(double energy) const;
    double pdf_integral() const;
public:
    double SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    ModifiedMoyalPlusExponentialEnergyDistribution(double energyMin, double energyMax, double mu, double sigma, double A, double l, double B, bool has_physical_normalization=false);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(::cereal::make_nvp("ParameterMu", mu));
            archive(::cereal::make_nvp("ParameterSigma", sigma));
            archive(::cereal::make_nvp("ParameterA", A));
            archive(::cereal::make_nvp("ParameterL", l));
            archive(::cereal::make_nvp("ParameterB", B));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("ModifiedMoyalPlusExponentialEnergyDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<ModifiedMoyalPlusExponentialEnergyDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double min, max, mu, s, a, l, b;
            archive(::cereal::make_nvp("EnergyMin", min));
            archive(::cereal::make_nvp("EnergyMax", max));
            archive(::cereal::make_nvp("ParameterMu", mu));
            archive(::cereal::make_nvp("ParameterSigma", s));
            archive(::cereal::make_nvp("ParameterA", a));
            archive(::cereal::make_nvp("ParameterL", l));
            archive(::cereal::make_nvp("ParameterB", b));
            construct(min, max, mu, s, a, l, b);
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("ModifiedMoyalPlusExponentialEnergyDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::ModifiedMoyalPlusExponentialEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::ModifiedMoyalPlusExponentialEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryEnergyDistribution, siren::distributions::ModifiedMoyalPlusExponentialEnergyDistribution);

#endif // SIREN_ModifiedMoyalPlusExponentialEnergyDistribution_H
