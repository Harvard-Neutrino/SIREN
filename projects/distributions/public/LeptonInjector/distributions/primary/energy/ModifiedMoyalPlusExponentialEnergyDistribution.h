#pragma once
#ifndef LI_ModifiedMoyalPlusExponentialEnergyDistribution_H
#define LI_ModifiedMoyalPlusExponentialEnergyDistribution_H

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"

namespace LI {
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
    double unnormed_pdf(double energy) const ;
    double pdf(double energy) const;
public:
    double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    ModifiedMoyalPlusExponentialEnergyDistribution(double energyMin, double energyMax, double mu, double sigma, double A, double l, double B, bool has_physical_normalization=false);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
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
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::ModifiedMoyalPlusExponentialEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::ModifiedMoyalPlusExponentialEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryEnergyDistribution, LI::distributions::ModifiedMoyalPlusExponentialEnergyDistribution);

#endif // LI_ModifiedMoyalPlusExponentialEnergyDistribution_H
