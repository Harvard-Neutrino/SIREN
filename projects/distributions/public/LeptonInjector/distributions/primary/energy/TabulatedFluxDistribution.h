#pragma once
#ifndef LI_TabulatedFluxDistribution_H
#define LI_TabulatedFluxDistribution_H

#include <string>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/utilities/Interpolator.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"

namespace LI {
namespace distributions {

class TabulatedFluxDistribution : virtual public PrimaryEnergyDistribution {
// Assumes table is in units of nu cm^-2 GeV^-1 Livetime^-1
friend cereal::access;
protected:
    TabulatedFluxDistribution();
    void ComputeIntegral();
private:
    double energyMin;
    double energyMax;
    bool bounds_set;
    std::string fluxTableFilename;
    LI::utilities::Interpolator1D<double> fluxTable;
    double integral;
    const size_t burnin = 40;
    double unnormed_pdf(double energy) const;
    double pdf(double energy) const;
    void LoadFluxTable();
public:
    double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override;
    void SetEnergyBounds(double energyMin, double energyMax);
    TabulatedFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization=false);
    TabulatedFluxDistribution(double energyMin, double energyMax, std::string fluxTableFilename, bool has_physical_normalization=false);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(::cereal::make_nvp("FluxTable", fluxTable));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("TabulatedFluxDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(::cereal::make_nvp("FluxTable", fluxTable));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
            ComputeIntegral();
        } else {
            throw std::runtime_error("TabulatedFluxDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::TabulatedFluxDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::TabulatedFluxDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryEnergyDistribution, LI::distributions::TabulatedFluxDistribution);

#endif // LI_TabulatedFluxDistribution_H

