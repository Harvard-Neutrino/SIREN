#pragma once
#ifndef LI_TabulatedFluxDistribution_H
#define LI_TabulatedFluxDistribution_H

#include <memory>
#include <string>
#include <cstdint>
#include <stddef.h>
#include <stdexcept>
#include <vector>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "LeptonInjector/utilities/Interpolator.h"

namespace LI { namespace crosssections { class CrossSectionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class EarthModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace utilities { class LI_random; } }

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
    LI::utilities::Interpolator1D<double> inverseCdfTable;
    double integral;
    std::vector<double> cdf;
    std::vector<double> energy_nodes;
    std::vector<double> cdf_energy_nodes;
    const size_t burnin = 40; //original burnin parameter for MH sampling
    double unnormed_pdf(double energy) const;
    double pdf(double energy) const;
    void LoadFluxTable();
public:
    double SamplePDF(double energy) const; 
    double SampleUnnormedPDF(double energy) const; 
    double GetIntegral() const; 
    void ComputeCDF();
    std::vector<double> GetCDF() const;
    std::vector<double> GetEnergyNodes() const;
    std::vector<double> GetCDFEnergyNodes() const;
    double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const override;
    //double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record, const std::string& sampleMethod) const;
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

