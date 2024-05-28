#pragma once
#ifndef SIREN_TabulatedFluxDistribution_H
#define SIREN_TabulatedFluxDistribution_H

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

#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "SIREN/utilities/Interpolator.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
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
    siren::utilities::Interpolator1D<double> fluxTable;
    siren::utilities::Interpolator1D<double> inverseCdfTable;
    double integral;
    std::vector<double> cdf;
    std::vector<double> energy_nodes;
    std::vector<double> cdf_energy_nodes;
    const size_t burnin = 40; //original burnin parameter for MH sampling
    double unnormed_pdf(double energy) const;
    double pdf(double energy) const;
    void LoadFluxTable();
    void LoadFluxTable(std::vector<double> & energies, std::vector<double> & flux);
public:
    double SamplePDF(double energy) const; 
    double SampleUnnormedPDF(double energy) const; 
    double GetIntegral() const; 
    void ComputeCDF();
    std::vector<double> GetCDF() const;
    std::vector<double> GetEnergyNodes() const;
    std::vector<double> GetCDFEnergyNodes() const;
    double SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    void SetEnergyBounds(double energyMin, double energyMax);
    TabulatedFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization=false);
    TabulatedFluxDistribution(double energyMin, double energyMax, std::string fluxTableFilename, bool has_physical_normalization=false);
    TabulatedFluxDistribution(std::vector<double> energies, std::vector<double> flux, bool has_physical_normalization=false);
    TabulatedFluxDistribution(double energyMin, double energyMax, std::vector<double> energies, std::vector<double> flux, bool has_physical_normalization=false);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
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
            bounds_set = true;
            ComputeIntegral();
            ComputeCDF();
        } else {
            throw std::runtime_error("TabulatedFluxDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::TabulatedFluxDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::TabulatedFluxDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryEnergyDistribution, siren::distributions::TabulatedFluxDistribution);

#endif // SIREN_TabulatedFluxDistribution_H

