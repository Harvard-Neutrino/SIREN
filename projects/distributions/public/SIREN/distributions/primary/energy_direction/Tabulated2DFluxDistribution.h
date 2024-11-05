#pragma once
#ifndef SIREN_Tabulated2DFluxDistribution_H
#define SIREN_Tabulated2DFluxDistribution_H

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

#include "SIREN/distributions/primary/energy_direction/PrimaryEnergyDirectionDistribution.h"
#include "SIREN/utilities/Interpolator.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class Tabulated2DFluxDistribution : virtual public PrimaryEnergyDirectionDistribution {
// Assumes table is in units of nu cm^-2 GeV^-1 Livetime^-1 sr^-1
friend cereal::access;
protected:
    Tabulated2DFluxDistribution();
    void ComputeIntegral();
    std::pair<double, double> SampleTablePoint(std::shared_ptr<siren::utilities::SIREN_random> rand) const;
private:
    double energyMin;
    double energyMax;
    double zenithMin;
    double zenithMax;
    bool energy_bounds_set;
    bool zenith_bounds_set;
    std::string fluxTableFilename;
    siren::utilities::Interpolator2D<double> fluxTable;
    double integral;
    std::vector<double> energy_nodes;
    std::vector<double> zenith_nodes;
    // metropolis hastings variables
    const size_t burnin = 40; // burnin parameter for MH sampling
    mutable size_t MH_sampled_points = 0; // to track the number of samples
    mutable double MH_density;

    double unnormed_pdf(double energy, double zenith) const;
    double pdf(double energy, double zenith) const;
    void LoadFluxTable();
    void LoadFluxTable(std::vector<double> & energies, std::vector<double> & zeniths, std::vector<double> & flux);
public:
    double SamplePDF(double energy, double zenith) const;
    double SampleUnnormedPDF(double energy, double zenith) const;
    double GetIntegral() const;
    std::vector<double> GetEnergyNodes() const;
    std::vector<double> GetZenithNodes() const;
    std::pair<double,siren::math::Vector3D> SampleEnergyAndDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    void SetEnergyBounds(double energyMin, double energyMax);
    void SetZenithBounds(double zenithMin, double zenithMax);
    Tabulated2DFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization=false);
    Tabulated2DFluxDistribution(double energyMin, double energyMax, std::string fluxTableFilename, bool has_physical_normalization=false);
    Tabulated2DFluxDistribution(double energyMin, double energyMax, double zenithMin, double zenithMax, std::string fluxTableFilename, bool has_physical_normalization=false);
    Tabulated2DFluxDistribution(std::vector<double> energies, std::vector<double> zeniths, std::vector<double> flux, bool has_physical_normalization=false);
    Tabulated2DFluxDistribution(double energyMin, double energyMax, std::vector<double> energies, std::vector<double> zeniths, std::vector<double> flux, bool has_physical_normalization=false);
    Tabulated2DFluxDistribution(double energyMin, double energyMax, double zenithMin, double zenithMax, std::vector<double> energies, std::vector<double> zeniths, std::vector<double> flux, bool has_physical_normalization=false);
    std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(::cereal::make_nvp("ZenithMin", zenithMin));
            archive(::cereal::make_nvp("ZenithMax", zenithMax));
            archive(::cereal::make_nvp("FluxTable", fluxTable));
            archive(cereal::virtual_base_class<PrimaryEnergyDirectionDistribution>(this));
        } else {
            throw std::runtime_error("Tabulated2DFluxDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(::cereal::make_nvp("ZenithMin", zenithMin));
            archive(::cereal::make_nvp("ZenithMax", zenithMax));
            archive(::cereal::make_nvp("FluxTable", fluxTable));
            archive(cereal::virtual_base_class<PrimaryEnergyDirectionDistribution>(this));
            energy_bounds_set = true;
            zenith_bounds_set = true;
            ComputeIntegral();
        } else {
            throw std::runtime_error("Tabulated2DFluxDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::Tabulated2DFluxDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::Tabulated2DFluxDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryEnergyDirectionDistribution, siren::distributions::Tabulated2DFluxDistribution);

#endif // SIREN_Tabulated2DFluxDistribution_H

