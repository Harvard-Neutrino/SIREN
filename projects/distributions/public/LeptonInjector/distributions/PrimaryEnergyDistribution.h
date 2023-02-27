#pragma once
#ifndef LI_PrimaryEnergyDistribution_H
#define LI_PrimaryEnergyDistribution_H

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

#include "LeptonInjector/utilities/Particle.h"
#include "LeptonInjector/utilities/Interpolator.h"

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace crosssections {
struct InteractionRecord;
struct InteractionSignature;
class CrossSectionCollection;
} // namespace crosssections
} // namespace LeptonInjector

namespace LI {
namespace distributions {

class PrimaryEnergyDistribution : virtual public InjectionDistribution, virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
private:
    virtual double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
            archive(cereal::virtual_base_class<PhysicallyNormalizedDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class PowerLaw : virtual public PrimaryEnergyDistribution {
friend cereal::access;
protected:
    PowerLaw() {};
private:
    double powerLawIndex;
    double energyMin;
    double energyMax;
public:
    PowerLaw(double powerLawIndex, double energyMin, double energyMax);
    double pdf(double energy) const;
    double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    std::string Name() const override;
    void SetNormalizationAtEnergy(double normalization, double energy);
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PowerLawIndex", powerLawIndex));
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("PowerLaw only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PowerLaw> & construct, std::uint32_t const version) {
        if(version == 0) {
            double gamma, min, max;
            archive(::cereal::make_nvp("PowerLawIndex", gamma));
            archive(::cereal::make_nvp("EnergyMin", min));
            archive(::cereal::make_nvp("EnergyMax", max));
            construct(gamma, min, max);
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PowerLaw only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

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
            archive(::cereal::make_nvp("ParameteB", B));
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
            archive(::cereal::make_nvp("ParameteB", b));
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
    double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
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

CEREAL_CLASS_VERSION(LI::distributions::PrimaryEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PhysicallyNormalizedDistribution, LI::distributions::PrimaryEnergyDistribution);

CEREAL_CLASS_VERSION(LI::distributions::PowerLaw, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PowerLaw);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryEnergyDistribution, LI::distributions::PowerLaw);

CEREAL_CLASS_VERSION(LI::distributions::ModifiedMoyalPlusExponentialEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::ModifiedMoyalPlusExponentialEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryEnergyDistribution, LI::distributions::ModifiedMoyalPlusExponentialEnergyDistribution);

CEREAL_CLASS_VERSION(LI::distributions::TabulatedFluxDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::TabulatedFluxDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryEnergyDistribution, LI::distributions::TabulatedFluxDistribution);

#endif // LI_PrimaryEnergyDistribution_H
