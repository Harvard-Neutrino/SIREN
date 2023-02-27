#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/detector/EarthModelCalculator.h"

#include "LeptonInjector/crosssections/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Particle.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/PrimaryEnergyDistribution.h"

#include "LeptonInjector/utilities/Errors.h"

namespace LI {
namespace distributions {
namespace {
    bool fexists(const std::string filename)
    {
            std::ifstream ifile(filename.c_str());
            return (bool)ifile;
    }
}

//---------------
// class PrimaryEnergyDistribution : InjectionDistribution
//---------------
void PrimaryEnergyDistribution::Sample(
        std::shared_ptr<LI::utilities::LI_random> rand,
        std::shared_ptr<LI::detector::EarthModel const> earth_model,
        std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections,
        LI::crosssections::InteractionRecord & record) const {
    record.primary_momentum[0] = SampleEnergy(rand, earth_model, cross_sections, record);
}

std::vector<std::string> PrimaryEnergyDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy"};
}

//---------------
// class PowerLaw : PrimaryEnergyDistribution
//---------------
PowerLaw::PowerLaw(double powerLawIndex, double energyMin, double energyMax)
    : powerLawIndex(powerLawIndex)
    , energyMin(energyMin)
    , energyMax(energyMax)
{}

double PowerLaw::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    if(energyMin == energyMax)
        return energyMin; //return the only allowed energy

    if(powerLawIndex == 1.0) //sample uniformly in log space
        return pow(10.0, rand->Uniform(log10(energyMin), log10(energyMax)));
    else {
        double u = rand->Uniform();
        double energyP = (1 - u) * pow(energyMin, 1 - powerLawIndex) + u * pow(energyMax, 1 - powerLawIndex);
        return pow(energyP, 1 / (1 - powerLawIndex));
    }
}

double PowerLaw::pdf(double energy) const {
    if(energyMin == energyMax)
        return 1.0; // only one allowed energy

    if(powerLawIndex == 1.0)
        return 1.0 / (energy * log(energyMax / energyMin));
    else {
        return pow(energy, -powerLawIndex) * (powerLawIndex - 1.0) * (pow(energyMin, powerLawIndex - 1.0) - pow(energyMax, powerLawIndex - 1.0));
    }
}

double PowerLaw::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    return pdf(record.primary_momentum[0]);
}

std::string PowerLaw::Name() const {
    return "PowerLaw";
}

std::shared_ptr<InjectionDistribution> PowerLaw::clone() const {
    return std::shared_ptr<InjectionDistribution>(new PowerLaw(*this));
}

bool PowerLaw::equal(WeightableDistribution const & other) const {
    const PowerLaw* x = dynamic_cast<const PowerLaw*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, powerLawIndex)
            ==
            std::tie(x->energyMin, x->energyMax, x->powerLawIndex);
}

bool PowerLaw::less(WeightableDistribution const & other) const {
    const PowerLaw* x = dynamic_cast<const PowerLaw*>(&other);
    return
        std::tie(energyMin, energyMax, powerLawIndex)
        <
        std::tie(x->energyMin, x->energyMax, x->powerLawIndex);
}

void PowerLaw::SetNormalizationAtEnergy(double norm, double energy) {
    SetNormalization(norm / pdf(energy));
}

//---------------
// class ModifiedMoyalPlusExponentialEnergyDistribution : PrimaryEnergyDistribution
//---------------

double ModifiedMoyalPlusExponentialEnergyDistribution::unnormed_pdf(double energy) const {
    double x = (energy - mu) / sigma;
    double moyal = (A / sigma) * std::exp(-(x + std::exp(-x))/2) / std::sqrt(2.0 * M_PI);
    double exponential = (B / l) * std::exp(-energy / l);
    return moyal + exponential;
}

double ModifiedMoyalPlusExponentialEnergyDistribution::pdf(double energy) const {
    return unnormed_pdf(energy) / integral;
}

ModifiedMoyalPlusExponentialEnergyDistribution::ModifiedMoyalPlusExponentialEnergyDistribution(double energyMin, double energyMax, double mu, double sigma, double A, double l, double B, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , mu(mu)
    , sigma(sigma)
    , A(A)
    , l(l)
    , B(B)
{
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    integral = LI::detector::Integration::rombergIntegrate(integrand, energyMin, energyMax);
    if(has_physical_normalization) {
        SetNormalization(integral);
    }
}

double ModifiedMoyalPlusExponentialEnergyDistribution::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    // Metropolis-Hastings algorithm to sample from PDF.
    // Pass in a function pointer for the PDF

    double energy, density, test_energy, test_density, odds;
    bool accept;

    // sample an initial point uniformly
    energy = rand->Uniform(energyMin, energyMax);
    density = pdf(energy);

    // Metropolis Hastings loop
    for (size_t j = 0; j <= burnin; ++j) {
        test_energy = rand->Uniform(energyMin, energyMax);
        test_density = pdf(test_energy);
        odds = test_density / density;
        accept = (odds > 1.) or (rand->Uniform(0,1) < odds);
        if(accept) {
            energy = test_energy;
            density = test_density;
        }
    }

    return energy;
}

double ModifiedMoyalPlusExponentialEnergyDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string ModifiedMoyalPlusExponentialEnergyDistribution::Name() const {
    return "ModifiedMoyalPlusExponentialEnergyDistribution";
}

std::shared_ptr<InjectionDistribution> ModifiedMoyalPlusExponentialEnergyDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ModifiedMoyalPlusExponentialEnergyDistribution(*this));
}

bool ModifiedMoyalPlusExponentialEnergyDistribution::equal(WeightableDistribution const & other) const {
    const ModifiedMoyalPlusExponentialEnergyDistribution* x = dynamic_cast<const ModifiedMoyalPlusExponentialEnergyDistribution*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, mu, sigma, A, l, B)
            ==
            std::tie(x->energyMin, x->energyMax, x->mu, x->sigma, x->A, x->l, x->B);
}

bool ModifiedMoyalPlusExponentialEnergyDistribution::less(WeightableDistribution const & other) const {
    const ModifiedMoyalPlusExponentialEnergyDistribution* x = dynamic_cast<const ModifiedMoyalPlusExponentialEnergyDistribution*>(&other);
    return
        std::tie(energyMin, energyMax, mu, sigma, A, l, B)
        <
        std::tie(x->energyMin, x->energyMax, x->mu, x->sigma, x->A, x->l, x->B);
}

//---------------
// class TabulatedFluxDistribution : PrimaryEnergyDistribution
//---------------
TabulatedFluxDistribution::TabulatedFluxDistribution() {}

void TabulatedFluxDistribution::ComputeIntegral() {
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    integral = LI::detector::Integration::rombergIntegrate(integrand, energyMin, energyMax);
}

void TabulatedFluxDistribution::LoadFluxTable() {
    if(fexists(fluxTableFilename)) {
        std::ifstream in(fluxTableFilename.c_str());
        std::string buf;
        std::string::size_type pos;
        LI::utilities::TableData1D<double> table_data;
        while(std::getline(in, buf)) {
            // Ignore comments and blank lines
            if((pos = buf.find('#')) != std::string::npos)
                buf.erase(pos);
            const char* whitespace=" \n\r\t\v";
            if((pos=buf.find_first_not_of(whitespace))!=0)
                buf.erase(0,pos);
            if(!buf.empty() && (pos=buf.find_last_not_of(whitespace))!=buf.size()-1)
                buf.erase(pos+1);
            if(buf.empty())
                continue;

            std::stringstream ss(buf);
            double x, f;
            ss >> x >> f;
            table_data.x.push_back(x);
            table_data.f.push_back(f);
        }
        // If no physical are manually set, use first/last entry of table
        if(not bounds_set) {
            energyMin = table_data.x[0];
            energyMax = table_data.x[table_data.x.size()-1];
        }
        fluxTable = LI::utilities::Interpolator1D<double>(table_data);
    } else {
        throw std::runtime_error("Failed to open flux table file!");
    }
}

double TabulatedFluxDistribution::unnormed_pdf(double energy) const {
    return fluxTable(energy);
}

double TabulatedFluxDistribution::pdf(double energy) const {
    return unnormed_pdf(energy) / integral;
}

void TabulatedFluxDistribution::SetEnergyBounds(double eMin, double eMax) {
    energyMin = eMin;
    energyMax = eMax;
    bounds_set = true;
    ComputeIntegral();
}

TabulatedFluxDistribution::TabulatedFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization)
    : bounds_set(false)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

TabulatedFluxDistribution::TabulatedFluxDistribution(double energyMin, double energyMax, std::string fluxTableFilename, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , bounds_set(true)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

double TabulatedFluxDistribution::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    // Metropolis-Hastings algorithm to sample from PDF.
    // Pass in a function pointer for the PDF

    double energy, density, test_energy, test_density, odds;
    bool accept;

    // sample an initial point uniformly
    energy = rand->Uniform(energyMin, energyMax);
    density = pdf(energy);

    // Metropolis Hastings loop
    for (size_t j = 0; j <= burnin; ++j) {
        test_energy = rand->Uniform(energyMin, energyMax);
        test_density = pdf(test_energy);
        odds = test_density / density;
        accept = (odds > 1.) or (rand->Uniform(0,1) < odds);
        if(accept) {
            energy = test_energy;
            density = test_density;
        }
    }

    return energy;
}

double TabulatedFluxDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string TabulatedFluxDistribution::Name() const {
    return "TabulatedFluxDistribution";
}

std::shared_ptr<InjectionDistribution> TabulatedFluxDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new TabulatedFluxDistribution(*this));
}

bool TabulatedFluxDistribution::equal(WeightableDistribution const & other) const {
    const TabulatedFluxDistribution* x = dynamic_cast<const TabulatedFluxDistribution*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, fluxTableFilename)
            ==
            std::tie(x->energyMin, x->energyMax, x->fluxTableFilename);
}

bool TabulatedFluxDistribution::less(WeightableDistribution const & other) const {
    const TabulatedFluxDistribution* x = dynamic_cast<const TabulatedFluxDistribution*>(&other);
    return
        std::tie(energyMin, energyMax, fluxTable)
        <
        std::tie(x->energyMin, x->energyMax, x->fluxTable);
}

} // namespace distributions
} // namespace LeptonInjector
