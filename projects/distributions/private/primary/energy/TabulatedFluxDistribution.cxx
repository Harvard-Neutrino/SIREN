#include "LeptonInjector/distributions/primary/energy/TabulatedFluxDistribution.h"
#include <array>                                           // for array
#include <tuple>                                           // for tie, opera...
#include <vector>                                          // for vector
#include <fstream>                                         // for basic_istream
#include <functional>                                      // for function

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/distributions/Distributions.h"    // for InjectionD...
#include "LeptonInjector/utilities/Integration.h"          // for rombergInt...
#include "LeptonInjector/utilities/Interpolator.h"         // for TableData1D
#include "LeptonInjector/utilities/Random.h"               // for LI_random

namespace LI { namespace crosssections { class CrossSectionCollection; } }
namespace LI { namespace detector { class EarthModel; } }

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
// class TabulatedFluxDistribution : PrimaryEnergyDistribution
//---------------
TabulatedFluxDistribution::TabulatedFluxDistribution() {}

void TabulatedFluxDistribution::ComputeIntegral() {
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    integral = LI::utilities::rombergIntegrate(integrand, energyMin, energyMax);
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

            energy_nodes.push_back(x);
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
    
double TabulatedFluxDistribution::SampleUnnormedPDF(double energy) const {
    return unnormed_pdf(energy);
}

double TabulatedFluxDistribution::GetIntegral() const {
    return integral;
}

std::vector<double> TabulatedFluxDistribution::GetEnergyNodes() const {
    return energy_nodes;
}

double TabulatedFluxDistribution::pdf(double energy) const {
    return unnormed_pdf(energy) / integral;
}

double TabulatedFluxDistribution::SamplePDF(double energy) const {
    return pdf(energy);
}

void TabulatedFluxDistribution::ComputeCDF() {

    // this part makes an energies vector that is within the specified energy range
    std::vector<double> energies;
    energies.push_back(energyMin);
    // Use the copy_if algorithm to copy only the elements that satisfy the condition
    std::copy_if(energy_nodes.begin(), energy_nodes.end(), std::back_inserter(energies),
       [this](double value) {
            return value > this->energyMin && value < this->energyMax;
       }
    );
    energies.push_back(energyMax);
    
    // assign the energies vector so it's accessible outside of the function
    cdf_energy_nodes = energies;

    // declare the cdf vector
    std::vector<double> cdfvector(energies.size());

    cdfvector[0] = 0; 
    for (std::size_t i = 1; i < energies.size(); ++i) {
       double y = pdf(energies[i]);
       double area = y * (energies[i]-energies[i-1]);
       cdfvector[i] = cdfvector[i-1] + area; 
    } 
    
    //std::cout << "ComputeCDF: computed cdf" << std::endl;

    // find the max of CDF (should be 1 since we computed from normalized PDF)
    auto max_it = std::max_element(cdfvector.begin(), cdfvector.end());
    double max_cdf = *max_it;

    //std::cout << "ComputeCDF: computed max of cdf" << std::endl;

    // should be normalized, but just to make sure in case energy nodes are too sparse
    for (double& value : cdfvector) {
       value *= 1/max_cdf;
    }

    //std::cout << "ComputeCDF: normalized cdf" << std::endl;

    // assign the cdf vector so it's accessible outside of the function
    cdf = cdfvector;

    LI::utilities::TableData1D<double> inverse_cdf_data;
    inverse_cdf_data.x = cdf; 
    inverse_cdf_data.f = energies; 

    //std::cout << "ComputeCDF: put cdf and energy nodes in interpolater table" << std::endl;

    inverseCdfTable = LI::utilities::Interpolator1D<double>(inverse_cdf_data);

}

std::vector<double> TabulatedFluxDistribution::GetCDF() const {
    return cdf;
}

std::vector<double> TabulatedFluxDistribution::GetCDFEnergyNodes() const {
    return cdf_energy_nodes;
}

void TabulatedFluxDistribution::SetEnergyBounds(double eMin, double eMax) {
    energyMin = eMin;
    energyMax = eMax;
    bounds_set = true;
    ComputeIntegral();
    ComputeCDF();
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
    
    ComputeCDF(); 
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

    ComputeCDF();
}

double TabulatedFluxDistribution::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    // inverse CDF algorithm to sample from PDF.
    double randomValue = rand->Uniform(0,1);

    return inverseCdfTable(randomValue);
}

// double TabulatedFluxDistribution::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
//     // Metropolis-Hastings algorithm to sample from PDF.
//     // Pass in a function pointer for the PDF
//     double energy; 
//     double density, test_energy, test_density, odds;
//     bool accept;

//     // sample an initial point uniformly
//     //energy = rand->Uniform(energyMin, energyMax);
//     energy = rand->Uniform(std::log10(energyMin), std::log10(energyMax));
//     energy = std::pow(10, energy);
    
//     density = pdf(energy);

//     // Metropolis Hastings loop
//     for (size_t j = 0; j <= burnin; ++j) {
//         //test_energy = rand->Uniform(energyMin, energyMax);
//         test_energy = rand->Uniform(std::log10(energyMin), std::log10(energyMax));
//         test_energy = std::pow(10, test_energy);
//         test_density = pdf(test_energy);
//         odds = test_density / density;
//         accept = (odds > 1.) or (rand->Uniform(0,1) < odds);
//         if(accept) {
//             energy = test_energy;
//             density = test_density;
//         }
//     }
//     return energy;
// }

double TabulatedFluxDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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

} // namespace LI

