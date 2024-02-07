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

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace detector { class DetectorModel; } }

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

void TabulatedFluxDistribution::LoadFluxTable(std::vector<double> & energies, std::vector<double> & flux) {
    
    assert(energies.size()==flux.size());

    LI::utilities::TableData1D<double> table_data;

    table_data.x = energies;
    table_data.f = flux;
    energy_nodes = energies;

    // If no physical are manually set, use first/last entry of table
    if(not bounds_set) {
        energyMin = table_data.x[0];
        energyMax = table_data.x[table_data.x.size()-1];
    }
    fluxTable = LI::utilities::Interpolator1D<double>(table_data);
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
    //cdf_energy_nodes = energies;

    // declare the cdf vectors
    std::vector<double> cdf_vector;
    std::vector<double> cdf_energy_nodes;

    cdf_vector.push_back(0); // start the CDF at zero
    cdf_energy_nodes.push_back(energies[0]);
    for(int i = 1; i < energies.size(); ++i) {
        double p1 = pdf(energies[i-1]);
        double p2 = pdf(energies[i]);

        // skip unsupported parts of pdf
        if((p1+p2)<=0) continue;

        // Check if we have skipped previous energy node
        if(cdf_energy_nodes.back() != energies[i-1]){
            // if so, we need to add small epsilon to CDF to
            // ignore unsupported part of PDF during interpolation
            cdf_energy_nodes.push_back(energies[i-1]);
            cdf_vector.push_back(cdf_vector.back()+1e-12);
        }

        // trapezoidal area
        double area = 0.5 * (p1 + p2) * (energies[i]-energies[i-1]);
        cdf_vector.push_back(cdf_vector.back() + area);
        cdf_energy_nodes.push_back(energies[i]);
    } 
    

    // find the max of CDF (should be 1 since we computed from normalized PDF)
    auto max_it = std::max_element(cdf_vector.begin(), cdf_vector.end());
    double max_cdf = *max_it;


    // should be normalized, but just to make sure in case energy nodes are too sparse
    for (double& value : cdf_vector) {
       value *= 1/max_cdf;
    }


    // assign the cdf vector so it's accessible outside of the function
    cdf = cdf_vector;

    LI::utilities::TableData1D<double> inverse_cdf_data;
    inverse_cdf_data.x = cdf; 
    inverse_cdf_data.f = cdf_energy_nodes;


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

TabulatedFluxDistribution::TabulatedFluxDistribution(std::vector<double> energies, std::vector<double> flux, bool has_physical_normalization)
    : bounds_set(false)
{
    LoadFluxTable(energies,flux);
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);

    ComputeCDF();
}

TabulatedFluxDistribution::TabulatedFluxDistribution(double energyMin, double energyMax, std::vector<double> energies, std::vector<double> flux, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , bounds_set(true)
{
    LoadFluxTable(energies,flux);
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);

    ComputeCDF();
}

TabulatedFluxDistribution::TabulatedFluxDistribution(std::vector<double> energies, std::vector<double> flux, bool has_physical_normalization)
    : bounds_set(false)
{
    LoadFluxTable(energies,flux);
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);

    ComputeCDF();
}

TabulatedFluxDistribution::TabulatedFluxDistribution(double energyMin, double energyMax, std::vector<double> energies, std::vector<double> flux, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , bounds_set(true)
{
    LoadFluxTable(energies,flux);
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);

    ComputeCDF();
}

double TabulatedFluxDistribution::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::PrimaryDistributionRecord & record) const {
    // inverse CDF algorithm to sample from PDF.
    double randomValue = rand->Uniform(0,1);

    return inverseCdfTable(randomValue);
}


double TabulatedFluxDistribution::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string TabulatedFluxDistribution::Name() const {
    return "TabulatedFluxDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> TabulatedFluxDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new TabulatedFluxDistribution(*this));
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

