#include "SIREN/distributions/primary/energy/Tabulated2DFluxDistribution.h"
#include <array>                                           // for array
#include <tuple>                                           // for tie, opera...
#include <vector>                                          // for vector
#include <fstream>                                         // for basic_istream
#include <functional>                                      // for function

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/utilities/Integration.h"          // for rombergInt...
#include "SIREN/utilities/Interpolator.h"         // for TableData1D
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace distributions {
namespace {
    bool fexists(const std::string filename)
    {
            std::ifstream ifile(filename.c_str());
            return (bool)ifile;
    }
}

//---------------
// class Tabulated2DFluxDistribution : PrimaryEnergyDirectionDistribution
//---------------
Tabulated2DFluxDistribution::Tabulated2DFluxDistribution() {}

void Tabulated2DFluxDistribution::ComputeIntegral() {
    std::function<double(double, double)> integrand = [&] (double x, double y) -> double {
        return unnormed_pdf(x,y);
    };
    integral = siren::utilities::rombergIntegrate2D(integrand, energyMin, energyMax, zenithMin, zenithMax);
}

void Tabulated2DFluxDistribution::LoadFluxTable() {
    if(fexists(fluxTableFilename)) {
        std::ifstream in(fluxTableFilename.c_str());
        std::string buf;
        std::string::size_type pos;
        siren::utilities::TableData2D<double> table_data;

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
            double x, y, f;
            ss >> x >> y >> f;
            table_data.x.push_back(x);
            table_data.y.push_back(y);
            table_data.f.push_back(f);

            energy_nodes.push_back(x);
            zenith_nodes.push_back(y);
        }
        // If no physical are manually set, use first/last entry of table
        if(not energy_bounds_set) {
            energyMin = table_data.x[0];
            energyMax = table_data.x[table_data.x.size()-1];
        }
        if(not zenith_bounds_set) {
            zenithMin = table_data.y[0];
            zenithMax = table_data.y[table_data.y.size()-1];
        }
        fluxTable = siren::utilities::Interpolator2D<double>(table_data);
    } else {
        throw std::runtime_error("Failed to open flux table file!");
    }
}

void Tabulated2DFluxDistribution::LoadFluxTable(std::vector<double> & energies, std::vector<double> & zeniths, std::vector<double> & flux) {

    assert(energies.size()==flux.size());
    assert(zeniths.size()==flux.size());

    siren::utilities::TableData2D<double> table_data;

    table_data.x = energies;
    table_data.y = zeniths;
    table_data.f = flux;
    energy_nodes = energies;
    zenith_nodes = energies;

    // If no physical are manually set, use first/last entry of table
    if(not energy_bounds_set) {
        energyMin = table_data.x[0];
        energyMax = table_data.x[table_data.x.size()-1];
    }
    if(not zenith_bounds_set) {
        zenithMin = table_data.y[0];
        zenithMax = table_data.y[table_data.y.size()-1];
    }
    fluxTable = siren::utilities::Interpolator2D<double>(table_data);
}

double Tabulated2DFluxDistribution::unnormed_pdf(double energy, double zenith) const {
    return fluxTable(energy, zenith);
}

double Tabulated2DFluxDistribution::SampleUnnormedPDF(double energy, double zenith) const {
    return unnormed_pdf(energy, zenith);
}

double Tabulated2DFluxDistribution::GetIntegral() const {
    return integral;
}

std::vector<double> Tabulated2DFluxDistribution::GetEnergyNodes() const {
    return energy_nodes;
}

double Tabulated2DFluxDistribution::pdf(double energy, double zenith) const {
    return unnormed_pdf(energy,zenith) / integral;
}

double Tabulated2DFluxDistribution::SamplePDF(double energy, double zenith) const {
    return pdf(energy, zenith);
}

void Tabulated2DFluxDistribution::SetEnergyBounds(double eMin, double eMax) {
    energyMin = eMin;
    energyMax = eMax;
    bounds_set = true;
    ComputeIntegral();
    ComputeCDF();
}

void Tabulated2DFluxDistribution::SetZenithBounds(double zMin, double zMax) {
    zenithMin = zMin;
    zenithMax = zMax;
    bounds_set = true;
    ComputeIntegral();
    ComputeCDF();
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization)
    : energy_bounds_set(false)
    , zenith_bounds_set(false)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, std::string fluxTableFilename, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , energy_bounds_set(true)
    , zenith_bounds_set(false)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, double zenithMin, double zenithMax, std::string fluxTableFilename, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , zenithMin(zenithMin)
    , zenithMax(zenithMax)
    , energy_bounds_set(true)
    , zenith_bounds_set(true)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(std::vector<double> energies, std::vector<double> zeniths, std::vector<double> flux, bool has_physical_normalization)
    : energy_bounds_set(false)
    , zenith_bounds_set(false)
{
    LoadFluxTable(energies,zeniths,flux);
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, std::vector<double> energies, std::vector<double> zeniths, std::vector<double> flux, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , energy_bounds_set(true)
    , zenith_bounds_set(false)
{
    LoadFluxTable(energies,zeniths,flux);
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, double zenithMin, double zenithMax, std::vector<double> energies, std::vector<double> zeniths, std::vector<double> flux, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    : zenithMin(zenithMin)
    , zenithMax(zenithMax)
    , energy_bounds_set(true)
    , zenith_bounds_set(true)
{
    LoadFluxTable(energies,zeniths,flux);
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

double Tabulated2DFluxDistribution::SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    // TODO: Use metropolis hastings
    return 0;
}


double Tabulated2DFluxDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    // TODO: check zenith
    else if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string Tabulated2DFluxDistribution::Name() const {
    return "Tabulated2DFluxDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> Tabulated2DFluxDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new Tabulated2DFluxDistribution(*this));
}

bool Tabulated2DFluxDistribution::equal(WeightableDistribution const & other) const {
    const Tabulated2DFluxDistribution* x = dynamic_cast<const Tabulated2DFluxDistribution*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, zenithMin, zenithMax, fluxTable)
            ==
            std::tie(x->energyMin, x->energyMax, x->zenithMin, x->zenithMax, x->fluxTable);
}

bool Tabulated2DFluxDistribution::less(WeightableDistribution const & other) const {
    const Tabulated2DFluxDistribution* x = dynamic_cast<const Tabulated2DFluxDistribution*>(&other);
    return
        std::tie(energyMin, energyMax, zenithMin, zenithMax, fluxTable)
        <
        std::tie(x->energyMin, x->energyMax, x->zenithMin, x->zenithMax, x->fluxTable);
}

} // namespace distributions

} // namespace siren

