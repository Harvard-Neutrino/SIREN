#include "SIREN/distributions/primary/energy_direction/Tabulated2DFluxDistribution.h"
#include <algorithm>                                       // for minmax_element
#include <array>                                           // for array
#include <cmath>                                           // for M_PI, cos, sin, sqrt, pow, log, log10, cbrt
#include <tuple>                                           // for tie, operator<
#include <vector>                                          // for vector
#include <fstream>                                         // for basic_istream
#include <sstream>                                         // for stringstream
#include <functional>                                      // for function
#include <stdexcept>                                       // for runtime_error

#include "SIREN/dataclasses/InteractionRecord.h"  // for InteractionRecord
#include "SIREN/distributions/Distributions.h"    // for InjectionDistribution
#include "SIREN/utilities/Integration.h"          // for simpsonIntegrate2D
#include "SIREN/utilities/Interpolator.h"         // for TableData2D
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace distributions {

//---------------
// class Tabulated2DFluxDistribution : PrimaryEnergyDirectionDistribution
//---------------
Tabulated2DFluxDistribution::Tabulated2DFluxDistribution() {}

void Tabulated2DFluxDistribution::ComputeIntegral() {
    // Integrate in log-energy space when the table is log-spaced in x.
    // Assumes the table is never log-spaced in y (cosZenith).
    if (fluxTable.IsLogY()) {
        throw std::runtime_error("Tabulated2DFluxDistribution does not support log-spaced cosZenith axis");
    }
    double eMin = energyMin;
    double eMax = energyMax;
    std::function<double(double, double)> integrand = [&] (double x, double y) -> double {
       return unnormed_pdf(x,y);
    };
    if (fluxTable.IsLogX()) {
        if(energyMin <= 0 || energyMax <= 0) {
            throw std::runtime_error("Energy bounds must be positive for log-spaced energy tables");
        }
        eMin = log10(energyMin);
        eMax = log10(energyMax);
        integrand =  [&] (double x, double y) -> double {
            return log(10)*pow(10,x)*unnormed_pdf(pow(10,x),y);
        };
    }

    integral = siren::utilities::simpsonIntegrate2D(integrand, eMin, eMax, cosZenithMin, cosZenithMax);
}

void Tabulated2DFluxDistribution::LoadFluxTable() {
    std::ifstream in(fluxTableFilename.c_str());
    if(!in.is_open()) {
        throw std::runtime_error("Failed to open flux table file: " + fluxTableFilename);
    }
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
        if(!(ss >> x >> y >> f)) {
            throw std::runtime_error("Malformed row in flux table: " + buf);
        }
        table_data.x.push_back(x);
        table_data.y.push_back(y);
        table_data.f.push_back(f);

        energy_nodes.push_back(x);
        cosZenith_nodes.push_back(y);
    }
    if(table_data.x.empty()) {
        throw std::runtime_error("No valid data rows in flux table: " + fluxTableFilename);
    }
    // If no bounds are manually set, derive from min/max of parsed values
    if(not energy_bounds_set) {
        auto mm = std::minmax_element(table_data.x.begin(), table_data.x.end());
        energyMin = *mm.first;
        energyMax = *mm.second;
    }
    if(not cosZenith_bounds_set) {
        auto mm = std::minmax_element(table_data.y.begin(), table_data.y.end());
        cosZenithMin = *mm.first;
        cosZenithMax = *mm.second;
    }
    fluxTable = siren::utilities::Interpolator2D<double>(table_data);
}

void Tabulated2DFluxDistribution::LoadFluxTable(std::vector<double> & energies, std::vector<double> & cosZeniths, std::vector<double> & flux) {

    if(energies.empty()) {
        throw std::runtime_error("Empty vectors passed to LoadFluxTable");
    }
    if(energies.size() != flux.size() || cosZeniths.size() != flux.size()) {
        throw std::runtime_error("LoadFluxTable: energies, cosZeniths, and flux vectors must have the same size");
    }

    siren::utilities::TableData2D<double> table_data;

    table_data.x = energies;
    table_data.y = cosZeniths;
    table_data.f = flux;
    energy_nodes = energies;
    cosZenith_nodes = cosZeniths;

    // If no bounds are manually set, derive from min/max of input values
    if(not energy_bounds_set) {
        auto mm = std::minmax_element(table_data.x.begin(), table_data.x.end());
        energyMin = *mm.first;
        energyMax = *mm.second;
    }
    if(not cosZenith_bounds_set) {
        auto mm = std::minmax_element(table_data.y.begin(), table_data.y.end());
        cosZenithMin = *mm.first;
        cosZenithMax = *mm.second;
    }
    fluxTable = siren::utilities::Interpolator2D<double>(table_data);
}

double Tabulated2DFluxDistribution::unnormed_pdf(double energy, double cosZenith) const {
    return fluxTable(energy, cosZenith);
}

double Tabulated2DFluxDistribution::SampleUnnormedPDF(double energy, double cosZenith) const {
    return unnormed_pdf(energy, cosZenith);
}

double Tabulated2DFluxDistribution::GetIntegral() const {
    return integral;
}

std::vector<double> Tabulated2DFluxDistribution::GetEnergyNodes() const {
    return energy_nodes;
}

std::vector<double> Tabulated2DFluxDistribution::GetCosZenithNodes() const {
    return cosZenith_nodes;
}

double Tabulated2DFluxDistribution::pdf(double energy, double cosZenith) const {
    return unnormed_pdf(energy,cosZenith) / integral;
}

double Tabulated2DFluxDistribution::SamplePDF(double energy, double cosZenith) const {
    return pdf(energy, cosZenith);
}

void Tabulated2DFluxDistribution::SetEnergyBounds(double eMin, double eMax) {
    energyMin = eMin;
    energyMax = eMax;
    energy_bounds_set = true;
    MH_sampled_points = 0;
    MH_density = 0;
    ComputeIntegral();
}

void Tabulated2DFluxDistribution::SetCosZenithBounds(double czMin, double czMax) {
    if(czMin < -1.0 || czMax > 1.0) {
        throw std::runtime_error("cosZenith bounds must be within [-1, 1]");
    }
    if(czMin > czMax) {
        throw std::runtime_error("cosZenithMin must be <= cosZenithMax");
    }
    cosZenithMin = czMin;
    cosZenithMax = czMax;
    cosZenith_bounds_set = true;
    MH_sampled_points = 0;
    MH_density = 0;
    ComputeIntegral();
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization)
    : energy_bounds_set(false)
    , cosZenith_bounds_set(false)
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
    , cosZenith_bounds_set(false)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, double cosZenithMin, double cosZenithMax, std::string fluxTableFilename, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , cosZenithMin(cosZenithMin)
    , cosZenithMax(cosZenithMax)
    , energy_bounds_set(true)
    , cosZenith_bounds_set(true)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(std::vector<double> energies, std::vector<double> cosZeniths, std::vector<double> flux, bool has_physical_normalization)
    : energy_bounds_set(false)
    , cosZenith_bounds_set(false)
{
    LoadFluxTable(energies,cosZeniths,flux);
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, std::vector<double> energies, std::vector<double> cosZeniths, std::vector<double> flux, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , energy_bounds_set(true)
    , cosZenith_bounds_set(false)
{
    LoadFluxTable(energies,cosZeniths,flux);
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

Tabulated2DFluxDistribution::Tabulated2DFluxDistribution(double energyMin, double energyMax, double cosZenithMin, double cosZenithMax, std::vector<double> energies, std::vector<double> cosZeniths, std::vector<double> flux, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , cosZenithMin(cosZenithMin)
    , cosZenithMax(cosZenithMax)
    , energy_bounds_set(true)
    , cosZenith_bounds_set(true)
{
    LoadFluxTable(energies,cosZeniths,flux);
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

std::pair<double, double> Tabulated2DFluxDistribution::SampleTablePoint(std::shared_ptr<siren::utilities::SIREN_random> rand) const {
    // Sample uniformly in linear or log space depending on the table
    double u1 = rand->Uniform();
    double u2 = rand->Uniform();
    double energy, cosZenith;
    if (fluxTable.IsLogX()) {
        double logEnergyMin = log10(energyMin);
        double logEnergyMax = log10(energyMax);
        energy = pow(10, logEnergyMin + u1 * (logEnergyMax - logEnergyMin));
    }
    else {
        energy = energyMin + u1 * (energyMax - energyMin);
    }
    // cosZenith axis is always linear (log-Y is rejected in ComputeIntegral)
    cosZenith = cosZenithMin + u2 * (cosZenithMax - cosZenithMin);
    return std::make_pair(energy, cosZenith);
}

std::pair<double,siren::math::Vector3D> Tabulated2DFluxDistribution::SampleEnergyAndDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {

    double energy, cosZenith, test_density;
    std::pair<double,double> energy_cosZenith;

    // Metropolis-Hastings burn-in
    while (MH_sampled_points <= burnin) {
        energy_cosZenith = SampleTablePoint(rand);
        energy = energy_cosZenith.first;
        cosZenith = energy_cosZenith.second;
        test_density = pdf(energy, cosZenith);
        bool accept = (MH_density <= 0) || (test_density >= MH_density) || (rand->Uniform() < test_density / MH_density);
        if(accept) {
            MH_density = test_density;
            MH_sampled_points++;
        }
    }
    // Sample one point after burn-in
    bool accept = false;
    while(!accept) {
        energy_cosZenith = SampleTablePoint(rand);
        energy = energy_cosZenith.first;
        cosZenith = energy_cosZenith.second;
        test_density = pdf(energy, cosZenith);
        accept = (MH_density <= 0) || (test_density >= MH_density) || (rand->Uniform() < test_density / MH_density);
        if(accept) {
            MH_density = test_density;
            MH_sampled_points++;
        }
    }
    // Convert cos(zenith) to a 3D direction vector
    double sin_theta = sqrt(1.0 - cosZenith * cosZenith);
    double phi = rand->Uniform(-M_PI, M_PI);
    double nx = sin_theta * cos(phi);
    double ny = sin_theta * sin(phi);
    siren::math::Vector3D dir(nx, ny, cosZenith);
    dir.normalize();
    return std::make_pair(energy, dir);
}


double Tabulated2DFluxDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    double p_mag = sqrt(pow(record.primary_momentum[1], 2) + pow(record.primary_momentum[2], 2) + pow(record.primary_momentum[3], 2));
    if(p_mag == 0)
        return 0.0;
    double cosZenith = record.primary_momentum[3] / p_mag;
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else if(cosZenith < cosZenithMin or cosZenith > cosZenithMax)
        return 0.0;
    else
        return pdf(energy, cosZenith);
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
            std::tie(energyMin, energyMax, cosZenithMin, cosZenithMax, fluxTable)
            ==
            std::tie(x->energyMin, x->energyMax, x->cosZenithMin, x->cosZenithMax, x->fluxTable);
}

bool Tabulated2DFluxDistribution::less(WeightableDistribution const & other) const {
    const Tabulated2DFluxDistribution* x = dynamic_cast<const Tabulated2DFluxDistribution*>(&other);
    if(!x)
        return false;
    return
        std::tie(energyMin, energyMax, cosZenithMin, cosZenithMax, fluxTable)
        <
        std::tie(x->energyMin, x->energyMax, x->cosZenithMin, x->cosZenithMax, x->fluxTable);
}

} // namespace distributions

} // namespace siren
