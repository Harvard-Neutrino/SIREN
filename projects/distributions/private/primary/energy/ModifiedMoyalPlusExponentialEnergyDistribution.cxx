#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/utilities/Integration.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "LeptonInjector/distributions/primary/energy/ModifiedMoyalPlusExponentialEnergyDistribution.h"

namespace LI {
namespace distributions {

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

double ModifiedMoyalPlusExponentialEnergyDistribution::pdf_integral() const {
    double exponential = B * (exp(-energyMin/l) - exp(-energyMax/l));
    double moyal = A * (std::erf(exp((mu - energyMin)/(2.0 * sigma)) / sqrt(2.0)) - std::erf(exp((mu - energyMax)/(2.0 * sigma)) / sqrt(2.0)));
    return exponential + moyal;
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
    integral = pdf_integral();
    std::function<double(double)> integrand = [&] (double x) -> double {
        return pdf(x);
    };
    double test_norm = LI::utilities::rombergIntegrate(integrand, energyMin, energyMax, 1e-8);
    if(std::abs(1.0 - test_norm) < 1e-6) {
        integral = 1.0;
        integral = LI::utilities::rombergIntegrate(integrand, energyMin, energyMax, 1e-8);
    }
    if(has_physical_normalization) {
        SetNormalization(integral);
    }
}

double ModifiedMoyalPlusExponentialEnergyDistribution::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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

double ModifiedMoyalPlusExponentialEnergyDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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

} // namespace distributions
} // namespace LI

