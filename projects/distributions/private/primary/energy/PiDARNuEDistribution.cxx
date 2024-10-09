#include "SIREN/distributions/primary/energy/PiDARNuEDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for exp, sqrt
#include <tuple>                                           // for tie, opera...
#include <string>                                          // for basic_string
#include <stdlib.h>                                        // for abs, size_t
#include <functional>                                      // for function
#include <cmath>

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/utilities/Integration.h"          // for rombergInt...
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace distributions {

//---------------
// class PiDARNuEDistribution : PrimaryEnergyDistribution
//---------------

//constant for muon mass
double const & siren::distributions::PiDARNuEDistribution::mu_mass = siren::utilities::Constants::muonMass;

//pdf for electron neutrinos from pion decay at rest
//this function is already normalized
double PiDARNuEDistribution::pdf(double energy) const {
    double phi_0 = 1.0;
    double term1 = (192 * std::pow(energy, 2)) / std::pow(mu_mass, 3);
    double term2 = (0.5 - (energy / mu_mass));
    double pdf = phi_0 * term1 * term2;

    return pdf;

}

//cdf function (integral of the pdf from 0 to energy)
double PiDARNuEDistribution::cdf(double energy) const {
    double phi_0 = 1.0;
    double A = (96.0 / std::pow(mu_mass, 3));
    double B = (192.0 / std::pow(mu_mass, 4));

    // CDF analytical expression
    return phi_0 * (A * (std::pow(energy, 3) /3.0) - B * (std::pow(energy, 4) /4.0));
}


//analytical inverse cdf function
double PiDARNuEDistribution::inv_cdf(double x) const {
    const double m = mu_mass;
    return (2*m - sqrt(4*(m * m) + (6*pow(m, 4)*x)/
        pow(pow(m, 6)*x - sqrt(-(pow(m, 12)*(-1 + x)*(x * x))),1./3.)
         + 6*pow(pow(m, 6)*x - sqrt(-(pow(m, 12)*(-1 + x)*(x * x))),
         1./3.)) - sqrt(2)*
      sqrt(4*(m * m) - (3*pow(m, 4)*x)/
         pow(pow(m, 6)*x - sqrt(-(pow(m, 12)*(-1 + x)*(x * x))),
          1./3.) - 3*pow(pow(m, 6)*x -
           sqrt(-(pow(m, 12)*(-1 + x)*(x * x))),1./3.) -
        (4*(m * m * m))/
         sqrt((m * m) + (3*pow(m, 4)*x)/
            (2.*pow(pow(m, 6)*x - sqrt(-(pow(m, 12)*(-1 + x)*(x * x))),
               1./3.)) +
           (3*pow(pow(m, 6)*x - sqrt(-(pow(m, 12)*(-1 + x)*(x * x))),
               1./3.))/2.)))/12.;
}

//Use inverse cdf
double PiDARNuEDistribution::SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {

    double u = rand->Uniform();
    double energy = inv_cdf(u);

    return energy;
}

double PiDARNuEDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    double const energyMin = 0.0;
    double const energyMax = mu_mass/2.0;

    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string PiDARNuEDistribution::Name() const {
    return "PiDARNuEDistribution";
}
//clone, equal , less???
std::shared_ptr<PrimaryInjectionDistribution> PiDARNuEDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new PiDARNuEDistribution(*this));
}

bool PiDARNuEDistribution::equal(WeightableDistribution const & other) const {
    const PiDARNuEDistribution* x = dynamic_cast<const PiDARNuEDistribution*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool PiDARNuEDistribution::less(WeightableDistribution const & other) const {
    const PiDARNuEDistribution* x = dynamic_cast<const PiDARNuEDistribution*>(&other);
    return false;
}

} // namespace distributions
} // namespace siren

