#include "SIREN/distributions/primary/energy/PowerLaw.h"

#include <array>                                           // for array
#include <algorithm>                                       // for min, max
#include <cmath>                                           // for isfinite, isnan
#include <tuple>                                           // for tie, opera...
#include <math.h>                                          // for pow, log10
#include <stdexcept>                                       // for invalid_ar...
#include <string>                                          // for basic_string

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace distributions {

//---------------
// class PowerLaw : PrimaryEnergyDistribution
//---------------
PowerLaw::PowerLaw(double powerLawIndex, double energyMin, double energyMax)
    : powerLawIndex(powerLawIndex)
    , energyMin(energyMin)
    , energyMax(energyMax)
{
    // Reject bounds that would sample with zero or undefined density.
    if(!std::isfinite(powerLawIndex) || !std::isfinite(energyMin) || std::isnan(energyMax)
            || energyMin < 0.0 || energyMin > energyMax)
        throw std::invalid_argument("PowerLaw: need a finite index and 0 <= energyMin <= energyMax with energyMin finite");
    bool divergent = (energyMin == 0.0 && powerLawIndex >= 1.0) || (std::isinf(energyMax) && powerLawIndex <= 1.0);
    if(energyMin < energyMax && divergent)
        throw std::invalid_argument("PowerLaw: E^-index is not normalizable on [energyMin, energyMax]");
}

double PowerLaw::SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    if(energyMin == energyMax)
        return energyMin; //return the only allowed energy

    double energy;
    if(powerLawIndex == 1.0) //sample uniformly in log space
        energy = pow(10.0, rand->Uniform(log10(energyMin), log10(energyMax)));
    else {
        double u = rand->Uniform();
        double energyP = (1 - u) * pow(energyMin, 1 - powerLawIndex) + u * pow(energyMax, 1 - powerLawIndex);
        energy = pow(energyP, 1 / (1 - powerLawIndex));
    }
    // Rounding can land endpoint samples a few ulp outside, where pdf is zero.
    return std::min(std::max(energy, energyMin), energyMax);
}

double PowerLaw::pdf(double energy) const {
    if(!std::isfinite(energy))
        return 0.0;
    if(energyMin == energyMax) // Monoenergetic's tolerance
        return std::abs(energy - energyMin) <= 1e-6 * energyMin ? 1.0 : 0.0;
    if(energy < energyMin || energy > energyMax)
        return 0.0;
    return UnboundedDensity(energy);
}

double PowerLaw::UnboundedDensity(double energy) const {
    if(energyMin == energyMax)
        return 1.0; // only one allowed energy

    if(powerLawIndex == 1.0)
        return 1.0 / (energy * log(energyMax / energyMin));
    else {
        return pow(energy, -powerLawIndex) * (-1.0 + powerLawIndex) / (pow(energyMin, 1.0 - powerLawIndex) - pow(energyMax, 1.0 - powerLawIndex));
    }
}

double PowerLaw::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    return pdf(record.primary_momentum[0]);
}

std::string PowerLaw::Name() const {
    return "PowerLaw";
}

std::shared_ptr<PrimaryInjectionDistribution> PowerLaw::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new PowerLaw(*this));
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
    SetNormalization(norm / UnboundedDensity(energy));
}

} // namespace distributions
} // namespace siren

