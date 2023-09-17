#include "LeptonInjector/distributions/primary/energy/PowerLaw.h"

#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"

namespace LI {
namespace distributions {

//---------------
// class PowerLaw : PrimaryEnergyDistribution
//---------------
PowerLaw::PowerLaw(double powerLawIndex, double energyMin, double energyMax)
    : powerLawIndex(powerLawIndex)
    , energyMin(energyMin)
    , energyMax(energyMax)
{}

double PowerLaw::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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
        return pow(energy, -powerLawIndex) * (-1.0 + powerLawIndex) / (pow(energyMin, 1.0 - powerLawIndex) - pow(energyMax, 1.0 - powerLawIndex));
    }
}

double PowerLaw::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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

} // namespace distributions
} // namespace LI

