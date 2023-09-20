#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"

#include <array>                                           // for array
#include <string>                                          // for basic_string

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...

namespace LI {
namespace distributions {

//---------------
// class PrimaryEnergyDistribution : InjectionDistribution
//---------------
void PrimaryEnergyDistribution::Sample(
        std::shared_ptr<LI::utilities::LI_random> rand,
        std::shared_ptr<LI::detector::EarthModel const> earth_model,
        std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections,
        LI::dataclasses::InteractionRecord & record) const {
    record.primary_momentum[0] = SampleEnergy(rand, earth_model, cross_sections, record);
}

std::vector<std::string> PrimaryEnergyDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy"};
}

} // namespace distributions
} // namespace LeptonInjector
