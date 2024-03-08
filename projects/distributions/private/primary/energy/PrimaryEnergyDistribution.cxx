#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"

#include <array>                                           // for array
#include <string>                                          // for basic_string

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...

namespace SI {
namespace distributions {

//---------------
// class PrimaryEnergyDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryEnergyDistribution::Sample(
        std::shared_ptr<SI::utilities::LI_random> rand,
        std::shared_ptr<SI::detector::DetectorModel const> detector_model,
        std::shared_ptr<SI::interactions::InteractionCollection const> interactions,
        SI::dataclasses::PrimaryDistributionRecord & record) const {

    double energy = SampleEnergy(rand, detector_model, interactions, record);
    record.SetEnergy(energy);
}

std::vector<std::string> PrimaryEnergyDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy"};
}

} // namespace distributions
} // namespace SIREN
