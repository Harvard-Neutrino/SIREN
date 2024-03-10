#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"

#include <array>                                           // for array
#include <string>                                          // for basic_string

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...

namespace siren {
namespace distributions {

//---------------
// class PrimaryEnergyDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryEnergyDistribution::Sample(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::PrimaryDistributionRecord & record) const {

    double energy = SampleEnergy(rand, detector_model, interactions, record);
    record.SetEnergy(energy);
}

std::vector<std::string> PrimaryEnergyDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy"};
}

} // namespace distributions
} // namespace sirenREN
