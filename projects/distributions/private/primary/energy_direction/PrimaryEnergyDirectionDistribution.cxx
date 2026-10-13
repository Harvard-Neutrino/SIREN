#include "SIREN/distributions/primary/energy_direction/PrimaryEnergyDirectionDistribution.h"

#include <set>                                             // for set
#include <array>                                           // for array
#include <string>                                          // for basic_string

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/DistributionVariable.h"  // for DistributionVariable

namespace siren {
namespace distributions {

//---------------
// class PrimaryEnergyDirectionDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryEnergyDirectionDistribution::Sample(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::PrimaryDistributionRecord & record) const {

    std::pair<double,siren::math::Vector3D> energy_and_direction = SampleEnergyAndDirection(rand, detector_model, interactions, record);
    record.SetEnergy(energy_and_direction.first);
    record.SetDirection(energy_and_direction.second);
}

std::set<DistributionVariable> PrimaryEnergyDirectionDistribution::SetVariables() const {
    return {DistributionVariable::PrimaryEnergy, DistributionVariable::PrimaryDirection};
}

std::vector<std::string> PrimaryEnergyDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy","PrimaryDirection"};
}

} // namespace distributions
} // namespace siren
