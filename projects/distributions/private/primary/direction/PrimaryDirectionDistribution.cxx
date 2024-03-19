#include "SIREN/distributions/primary/direction/PrimaryDirectionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt
#include <string>                                          // for basic_string

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/math/Vector3D.h"                  // for Vector3D

namespace siren {
namespace distributions {

//---------------
// class PrimaryDirectionDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryDirectionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    siren::math::Vector3D dir = SampleDirection(rand, detector_model, interactions, record);
    record.SetDirection(dir);
}

std::vector<std::string> PrimaryDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryDirection"};
}

} // namespace distributions
} // namespace siren

