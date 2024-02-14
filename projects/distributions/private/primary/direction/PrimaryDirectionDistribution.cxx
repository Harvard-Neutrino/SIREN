#include "LeptonInjector/distributions/primary/direction/PrimaryDirectionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt
#include <string>                                          // for basic_string

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/math/Vector3D.h"                  // for Vector3D

namespace LI {
namespace distributions {

//---------------
// class PrimaryDirectionDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryDirectionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::PrimaryDistributionRecord & record) const {
    LI::math::Vector3D dir = SampleDirection(rand, detector_model, interactions, record);
    record.SetDirection(dir);
}

std::vector<std::string> PrimaryDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryDirection"};
}

} // namespace distributions
} // namespace LI

