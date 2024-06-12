#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"

#include <array>                                                  // for array
#include <string>                                                 // for bas...

#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...
#include "SIREN/math/Vector3D.h"                         // for Vec...

namespace siren {
namespace distributions {

//---------------
// class SecondaryVertexPositionDistribution : InjectionDistribution
//---------------
void SecondaryVertexPositionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const {
    // std::cout << "sampling vertex" << std::endl;
    SampleVertex(rand, detector_model, interactions, record);
}

std::vector<std::string> SecondaryVertexPositionDistribution::DensityVariables() const {
    return std::vector<std::string>{"Length"};
}

bool SecondaryVertexPositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution) and detector_model->operator==(*second_detector_model) and interactions->operator==(*second_interactions);
}

} // namespace distributions
} // namespace sirenREN
