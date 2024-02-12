#include "LeptonInjector/distributions/secondary/SecondaryInjectionDistribution.h"

#include <array>                                                  // for array
#include <string>                                                 // for bas...

#include "LeptonInjector/interactions/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/detector/DetectorModel.h"                   // for Ear...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...

namespace LI {
namespace distributions {

//---------------
// class SecondaryInjectionDistribution : InjectionDistribution
//---------------

void SecondaryInjectionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const {
    LI::dataclasses::SecondaryDistributionRecord secondary_record(record);
    Sample(rand, detector_model, interactions, secondary_record);
    secondary_record.Finalize(record);
}

std::vector<std::string> SecondaryInjectionDistribution::DensityVariables() const {
    return {};
}

bool SecondaryInjectionDistribution::AreEquivalent(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::DetectorModel const> second_detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution) and detector_model->operator==(*second_detector_model) and interactions->operator==(*second_interactions);
}

} // namespace distributions
} // namespace LeptonInjector
