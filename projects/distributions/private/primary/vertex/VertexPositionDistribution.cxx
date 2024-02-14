#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"

#include <array>                                                  // for array
#include <string>                                                 // for bas...

#include "LeptonInjector/interactions/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/detector/DetectorModel.h"                   // for Ear...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...

namespace LI {
namespace distributions {

//---------------
// class VertexPositionDistribution : PrimaryInjectionDistribution
//---------------
void VertexPositionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::PrimaryDistributionRecord & record) const {
    std::tuple<LI::math::Vector3D, LI::math::Vector3D> init_and_pos = SamplePosition(rand, detector_model, interactions, record);
    LI::math::Vector3D const & init = std::get<0>(init_and_pos);
    LI::math::Vector3D const & pos = std::get<1>(init_and_pos);
    record.SetInitialPosition((std::array<double, 3>)init);
    record.SetInteractionVertex((std::array<double, 3>)pos);
}

std::vector<std::string> VertexPositionDistribution::DensityVariables() const {
    return std::vector<std::string>{"InteractionVertexPosition"};
}

std::tuple<LI::math::Vector3D, LI::math::Vector3D> VertexPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionTreeDatum const & datum) const {
  return InjectionBounds(detector_model, interactions, datum.record);
}

bool VertexPositionDistribution::AreEquivalent(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::DetectorModel const> second_detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution) and detector_model->operator==(*second_detector_model) and interactions->operator==(*second_interactions);
}

} // namespace distributions
} // namespace LeptonInjector
