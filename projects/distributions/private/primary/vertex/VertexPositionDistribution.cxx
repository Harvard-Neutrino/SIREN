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
// class VertexPositionDistribution : InjectionDistribution
//---------------
void VertexPositionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D pos = SamplePosition(rand, detector_model, interactions, record);
    record.interaction_vertex[0] = pos.GetX();
    record.interaction_vertex[1] = pos.GetY();
    record.interaction_vertex[2] = pos.GetZ();
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
