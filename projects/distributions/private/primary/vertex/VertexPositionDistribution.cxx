#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"

#include <array>                                                  // for array
#include <string>                                                 // for bas...

#include "LeptonInjector/crosssections/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/detector/EarthModel.h"                   // for Ear...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...

namespace LI {
namespace distributions {

//---------------
// class VertexPositionDistribution : InjectionDistribution
//---------------
LI::math::Vector3D VertexPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum & datum) const {
  return SamplePosition(rand, earth_model, cross_sections, datum.record);
}

void VertexPositionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, record);
    record.interaction_vertex[0] = pos.GetX();
    record.interaction_vertex[1] = pos.GetY();
    record.interaction_vertex[2] = pos.GetZ();
}

std::vector<std::string> VertexPositionDistribution::DensityVariables() const {
    return std::vector<std::string>{"InteractionVertexPosition"};
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> VertexPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum const & datum) const {
  return InjectionBounds(earth_model, cross_sections, datum.record);
}

bool VertexPositionDistribution::AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> second_cross_sections) const {
    return this->operator==(*distribution) and earth_model->operator==(*second_earth_model) and cross_sections->operator==(*second_cross_sections);
}

} // namespace distributions
} // namespace LeptonInjector
