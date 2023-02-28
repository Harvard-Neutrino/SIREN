#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"

#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"

namespace LI {
namespace distributions {

//---------------
// class VertexPositionDistribution : InjectionDistribution
//---------------
void VertexPositionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, record);
    record.interaction_vertex[0] = pos.GetX();
    record.interaction_vertex[1] = pos.GetY();
    record.interaction_vertex[2] = pos.GetZ();
}

std::vector<std::string> VertexPositionDistribution::DensityVariables() const {
    return std::vector<std::string>{"InteractionVertexPosition"};
}

bool VertexPositionDistribution::AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> second_cross_sections) const {
    return this->operator==(*distribution) and earth_model->operator==(*second_earth_model) and cross_sections->operator==(*second_cross_sections);
}

} // namespace distributions
} // namespace LeptonInjector
