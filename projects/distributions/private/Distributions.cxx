#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace distributions {

//---------------
// PhysicallyNormalizedDistribution
//---------------
PhysicallyNormalizedDistribution::PhysicallyNormalizedDistribution() {
    SetNormalization(1.0);
}

PhysicallyNormalizedDistribution::PhysicallyNormalizedDistribution(double norm) {
    SetNormalization(norm);
}

void PhysicallyNormalizedDistribution::SetNormalization(double norm) {
    normalization = norm;
}

double PhysicallyNormalizedDistribution::GetNormalization() const {
    return normalization;
}

bool PhysicallyNormalizedDistribution::IsNormalizationSet() const {
    return normalization != 1.0;
}


//---------------
// class WeightableDistribution
//---------------
std::vector<std::string> WeightableDistribution::DensityVariables() const {
    return {};
}

bool WeightableDistribution::operator==(WeightableDistribution const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool WeightableDistribution::operator<(WeightableDistribution const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

bool WeightableDistribution::AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> second_cross_sections) const {
    return this->operator==(*distribution);
}

double WeightableDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum const & datum) const {
    return GenerationProbability(earth_model, cross_sections, datum.record);
}

//---------------
// class NormalizationConstant : WeightableDistribution, PhysicallyNormalizedDistribution
//---------------
//
NormalizationConstant::NormalizationConstant() {}

NormalizationConstant::NormalizationConstant(double norm) {
    SetNormalization(norm);
}

double NormalizationConstant::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    return 1.0;
}

std::string NormalizationConstant::Name() const {
    return "NormalizationConstant";
}

bool NormalizationConstant::equal(WeightableDistribution const & distribution) const {
    return false;
}

bool NormalizationConstant::less(WeightableDistribution const & distribution) const {
    return false;
}

//---------------
// class InjectionDistribution
//---------------

void InjectionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
}

void InjectionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum & datum) const {
  Sample(rand, earth_model, cross_sections, datum.record);
}

} // namespace distributions
} // namespace LeptonInjector
