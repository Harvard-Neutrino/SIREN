#include "SIREN/distributions/Distributions.h"

#include <string>     // for basic_string
#include <vector>     // for vector
#include <typeinfo>   // for type_info
#include <typeindex>  // for type_index

namespace siren {
namespace distributions {

//---------------
// PhysicallyNormalizedDistribution
//---------------
PhysicallyNormalizedDistribution::PhysicallyNormalizedDistribution()
    : normalization_set(false), normalization(1.0) {}

PhysicallyNormalizedDistribution::PhysicallyNormalizedDistribution(double norm) {
    SetNormalization(norm);
}

void PhysicallyNormalizedDistribution::SetNormalization(double norm) {
    normalization = norm;
    normalization_set = true;
}

double PhysicallyNormalizedDistribution::GetNormalization() const {
    return normalization;
}

bool PhysicallyNormalizedDistribution::IsNormalizationSet() const {
    return normalization_set;
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

bool WeightableDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution);
}

//---------------
// class NormalizationConstant : WeightableDistribution, PhysicallyNormalizedDistribution
//---------------
//
NormalizationConstant::NormalizationConstant() {}

NormalizationConstant::NormalizationConstant(double norm) {
    SetNormalization(norm);
}

double NormalizationConstant::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    return 1.0;
}

std::string NormalizationConstant::Name() const {
    return "NormalizationConstant";
}

bool NormalizationConstant::equal(WeightableDistribution const & distribution) const {
    const PhysicallyNormalizedDistribution* dist = dynamic_cast<const PhysicallyNormalizedDistribution*>(&distribution);
    if(!dist) return false;
    return normalization==dist->GetNormalization();
}

bool NormalizationConstant::less(WeightableDistribution const & distribution) const {
    const PhysicallyNormalizedDistribution* dist = dynamic_cast<const PhysicallyNormalizedDistribution*>(&distribution);
    if(!dist) return false;
    return normalization<dist->GetNormalization();
}

} // namespace distributions
} // namespace sirenREN
