#include "LeptonInjector/distributions/target/momentum/TargetMomentumDistribution.h"

#include <string>                                          // for basic_string

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/distributions/Distributions.h"    // for InjectionD...

namespace LI {
namespace distributions {

//---------------
// class TargetMomentumDistribution : InjectionDistribution
//---------------

void TargetMomentumDistribution::Sample(
        std::shared_ptr<LI::utilities::LI_random> rand,
        std::shared_ptr<LI::detector::EarthModel const> earth_model,
        std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections,
        LI::dataclasses::InteractionRecord & record) const {
    record.target_momentum = SampleMomentum(rand, earth_model, cross_sections, record);
}

std::vector<std::string> TargetMomentumDistribution::DensityVariables() const {
    return std::vector<std::string>{"TargetMomentum"};
}

//---------------
// class TargetAtRest : TargetMomentumDistribution : InjectionDistribution
//---------------
std::array<double, 4> TargetAtRest::SampleMomentum(
        std::shared_ptr<LI::utilities::LI_random> rand,
        std::shared_ptr<LI::detector::EarthModel const> earth_model,
        std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections,
        LI::dataclasses::InteractionRecord const & record) const {
    return std::array<double, 4>{record.target_mass, 0, 0, 0};
}

double TargetAtRest::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    return 1.0;
}

std::vector<std::string> TargetAtRest::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<InjectionDistribution> TargetAtRest::clone() const {
    return std::shared_ptr<InjectionDistribution>(new TargetAtRest(*this));
}

std::string TargetAtRest::Name() const {
    return "TargetAtRest";
}

bool TargetAtRest::equal(WeightableDistribution const & other) const {
    const TargetAtRest* x = dynamic_cast<const TargetAtRest*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool TargetAtRest::less(WeightableDistribution const & other) const {
    return false;
}

} // namespace distributions
} // namespace LeptonInjector
