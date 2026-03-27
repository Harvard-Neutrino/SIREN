#include "SIREN/distributions/primary/vertex/SamplePrimaryPositionDistribution.h"

#include <array>                                                  // for array
#include <tuple>                                                  // for tie, tuple
#include <string>                                                 // for string
#include <vector>                                                 // for vector
#include <algorithm>                                              // for lower_bound
#include <stdexcept>                                              // for runtime_error

#include "SIREN/dataclasses/InteractionRecord.h"         // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"                // for DetectorModel
#include "SIREN/math/Vector3D.h"                         // for Vector3D
#include "SIREN/utilities/Random.h"                      // for SIREN_random

namespace siren {
namespace distributions {

SamplePrimaryPositionDistribution::SamplePrimaryPositionDistribution(std::vector<siren::math::Vector3D> positions, std::vector<double> weights)
    : positions_(positions), weights_(weights) {
    if (weights_.empty()) {
        std::cout << "Caution: No weights are used in association with the position data. This assumes you're data is unbiased and it will be sampled uniformly!" << std::endl;
        weights_.assign(positions_.size(), 1.0); 
    }
    if (positions_.size() != weights_.size()) {
        throw std::runtime_error("Sizes of positions and weights must match.");
    }
    double sum = 0.0;
    sum_weights_.resize(weights_.size());
    for (size_t i = 0; i < weights_.size(); ++i) {
        sum += weights_[i];
        sum_weights_[i] = sum;
    }
    total_sum_weights_ = sum;
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePrimaryPositionDistribution::SamplePosition(
    std::shared_ptr<siren::utilities::SIREN_random> rand,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::PrimaryDistributionRecord & record) const {
    if (positions_.empty()) {
        throw std::runtime_error("No positions available for sampling.");
    }
    if (total_sum_weights_ == 0.0) {
        throw std::runtime_error("Sum of weights is zero; cannot sample.");
    }
    double u = rand->Uniform(0, total_sum_weights_);
    auto it = std::lower_bound(sum_weights_.begin(), sum_weights_.end(), u);
    size_t idx = std::distance(sum_weights_.begin(), it);
    siren::math::Vector3D pos = positions_[idx];
    return {pos, pos}; // first pos is where primary particle is generated. Second pos is where interaction occurs. We are assuming the interactions are occuring where we are spawning our primary particle.
}

double SamplePrimaryPositionDistribution::GenerationProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record) const {
    return 1.0; // Since primary particle position input data is already modeled with Geant4 this can just be set to 1.0
}

std::string SamplePrimaryPositionDistribution::Name() const {
    return "SamplePrimaryPositionDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> SamplePrimaryPositionDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new SamplePrimaryPositionDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePrimaryPositionDistribution::InjectionBounds(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & interaction) const {
    siren::math::Vector3D pos(interaction.interaction_vertex[0], interaction.interaction_vertex[1], interaction.interaction_vertex[2]);
    return {pos, pos};
}

bool SamplePrimaryPositionDistribution::equal(WeightableDistribution const & distribution) const {
    const SamplePrimaryPositionDistribution* x = dynamic_cast<const SamplePrimaryPositionDistribution*>(&distribution);
    if (!x) return false;
    return (positions_ == x->positions_) && (weights_ == x->weights_);
}

bool SamplePrimaryPositionDistribution::less(WeightableDistribution const & distribution) const {
    const SamplePrimaryPositionDistribution* x = dynamic_cast<const SamplePrimaryPositionDistribution*>(&distribution);
    if (!x) return false;
    if (positions_ < x->positions_) return true;
    if (positions_ > x->positions_) return false;
    return weights_ < x->weights_;
}

} // namespace distributions
} // namespace siren