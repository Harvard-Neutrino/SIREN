#include "SIREN/distributions/primary/energy/SamplePrimaryEnergyDistribution.h"

#include <array>                                                  // for array
#include <tuple>                                                  // for tie, tuple
#include <string>                                                 // for string
#include <vector>                                                 // for vector
#include <algorithm>                                              // for lower_bound
#include <stdexcept>                                              // for runtime_error

#include "SIREN/dataclasses/InteractionRecord.h"         // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"                // for DetectorModel
#include "SIREN/utilities/Random.h"                      // for SIREN_random

namespace siren {
namespace distributions {

SamplePrimaryEnergyDistribution::SamplePrimaryEnergyDistribution(std::vector<double> energies, std::vector<double> weights)
    : energies_(energies), weights_(weights) {
    if (weights_.empty()) {
        std::cout << "Caution: No weights are used in association with the position data. This assumes you're data is unbiased and it will be sampled uniformly!" << std::endl;
        weights_.assign(energies_.size(), 1.0);
    }
    if (energies_.size() != weights_.size()) {
        throw std::runtime_error("Sizes of energies and weights must match.");
    }
    double sum = 0.0;
    sum_weights_.resize(weights_.size());
    for (size_t i = 0; i < weights_.size(); ++i) {
        sum += weights_[i];
        sum_weights_[i] = sum;
    }
    total_sum_weights_ = sum;
}

double SamplePrimaryEnergyDistribution::pdf(double energy) const {
    return 1.0;  // Constant density, since we're ignoring in generation probability
}

double SamplePrimaryEnergyDistribution::SampleEnergy(
    std::shared_ptr<siren::utilities::SIREN_random> rand,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::PrimaryDistributionRecord & record) const {
    if (energies_.empty()) {
        throw std::runtime_error("No energies available for sampling.");
    }
    if (total_sum_weights_ == 0.0) {
        throw std::runtime_error("Sum of weights is zero; cannot sample.");
    }
    double u = rand->Uniform(0, total_sum_weights_);
    auto it = std::lower_bound(sum_weights_.begin(), sum_weights_.end(), u);
    size_t idx = std::distance(sum_weights_.begin(), it);
    return energies_[idx];
}

double SamplePrimaryEnergyDistribution::GenerationProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record) const {
    return 1.0;
}

std::string SamplePrimaryEnergyDistribution::Name() const {
    return "SamplePrimaryEnergyDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> SamplePrimaryEnergyDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new SamplePrimaryEnergyDistribution(*this));
}

bool SamplePrimaryEnergyDistribution::equal(WeightableDistribution const & distribution) const {
    const SamplePrimaryEnergyDistribution* x = dynamic_cast<const SamplePrimaryEnergyDistribution*>(&distribution);
    if (!x) return false;
    return (energies_ == x->energies_) && (weights_ == x->weights_);
}

bool SamplePrimaryEnergyDistribution::less(WeightableDistribution const & distribution) const {
    const SamplePrimaryEnergyDistribution* x = dynamic_cast<const SamplePrimaryEnergyDistribution*>(&distribution);
    if (!x) return false;
    if (energies_ < x->energies_) return true;
    if (energies_ > x->energies_) return false;
    return weights_ < x->weights_;
}

} // namespace distributions
} // namespace siren