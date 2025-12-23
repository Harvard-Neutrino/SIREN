#include "SIREN/distributions/primary/sampled_kinematics/SamplePrimaryKinematicsDistribution.h"

#include <array>                                                  // for array
#include <tuple>                                                  // for tuple
#include <string>                                                 // for string
#include <vector>                                                 // for vector
#include <algorithm>                                              // for lower_bound, find_if
#include <stdexcept>                                              // for runtime_error
#include <cmath>                                                  // for sqrt

#include "SIREN/dataclasses/InteractionRecord.h"         // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"                // for DetectorModel
#include "SIREN/utilities/Random.h"                      // for SIREN_random
#include "SIREN/math/Vector3D.h"                         // for Vector3D

namespace siren {
namespace distributions {

SamplePrimaryKinematicsDistribution::SamplePrimaryKinematicsDistribution(std::vector<double> energies, std::vector<siren::math::Vector3D> positions, std::vector<siren::math::Vector3D> directions, double primary_mass, std::vector<double> weights)
    : energies_(energies), positions_(positions), directions_(directions), primary_mass_(primary_mass), weights_(weights) {
    if (weights_.empty()) {
        std::cout << "Caution: No weights are used in association with the data. This assumes your data is unbiased and it will be sampled uniformly!" << std::endl;
        weights_.assign(energies_.size(), 1.0);
    }
    if (energies_.size() != positions_.size() || energies_.size() != directions_.size() || energies_.size() != weights_.size()) {
        throw std::runtime_error("Sizes of energies, positions, directions, and weights must match.");
    }
    double sum = 0.0;
    sum_weights_.resize(weights_.size());
    for (size_t i = 0; i < weights_.size(); ++i) {
        sum += weights_[i];
        sum_weights_[i] = sum;
    }
    total_sum_weights_ = sum;
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePrimaryKinematicsDistribution::SamplePosition(
    std::shared_ptr<siren::utilities::SIREN_random> rand,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::PrimaryDistributionRecord & record) const {
    throw std::runtime_error("SamplePrimaryKinematicsDistribution::SamplePosition is not implemented!");
    return {siren::math::Vector3D(), siren::math::Vector3D()};
}

void SamplePrimaryKinematicsDistribution::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> rand,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::PrimaryDistributionRecord & record) const {
    if (energies_.empty()) {
        throw std::runtime_error("No data available for sampling.");
    }
    if (total_sum_weights_ == 0.0) {
        throw std::runtime_error("Sum of weights is zero; cannot sample.");
    }
    double u = rand->Uniform(0, total_sum_weights_);
    auto it = std::lower_bound(sum_weights_.begin(), sum_weights_.end(), u);
    size_t idx = std::distance(sum_weights_.begin(), it);
    double step_length = weights_[idx];
    double energy = energies_[idx];
    siren::math::Vector3D pos = positions_[idx];
    siren::math::Vector3D dir = directions_[idx];
    // Do not normalize dir as per user instruction; assume it's already a unit vector
    std::array<double, 3> position_array = {pos.GetX(), pos.GetY(), pos.GetZ()};
    std::array<double, 3> direction_array = {dir.GetX(), dir.GetY(), dir.GetZ()};
    record.SetEnergy(energy);
    record.SetInitialPosition(position_array);
    record.SetInteractionVertex(position_array); // Interaction at spawn point
    record.SetDirection(direction_array);
    record.SetInjectionStepLength(step_length);
}

double SamplePrimaryKinematicsDistribution::GenerationProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record) const {
    return 1.0 / (total_sum_weights_ / simulated_pot_); // Constant probability, as data is pre-modeled
}

std::vector<std::string> SamplePrimaryKinematicsDistribution::DensityVariables() const {
    return {"PrimaryEnergy", "InteractionVertexPosition", "PrimaryDirection"};
}

std::string SamplePrimaryKinematicsDistribution::Name() const {
    return "SamplePrimaryKinematicsDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> SamplePrimaryKinematicsDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new SamplePrimaryKinematicsDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePrimaryKinematicsDistribution::InjectionBounds(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & interaction) const {
    siren::math::Vector3D vertex(interaction.interaction_vertex[0], interaction.interaction_vertex[1], interaction.interaction_vertex[2]);
    siren::math::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D end = vertex + interaction.primary_injection_step_length*dir;
    return {vertex, end};
}

bool SamplePrimaryKinematicsDistribution::equal(WeightableDistribution const & distribution) const {
    const SamplePrimaryKinematicsDistribution* x = dynamic_cast<const SamplePrimaryKinematicsDistribution*>(&distribution);
    if (!x) return false;
    return (energies_ == x->energies_) &&
           (positions_ == x->positions_) &&
           (directions_ == x->directions_) &&
           (primary_mass_ == x->primary_mass_) &&
           (weights_ == x->weights_);
}

bool SamplePrimaryKinematicsDistribution::less(WeightableDistribution const & distribution) const {
    const SamplePrimaryKinematicsDistribution* x = dynamic_cast<const SamplePrimaryKinematicsDistribution*>(&distribution);
    if (!x) return false;
    if (energies_ < x->energies_) return true;
    if (energies_ > x->energies_) return false;
    if (positions_ < x->positions_) return true;
    if (positions_ > x->positions_) return false;
    if (directions_ < x->directions_) return true;
    if (directions_ > x->directions_) return false;
    if (primary_mass_ < x->primary_mass_) return true;
    if (primary_mass_ > x->primary_mass_) return false;
    return weights_ < x->weights_;
}

} // namespace distributions
} // namespace siren
