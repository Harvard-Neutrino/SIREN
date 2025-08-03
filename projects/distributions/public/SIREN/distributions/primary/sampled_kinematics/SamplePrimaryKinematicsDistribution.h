#pragma once
#ifndef SIREN_SamplePrimaryKinematicsDistribution_H
#define SIREN_SamplePrimaryKinematicsDistribution_H

#include <tuple>                                         // for tuple
#include <memory>                                        // for shared_ptr
#include <string>                                        // for string
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#include "SIREN/dataclasses/InteractionRecord.h"  // for InteractionRecord
#include "SIREN/distributions/Distributions.h"   // for WeightableDistribution
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h" // for VertexPositionDistribution
#include "SIREN/math/Vector3D.h"                 // for Vector3D

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class SamplePrimaryKinematicsDistribution : public VertexPositionDistribution {
friend cereal::access;
public:
    SamplePrimaryKinematicsDistribution() = default;
    SamplePrimaryKinematicsDistribution(std::vector<double> energies, std::vector<siren::math::Vector3D> positions, std::vector<siren::math::Vector3D> directions, double primary_mass, std::vector<double> weights = {});
    virtual ~SamplePrimaryKinematicsDistribution() {};
private:
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("energies", energies_));
            archive(::cereal::make_nvp("positions", positions_));
            archive(::cereal::make_nvp("directions", directions_));
            archive(::cereal::make_nvp("primary_mass", primary_mass_));
            archive(::cereal::make_nvp("weights", weights_));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("SamplePrimaryKinematicsDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("energies", energies_));
            archive(::cereal::make_nvp("positions", positions_));
            archive(::cereal::make_nvp("directions", directions_));
            archive(::cereal::make_nvp("primary_mass", primary_mass_));
            archive(::cereal::make_nvp("weights", weights_));
            if (weights_.empty()) {
                weights_.assign(energies_.size(), 1.0);
            }
            if (energies_.size() != positions_.size() || energies_.size() != directions_.size() || energies_.size() != weights_.size()) {
                throw std::runtime_error("Sizes of energies, positions, directions, and weights must match.");
            }
            double sum = 0.0;
            sum_weights_.resize(weights_.size());
            for(size_t i = 0; i < weights_.size(); ++i) {
                sum += weights_[i];
                sum_weights_[i] = sum;
            }
            total_sum_weights_ = sum;
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("SamplePrimaryKinematicsDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
private:
    std::vector<double> energies_;
    std::vector<siren::math::Vector3D> positions_;
    std::vector<siren::math::Vector3D> directions_;
    double primary_mass_;
    std::vector<double> weights_;
    mutable std::vector<double> sum_weights_; // mutable to allow modification in const methods
    mutable double total_sum_weights_; // mutable to allow modification in const methods
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::SamplePrimaryKinematicsDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::SamplePrimaryKinematicsDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::SamplePrimaryKinematicsDistribution);

#endif // SIREN_SamplePrimaryKinematicsDistribution_H