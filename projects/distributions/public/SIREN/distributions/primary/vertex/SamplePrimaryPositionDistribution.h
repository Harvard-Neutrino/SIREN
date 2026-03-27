#pragma once
#ifndef SIREN_SamplePrimaryPositionDistribution_H
#define SIREN_SamplePrimaryPositionDistribution_H

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

#include "SIREN/dataclasses/InteractionRecord.h"  // for InteractionTree...
#include "SIREN/distributions/Distributions.h"  // for WeightableDistri...
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/math/Vector3D.h"                // for Vector3D

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

class SamplePrimaryPositionDistribution : public VertexPositionDistribution {
friend cereal::access;
public:
    SamplePrimaryPositionDistribution() = default;
    SamplePrimaryPositionDistribution(std::vector<siren::math::Vector3D> positions, std::vector<double> weights = {});
    virtual ~SamplePrimaryPositionDistribution() {};
private:
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("positions", positions_));
            archive(::cereal::make_nvp("weights", weights_));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("SamplePrimaryPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("positions", positions_));
            archive(::cereal::make_nvp("weights", weights_));
            if (weights_.empty()) {
                weights_.assign(positions_.size(), 1.0);
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
            throw std::runtime_error("SamplePrimaryPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
private:
    std::vector<siren::math::Vector3D> positions_;
    std::vector<double> weights_;
    std::vector<double> sum_weights_;
    double total_sum_weights_;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::SamplePrimaryPositionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::SamplePrimaryPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::SamplePrimaryPositionDistribution);

#endif // SIREN_SamplePrimaryPositionDistribution_H