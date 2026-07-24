#pragma once
#ifndef SIREN_PrimaryExternalDistribution_H
#define SIREN_PrimaryExternalDistribution_H

#include <memory>                                        // for shared_ptr
#include <array>                                         // for array
#include <set>                                           // for set
#include <string>                                        // for string
#include <tuple>                                         // for tuple
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>

#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/DistributionVariable.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }
namespace cereal { class access; }

namespace siren {
namespace distributions {

// Table-driven primary distribution: each row fixes the primary's kinematics
// (and optionally its vertex, initial position, time, and extra interaction
// parameters).
//
// Rows are sampled UNIFORMLY unless explicit sampling weights are supplied:
// the supplied parent list is taken to be the intended injection ensemble
// (a dk2nu production's importance reweighting deliberately oversamples the
// phase space it cares about, and injection respects that choice).
//
// A column named "weight" declares per-row PHYSICAL weights: row i represents
// weight_i physical primaries per unit exposure (for dk2nu tables,
// nimpwt/POT). The weights encode how to return from the supplied ensemble to
// the physical one: on the physical side of a process the distribution
// reports the physical row density through PhysicalDensity (never cancelled
// against the injection side, whose GenerationProbability stays the sampling
// density) and the weight total through the PhysicallyNormalizedDistribution
// interface. One instance shared between the injection and physical sides is
// therefore exact for uniform sampling and for any explicit sampling bias.
// Tables without a "weight" column keep the legacy semantics: flat physical
// density, no normalization, injection/physical cancellation as usual.
class PrimaryExternalDistribution : virtual public VertexPositionDistribution,
                                    virtual public PhysicallyNormalizedDistribution {
friend cereal::access;
protected:
    PrimaryExternalDistribution() {};
    void LoadInputFile(std::string const & _filename);
private:
    std::string filename;
    std::vector<std::vector<double>> input_data;
    std::vector<std::string> keys;
    std::vector<double> sampling_weights_;
    double sampling_weights_sum_ = 0;
    std::vector<double> sampling_cdf_;
    // Physical row weights from the "weight" column (empty without one).
    // Derived state: rebuilt from input_data on construction, refiltering,
    // and deserialization; never archived.
    std::vector<double> physical_weights_;
    double physical_weights_sum_ = 0;
    bool init_pos_set = false;
    bool vertex_set = false;
    bool mom_set = false;
    double emin = 0;
    std::set<DistributionVariable> set_variables_;
    mutable std::array<double, 3> _cached_position = {0.0, 0.0, 0.0};
    void ComputeSetVariables();
    void BuildSamplingCDF();
    void DerivePhysicalRowWeights();
public:
    PrimaryExternalDistribution(std::string _filename);
    PrimaryExternalDistribution(std::string _filename, double emin);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, double emin);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, std::vector<double> _sampling_weights);
    PrimaryExternalDistribution(std::vector<std::string> _keys, std::vector<std::vector<double>> _data, std::vector<double> _sampling_weights, double emin);
    PrimaryExternalDistribution(PrimaryExternalDistribution const & other) = default;
    size_t GetPhysicalNumEvents() const;
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual double PhysicalDensity(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    virtual bool PhysicalDensityDiffers() const override;
    virtual std::set<DistributionVariable> SetVariables() const override;
    virtual std::set<DistributionVariable> RequiredVariables() const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
        } else if(version == 1) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            archive(::cereal::make_nvp("SamplingWeights", sampling_weights_));
        } else {
            throw std::runtime_error("PrimaryExternalDistribution only supports version <= 1!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            ComputeSetVariables();
            DerivePhysicalRowWeights();
        } else if(version == 1) {
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
            archive(::cereal::make_nvp("Emin", emin));
            archive(::cereal::make_nvp("Keys", keys));
            archive(::cereal::make_nvp("InputData", input_data));
            archive(::cereal::make_nvp("InitPosSet", init_pos_set));
            archive(::cereal::make_nvp("VertexSet", vertex_set));
            archive(::cereal::make_nvp("MomSet", mom_set));
            archive(::cereal::make_nvp("SamplingWeights", sampling_weights_));
            ComputeSetVariables();
            DerivePhysicalRowWeights();
            BuildSamplingCDF();
        } else {
            throw std::runtime_error("PrimaryExternalDistribution only supports version <= 1!");
        }
    }
private:
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::PrimaryExternalDistribution, 1);
CEREAL_REGISTER_TYPE(siren::distributions::PrimaryExternalDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::VertexPositionDistribution, siren::distributions::PrimaryExternalDistribution);

#endif // SIREN_PrimaryExternalDistribution_H
